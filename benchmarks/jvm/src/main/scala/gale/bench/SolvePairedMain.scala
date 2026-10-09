package gale.bench

import breeze.linalg.{DenseMatrix as BDM, inv}
import dev.ludovic.netlib.lapack.LAPACK
import gale.backend.Backend
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import org.netlib.util.intW

/** Interleaved paired timing for the dense solve paths (multi-RHS LU and Cholesky
  * solves, inverse, and the small-n factorizations), gale vs Breeze, with the
  * same alternating-batch protocol as [[PairedRatioMain]]: the figure is the
  * median per-round speed ratio `time(breeze) / time(gale)` (above 1 means gale
  * is faster) with its interquartile range.
  *
  * Run with a plain `java -cp <benchmarksJVM runtime classpath>` (no
  * `--add-modules=jdk.incubator.vector`) so Breeze resolves its scalar lane-A
  * BLAS. `rounds` defaults to 41; the first argument filters case names.
  */
object SolvePairedMain:
  final case class Case(name: String, gale: () => Double, breeze: () => Double)

  @volatile private var sink = 0.0
  private given Backend = Backend.pure

  def main(args: Array[String]): Unit =
    val filter = args.headOption.getOrElse("")
    val rounds = sys.props.getOrElse("rounds", "41").toInt
    println(s"blas=${dev.ludovic.netlib.blas.BLAS.getInstance().getClass.getName}")
    for c <- cases if c.name.contains(filter) do run(c, rounds)

  private def first(m: Either[LinAlgError, DMat]): Double = m.fold(e => throw e, _(0, 0))

  private def cases: Seq[Case] =
    val out = Seq.newBuilder[Case]
    for (n, k) <- Seq((256, 64), (64, 16)) do
      val aData = diagonallyDominant(n, 1700L)
      val sData = spd(n, 1800L)
      val bData = matrixData(n, k, 1900L)
      val gA = galeMatrix(aData)
      val gS = galeMatrix(sData)
      val gB = galeMatrix(bData)
      val bA = breezeMatrix(aData)
      val bS = breezeMatrix(sData)
      val bB = breezeMatrix(bData)
      val lapack = LAPACK.getInstance()
      out += Case(s"luSolve n=$n k=$k", () => first(gA.solve(gB)), () => { val r: BDM[Double] = bA \ bB; r(0, 0) })
      out += Case(
        s"cholSolve n=$n k=$k",
        () => first(gS.cholesky.flatMap(_.solve(gB))),
        () =>
          val a = bS.data.clone()
          val x = bB.data.clone()
          val info = new intW(0)
          lapack.dpotrf("L", n, a, n, info)
          lapack.dpotrs("L", n, k, a, n, x, n, info)
          x(0)
      )
      // Solve phase alone: factors precomputed on both sides.
      val gLu = gA.lu.fold(e => throw e, identity)
      val gCh = gS.cholesky.fold(e => throw e, identity)
      val luA = bA.data.clone()
      val ipiv = new Array[Int](n)
      lapack.dgetrf(n, n, luA, n, ipiv, new intW(0))
      val chA = bS.data.clone()
      lapack.dpotrf("L", n, chA, n, new intW(0))
      out += Case(
        s"luSolveOnly n=$n k=$k",
        () => first(gLu.solve(gB)),
        () =>
          val x = bB.data.clone()
          lapack.dgetrs("N", n, k, luA, n, ipiv, x, n, new intW(0))
          x(0)
      )
      out += Case(
        s"cholSolveOnly n=$n k=$k",
        () => first(gCh.solve(gB)),
        () =>
          val x = bB.data.clone()
          lapack.dpotrs("L", n, k, chA, n, x, n, new intW(0))
          x(0)
      )
    // Live solve vs the pre-dtrsmLeft loops ([[SolveRef]]); "gale" is the live side.
    for (n, k) <- Seq((8, 1), (16, 16), (64, 16), (256, 64)) do
      val rng = new scala.util.Random(11L)
      val gA = Matrix.tabulate(n, n)((i, j) => if i == j then n.toDouble else rng.nextDouble() * 2.0 - 1.0)
      val gS = (gA * gA.t) + Matrix.eye(n) * n.toDouble
      val gB = Matrix.tabulate(n, k)((_, _) => rng.nextDouble() * 2.0 - 1.0)
      val lu = gA.lu.fold(e => throw e, identity)
      val ch = gS.cholesky.fold(e => throw e, identity)
      out += Case(s"cur/ref luSolveOnly n=$n k=$k", () => first(lu.solve(gB)), () => SolveRef.luSolve(lu, gB)(0, 0))
      out += Case(s"cur/ref cholSolveOnly n=$n k=$k", () => first(ch.solve(gB)), () => SolveRef.choleskySolve(ch, gB)(0, 0))
    for n <- Seq(4, 16, 64, 256, 512) do
      val aData = diagonallyDominant(n, 2100L)
      val gA = galeMatrix(aData)
      val bA = breezeMatrix(aData)
      val gEye = Matrix.eye(n)
      out += Case(s"inv n=$n", () => first(gA.inverse), () => { val r: BDM[Double] = inv(bA); r(0, 0) })
      out += Case(s"invSolveI n=$n", () => first(gA.solve(gEye)), () => { val r: BDM[Double] = inv(bA); r(0, 0) })
    out.result()

  private def time(f: () => Double, calls: Int): Long =
    val t0 = System.nanoTime()
    var k = 0
    var s = 0.0
    while k < calls do
      s += f()
      k += 1
    val t = System.nanoTime() - t0
    sink += s
    t

  private def run(c: Case, rounds: Int): Unit =
    var calls = 1
    while time(c.gale, calls) < 20_000_000L do calls *= 2
    val warmEnd = System.nanoTime() + 2_000_000_000L
    while System.nanoTime() < warmEnd do
      time(c.gale, calls)
      time(c.breeze, calls)
    val ratios = Array.tabulate(rounds) { r =>
      if r % 2 == 0 then
        val tg = time(c.gale, calls)
        val tb = time(c.breeze, calls)
        tb.toDouble / tg
      else
        val tb = time(c.breeze, calls)
        val tg = time(c.gale, calls)
        tb.toDouble / tg
    }
    java.util.Arrays.sort(ratios)
    val med = ratios(rounds / 2)
    val q1 = ratios(rounds / 4)
    val q3 = ratios(3 * rounds / 4)
    val verdict = if q1 > 1.0 then "gale faster" else if q3 < 1.0 then "breeze faster" else "tie"
    println(f"${c.name}%-28s ${med}%5.2fx  IQR [${q1}%5.2f, ${q3}%5.2f]  $verdict")
