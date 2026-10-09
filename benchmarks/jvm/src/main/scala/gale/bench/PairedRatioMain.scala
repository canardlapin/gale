package gale.bench

import breeze.linalg.{SparseVector as BSV, *}
import gale.bench.BreezeBenchData.*
import gale.kernel.DoubleKernels
import gale.linalg.{Axis, Numerics}
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.sparse.SparseVector
import gale.syntax.all.*

/** Interleaved paired timing for a loaded, heterogeneous (P/E-core) host where
  * independent JMH trials drift by +-50%. Each case alternates short timed
  * batches of its two sides (order flipped every round) in one JVM, so both
  * sides see the same load window; the reported figure is the median of the
  * per-round speed ratios `time(b) / time(a)` with its interquartile range.
  *
  * `runMain gale.bench.PairedRatioMain [filter]`; lane A when the JVM lacks
  * `--add-modules=jdk.incubator.vector` (Breeze then resolves Java11BLAS).
  */
object PairedRatioMain:
  final case class Case(name: String, a: () => Double, b: () => Double)

  @volatile private var sink = 0.0

  def main(args: Array[String]): Unit =
    val filter = args.headOption.getOrElse("")
    val rounds = sys.props.getOrElse("rounds", "61").toInt
    println(s"blas=${dev.ludovic.netlib.blas.BLAS.getInstance().getClass.getName}")
    for c <- cases if c.name.contains(filter) do run(c, rounds)

  private def cases: Seq[Case] =
    val out = Seq.newBuilder[Case]
    for n <- Seq(1024, 65536, 1048576) do
      val xd = vectorData(n, 1L)
      val yd = vectorData(n, 2L)
      val gx = galeVector(xd)
      val gy = galeVector(yd)
      val bx = breezeVector(xd)
      val by = breezeVector(yd)
      val gw = gy.mutableCopy
      val bw = by.copy
      val rx = DoubleArray.adopt(xd.clone())
      val ry = DoubleArray.adopt(yd.clone())
      val ro = DoubleArray.alloc(n)
      out += Case(s"axpy gale/breeze n=$n", () => { gw.axpyInPlace(1.5, gx); gw.axpyInPlace(-1.5, gx); gw(0) }, () => {
        axpy(1.5, bx, bw); axpy(-1.5, bx, bw); bw(0)
      })
      out += Case(s"axpy cur/ref n=$n", () => { DoubleKernels.daxpy(n, 1.5, rx, 0, 1, ry, 0, 1); DoubleKernels.daxpy(n, -1.5, rx, 0, 1, ry, 0, 1); ry(0) }, () => {
        L1HillClimbRef.refDaxpy(n, 1.5, rx, 0, 1, ry, 0, 1); L1HillClimbRef.refDaxpy(n, -1.5, rx, 0, 1, ry, 0, 1); ry(0)
      })
      out += Case(s"dot gale/breeze n=$n", () => gx.dot(gy), () => bx.dot(by))
      out += Case(s"dot cur/ref n=$n", () => DoubleKernels.ddot(n, rx, 0, 1, ry, 0, 1), () => L1HillClimbRef.refDdot(n, rx, 0, 1, ry, 0, 1))
      out += Case(s"norm gale/breeze n=$n", () => gy.norm2, () => breeze.linalg.norm(by))
      out += Case(s"max gale/breeze n=$n", () => gx.max, () => max(bx))
      out += Case(s"max cur/ref n=$n", () => DoubleKernels.dmax(n, rx, 0, 1), () => rx(L1HillClimbRef.refMaxIndex(n, rx, 0, 1)))
      out += Case(s"exp gale/breeze n=$n", () => Numerics.exp(gx)(0), () => { val r: breeze.linalg.DenseVector[Double] = breeze.numerics.exp(bx); r(0) })
      out += Case(s"sigmoid gale/breeze n=$n", () => Numerics.sigmoid(gx)(0), () => { val r: breeze.linalg.DenseVector[Double] = breeze.numerics.sigmoid(bx); r(0) })
      out += Case(s"sigmoid cur/ref n=$n", () => { DoubleKernels.dsigmoidInto(n, rx, 0, 1, ro, 0, 1); ro(1) }, () => {
        L1HillClimbRef.refSigmoidInto(n, rx, ro); ro(1)
      })
      out += Case(s"logSumExp gale/breeze n=$n", () => Numerics.logSumExp(gx), () => softmax(bx))
    for n <- Seq(256, 1024) do
      val ad = matrixData(n, n, 5100L)
      val bd = matrixData(n, n, 5200L)
      val gA = galeMatrix(ad)
      val gB = galeMatrix(bd)
      val bA = breezeMatrix(ad)
      val bB = breezeMatrix(bd)
      out += Case(s"add gale/breeze n=$n", () => (gA + gB)(0, 0), () => { val r: breeze.linalg.DenseMatrix[Double] = bA + bB; r(0, 0) })
      out += Case(s"sub gale/breeze n=$n", () => (gA - gB)(0, 0), () => { val r: breeze.linalg.DenseMatrix[Double] = bA - bB; r(0, 0) })
      out += Case(s"hadamard gale/breeze n=$n", () => (gA.pointwise * gB)(0, 0), () => { val r: breeze.linalg.DenseMatrix[Double] = bA *:* bB; r(0, 0) })
      out += Case(s"maxCols gale/breeze n=$n", () => gA.max(Axis.Cols)(0), () => { val r: breeze.linalg.Transpose[breeze.linalg.DenseVector[Double]] = max(bA(::, breeze.linalg.*)); r.inner(0) })
    for nnz <- Seq(1000, 10000) do
      val (xi, xv) = sparseEntries(100000, nnz, 4100L)
      val (yi, yv) = sparseEntries(100000, nnz, 4200L)
      val gsx = SparseVector.tryFromArrays(100000, xi, xv).fold(e => throw e, identity)
      val gsy = SparseVector.tryFromArrays(100000, yi, yv).fold(e => throw e, identity)
      val bsx = new BSV[Double](xi.clone(), xv.clone(), 100000)
      val bsy = new BSV[Double](yi.clone(), yv.clone(), 100000)
      out += Case(s"svAdd gale/breeze nnz=$nnz", () => (gsx + gsy).activeSize.toDouble, () => { val r: BSV[Double] = bsx + bsy; r.activeSize.toDouble })
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
    // Calibrate a batch to ~20 ms of side `a`, then warm both sides ~1 s each.
    var calls = 1
    while time(c.a, calls) < 20_000_000L do calls *= 2
    val warmEnd = System.nanoTime() + 2_000_000_000L
    while System.nanoTime() < warmEnd do
      time(c.a, calls)
      time(c.b, calls)
    val ratios = Array.tabulate(rounds) { r =>
      if r % 2 == 0 then
        val ta = time(c.a, calls)
        val tb = time(c.b, calls)
        tb.toDouble / ta
      else
        val tb = time(c.b, calls)
        val ta = time(c.a, calls)
        tb.toDouble / ta
    }
    java.util.Arrays.sort(ratios)
    val med = ratios(rounds / 2)
    val q1 = ratios(rounds / 4)
    val q3 = ratios(3 * rounds / 4)
    val verdict = if q1 > 1.0 then "a faster" else if q3 < 1.0 then "b faster" else "tie"
    println(f"${c.name}%-36s ${med}%5.2fx  IQR [${q1}%5.2f, ${q3}%5.2f]  $verdict")
