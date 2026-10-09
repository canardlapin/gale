package gale.parity

import breeze.linalg.Axis
import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.SparseVector as BSV
import breeze.linalg.argmax as bArgmax
import breeze.linalg.argmin as bArgmin
import breeze.linalg.cholesky as bCholesky
import breeze.linalg.cond as bCond
import breeze.linalg.det as bDet
import breeze.linalg.eigSym as bEigSym
import breeze.linalg.inv as bInv
import breeze.linalg.kron as bKron
import breeze.linalg.max as bMax
import breeze.linalg.min as bMin
import breeze.linalg.norm as bNorm
import breeze.linalg.pinv as bPinv
import breeze.linalg.qr as bQr
import breeze.linalg.softmax as bLogSumExp
import breeze.linalg.sum as bSum
import breeze.linalg.svd as bSvd
import breeze.numerics.abs as bAbs
import breeze.numerics.exp as bExp
import breeze.numerics.log as bLog
import breeze.numerics.sigmoid as bSigmoid
import breeze.stats.mean as bMean

import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import gale.golden.GoldenArray
import gale.golden.GoldenCase
import gale.golden.GoldenTolerance

import scala.collection.mutable.ArrayBuffer

/** Generates `core/shared/src/test/scala/gale/golden/BreezeGoldens.scala`: Breeze 2.1.0 reference results over a
  * fixed seeded corpus, replayed by the platform-independent `BreezeGoldenSuite` on the JVM and on Scala.js.
  *
  * {{{
  * sbt breezeGoldens        # --write: regenerate the checked-in source
  * sbt breezeGoldensCheck   # --check: regenerate in memory and compare
  * }}}
  *
  * Inputs are built with plain IEEE arithmetic, `sqrt` and `StrictMath` only, so they are bit-identical on every
  * host. Breeze's references are not: netlib picks native, SIMD or scalar Java BLAS per host, and the libm
  * intrinsics differ by architecture. `--check` therefore requires an identical corpus (families, names, shapes and
  * input bits) and accepts a reference that moved by at most 64 units of the replay's forward unit. See
  * `parity/README.md`.
  */
object GenerateBreezeGoldens:

  val DefaultTarget = "core/shared/src/test/scala/gale/golden/BreezeGoldens.scala"

  /** One row-major array; `input` arrays are generator-built (host-independent), the rest are Breeze results. */
  final case class Arr(key: String, rows: Int, cols: Int, values: Array[Double], input: Boolean)
  final case class Case(family: String, name: String, arrays: Vector[Arr])

  // ---------------------------------------------------------------------------
  // Entry point
  // ---------------------------------------------------------------------------

  def main(args: Array[String]): Unit =
    val mode = args.headOption.getOrElse("--write")
    val target = resolve(args.lift(1).getOrElse(DefaultTarget))
    val cases = corpus()
    val source = render(cases)
    mode match
      case "--write" =>
        Files.createDirectories(target.getParent)
        Files.write(target, source.getBytes(StandardCharsets.UTF_8))
        val doubles = cases.iterator.flatMap(_.arrays).map(_.values.length).sum
        println(s"wrote ${cases.size} cases ($doubles doubles, ${source.length} bytes) to $target")
        println(provenance)
      case "--check" =>
        if !Files.exists(target) then fail(s"$target does not exist; run `sbt breezeGoldens`")
        val committed = new String(Files.readAllBytes(target), StandardCharsets.UTF_8)
        if committed == source then println(s"Breeze goldens are fresh (byte-identical): $target")
        else
          val problems = compare(parse(committed), cases)
          if problems.nonEmpty then
            fail(
              (s"Breeze goldens are stale ($target); run `sbt breezeGoldens` and commit the result." +:
                problems.take(25)).mkString("\n  ")
            )
          println(
            s"Breeze goldens are fresh: corpus and inputs identical; references moved only within the freshness " +
              s"band (expected when this host's BLAS/libm differs from the generating host; here: $provenance)"
          )
      case other => fail(s"unknown mode $other (expected --write or --check)")

  private def fail(message: String): Nothing =
    System.err.println(message)
    sys.exit(1)

  /** Resolves a relative path against the repository root (the forked JVM runs in `parity/`). */
  private def resolve(path: String): Path =
    val p = Paths.get(path)
    if p.isAbsolute then p
    else
      var dir = Paths.get("").toAbsolutePath
      while dir != null && !Files.exists(dir.resolve("build.sbt")) do dir = dir.getParent
      if dir == null then fail("cannot locate the repository root (no build.sbt above the working directory)")
      dir.resolve(p)

  /** Breeze's version from its jar manifest (the build pins 2.1.0). */
  private def breezeVersion: String =
    Option(classOf[BDV[?]].getPackage).flatMap(p => Option(p.getImplementationVersion)).getOrElse("2.1.0 (no manifest)")

  private def provenance: String =
    val blas = dev.ludovic.netlib.blas.BLAS.getInstance().getClass.getName
    val lapack = dev.ludovic.netlib.lapack.LAPACK.getInstance().getClass.getName
    s"breeze=$breezeVersion blas=$blas lapack=$lapack java=${System.getProperty("java.version")} " +
      s"arch=${System.getProperty("os.arch")}"

  // ---------------------------------------------------------------------------
  // Deterministic, host-independent inputs
  // ---------------------------------------------------------------------------

  private def uniform(n: Int, seed: Long, lo: Double = -1.0, hi: Double = 1.0): Array[Double] =
    val rng = new scala.util.Random(seed)
    Array.fill(n)(lo + (hi - lo) * rng.nextDouble())

  /** Row-major `rows × cols` with entries uniform in `[-1, 1)`. */
  private def general(rows: Int, cols: Int, seed: Long): Array[Double] = uniform(rows * cols, seed)

  private def diagonallyDominant(n: Int, seed: Long): Array[Double] =
    val a = general(n, n, seed)
    for i <- 0 until n do
      var s = 0.0
      for j <- 0 until n if j != i do s += math.abs(a(i * n + j))
      a(i * n + i) = s + 1.0
    a

  /** `B Bᵀ + n I`, exactly symmetric. */
  private def spd(n: Int, seed: Long): Array[Double] =
    val b = general(n, n, seed)
    val a = new Array[Double](n * n)
    for i <- 0 until n; j <- 0 to i do
      var s = 0.0
      for k <- 0 until n do s += b(i * n + k) * b(j * n + k)
      if i == j then s += n.toDouble
      a(i * n + j) = s
      a(j * n + i) = s
    a

  private def symmetric(n: Int, seed: Long): Array[Double] =
    val src = general(n, n, seed)
    val a = new Array[Double](n * n)
    for i <- 0 until n; j <- 0 to i do
      a(i * n + j) = src(i * n + j)
      a(j * n + i) = src(i * n + j)
    a

  private def hilbert(n: Int): Array[Double] =
    Array.tabulate(n * n)(k => 1.0 / (k / n + k % n + 1).toDouble)

  /** Columns of an `n × n` orthogonal matrix (row-major) from twice-applied modified Gram–Schmidt. */
  private def orthogonal(n: Int, seed: Long): Array[Double] =
    val cols = Array.tabulate(n)(j => uniform(n, seed * 1009 + j))
    for _ <- 0 until 2; j <- 0 until n do
      val v = cols(j)
      for p <- 0 until j do
        val q = cols(p)
        var d = 0.0
        for i <- 0 until n do d += q(i) * v(i)
        for i <- 0 until n do v(i) -= d * q(i)
      var s = 0.0
      for i <- 0 until n do s += v(i) * v(i)
      val inv = 1.0 / math.sqrt(s)
      for i <- 0 until n do v(i) *= inv
    Array.tabulate(n * n)(k => cols(k % n)(k / n))

  /** `U[:, :p] diag(sigma) V[:, :p]ᵀ`, `p = sigma.length ≤ min(rows, cols)`. */
  private def withSingularValues(rows: Int, cols: Int, sigma: Array[Double], seed: Long): Array[Double] =
    val u = orthogonal(rows, seed)
    val v = orthogonal(cols, seed + 7)
    Array.tabulate(rows * cols) { k =>
      val i = k / cols
      val j = k % cols
      var s = 0.0
      for p <- sigma.indices do s += u(i * rows + p) * sigma(p) * v(j * cols + p)
      s
    }

  /** `Q diag(spectrum) Qᵀ`, built exactly symmetric. */
  private def symmetricWithSpectrum(spectrum: Array[Double], seed: Long): Array[Double] =
    val n = spectrum.length
    val q = orthogonal(n, seed)
    val a = new Array[Double](n * n)
    for i <- 0 until n; j <- 0 to i do
      var s = 0.0
      for k <- 0 until n do s += q(i * n + k) * spectrum(k) * q(j * n + k)
      a(i * n + j) = s
      a(j * n + i) = s
    a

  /** `n` values from 1 down to `1/kappa`, geometrically spaced (StrictMath for host independence). */
  private def geometric(n: Int, kappa: Double): Array[Double] =
    if n == 1 then Array(1.0) else Array.tabulate(n)(i => StrictMath.pow(kappa, -i.toDouble / (n - 1)))

  /** `count` distinct sorted indices in `[0, length)`. */
  private def sparseIndices(length: Int, count: Int, seed: Long): Array[Int] =
    new scala.util.Random(seed).shuffle((0 until length).toVector).take(count).sorted.toArray

  // ---------------------------------------------------------------------------
  // Breeze bridges
  // ---------------------------------------------------------------------------

  private def bm(rows: Int, cols: Int, a: Array[Double]): BDM[Double] =
    BDM.tabulate(rows, cols)((i, j) => a(i * cols + j))

  private def bv(a: Array[Double]): BDV[Double] = BDV(a.clone())

  private def fromBm(m: BDM[Double]): Array[Double] =
    Array.tabulate(m.rows * m.cols)(k => m(k / m.cols, k % m.cols))

  private def fromBv(v: BDV[Double]): Array[Double] = Array.tabulate(v.length)(v(_))

  /** The input arrays of each family; every other array of a case is a Breeze reference. */
  private val InputKeys: Map[String, Set[String]] = Map(
    "matvec" -> Set("A", "x", "xt"),
    "matmul" -> Set("A", "B", "D"),
    "vector" -> Set("x", "y", "alpha"),
    "lu" -> Set("A", "b"),
    "cholesky" -> Set("A", "b"),
    "qr" -> Set("A", "b"),
    "eigsym" -> Set("A"),
    "svd" -> Set("A"),
    "pinv" -> Set("A"),
    "kron" -> Set("A", "B"),
    "reduce-vector" -> Set("x"),
    "reduce-matrix" -> Set("A"),
    "numerics" -> Set("x", "p"),
    "sparse" -> Set("xi", "xv", "yi", "yv", "d")
  )

  private final class CaseBuilder(family: String, name: String):
    private val out = Vector.newBuilder[Arr]
    def mat(key: String, rows: Int, cols: Int, a: Array[Double]): this.type =
      require(a.length == rows * cols, s"$name/$key")
      out += Arr(key, rows, cols, a.clone(), InputKeys(family)(key)); this
    def mat(key: String, m: BDM[Double]): this.type = mat(key, m.rows, m.cols, fromBm(m))
    def vec(key: String, a: Array[Double]): this.type = mat(key, a.length, 1, a)
    def vec(key: String, v: BDV[Double]): this.type = vec(key, fromBv(v))
    def scalar(key: String, x: Double): this.type = mat(key, 1, 1, Array(x))
    def build: Case = Case(family, name, out.result())

  private def kappa2(m: BDM[Double]): Double =
    val s = bSvd(m).singularValues
    s(0) / s(s.length - 1)

  // ---------------------------------------------------------------------------
  // Corpus
  // ---------------------------------------------------------------------------

  def corpus(): Vector[Case] =
    val cases = ArrayBuffer.empty[Case]
    def add(family: String, name: String)(fill: CaseBuilder => Unit): Unit =
      val b = new CaseBuilder(family, s"$family/$name")
      fill(b)
      cases += b.build

    // Dense level-2 / level-3 products.
    for (m, n, seed) <- Seq((1, 1, 11L), (7, 5, 12L), (33, 33, 13L), (64, 64, 14L), (100, 100, 15L)) do
      add("matvec", s"${m}x$n") { c =>
        val a = general(m, n, seed); val x = uniform(n, seed + 1); val xt = uniform(m, seed + 2)
        val ba = bm(m, n, a)
        c.mat("A", m, n, a).vec("x", x).vec("xt", xt).vec("Ax", ba * bv(x)).vec("Atxt", ba.t * bv(xt))
      }
    for (m, k, n, seed) <- Seq((1, 1, 1, 21L), (5, 3, 4, 22L), (17, 31, 13, 23L), (48, 48, 48, 24L)) do
      add("matmul", s"${m}x${k}x$n") { c =>
        val a = general(m, k, seed); val b = general(k, n, seed + 1); val d = general(m, n, seed + 2)
        val ba = bm(m, k, a)
        c.mat("A", m, k, a).mat("B", k, n, b).mat("D", m, n, d)
        c.mat("AB", ba * bm(k, n, b)).mat("AtD", ba.t * bm(m, n, d))
      }
    for (n, seed) <- Seq((1, 31L), (10, 32L), (257, 33L), (1000, 34L)) do
      add("vector", s"n$n") { c =>
        val x = uniform(n, seed); val y = uniform(n, seed + 1); val alpha = -1.75
        c.vec("x", x).vec("y", y).scalar("alpha", alpha)
        c.scalar("dot", bv(x).dot(bv(y))).vec("axpy", bv(x) * alpha + bv(y))
      }

    // LU: solve, det, inverse.
    val luInputs: Seq[(String, Int, Array[Double])] =
      Seq(1, 2, 8, 33, 64).map(n => (s"general-n$n", n, general(n, n, 40L + n))) ++
        Seq(("dominant-n16", 16, diagonallyDominant(16, 41L))) ++
        Seq(1e2, 1e6, 1e10).map(k => (f"kappa$k%.0e-n24", 24, withSingularValues(24, 24, geometric(24, k), 42L)))
    for (name, n, a) <- luInputs do
      add("lu", name) { c =>
        val ba = bm(n, n, a); val b = uniform(n, 43L + n)
        c.mat("A", n, n, a).vec("b", b).scalar("kappa", kappa2(ba))
        c.vec("x", ba \ bv(b)).scalar("det", bDet(ba))
        if n <= 33 then c.mat("inv", bInv(ba))
      }

    // Cholesky: lower factor and solve.
    val cholInputs: Seq[(String, Int, Array[Double])] =
      Seq(1, 5, 16, 64).map(n => (s"spd-n$n", n, spd(n, 50L + n))) ++
        Seq(("hilbert-n8", 8, hilbert(8)), ("kappa1e8-n20", 20, symmetricWithSpectrum(geometric(20, 1e8), 51L)))
    for (name, n, a) <- cholInputs do
      add("cholesky", name) { c =>
        val ba = bm(n, n, a); val b = uniform(n, 52L + n)
        c.mat("A", n, n, a).vec("b", b).scalar("kappa", kappa2(ba))
        c.mat("L", bCholesky(ba)).vec("x", ba \ bv(b))
      }

    // Householder QR: sign-normalized R and least squares.
    val qrInputs: Seq[(String, Int, Int, Array[Double])] =
      Seq((1, 1), (6, 4), (40, 17), (64, 32)).map((m, n) => (s"${m}x$n", m, n, general(m, n, 60L + m + n))) ++
        Seq(("kappa1e6-30x10", 30, 10, withSingularValues(30, 10, geometric(10, 1e6), 61L)))
    for (name, m, n, a) <- qrInputs do
      add("qr", name) { c =>
        val ba = bm(m, n, a); val b = uniform(m, 62L + m)
        val r = bQr(ba).r(0 until n, ::).copy
        for i <- 0 until n if r(i, i) < 0.0 do r(i, ::) := r(i, ::) * -1.0
        val x = ba \ bv(b)
        c.mat("A", m, n, a).vec("b", b).scalar("kappa", kappa2(ba))
        c.mat("R", r).vec("x", x).scalar("residual", bNorm(bv(b) - ba * x))
      }

    // Symmetric eigen: ascending values and (optionally) vectors.
    val eigInputs: Seq[(String, Int, Array[Double], Boolean)] =
      Seq(1, 4, 16, 48).map(n => (s"random-n$n", n, symmetric(n, 70L + n), true)) ++
        Seq(
          ("random-n96-values", 96, symmetric(96, 71L), false),
          (
            "clustered-n12",
            12,
            symmetricWithSpectrum(Array(-2.0, -2.0, -2.0, 0.5, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 7.0, 9.0), 72L),
            true
          ),
          ("graded1e10-n20", 20, symmetricWithSpectrum(geometric(20, 1e10), 73L), true)
        )
    for (name, n, a, withVectors) <- eigInputs do
      add("eigsym", name) { c =>
        val es = bEigSym(bm(n, n, a))
        c.mat("A", n, n, a).vec("values", es.eigenvalues)
        if withVectors then c.mat("vectors", es.eigenvectors)
      }

    // Singular values.
    val svdInputs: Seq[(String, Int, Int, Array[Double])] =
      Seq((1, 1), (8, 5), (5, 8), (30, 30), (100, 60)).map((m, n) => (s"${m}x$n", m, n, general(m, n, 80L + m * n))) ++
        Seq(("kappa1e8-20x12", 20, 12, withSingularValues(20, 12, geometric(12, 1e8), 81L)))
    for (name, m, n, a) <- svdInputs do
      add("svd", name) { c =>
        c.mat("A", m, n, a).vec("sigma", bSvd(bm(m, n, a)).singularValues)
      }

    // Full-rank pseudoinverse.
    val pinvInputs: Seq[(String, Int, Int, Array[Double])] =
      Seq((12, 5), (40, 20), (10, 10)).map((m, n) => (s"${m}x$n", m, n, general(m, n, 90L + m + n))) ++
        Seq(("kappa1e4-25x8", 25, 8, withSingularValues(25, 8, geometric(8, 1e4), 91L)))
    for (name, m, n, a) <- pinvInputs do
      add("pinv", name) { c =>
        val ba = bm(m, n, a)
        c.mat("A", m, n, a).scalar("kappa", kappa2(ba)).mat("P", bPinv(ba))
      }

    for (am, an, cm, cn, seed) <- Seq((1, 1, 3, 3, 100L), (2, 3, 3, 2, 101L), (4, 4, 5, 3, 102L)) do
      add("kron", s"${am}x$an-${cm}x$cn") { c =>
        val a = general(am, an, seed); val b = general(cm, cn, seed + 1)
        c.mat("A", am, an, a).mat("B", cm, cn, b).mat("K", bKron(bm(am, an, a), bm(cm, cn, b)))
      }

    // Reductions and norms.
    val reduceVectors: Seq[(String, Array[Double])] =
      Seq(1, 7, 64, 1000).map(n => (s"n$n", uniform(n, 110L + n))) ++
        Seq(("wide-range-n200", {
          val mag = uniform(200, 111L, -8.0, 8.0); val sign = uniform(200, 112L)
          Array.tabulate(200)(i => math.signum(sign(i)) * StrictMath.pow(10.0, mag(i)))
        }))
    for (name, x) <- reduceVectors do
      add("reduce-vector", name) { c =>
        val v = bv(x)
        c.vec("x", x).scalar("sum", bSum(v)).scalar("mean", bMean(v))
        c.scalar("max", bMax(v)).scalar("min", bMin(v))
        c.scalar("argmax", bArgmax(v).toDouble).scalar("argmin", bArgmin(v).toDouble)
        c.scalar("norm1", bNorm(v, 1.0)).scalar("norm2", bNorm(v)).scalar("normInf", bNorm(v, Double.PositiveInfinity))
      }
    for (m, n, seed) <- Seq((1, 1, 120L), (5, 7, 121L), (40, 30, 122L)) do
      add("reduce-matrix", s"${m}x$n") { c =>
        val a = general(m, n, seed); val ba = bm(m, n, a)
        val (ai, aj) = bArgmax(ba)
        c.mat("A", m, n, a).scalar("sum", bSum(ba)).scalar("mean", bMean(ba))
        c.scalar("max", bMax(ba)).scalar("min", bMin(ba))
        c.scalar("argmaxRow", ai.toDouble).scalar("argmaxCol", aj.toDouble)
        c.scalar("frobenius", bNorm(ba.toDenseVector))
        c.scalar("norm1", bMax(bSum(bAbs(ba), Axis._0))).scalar("normInf", bMax(bSum(bAbs(ba), Axis._1)))
      }

    // Elementwise numerics; Breeze's `softmax(v)` is log-sum-exp.
    for (name, n, lo, hi, seed) <- Seq(
        ("n1", 1, -30.0, 30.0, 130L),
        ("n16", 16, -30.0, 30.0, 131L),
        ("n257", 257, -30.0, 30.0, 132L),
        ("large-n64", 64, -700.0, 700.0, 133L)
      )
    do
      add("numerics", name) { c =>
        val x = uniform(n, seed, lo, hi)
        val p = uniform(n, seed + 1, -20.0, 20.0).map(StrictMath.exp)
        val lse = bLogSumExp(bv(x))
        c.vec("x", x).vec("p", p)
        c.vec("exp", bExp(bv(x))).vec("log", bLog(bv(p))).vec("sigmoid", bSigmoid(bv(x)))
        c.scalar("logSumExp", lse).vec("softmax", bExp(bv(x) - lse))
      }

    // SparseVector arithmetic and norms.
    for (n, nx, ny, seed) <- Seq((1, 1, 1, 140L), (50, 10, 20, 141L), (1000, 40, 300, 142L), (64, 64, 32, 143L)) do
      add("sparse", s"n$n-nnz$nx-$ny") { c =>
        val xi = sparseIndices(n, nx, seed); val yi = sparseIndices(n, ny, seed + 1)
        val xv = uniform(nx, seed + 2); val yv = uniform(ny, seed + 3); val d = uniform(n, seed + 4)
        val sx = new BSV(xi.clone(), xv.clone(), n); val sy = new BSV(yi.clone(), yv.clone(), n)
        c.vec("xi", xi.map(_.toDouble)).vec("xv", xv).vec("yi", yi.map(_.toDouble)).vec("yv", yv).vec("d", d)
        c.scalar("dot", sx.dot(sy)).scalar("dotDense", sx.dot(bv(d)))
        c.vec("plus", (sx + sy).toDenseVector).vec("minus", (sx - sy).toDenseVector)
        c.scalar("sum", bSum(sx)).scalar("max", bMax(sx)).scalar("min", bMin(sx))
        c.scalar("norm1", bNorm(sx, 1.0)).scalar("norm2", bNorm(sx))
        c.scalar("normInf", bNorm(sx, Double.PositiveInfinity))
      }

    cases.toVector

  // ---------------------------------------------------------------------------
  // Rendering
  // ---------------------------------------------------------------------------

  private val DoublesPerLiteral = 64
  private val DoublesPerPart = 12000

  private def hex(x: Double): String =
    val s = java.lang.Long.toHexString(java.lang.Double.doubleToRawLongBits(x))
    ("0" * (16 - s.length)) + s

  def render(cases: Vector[Case]): String =
    val parts = ArrayBuffer(ArrayBuffer.empty[Case])
    var inPart = 0
    for c <- cases do
      val size = c.arrays.map(_.values.length).sum
      if inPart > 0 && inPart + size > DoublesPerPart then
        parts += ArrayBuffer.empty[Case]
        inPart = 0
      parts.last += c
      inPart += size
    val sb = new StringBuilder
    sb ++= "package gale.golden\n\n"
    sb ++= "import gale.golden.GoldenSyntax.in\n"
    sb ++= "import gale.golden.GoldenSyntax.ref\n\n"
    sb ++= "// GENERATED by parity/src/test/scala/gale/parity/GenerateBreezeGoldens.scala (`sbt breezeGoldens`).\n"
    sb ++= "// Do not edit by hand. Values are IEEE-754 bit patterns (16 hex digits per double, row-major).\n"
    sb ++= s"// Generated with: $provenance\n\n"
    sb ++= "/** Breeze 2.1.0 reference results over a fixed seeded corpus; replayed by `BreezeGoldenSuite`. */\n"
    sb ++= "object BreezeGoldens:\n"
    sb ++= "  lazy val cases: IndexedSeq[GoldenCase] =\n"
    sb ++= parts.indices.map(p => s"BreezeGoldensPart$p.cases").mkString("    ", " ++\n      ", "\n")
    for (part, p) <- parts.zipWithIndex do
      sb ++= s"\nprivate object BreezeGoldensPart$p:\n"
      sb ++= "  val cases: IndexedSeq[GoldenCase] = IndexedSeq(\n"
      for c <- part do
        sb ++= s"""    GoldenCase("${c.family}", "${c.name}", IndexedSeq(\n"""
        for a <- c.arrays do
          sb ++= s"""      ${if a.input then "in" else "ref"}("${a.key}", ${a.rows}, ${a.cols}"""
          for chunk <- a.values.grouped(DoublesPerLiteral) do
            sb ++= ",\n        \""
            chunk.foreach(x => sb ++= hex(x))
            sb += '"'
          sb ++= "),\n"
        sb ++= "    )),\n"
      sb ++= "  )\n"
    sb.result()

  // ---------------------------------------------------------------------------
  // Freshness check
  // ---------------------------------------------------------------------------

  private val Token = """GoldenCase\("([^"]+)", "([^"]+)"|(in|ref)\("([^"]+)", (\d+), (\d+)|"([0-9a-f]+)"""".r

  /** Parses a rendered golden source back into cases. */
  def parse(source: String): Vector[Case] =
    val cases = ArrayBuffer.empty[Case]
    var family = ""
    var name = ""
    val arrays = ArrayBuffer.empty[Arr]
    var key: String = null
    var input = false
    var rows = 0
    var cols = 0
    val hexes = new StringBuilder
    def flushArray(): Unit =
      if key != null then
        val h = hexes.result()
        arrays += Arr(
          key,
          rows,
          cols,
          Array.tabulate(h.length / 16)(i =>
            java.lang.Double.longBitsToDouble(java.lang.Long.parseUnsignedLong(h.substring(16 * i, 16 * i + 16), 16))
          ),
          input
        )
        key = null
        hexes.clear()
    def flushCase(): Unit =
      flushArray()
      if name.nonEmpty then cases += Case(family, name, arrays.toVector)
      arrays.clear()
    for m <- Token.findAllMatchIn(source) do
      if m.group(1) != null then
        flushCase(); family = m.group(1); name = m.group(2)
      else if m.group(3) != null then
        flushArray(); input = m.group(3) == "in"; key = m.group(4); rows = m.group(5).toInt; cols = m.group(6).toInt
      else hexes ++= m.group(7)
    flushCase()
    cases.toVector

  /** Freshness band, in units of the replay's own forward unit (`gale.golden.GoldenTolerance`, where the replay
    * allows `C = 8`): a regenerated reference may move by `FreshnessUnits` units across BLAS/libm hosts. Exact
    * operations have unit `0` and must not move at all.
    */
  private val FreshnessUnits = 64.0

  private def toGolden(c: Case): GoldenCase =
    GoldenCase(c.family, c.name, c.arrays.map(a => new GoldenArray(a.key, a.rows, a.cols, a.input, a.values.map(hex).mkString)))

  def compare(committed: Vector[Case], fresh: Vector[Case]): Vector[String] =
    val problems = Vector.newBuilder[String]
    def shape(c: Case) = (c.family, c.arrays.map(a => (a.key, a.rows, a.cols, a.input)))
    if committed.map(_.name) != fresh.map(_.name) then
      problems += s"case list differs: committed ${committed.size} cases, generator ${fresh.size}"
    else
      for (old, now) <- committed.zip(fresh) do
        if shape(old) != shape(now) then problems += s"${now.name}: structure differs"
        else
          val golden = toGolden(old)
          for (o, f) <- old.arrays.zip(now.arrays) do
            if o.input then
              if !o.values.indices.forall(i => java.lang.Double.compare(o.values(i), f.values(i)) == 0) then
                problems += s"${now.name}/${o.key}: input bits differ"
            else if o.key == "vectors" then
              val n = o.rows
              val ratio = GoldenTolerance.projectorRatio(
                (i, j) => o.values(i * n + j),
                (i, j) => f.values(i * n + j),
                old.arrays.find(_.key == "values").get.values
              )
              if !(ratio <= FreshnessUnits) then problems += s"${now.name}/vectors: subspaces moved ($ratio units)"
            else
              val unit = GoldenTolerance.forwardUnit(golden, o.key)
              for i <- o.values.indices do
                val (u, v) = (o.values(i), f.values(i))
                if java.lang.Double.compare(u, v) != 0 && !(math.abs(u - v) <= FreshnessUnits * unit(i)) then
                  problems += s"${now.name}/${o.key}($i): $u became $v (band ${FreshnessUnits * unit(i)})"
    problems.result()
