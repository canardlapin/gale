package gale.bench

import scala.compiletime.uninitialized

import gale.backend.PureBackend
import gale.linalg.*
import gale.platform.DoubleArray
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

private object FfmPivotedQrFixture:
  val options: QROptions = QROptions(pivoting = QRPivoting.Column)

  def matrix(rows: Int, cols: Int): DMat =
    Matrix.tabulate(rows, cols): (i, j) =>
      var state = (0x54535152L + rows * 31L + cols + i * 0x9e3779b9L) *
        6364136223846793005L + 1442695040888963407L
      var k = 0
      while k <= j do
        state = state * 6364136223846793005L + 1442695040888963407L
        k += 1
      ((state >>> 11).toDouble / (1L << 53).toDouble) * 2.0 - 1.0

private final class FfmPivotedQr(prototype: Dgeqp3Prototype):
  import FfmPivotedQrFixture.options

  private def copyColumnMajor(a: DMat): Array[Double] =
    val result = new Array[Double](a.rows * a.cols)
    var j = 0
    while j < a.cols do
      var i = 0
      while i < a.rows do
        result(j * a.rows + i) = a(i, j)
        i += 1
      j += 1
    result

  def factor(a: DMat): QR =
    val native = prototype.factorOwnedColumnMajor(copyColumnMajor(a), a.rows, a.cols)
    val packed = native.packedColumnMajor()
    val limit = math.min(a.rows, a.cols)
    val reflectors = Matrix.tabulate(a.rows, limit): (i, j) =>
      if i < j then 0.0
      else if i == j then 1.0
      else packed(j * a.rows + i)
    val r = Matrix.tabulate(a.rows, a.cols): (i, j) =>
      if i <= j then packed(j * a.rows + i) else 0.0
    val tolerance = options.rankTolerance.getOrElse(DenseDecompositions.rankToleranceFromMatrix(r))
    QR(
      reflectors,
      DoubleArray.adopt(native.tau()),
      r,
      FactorizationDiagnostics(
        rank = Some(DenseDecompositions.rankFromMatrix(r)),
        rankTolerance = Some(tolerance)
      ),
      ColumnPermutation.fromArray(native.permutation()),
      options
    )

  def nativeFootprint(a: DMat): (Int, Long) =
    val result = prototype.factorOwnedColumnMajor(copyColumnMajor(a), a.rows, a.cols)
    (result.lwork(), result.peakNativeBytes())

/** Copy-inclusive comparison of Gale's final portable column-pivoted QR and an
  * isolated FFM dgeqp3 prototype. Every native invocation includes all heap and
  * native copies, workspace query/allocation, permutation conversion, rank
  * decision, and owned Gale factor construction.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@Threads(1)
@State(Scope.Thread)
class FfmPivotedQrJmh:
  @Param(Array("512", "1024", "2048", "4096", "10000"))
  var n: Int = 0

  @Param(Array("3", "5", "6", "8", "16", "24"))
  var p: Int = 0

  private var design: DMat = uninitialized
  private var prototype: Dgeqp3Prototype = uninitialized
  private var nativeQr: FfmPivotedQr = uninitialized

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    require(p <= n, s"tall pivoted-QR court requires n >= p, got n=$n p=$p")
    require(sys.env.get("VECLIB_MAXIMUM_THREADS").contains("1"),
      "set VECLIB_MAXIMUM_THREADS=1 before starting the JMH JVM")
    design = FfmPivotedQrFixture.matrix(n, p)
    prototype = Dgeqp3Prototype.loadDefault()
    nativeQr = new FfmPivotedQr(prototype)
    require(design.qr(FfmPivotedQrFixture.options)(using PureBackend).diagnostics.rank.contains(p))
    require(nativeQr.factor(design).diagnostics.rank.contains(p))

  @TearDown(Level.Trial)
  def tearDownTrial(): Unit = prototype.close()

  @Benchmark def portablePivotedQr(blackhole: Blackhole): Unit =
    blackhole.consume(design.qr(FfmPivotedQrFixture.options)(using PureBackend))

  @Benchmark def ffmPivotedQr(blackhole: Blackhole): Unit =
    blackhole.consume(nativeQr.factor(design))

object FfmPivotedQrAssessmentMain:
  private def relativeError(actual: DMat, expected: DMat): Double =
    var differenceScale = 0.0
    var expectedScale = 0.0
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        differenceScale = math.max(differenceScale, math.abs(actual(i, j) - expected(i, j)))
        expectedScale = math.max(expectedScale, math.abs(expected(i, j)))
        j += 1
      i += 1
    if differenceScale == 0.0 then 0.0
    else if expectedScale == 0.0 then Double.PositiveInfinity
    else
      var differenceSquares = 0.0
      var expectedSquares = 0.0
      i = 0
      while i < actual.rows do
        var j = 0
        while j < actual.cols do
          val scaledDifference = (actual(i, j) - expected(i, j)) / differenceScale
          val scaledExpected = expected(i, j) / expectedScale
          differenceSquares += scaledDifference * scaledDifference
          expectedSquares += scaledExpected * scaledExpected
          j += 1
        i += 1
      (differenceScale / expectedScale) * math.sqrt(differenceSquares / expectedSquares)

  private def permuted(a: DMat, permutation: ColumnPermutation): DMat =
    Matrix.tabulate(a.rows, a.cols)((i, j) => a(i, permutation(j)))

  private def requireClose(actual: DMat, expected: DMat, tolerance: Double, label: String): Unit =
    val error = relativeError(actual, expected)
    require(error <= tolerance, s"$label relative error $error exceeds $tolerance")

  private def validateFactor(a: DMat, qr: QR, label: String): Unit =
    val reconstructed = qr.applyQ(qr.r).orThrow
    requireClose(reconstructed, permuted(a, qr.columnPermutation), 2e-11, s"$label reconstruction")
    val width = math.min(4, math.min(a.rows, math.max(1, a.cols)))
    val probe = Matrix.tabulate(a.rows, width)((i, j) => math.sin((i + 1.0) * (j + 2.0) * 0.017))
    val transformed = qr.applyQ(probe).orThrow
    requireClose(
      transformed.t.*(transformed)(using PureBackend),
      probe.t.*(probe)(using PureBackend),
      3e-11,
      s"$label orthogonality"
    )

  def main(args: Array[String]): Unit =
    require(sys.env.get("VECLIB_MAXIMUM_THREADS").contains("1"),
      "set VECLIB_MAXIMUM_THREADS=1 before starting the JVM")
    val prototype = Dgeqp3Prototype.loadDefault()
    try
      val nativeQr = new FfmPivotedQr(prototype)
      var permutationMismatches = 0
      for
        n <- Seq(512, 1024, 2048, 4096, 10000)
        p <- Seq(3, 5, 6, 8, 16, 24)
      do
        val a = FfmPivotedQrFixture.matrix(n, p)
        val pure = a.qr(FfmPivotedQrFixture.options)(using PureBackend)
        val ffm = nativeQr.factor(a)
        validateFactor(a, ffm, s"$n x $p")
        requireEquals(ffm.diagnostics.rank, pure.diagnostics.rank, s"$n x $p rank")
        if ffm.columnPermutation.toIndexSeq != pure.columnPermutation.toIndexSeq then
          permutationMismatches += 1
        val (lwork, bytes) = nativeQr.nativeFootprint(a)
        println(s"shape=$n x $p lwork=$lwork peakNativeBytes=$bytes permutationEqual=${ffm.columnPermutation.toIndexSeq == pure.columnPermutation.toIndexSeq}")

      for scale <- Seq(1e-150, 1e150) do
        val base = FfmPivotedQrFixture.matrix(128, 6)
        val scaled = Matrix.tabulate(base.rows, base.cols)((i, j) => scale * base(i, j))
        val qr = nativeQr.factor(scaled)
        validateFactor(scaled, qr, s"extreme scale $scale")
        val pureScaled = scaled.qr(FfmPivotedQrFixture.options)(using PureBackend)
        requireEquals(qr.diagnostics.rank, pureScaled.diagnostics.rank, s"extreme scale $scale rank")

      val full = FfmPivotedQrFixture.matrix(128, 6)
      val deficient = Matrix.tabulate(128, 6)((i, j) => if j == 5 then full(i, 2) else full(i, j))
      val pureDeficient = deficient.qr(FfmPivotedQrFixture.options)(using PureBackend)
      val ffmDeficient = nativeQr.factor(deficient)
      requireEquals(ffmDeficient.diagnostics.rank, pureDeficient.diagnostics.rank, "rank-deficient rank")

      val first = nativeQr.factor(full).columnPermutation.toIndexSeq
      val second = nativeQr.factor(full).columnPermutation.toIndexSeq
      requireEquals(first, second, "repeated native permutation")

      val tied = Matrix.tabulate(64, 6): (i, j) =>
        if j == 1 then math.sin((i + 1.0) * 2.0 * 0.031)
        else math.sin((i + 1.0) * (j + 2.0) * 0.031)
      val tiedPure = tied.qr(FfmPivotedQrFixture.options)(using PureBackend)
      val tiedFfm = nativeQr.factor(tied)
      validateFactor(tied, tiedFfm, "exact tied columns")
      val tiedPermutationEqual = tiedFfm.columnPermutation.toIndexSeq == tiedPure.columnPermutation.toIndexSeq

      var adversarialMismatches = 0
      var seed = 0
      while seed < 256 do
        val epsilon = math.pow(10.0, -8.0 - seed % 8)
        val nearTied = Matrix.tabulate(64, 6): (i, j) =>
          val shared = math.sin((i + 1.0) * (seed + 3.0) * 0.019)
          val perturbation = math.cos((i + 2.0) * (j + 5.0) * (seed + 1.0) * 0.007)
          if j < 5 then shared + epsilon * perturbation
          else math.sin((i + 1.0) * (seed + 11.0) * 0.023)
        val nearPure = nearTied.qr(FfmPivotedQrFixture.options)(using PureBackend)
        val nearFfm = nativeQr.factor(nearTied)
        validateFactor(nearTied, nearFfm, s"near-tied seed $seed")
        if nearFfm.columnPermutation.toIndexSeq != nearPure.columnPermutation.toIndexSeq then
          adversarialMismatches += 1
        seed += 1

      val truth = Vec(1.0, -2.0, 0.5, 3.0, -1.5, 0.25)
      val rhs = full.*(truth)(using PureBackend)
      val solved = nativeQr.factor(full).solveLeastSquares(rhs).orThrow
      require((solved - truth).norm2 <= 2e-10, "native least-squares solution")
      require(nativeQr.factor(full).solveLeastSquares(Vec.zeros(full.rows - 1)).left.exists(
        _.isInstanceOf[LinAlgError.DimensionMismatch]
      ), "native QR must preserve typed RHS dimension failure")

      println(s"library=${prototype.libraryName()}")
      println(s"permutationMismatches=$permutationMismatches")
      println(s"exactTiePermutationEqual=$tiedPermutationEqual")
      println(s"adversarialPermutationMismatches=$adversarialMismatches")
      println("status=ok")
    finally prototype.close()

  private def requireEquals[A](actual: A, expected: A, label: String): Unit =
    require(actual == expected, s"$label: $actual != $expected")
