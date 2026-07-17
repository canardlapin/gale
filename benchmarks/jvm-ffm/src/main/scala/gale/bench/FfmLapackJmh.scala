package gale.bench

import gale.backend.PureBackend
import gale.backend.jvm.blas.{FfmBlasBackend, FfmBlasThresholds}
import gale.linalg.{DMat, DVec, Matrix, Vec}
import gale.spectral.{Eigen, EigenSelection, EigenVectors}

import org.openjdk.jmh.annotations.*
import java.util.concurrent.TimeUnit

/** Heap-copy-inclusive LAPACK crossover through Gale's typed factorization
  * facades. Setup constructs deterministic inputs outside the timed region;
  * every benchmark includes result allocation and native copy-in/out.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class FfmLapackJmh:
  @Param(Array("32", "64", "128", "256", "512"))
  var n: Int = _

  private var square: DMat = _
  private var spd: DMat = _
  private var tall: DMat = _
  private var nativeBackend: FfmBlasBackend = _

  @Setup(Level.Trial)
  def setup(): Unit =
    square = Matrix.tabulate(n, n)((i, j) =>
      val base = (((i.toLong * 37 + j.toLong * 19 + 11) % 101).toDouble - 50.0) / 51.0
      if i == j then base + n.toDouble else base
    )
    spd = Matrix.tabulate(n, n)((i, j) =>
      if i == j then n.toDouble + 2.0
      else 1.0 / (1.0 + math.abs(i - j).toDouble)
    )
    tall = Matrix.tabulate(2 * n, n)((i, j) =>
      val base = math.sin((i + 1).toDouble * 0.013 + (j + 2).toDouble * 0.021)
      if i == j then base + 2.0 else base
    )
    nativeBackend = FfmBlasBackend.load(thresholds = Some(FfmBlasThresholds(
      nativeLuMinSize = 0,
      nativeCholeskyMinSize = 0,
      nativeQrMinSize = 0
    ))).fold(throw _, identity)

  @TearDown(Level.Trial)
  def tearDown(): Unit = nativeBackend.close()

  @Benchmark def pureLu(): Any = square.lu(using PureBackend)
  @Benchmark def ffmLu(): Any = square.lu(using nativeBackend)

  @Benchmark def pureCholesky(): Any = spd.cholesky(using PureBackend)
  @Benchmark def ffmCholesky(): Any = spd.cholesky(using nativeBackend)

  @Benchmark def pureQr(): Any = tall.qr(using PureBackend)
  @Benchmark def ffmQr(): Any = tall.qr(using nativeBackend)

  @Benchmark def pureEigValues(): Any =
    Eigen.eigSymmetric(spd, EigenSelection.All, EigenVectors.ValuesOnly)

  @Benchmark def ffmEigValues(): Any =
    nativeBackend.spectral.get.denseSymmetricEigen(spd, wantVectors = false)

/** End-to-end solver scenarios, including factorization plus the typed solve. */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class FfmSolverScenarioJmh:
  @Param(Array("64", "128", "256", "512"))
  var n: Int = _

  private var square: DMat = _
  private var tall: DMat = _
  private var squareRhs: DVec = _
  private var tallRhs: DVec = _
  private var nativeBackend: FfmBlasBackend = _

  @Setup(Level.Trial)
  def setup(): Unit =
    square = Matrix.tabulate(n, n)((i, j) =>
      val base = (((i.toLong * 29 + j.toLong * 17 + 5) % 83).toDouble - 41.0) / 43.0
      if i == j then base + n.toDouble else base
    )
    tall = Matrix.tabulate(2 * n, n)((i, j) =>
      val base = math.cos((i + 1).toDouble * 0.017 + (j + 3).toDouble * 0.019)
      if i == j then base + 2.0 else base
    )
    squareRhs = Vec.tabulate(n)(i => ((i * 7) % 19 - 9).toDouble / 10.0)
    val truth = Vec.tabulate(n)(i => ((i * 11) % 23 - 11).toDouble / 12.0)
    tallRhs = tall.*(truth)(using PureBackend)
    nativeBackend = FfmBlasBackend.load(thresholds = Some(FfmBlasThresholds(
      nativeLuMinSize = 0,
      nativeCholeskyMinSize = 0,
      nativeQrMinSize = 0
    ))).fold(throw _, identity)

  @TearDown(Level.Trial)
  def tearDown(): Unit = nativeBackend.close()

  @Benchmark def pureSolve(): Any = square.solve(squareRhs)(using PureBackend)
  @Benchmark def ffmSolve(): Any = square.solve(squareRhs)(using nativeBackend)

  @Benchmark def pureLeastSquares(): Any = tall.leastSquares(tallRhs)(using PureBackend)
  @Benchmark def ffmLeastSquares(): Any = tall.leastSquares(tallRhs)(using nativeBackend)

/** Measures the cost of selecting but declining a native provider. At n=32 the
  * default Accelerate LU threshold is not met and QR is default-disabled, so all
  * four methods execute Gale's pure algorithms after their normal public dispatch.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class FfmDispatchJmh:
  private var square: DMat = _
  private var tall: DMat = _
  private var backend: FfmBlasBackend = _

  @Setup(Level.Trial)
  def setup(): Unit =
    val n = 32
    square = Matrix.tabulate(n, n)((i, j) =>
      val base = (((i * 29 + j * 17 + 5) % 83).toDouble - 41.0) / 43.0
      if i == j then base + n.toDouble else base
    )
    tall = Matrix.tabulate(2 * n, n)((i, j) =>
      val base = math.cos((i + 1).toDouble * 0.017 + (j + 3).toDouble * 0.019)
      if i == j then base + 2.0 else base
    )
    backend = FfmBlasBackend.load().fold(throw _, identity)

  @TearDown(Level.Trial)
  def tearDown(): Unit = backend.close()

  @Benchmark def pureLu(): Any = square.lu(using PureBackend)
  @Benchmark def declinedNativeLu(): Any = square.lu(using backend)
  @Benchmark def pureQr(): Any = tall.qr(using PureBackend)
  @Benchmark def declinedNativeQr(): Any = tall.qr(using backend)
