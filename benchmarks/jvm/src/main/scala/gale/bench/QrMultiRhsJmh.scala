package gale.bench

import scala.compiletime.uninitialized

import gale.bench.BreezeBenchData.*
import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Copy-inclusive court for applying a stored pivoted QR factor to dense matrix right-hand sides.
  *
  * The QR factor and immutable RHS are prepared outside the timed boundary. Both timed factor-application methods
  * include Gale's mandatory RHS copy and owned result construction. `solveLeastSquares` additionally includes the
  * strided triangular solves, pivot permutation, and coefficient ownership. `factorPivotedQr` is a protected control
  * for the factorization path; it does not materialize `Q`.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@Threads(1)
@State(Scope.Thread)
class QrMultiRhsJmh:
  @Param(Array("512", "2048", "10000"))
  var n: Int = 0

  @Param(Array("6", "24"))
  var p: Int = 0

  @Param(Array("1", "8", "16", "32", "100"))
  var q: Int = 0

  private var design: DMat = uninitialized
  private var responses: DMat = uninitialized
  private var response: DVec = uninitialized
  private var factor: QR = uninitialized
  private var matrixWorkspace: DenseWorkspace = uninitialized
  private var vectorWorkspace: DenseWorkspace = uninitialized

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    require(p <= n, s"QR multi-RHS court requires a tall design, got n=$n p=$p")
    design = galeMatrix(matrixData(n, p, seed = 0x514d5253L + n * 31L + p))
    responses = galeMatrix(matrixData(n, q, seed = 0x52534853L + n * 37L + q))
    response = responses.col(0)
    factor = design.qr(QROptions(pivoting = QRPivoting.Column))
    matrixWorkspace = DenseWorkspace.forQRSolve(n, q)
    vectorWorkspace = DenseWorkspace.forQRSolve(n)
    require(factor.diagnostics.rank.contains(p), s"fixture is not full rank: ${factor.diagnostics}")

  @Benchmark
  def applyQtOwned(blackhole: Blackhole): Unit =
    blackhole.consume(factor.applyQT(responses))

  @Benchmark
  def solveLeastSquaresOwned(blackhole: Blackhole): Unit =
    blackhole.consume(factor.solveLeastSquares(responses))

  @Benchmark
  def solveLeastSquaresVectorOwned(blackhole: Blackhole): Unit =
    blackhole.consume(factor.solveLeastSquares(response))

  @Benchmark
  def solveLeastSquaresWithWorkspace(blackhole: Blackhole): Unit =
    blackhole.consume(factor.solveLeastSquaresWith(responses, matrixWorkspace))

  @Benchmark
  def solveLeastSquaresVectorWithWorkspace(blackhole: Blackhole): Unit =
    blackhole.consume(factor.solveLeastSquaresWith(response, vectorWorkspace))

  @Benchmark
  def factorPivotedQr(blackhole: Blackhole): Unit =
    blackhole.consume(design.qr(QROptions(pivoting = QRPivoting.Column)))
