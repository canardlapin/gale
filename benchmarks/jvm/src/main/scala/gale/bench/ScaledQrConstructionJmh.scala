package gale.bench

import scala.compiletime.uninitialized

import gale.backend.PureBackend
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Copy-inclusive court for transient QR construction and algebraic row
  * scaling. Immutable source data is prepared outside the timed boundary;
  * construction, scaling, factorization, solve work, and owned results are
  * timed.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 10, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@Threads(1)
@State(Scope.Thread)
class ScaledQrConstructionJmh:
  @Param(Array("1024"))
  var n: Int = 0

  @Param(Array("5"))
  var p: Int = 0

  @Param(Array("8"))
  var q: Int = 0

  private var design: DMat = uninitialized
  private var scales: DVec = uninitialized
  private var responses: DMat = uninitialized
  private var workspace: DenseWorkspace = uninitialized
  private val options = QROptions(pivoting = QRPivoting.Column)

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    design = galeMatrix(matrixData(n, p, seed = 0x53515243L + n * 31L + p))
    scales = Vec.tabulate(n)(row => 0.5 + (row % 17).toDouble / 16.0)
    responses = galeMatrix(matrixData(n, q, seed = 0x53525248L + n * 37L + q))
    workspace = DenseWorkspace.empty
    workspace.reserve(DenseWorkspace.qrRequirement(n, p, options).orThrow)
    workspace.reserve(DenseWorkspace.qrSolveRequirement(n, q).orThrow)

  @Benchmark
  def builderResultThenQr(blackhole: Blackhole): Unit =
    val builder = copyBuilder(scaleRows = false)
    blackhole.consume(builder.result().qr(options)(using PureBackend))

  @Benchmark
  def materializedScaledQr(blackhole: Blackhole): Unit =
    val builder = copyBuilder(scaleRows = true)
    blackhole.consume(builder.result().qr(options)(using PureBackend))

  @Benchmark
  def materializedScaledSolve(blackhole: Blackhole): Unit =
    val designBuilder = copyBuilder(scaleRows = true)
    val responseBuilder = DMatBuilder.zeros(n, q)
    var row = 0
    while row < n do
      val scale = scales(row)
      var col = 0
      while col < q do
        responseBuilder(row, col) = scale * responses(row, col)
        col += 1
      row += 1
    val factor = designBuilder.result().qr(options)(using PureBackend)
    blackhole.consume(factor.solveLeastSquares(responseBuilder.result()))

  @Benchmark
  def builderConsumeQr(blackhole: Blackhole): Unit =
    val builder = copyBuilder(scaleRows = false)
    blackhole.consume(builder.consumeQR(options, workspace)(using PureBackend))

  @Benchmark
  def directScaledQr(blackhole: Blackhole): Unit =
    blackhole.consume(design.qrScaledRows(scales, options, workspace)(using PureBackend))

  @Benchmark
  def directScaledSolve(blackhole: Blackhole): Unit =
    val factor = design.qrScaledRows(scales, options, workspace)(using PureBackend).orThrow
    blackhole.consume(factor.solveLeastSquaresScaledRowsWith(responses, scales, workspace))

  private def copyBuilder(scaleRows: Boolean): DMatBuilder =
    val builder = DMatBuilder.zeros(n, p)
    var row = 0
    while row < n do
      val scale = if scaleRows then scales(row) else 1.0
      var col = 0
      while col < p do
        builder(row, col) = scale * design(row, col)
        col += 1
      row += 1
    builder
