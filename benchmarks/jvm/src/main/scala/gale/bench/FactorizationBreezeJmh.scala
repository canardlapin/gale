package gale.bench

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.LU as BreezeLU
import breeze.linalg.cholesky
import breeze.linalg.qr
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Dense factorization / solve paired benchmarks (`O(n³)`), gale vs Breeze.
  *
  *   - `solve`: gale `A.solve(b)` (LU) vs breeze `A \ b`.
  *   - `lu`: the factorization only — gale `A.lu` vs breeze `LU.primitive` (raw
  *     `dgetrf`, no `P`/`L`/`U` assembly), the honest factorization-cost match.
  *   - `chol`: gale `S.cholesky` vs breeze `cholesky(S)` on an SPD matrix.
  *   - `qr`: the reflector factorization without materialising `Q` — gale `A.qr`
  *     (Q is lazy) vs breeze `qr.justR(A)`.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class FactorizationBreezeJmh:
  @Param(Array("16", "64", "256"))
  var n: Int = 0

  private var gA: DMat       = _
  private var gS: DMat       = _
  private var gb: DVec       = _
  private var bA: BDM[Double] = _
  private var bS: BDM[Double] = _
  private var bb: BDV[Double] = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    val aData = diagonallyDominant(n, 100L)
    val sData = spd(n, 200L)
    val bData = vectorData(n, 300L)
    gA = galeMatrix(aData)
    gS = galeMatrix(sData)
    gb = galeVector(bData)
    bA = breezeMatrix(aData)
    bS = breezeMatrix(sData)
    bb = breezeVector(bData)

  @Benchmark def galeSolve(bh: Blackhole): Unit   = bh.consume(gA.solve(gb))
  @Benchmark def breezeSolve(bh: Blackhole): Unit = bh.consume(bA \ bb)

  @Benchmark def galeLu(bh: Blackhole): Unit   = bh.consume(gA.lu)
  @Benchmark def breezeLu(bh: Blackhole): Unit = bh.consume(BreezeLU.primitive(bA))

  @Benchmark def galeChol(bh: Blackhole): Unit   = bh.consume(gS.cholesky)
  @Benchmark def breezeChol(bh: Blackhole): Unit = bh.consume(cholesky(bS))

  @Benchmark def galeQr(bh: Blackhole): Unit   = bh.consume(gA.qr)
  @Benchmark def breezeQr(bh: Blackhole): Unit = bh.consume(qr.justR(bA))

/** Tall least-squares paired benchmark: an overdetermined `m × n` system with
  * `m = 4n`, solved via QR both ways — gale `A.leastSquares(b)` vs breeze `A \ b`
  * (which dispatches to LAPACK `dgels` for a non-square operand).
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class LeastSquaresBreezeJmh:
  @Param(Array("16", "64", "256"))
  var n: Int = 0

  private var gA: DMat       = _
  private var gb: DVec       = _
  private var bA: BDM[Double] = _
  private var bb: BDV[Double] = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    val m     = 4 * n
    val aData = matrixData(m, n, 400L)
    val bData = vectorData(m, 500L)
    gA = galeMatrix(aData)
    gb = galeVector(bData)
    bA = breezeMatrix(aData)
    bb = breezeVector(bData)

  @Benchmark def galeLstsq(bh: Blackhole): Unit   = bh.consume(gA.leastSquares(gb))
  @Benchmark def breezeLstsq(bh: Blackhole): Unit = bh.consume(bA \ bb)

/** Representative Scalafim migration workloads. The shapes model a modest fMRI
  * design (`240 x 24`) applied to a batch of 128 voxel time series. Keep these
  * benchmarks allocation-profileable with `-prof gc`: each timed path consumes
  * a matrix result, so JMH reports the owned result cost rather than allowing it
  * to be optimized away.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 3, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Fork(1)
@State(Scope.Thread)
class ScalafimMigrationJmh:
  private val observations = 240
  private val regressors   = 24
  private val responses    = 128

  private var design: DMat           = _
  private var responseBatch: DMat    = _
  private var coefficientBatch: DMat = _
  private var qr: QR                  = _
  private var cholesky: Cholesky      = _
  private var designTranspose: DMat   = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    design = galeMatrix(matrixData(observations, regressors, 600L))
    responseBatch = galeMatrix(matrixData(observations, responses, 700L))
    coefficientBatch = galeMatrix(matrixData(regressors, responses, 800L))
    qr = design.qr(QROptions(pivoting = QRPivoting.Column))
    cholesky = galeMatrix(spd(regressors, 900L)).cholesky.orThrow
    designTranspose = design.t

  @Benchmark def pivotedQrBatchSolve(bh: Blackhole): Unit =
    bh.consume(qr.solveLeastSquares(responseBatch))

  @Benchmark def qrBatchTransform(bh: Blackhole): Unit =
    bh.consume(qr.applyQT(responseBatch))

  @Benchmark def choleskyBatchSolve(bh: Blackhole): Unit =
    bh.consume(cholesky.solve(coefficientBatch))

  @Benchmark def operatorBatchApply(bh: Blackhole): Unit =
    bh.consume(designTranspose.applyTo(responseBatch))
