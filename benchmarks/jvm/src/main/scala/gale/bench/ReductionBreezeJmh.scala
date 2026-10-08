package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.*
import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.bench.BreezeBenchData.*
import gale.linalg.Axis
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Numerics
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams
import org.openjdk.jmh.infra.Blackhole

/** Vector reductions, norms and elementwise numerics, gale vs Breeze, at length `n`.
  * Scalar results are returned (JMH consumes them); vector results go to the
  * `Blackhole`. Pairings, each semantically equal:
  *
  *   - `sum`, `max`, `argmax`: gale `x.sum` (fast, reassociating), `x.max`, `x.argmax`
  *     vs `breeze.linalg.sum`/`max`/`argmax`.
  *   - `mean`: gale `x.mean` vs `breeze.stats.mean`.
  *   - `norm1`, `normInf`: gale `x.norm1`, `x.normInf` vs breeze `norm(x, 1.0)`,
  *     `norm(x, Double.PositiveInfinity)`.
  *   - `exp`, `sigmoid`: gale `Numerics.exp`/`sigmoid` vs `breeze.numerics.exp`/`sigmoid`
  *     (one allocated result each).
  *   - `logSumExp`: gale `Numerics.logSumExp(x)` vs breeze `softmax(x)` — '''Breeze's
  *     `softmax` is the scalar log-sum-exp''', not the normalized vector.
  *   - `softmax` (probabilities): gale `Numerics.softmax(x)` vs the Breeze idiom
  *     `exp(x - softmax(x))`. Asymmetric by construction: Breeze has no fused
  *     normalized softmax, so its side allocates two vectors (`x - lse`, then `exp`)
  *     and makes an extra pass; gale allocates one.
  *
  * Backend-insensitive: these gale ops take no `Backend` (W2 may route them).
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class ReductionBreezeJmh:
  @Param(Array("1024", "65536", "1048576"))
  var n: Int = 0

  private var gx: DVec        = uninitialized
  private var bx: BDV[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val data = vectorData(n, 5300L)
    gx = galeVector(data)
    bx = breezeVector(data)

  @Benchmark def galeSum(): Double   = gx.sum
  @Benchmark def breezeSum(): Double = sum(bx)

  @Benchmark def galeMax(): Double   = gx.max
  @Benchmark def breezeMax(): Double = max(bx)

  @Benchmark def galeArgmax(): Int   = gx.argmax
  @Benchmark def breezeArgmax(): Int = argmax(bx)

  @Benchmark def galeMean(): Double   = gx.mean
  @Benchmark def breezeMean(): Double = breeze.stats.mean(bx)

  @Benchmark def galeNorm1(): Double   = gx.norm1
  @Benchmark def breezeNorm1(): Double = norm(bx, 1.0)

  @Benchmark def galeNormInf(): Double   = gx.normInf
  @Benchmark def breezeNormInf(): Double = norm(bx, Double.PositiveInfinity)

  @Benchmark def galeExp(bh: Blackhole): Unit   = bh.consume(Numerics.exp(gx))
  @Benchmark def breezeExp(bh: Blackhole): Unit = bh.consume(breeze.numerics.exp(bx))

  @Benchmark def galeSigmoid(bh: Blackhole): Unit   = bh.consume(Numerics.sigmoid(gx))
  @Benchmark def breezeSigmoid(bh: Blackhole): Unit = bh.consume(breeze.numerics.sigmoid(bx))

  @Benchmark def galeLogSumExp(): Double   = Numerics.logSumExp(gx)
  @Benchmark def breezeLogSumExp(): Double = softmax(bx)

  @Benchmark def galeSoftmax(bh: Blackhole): Unit = bh.consume(Numerics.softmax(gx))
  @Benchmark def breezeSoftmax(bh: Blackhole): Unit =
    bh.consume(breeze.numerics.exp(bx - softmax(bx)))

/** Matrix reductions on a square `n × n` matrix (default `1024 × 1024`), gale vs
  * Breeze. gale's [[gale.linalg.Axis]] names what the result is indexed by:
  * `Axis.Rows` gives one value per row (Breeze `A(*, ::)`), `Axis.Cols` one value
  * per column (Breeze `A(::, *)`, whose result Breeze returns transposed).
  *
  *   - `sumRows` / `sumCols`: `A.sum(Axis.Rows|Cols)` vs `sum(A(*, ::))` / `sum(A(::, *))`.
  *   - `maxRows` / `maxCols`: `A.max(Axis.Rows|Cols)` vs `max(A(*, ::))` / `max(A(::, *))`.
  *   - `sum`: whole-matrix `A.sum` vs `sum(A)`.
  *   - `normFrobenius`: gale `A.normFrobenius` (overflow-safe scaled recurrence) vs
  *     breeze `norm(A.flatten())` (a zero-copy view of the contiguous column-major
  *     data; plain sum of squares, no scaling).
  *   - `softmaxRows`: per-row normalized softmax — gale `Numerics.softmax(A, Axis.Rows)`
  *     vs the Breeze idiom `exp(A(::, *) - softmax(A(*, ::)))`: per-row log-sum-exp
  *     broadcast down the columns, then `exp`. Breeze allocates the shifted matrix and
  *     the result; gale allocates the result only.
  *   - `logSumExpRows`: gale `Numerics.logSumExp(A, Axis.Rows)` vs breeze
  *     `softmax(A(*, ::))` (Breeze's `softmax` is log-sum-exp).
  *
  * Backend-insensitive.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class MatrixReductionBreezeJmh:
  @Param(Array("1024"))
  var n: Int = 0

  private var gA: DMat        = uninitialized
  private var bA: BDM[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val data = matrixData(n, n, 5400L)
    gA = galeMatrix(data)
    bA = breezeMatrix(data)

  @Benchmark def galeSumRows(bh: Blackhole): Unit   = bh.consume(gA.sum(Axis.Rows))
  @Benchmark def breezeSumRows(bh: Blackhole): Unit = bh.consume(sum(bA(*, ::)))

  @Benchmark def galeSumCols(bh: Blackhole): Unit   = bh.consume(gA.sum(Axis.Cols))
  @Benchmark def breezeSumCols(bh: Blackhole): Unit = bh.consume(sum(bA(::, *)))

  @Benchmark def galeMaxRows(bh: Blackhole): Unit   = bh.consume(gA.max(Axis.Rows))
  @Benchmark def breezeMaxRows(bh: Blackhole): Unit = bh.consume(max(bA(*, ::)))

  @Benchmark def galeMaxCols(bh: Blackhole): Unit   = bh.consume(gA.max(Axis.Cols))
  @Benchmark def breezeMaxCols(bh: Blackhole): Unit = bh.consume(max(bA(::, *)))

  @Benchmark def galeSum(): Double   = gA.sum
  @Benchmark def breezeSum(): Double = sum(bA)

  @Benchmark def galeNormFrobenius(): Double   = gA.normFrobenius
  @Benchmark def breezeNormFrobenius(): Double = norm(bA.flatten())

  @Benchmark def galeSoftmaxRows(bh: Blackhole): Unit = bh.consume(Numerics.softmax(gA, Axis.Rows))
  @Benchmark def breezeSoftmaxRows(bh: Blackhole): Unit =
    bh.consume(breeze.numerics.exp(bA(::, *) - softmax(bA(*, ::))))

  @Benchmark def galeLogSumExpRows(bh: Blackhole): Unit = bh.consume(Numerics.logSumExp(gA, Axis.Rows))
  @Benchmark def breezeLogSumExpRows(bh: Blackhole): Unit = bh.consume(softmax(bA(*, ::)))
