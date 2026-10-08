package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.CSCMatrix as BCSC
import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.SparseVector as BSV
import breeze.linalg.axpy
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import gale.sparse.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams
import org.openjdk.jmh.infra.Blackhole

/** Sparse matrix products, gale vs Breeze `CSCMatrix`, on one seeded `n × n` matrix
  * (each entry stored with probability `density`) given to both libraries:
  *
  *   - `cscMatvec`: gale `CSC * x` vs breeze `CSC * x`, both allocating the result.
  *     gale's CSC product runs the generic operator path (column scatter into a
  *     fresh vector, then an immutable snapshot copy) — a real API cost, kept.
  *   - `csrMatvec`: gale `CSR * x` (row dot products) vs breeze `CSC * x`. Breeze has
  *     no CSR type, so its twin repeats the CSC product on the same matrix.
  *   - `cscMatmul`: gale `CSC.applyTo(B)` (the operator's column-by-column product;
  *     gale has no dedicated CSC × dense kernel) vs breeze `CSC * B`, with `B` dense
  *     `n × 32`.
  *   - `csrMatmul`: gale `CSR * B` vs breeze `CSC * B` (again the only Breeze format).
  *
  * Backend-insensitive: gale's sparse products take no `Backend`. The gale CSR is
  * converted from the CSC once in setup; neither side converts in a timed method.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SparseMatrixBreezeJmh:
  @Param(Array("1000", "10000"))
  var n: Int = 0

  @Param(Array("0.01", "0.1"))
  var density: Double = 0.0

  private val denseCols = 32

  private var gCsc: CSC        = uninitialized
  private var gCsr: CSR        = uninitialized
  private var gx: DVec         = uninitialized
  private var gB: DMat         = uninitialized
  private var bCsc: BCSC[Double] = uninitialized
  private var bx: BDV[Double]  = uninitialized
  private var bB: BDM[Double]  = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val data  = sparseCsc(n, density, 3100L)
    val xData = vectorData(n, 3200L)
    val bData = matrixData(n, denseCols, 3300L)
    gCsc = CSCPattern.checked(n, n, data.colPtr, data.rowIdx).flatMap(_.bind(data.values)).orThrow
    gCsr = gCsc.toCSR
    gx = galeVector(xData)
    gB = galeMatrix(bData)
    bCsc = new BCSC[Double](data.values.clone(), n, n, data.colPtr.clone(), data.rowIdx.clone())
    bx = breezeVector(xData)
    bB = breezeMatrix(bData)

  @Benchmark def galeCscMatvec(bh: Blackhole): Unit   = bh.consume(gCsc * gx)
  @Benchmark def breezeCscMatvec(bh: Blackhole): Unit = bh.consume(bCsc * bx)

  @Benchmark def galeCsrMatvec(bh: Blackhole): Unit   = bh.consume(gCsr * gx)
  @Benchmark def breezeCsrMatvec(bh: Blackhole): Unit = bh.consume(bCsc * bx)

  @Benchmark def galeCscMatmul(bh: Blackhole): Unit   = bh.consume(gCsc.applyTo(gB))
  @Benchmark def breezeCscMatmul(bh: Blackhole): Unit = bh.consume(bCsc * bB)

  @Benchmark def galeCsrMatmul(bh: Blackhole): Unit   = bh.consume(gCsr * gB)
  @Benchmark def breezeCsrMatmul(bh: Blackhole): Unit = bh.consume(bCsc * bB)

/** Sparse-vector kernels, gale `SparseVector` vs Breeze `SparseVector`, length
  * `100000` with `nnz` stored entries; `x` and `y` have independent seeded patterns.
  *
  *   - `dot`: sparse · sparse (gale two-pointer merge; Breeze's sparse dot).
  *   - `dotDense`: sparse · dense vector.
  *   - `axpy`: `y += 1.5 x` into a dense vector — gale `x.axpyInto(1.5, y)` vs breeze
  *     `axpy(1.5, x, y)`; both mutate a work vector reset each iteration, so neither
  *     allocates in the timed method.
  *   - `add`: sparse `x + y` (union pattern), both allocating a new sparse vector.
  *
  * Backend-insensitive: gale's sparse-vector ops take no `Backend`.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SparseVectorBreezeJmh:
  @Param(Array("100000"))
  var length: Int = 0

  @Param(Array("1000", "10000"))
  var nnz: Int = 0

  private val alpha = 1.5

  private var gx: SparseVector    = uninitialized
  private var gy: SparseVector    = uninitialized
  private var gDense: DVec        = uninitialized
  private var bx: BSV[Double]     = uninitialized
  private var by: BSV[Double]     = uninitialized
  private var bDense: BDV[Double] = uninitialized
  private var gWork: MutableDVec  = uninitialized
  private var bWork: BDV[Double]  = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val (xIdx, xVal) = sparseEntries(length, nnz, 4100L)
    val (yIdx, yVal) = sparseEntries(length, nnz, 4200L)
    val denseData    = vectorData(length, 4300L)
    gx = SparseVector.tryFromArrays(length, xIdx, xVal).orThrow
    gy = SparseVector.tryFromArrays(length, yIdx, yVal).orThrow
    gDense = galeVector(denseData)
    bx = new BSV[Double](xIdx.clone(), xVal.clone(), length)
    by = new BSV[Double](yIdx.clone(), yVal.clone(), length)
    bDense = breezeVector(denseData)

  @Setup(Level.Iteration)
  def setupIteration(): Unit =
    gWork = gDense.mutableCopy
    bWork = bDense.copy

  @Benchmark def galeDot(): Double   = gx.dot(gy)
  @Benchmark def breezeDot(): Double = bx.dot(by)

  @Benchmark def galeDotDense(): Double   = gx.dot(gDense)
  @Benchmark def breezeDotDense(): Double = bx.dot(bDense)

  @Benchmark def galeAxpy(): Double =
    gx.axpyInto(alpha, gWork)
    gWork(0)
  @Benchmark def breezeAxpy(): Double =
    axpy(alpha, bx, bWork)
    bWork(0)

  @Benchmark def galeAdd(bh: Blackhole): Unit   = bh.consume(gx + gy)
  @Benchmark def breezeAdd(bh: Blackhole): Unit = bh.consume(bx + by)
