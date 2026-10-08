package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.backend.Backend
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams
import org.openjdk.jmh.infra.Blackhole

/** BLAS-1 (`O(n)`) paired benchmarks: dot, in-place axpy, and 2-norm, gale vs
  * Breeze on identical length-`n` vectors. Both axpy variants mutate a preallocated
  * work vector reset each iteration, so neither allocates in the timed method.
  *
  * Backend-insensitive: gale's `dot`, `norm2` and `axpyInPlace` take no `Backend`,
  * so these gale methods have no `GaleBackendState` and run once per size.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class BlasL1BreezeJmh:
  @Param(Array("65536", "262144", "1048576"))
  var n: Int = 0

  private val alpha = 1.5

  private var gx: DVec       = uninitialized
  private var gy: DVec       = uninitialized
  private var bx: BDV[Double] = uninitialized
  private var by: BDV[Double] = uninitialized
  private var gWork: MutableDVec = uninitialized
  private var bWork: BDV[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val xData = vectorData(n, 1L)
    val yData = vectorData(n, 2L)
    gx = galeVector(xData)
    gy = galeVector(yData)
    bx = breezeVector(xData)
    by = breezeVector(yData)

  @Setup(Level.Iteration)
  def setupIteration(): Unit =
    gWork = gy.mutableCopy
    bWork = by.copy

  @Benchmark def galeDot(): Double       = gx.dot(gy)
  @Benchmark def breezeDot(): Double      = bx.dot(by)

  @Benchmark def galeAxpy(): Double =
    gWork.axpyInPlace(alpha, gx)
    gWork(0)
  @Benchmark def breezeAxpy(): Double =
    breeze.linalg.axpy(alpha, bx, bWork)
    bWork(0)

  @Benchmark def galeNorm(): Double      = gy.norm2
  @Benchmark def breezeNorm(): Double     = breeze.linalg.norm(by)

/** BLAS-2 (`O(n²)`) paired benchmarks: matrix–vector product `A·x` and its
  * transpose `Aᵀ·x`, both allocating their result (the directly comparable form;
  * gale's zero-allocation `mulInto` path is covered by `DenseKernelJmh`).
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class BlasL2BreezeJmh:
  @Param(Array("256", "1024", "2048"))
  var n: Int = 0

  private var ga: DMat       = uninitialized
  private var gx: DVec       = uninitialized
  private var ba: BDM[Double] = uninitialized
  private var bx: BDV[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val aData = matrixData(n, n, 10L)
    val xData = vectorData(n, 20L)
    ga = galeMatrix(aData)
    gx = galeVector(xData)
    ba = breezeMatrix(aData)
    bx = breezeVector(xData)

  @Benchmark def galeGemv(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(ga * gx)
  @Benchmark def breezeGemv(bh: Blackhole): Unit = bh.consume(ba * bx)

  @Benchmark def galeGemvT(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(ga.t * gx)
  @Benchmark def breezeGemvT(bh: Blackhole): Unit = bh.consume(ba.t * bx)

/** BLAS-3 (`O(n³)`) paired benchmarks: square `A·B`, tall-skinny `T·B`
  * (`T` is `4n × n`), and the transpose product `TᵀT`.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class BlasL3BreezeJmh:
  @Param(Array("16", "64", "256"))
  var n: Int = 0

  private var ga: DMat       = uninitialized
  private var gb: DMat       = uninitialized
  private var gt: DMat       = uninitialized
  private var ba: BDM[Double] = uninitialized
  private var bb: BDM[Double] = uninitialized
  private var bt: BDM[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val aData = matrixData(n, n, 30L)
    val bData = matrixData(n, n, 40L)
    val tData = matrixData(4 * n, n, 50L)
    ga = galeMatrix(aData)
    gb = galeMatrix(bData)
    gt = galeMatrix(tData)
    ba = breezeMatrix(aData)
    bb = breezeMatrix(bData)
    bt = breezeMatrix(tData)

  @Benchmark def galeGemm(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(ga * gb)
  @Benchmark def breezeGemm(bh: Blackhole): Unit = bh.consume(ba * bb)

  @Benchmark def galeGemmTall(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gt * gb)
  @Benchmark def breezeGemmTall(bh: Blackhole): Unit = bh.consume(bt * bb)

  @Benchmark def galeAtA(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gt.t * gt)
  @Benchmark def breezeAtA(bh: Blackhole): Unit = bh.consume(bt.t * bt)
