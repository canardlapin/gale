package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseVector as BDV
import gale.backend.Backend
import gale.backend.BackendConfig
import gale.backend.BackendThresholds
import gale.backend.Capability
import gale.backend.DenseDoubleKernel
import gale.backend.PureBackend
import gale.backend.PureDenseDoubleKernel
import gale.backend.PureThresholds
import gale.backend.jvm.vector.VectorBackend
import gale.backend.jvm.vector.VectorL1Kernels
import gale.bench.BreezeBenchData.*
import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams

/** W2.1 spike (ADR A-2b gate): explicit Vector API level-1 and reduction kernels
  * ([[gale.backend.jvm.vector.VectorL1Kernels]]) against the pure static kernels
  * ([[gale.kernel.DoubleKernels]]) and Breeze. Lane B only: it needs
  * `--add-modules=jdk.incubator.vector`, under which Breeze's netlib resolves
  * VectorBLAS (logged by [[BreezeBenchData.recordNetlib]]).
  *
  * The Breeze `dot`/`axpy`/`nrm2` twins call the resolved netlib `BLAS` instance
  * directly, so they measure VectorBLAS itself rather than Breeze's dispatch.
  * `sum`/`max`/`argmax`/`exp` go through `breeze.linalg`/`breeze.numerics`;
  * `breezeExp` allocates its result while the gale twins write a preallocated
  * output, so it carries one `n`-double allocation the others do not.
  *
  * Every in-place `axpy` mutates a work vector reset each iteration; with
  * `alpha = 1e-3` it stays finite for any realistic invocation count.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SimdL1SpikeJmh:
  @Param(Array("1024", "4096", "65536", "1048576"))
  var n: Int = 0

  private val alpha = 1e-3

  private var x: DoubleArray = uninitialized
  private var y: DoubleArray = uninitialized
  private var out: DoubleArray = uninitialized
  private var xs: Array[Double] = uninitialized
  private var ys: Array[Double] = uninitialized
  private var bx: BDV[Double] = uninitialized
  private var blas: dev.ludovic.netlib.blas.BLAS = uninitialized

  private var simdWork: DoubleArray = uninitialized
  private var pureWork: DoubleArray = uninitialized
  private var breezeWork: Array[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    blas = dev.ludovic.netlib.blas.BLAS.getInstance()
    xs = vectorData(n, 1L)
    ys = vectorData(n, 2L)
    x = DoubleArray.fromArray(xs)
    y = DoubleArray.fromArray(ys)
    out = DoubleArray.alloc(n)
    bx = breezeVector(xs)

  @Setup(Level.Iteration)
  def setupIteration(): Unit =
    simdWork = DoubleArray.fromArray(ys)
    pureWork = DoubleArray.fromArray(ys)
    breezeWork = ys.clone()

  // --- dot -----------------------------------------------------------------------
  @Benchmark def simdDot(): Double = VectorL1Kernels.ddot(n, x, 0, 1, y, 0, 1)
  @Benchmark def pureDot(): Double = DoubleKernels.ddot(n, x, 0, 1, y, 0, 1)
  @Benchmark def breezeDot(): Double = blas.ddot(n, xs, 1, ys, 1)

  // --- axpy (in place) -------------------------------------------------------------
  @Benchmark def simdAxpy(): Double =
    VectorL1Kernels.daxpy(n, alpha, x, 0, 1, simdWork, 0, 1)
    DoubleArray.asArray(simdWork)(0)
  @Benchmark def pureAxpy(): Double =
    DoubleKernels.daxpy(n, alpha, x, 0, 1, pureWork, 0, 1)
    DoubleArray.asArray(pureWork)(0)
  @Benchmark def breezeAxpy(): Double =
    blas.daxpy(n, alpha, xs, 1, breezeWork, 1)
    breezeWork(0)

  // --- nrm2 ------------------------------------------------------------------------
  @Benchmark def simdNrm2(): Double = VectorL1Kernels.dnrm2(n, y, 0, 1)
  @Benchmark def pureNrm2(): Double = DoubleKernels.dnrm2(n, y, 0, 1)
  @Benchmark def breezeNrm2(): Double = blas.dnrm2(n, ys, 1)

  // --- sum -------------------------------------------------------------------------
  @Benchmark def simdSum(): Double = VectorL1Kernels.dsum(n, x, 0, 1)
  @Benchmark def pureSum(): Double = DoubleKernels.dsum(n, x, 0, 1)
  @Benchmark def breezeSum(): Double = breeze.linalg.sum(bx)

  // --- argmax (and Breeze's value max, for reference only) -------------------------
  @Benchmark def simdArgmax(): Int = VectorL1Kernels.dmaxIndex(n, x, 0, 1)
  @Benchmark def pureArgmax(): Int = DoubleKernels.dmaxIndex(n, x, 0, 1)
  @Benchmark def breezeArgmax(): Int = breeze.linalg.argmax(bx)
  @Benchmark def breezeMax(): Double = breeze.linalg.max(bx)

  // --- exp map ---------------------------------------------------------------------
  @Benchmark def simdExp(): Double =
    VectorL1Kernels.dexpInto(n, x, 0, 1, out, 0, 1)
    DoubleArray.asArray(out)(0)
  @Benchmark def pureExp(): Double =
    DoubleKernels.dexpInto(n, x, 0, 1, out, 0, 1)
    DoubleArray.asArray(out)(0)
  @Benchmark def breezeExp(): Double =
    val result: BDV[Double] = breeze.numerics.exp(bx)
    result(0)

end SimdL1SpikeJmh

/** The `(using Backend)` call shape a W2.2 facade would compile to: two interface
  * calls (`Backend.denseDouble`, then `DenseDoubleKernel.dot`). Kept out of line in
  * one method so its call-site type profile is shared by every caller.
  */
object BackendDispatchSite:
  def dot(n: Int, x: DoubleArray, y: DoubleArray)(using backend: Backend): Double =
    backend.denseDouble.dot(n, x, 0, 1, y, 0, 1)

/** A third `Backend`/`DenseDoubleKernel` class, used only to push the dispatch call
  * site past HotSpot's bimorphic inline cache (pure + vector + this one).
  */
object ForwardingDenseDoubleKernel extends DenseDoubleKernel:
  private val p = PureDenseDoubleKernel
  def dot(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Double =
    p.dot(n, x, xOffset, xStride, y, yOffset, yStride)
  def nrm2(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double = p.nrm2(n, x, xOffset, xStride)
  def copy(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    p.copy(n, x, xOffset, xStride, y, yOffset, yStride)
  def axpy(
      n: Int,
      alpha: Double,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit = p.axpy(n, alpha, x, xOffset, xStride, y, yOffset, yStride)
  def scal(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, xStride: Int): Unit =
    p.scal(n, alpha, x, xOffset, xStride)
  def gemv(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      rowStride: Int,
      colStride: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      beta: Double,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit = p.gemv(rows, cols, alpha, a, aOffset, rowStride, colStride, x, xOffset, xStride, beta, y, yOffset, yStride)
  def gemm(
      rows: Int,
      cols: Int,
      shared: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      aColStride: Int,
      b: DoubleArray,
      bOffset: Int,
      bRowStride: Int,
      bColStride: Int,
      beta: Double,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int,
      cColStride: Int
  ): Unit =
    p.gemm(
      rows, cols, shared, alpha, a, aOffset, aRowStride, aColStride,
      b, bOffset, bRowStride, bColStride, beta, c, cOffset, cRowStride, cColStride
    )
  def syrk(
      m: Int,
      k: Int,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int
  ): Unit = p.syrk(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)

object ForwardingBackend extends Backend:
  val name: String = "bench-forwarding"
  val capabilities: Set[Capability] = Set.empty
  val denseDouble: DenseDoubleKernel = ForwardingDenseDoubleKernel
  val thresholds: BackendThresholds = PureThresholds
  val config: BackendConfig = BackendConfig.singleThreaded

/** Requirement A-R1: the no-import `(using Backend)` path must be allocation-free and
  * within JMH noise of the direct static call. `implementations = 1` keeps the
  * dispatch site monomorphic (only `PureBackend` ever reaches it);
  * `implementations = 3` first drives it with the pure, vector and a forwarding
  * backend, so the measured `PureBackend` call goes through a megamorphic site.
  * Run with `-prof gc` to check the allocation half of A-R1.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class BackendDispatchJmh:
  @Param(Array("16", "32", "64", "128", "256"))
  var n: Int = 0

  @Param(Array("1", "3"))
  var implementations: Int = 0

  private var x: DoubleArray = uninitialized
  private var y: DoubleArray = uninitialized
  private var backend: Backend = PureBackend

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    x = DoubleArray.fromArray(vectorData(n, 3L))
    y = DoubleArray.fromArray(vectorData(n, 4L))
    if implementations == 3 then
      val all = Array[Backend](PureBackend, VectorBackend, ForwardingBackend)
      var sink = 0.0
      var i = 0
      while i < 60000 do
        sink += BackendDispatchSite.dot(n, x, y)(using all(i % 3))
        i += 1
      if sink.isNaN then println(sink)
    backend = PureBackend

  @Benchmark def directStatic(): Double = DoubleKernels.ddot(n, x, 0, 1, y, 0, 1)

  @Benchmark def viaUsingBackend(): Double = BackendDispatchSite.dot(n, x, y)(using backend)

end BackendDispatchJmh
