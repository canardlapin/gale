package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.backend.PureBackend
import gale.backend.jvm.vector.VectorBackend
import gale.linalg.{DMat, DVec, Matrix, Vec}
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

import java.util.concurrent.TimeUnit

/** Copy-inclusive public-facade comparison used to derive Vector backend thresholds. */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class VectorGemmJmh:
  @Param(Array("64", "128", "256", "512"))
  var n: Int = uninitialized

  private var a: DMat = uninitialized
  private var b: DMat = uninitialized
  private var breezeA: BDM[Double] = uninitialized
  private var breezeB: BDM[Double] = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    a = Matrix.tabulate(n, n)((i, j) => (((i.toLong * 131 + j.toLong * 17 + 3) % 97).toDouble - 48) / 49)
    b = Matrix.tabulate(n, n)((i, j) => (((i.toLong * 113 + j.toLong * 19 + 11) % 89).toDouble - 44) / 47)
    breezeA = BDM.tabulate(n, n)((i, j) => a(i, j))
    breezeB = BDM.tabulate(n, n)((i, j) => b(i, j))

  @Benchmark
  def pure(bh: Blackhole): Unit =
    bh.consume(a.*(b)(using PureBackend))

  @Benchmark
  def vector(bh: Blackhole): Unit =
    bh.consume(a.*(b)(using VectorBackend))

  /** Same-process control: with the Vector module enabled, Breeze 2.1 delegates
    * this operation to dev.ludovic's pure-Java VectorBLAS on JDK 22.
    */
  @Benchmark
  def breezeVectorBlas(bh: Blackhole): Unit =
    bh.consume(breezeA * breezeB)

/** Public-facade GEMV crossover for the optional Vector backend. */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class VectorGemvJmh:
  @Param(Array("64", "128", "256", "1024", "2048"))
  var n: Int = uninitialized

  private var a: DMat = uninitialized
  private var x: DVec = uninitialized
  private var breezeA: BDM[Double] = uninitialized
  private var breezeX: BDV[Double] = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    a = Matrix.tabulate(n, n)((i, j) => (((i.toLong * 131 + j.toLong * 17 + 3) % 97).toDouble - 48) / 49)
    x = Vec.tabulate(n)(i => (((i.toLong * 29 + 5) % 53).toDouble - 26) / 27)
    breezeA = BDM.tabulate(n, n)((i, j) => a(i, j))
    breezeX = BDV.tabulate(n)(i => x(i))

  @Benchmark
  def pure(bh: Blackhole): Unit =
    bh.consume(a.*(x)(using PureBackend))

  @Benchmark
  def vector(bh: Blackhole): Unit =
    bh.consume(a.*(x)(using VectorBackend))

  @Benchmark
  def breezeVectorBlas(bh: Blackhole): Unit =
    bh.consume(breezeA * breezeX)

/** `AᵀA` (`syrk`) crossover for the optional Vector backend across tall, square
  * and short-wide shapes (`m × k`, `A` row-major). `kernel` calls the selected
  * backend's `syrk` directly into a preallocated `C`, bypassing the facade's
  * routing gate but not the Vector kernel's own shape guard (`backend=pure` is
  * the pure kernel, `backend=vector` the guarded Vector kernel). `simdForced`
  * times the SIMD path even where that guard would choose pure, so the guard
  * itself can be validated; it ignores the backend parameter (needs lane B).
  * `facade` is the routed `A.t * A`.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class VectorSyrkJmh:
  @Param(Array("64x64", "128x64", "256x64", "128x128", "1024x256", "4x1024", "8x1024", "16x1024", "32x512", "8x4096"))
  var shape: String = uninitialized

  private var a: DMat = uninitialized
  private var aRaw: Array[Double] = uninitialized
  private var c: Array[Double] = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    val Array(m, k) = shape.split('x').map(_.toInt)
    a = Matrix.tabulate(m, k)((i, j) => (((i.toLong * 131 + j.toLong * 17 + 3) % 97).toDouble - 48) / 49)
    aRaw = Array.tabulate(m * k)(idx => a(idx / k, idx % k))
    c = new Array[Double](k * k)

  @Benchmark
  def kernel(g: GaleBackendState): Double =
    g.selected.denseDouble.syrk(
      a.rows, a.cols, a.data, 0, a.cols, // tabulate: offset 0, row stride k
      gale.platform.DoubleArray.adopt(c), 0, a.cols
    )
    c(0)

  @Benchmark
  def simdForced(): Double =
    gale.backend.jvm.vector.VectorSyrkBenchHook.syrkSimdForced(a.rows, a.cols, aRaw, 0, a.cols, c, 0, a.cols)
    c(0)

  @Benchmark
  def facade(g: GaleBackendState, bh: Blackhole): Unit =
    given gale.backend.Backend = g.selected
    bh.consume(a.t * a)
