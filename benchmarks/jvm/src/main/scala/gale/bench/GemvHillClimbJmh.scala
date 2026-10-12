package gale.bench

import scala.compiletime.uninitialized

import dev.ludovic.netlib.blas.JavaBLAS
import gale.backend.jvm.vector.VectorDenseDoubleKernel
import gale.backend.jvm.vector.VectorGemvBenchHook
import gale.platform.DoubleArray
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*

/** Paired kernel-level GEMV hill climb (HC-4) for the Vector backend: in one fork,
  *
  *   - `refGemv` runs the frozen pre-climb row-major kernel (`VectorGemvBaseline`),
  *   - `curGemv` the live `VectorDenseDoubleKernel.gemv` on the same row-major `A`,
  *   - `noTreeGemv` the live tile without the eight-lane row-sum tree (identical to
  *     `curGemv` on a species narrower than eight doubles),
  *   - `vblasGemv` Breeze's backend for `A * x`, VectorBLAS `dgemv("N")` on a
  *     column-major copy of the same `A`,
  *   - `curGemvT` / `vblasGemvT` the transpose product `Aᵀx`: gale's column-contiguous
  *     sweep over the row-major `A` against VectorBLAS `dgemv("T")` on the column-major
  *     copy (Breeze's `A.t * x`).
  *
  * Shapes are `rows x cols` of `A`: the square `BlasL2BreezeJmh` sizes and the
  * `1000 x 50` logistic-regression design of `LbfgsBreezeJmh`. Every method writes a
  * preallocated `y` with `beta = 0`, so nothing allocates and the facade is out of the
  * picture. Needs `--add-modules=jdk.incubator.vector` (inherited by the JMH forks);
  * setup fails unless VectorBLAS is the resolved Java BLAS.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(3)
@State(Scope.Thread)
class GemvHillClimbJmh:
  @Param(Array("256x256", "1024x1024", "2048x2048", "1000x50"))
  var shape: String = uninitialized

  private var rows: Int = 0
  private var cols: Int = 0
  private var aRow: Array[Double] = uninitialized
  private var aCol: Array[Double] = uninitialized
  private var x: Array[Double] = uninitialized
  private var xt: Array[Double] = uninitialized
  private var y: Array[Double] = uninitialized
  private var yt: Array[Double] = uninitialized
  private var blas: JavaBLAS = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    val Array(r, c) = shape.split('x').map(_.toInt)
    rows = r
    cols = c
    val rng = new java.util.Random(4041L)
    aRow = Array.fill(r * c)(rng.nextDouble() * 2.0 - 1.0)
    aCol = new Array[Double](r * c)
    for i <- 0 until r; j <- 0 until c do aCol(j * r + i) = aRow(i * c + j)
    x = Array.fill(c)(rng.nextDouble() * 2.0 - 1.0)
    xt = Array.fill(r)(rng.nextDouble() * 2.0 - 1.0)
    y = new Array[Double](r)
    yt = new Array[Double](c)
    blas = JavaBLAS.getInstance()
    val impl = blas.getClass.getName
    System.err.println(s"[gemv-hill-climb] shape=$shape blas=$impl")
    if impl != "dev.ludovic.netlib.blas.VectorBLAS" then
      throw new IllegalStateException(s"expected VectorBLAS, got $impl (run with --add-modules=jdk.incubator.vector)")

  @Benchmark def refGemv(): Double =
    VectorGemvBenchHook.gemvRowMajorBaseline(rows, cols, 1.0, aRow, 0, cols, x, 0, 0.0, y, 0)
    y(0)

  @Benchmark def curGemv(): Double =
    VectorDenseDoubleKernel.gemv(
      rows, cols, 1.0, DoubleArray.adopt(aRow), 0, cols, 1,
      DoubleArray.adopt(x), 0, 1, 0.0, DoubleArray.adopt(y), 0, 1
    )
    y(0)

  @Benchmark def noTreeGemv(): Double =
    VectorGemvBenchHook.gemvRowMajorNoTree(rows, cols, 1.0, aRow, 0, cols, x, 0, 0.0, y, 0)
    y(0)

  @Benchmark def vblasGemv(): Double =
    blas.dgemv("N", rows, cols, 1.0, aCol, rows, x, 1, 0.0, y, 1)
    y(0)

  @Benchmark def curGemvT(): Double =
    VectorDenseDoubleKernel.gemv(
      cols, rows, 1.0, DoubleArray.adopt(aRow), 0, 1, cols,
      DoubleArray.adopt(xt), 0, 1, 0.0, DoubleArray.adopt(yt), 0, 1
    )
    yt(0)

  @Benchmark def vblasGemvT(): Double =
    blas.dgemv("T", rows, cols, 1.0, aCol, rows, xt, 1, 0.0, yt, 1)
    yt(0)
