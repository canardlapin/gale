package gale.bench

import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*

/** Dense BLAS-1/2 kernels parameterised over problem size.
  *
  * `dot` and `axpy` run on length-`n` vectors; `gemv` multiplies a
  * `gemvRows x n` matrix (rows capped so the largest `n` stays in memory).
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class DenseKernelJmh:
  @Param(Array("256", "4096", "65536"))
  var n: Int = 0

  private val gemvRows = 64

  private var x: DVec = _
  private var y: DVec = _
  private var yWork: MutableDVec = _
  private var a: DMat = _
  private var gemvY: MutableDVec = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    x = Vec.tabulate(n)(i => ((i % 17) - 8).toDouble * 0.125)
    y = Vec.tabulate(n)(i => ((i % 11) - 5).toDouble * 0.25)
    a = Matrix.tabulate(gemvRows, n)((i, j) => ((i * 3 + j * 5) % 19).toDouble * 0.03125)
    gemvY = MutableVec.zeros(gemvRows)

  @Setup(Level.Iteration)
  def setupIteration(): Unit =
    yWork = y.mutableCopy

  @Benchmark
  def dot(): Double =
    x.dot(y)

  @Benchmark
  def axpy(): Double =
    yWork += x
    yWork(0)

  @Benchmark
  def gemv(): Double =
    a.mulInto(x, gemvY)
    gemvY(0)
