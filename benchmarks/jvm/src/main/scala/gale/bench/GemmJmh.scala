package gale.bench

import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*

/** Dense matrix-matrix product over `n x n` operands. The `n = 512` point
  * exercises the blocked row-major path; smaller `n` use the unblocked i-k-j
  * path.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class GemmJmh:
  @Param(Array("32", "128", "512"))
  var n: Int = 0

  private var a: DMat = _
  private var b: DMat = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    a = Matrix.tabulate(n, n)((i, j) => ((i * 3 + j * 5) % 19).toDouble * 0.03125)
    b = Matrix.tabulate(n, n)((i, j) => ((i * 7 - j * 2) % 23).toDouble * 0.015625)

  @Benchmark
  def gemm(): DMat =
    a * b
