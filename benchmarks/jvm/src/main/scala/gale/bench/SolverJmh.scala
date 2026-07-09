package gale.bench

import gale.linalg.*
import gale.solvers.*
import gale.sparse.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*

/** Solver-scenario benchmark: conjugate gradient on a 2D 5-point Laplacian over
  * a 64x64 grid (a 4096 x 4096 SPD system stored as CSR).
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SolverJmh:
  private val grid = 64

  private var a: CSR = _
  private var b: DVec = _
  private var config: SolverConfig = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    val n = grid * grid
    val builder = Sparse.coo(n, n)
    var idx = 0
    while idx < n do
      val r = idx / grid
      val c = idx % grid
      builder.add(idx, idx, 4.0)
      if r > 0 then builder.add(idx, idx - grid, -1.0)
      if r < grid - 1 then builder.add(idx, idx + grid, -1.0)
      if c > 0 then builder.add(idx, idx - 1, -1.0)
      if c < grid - 1 then builder.add(idx, idx + 1, -1.0)
      idx += 1
    a = builder.toCSR()
    b = Vec.fill(n)(1.0)
    config = SolverConfig(tolerance = 1e-8, maxIterations = 500)

  @Benchmark
  def cgLaplacian(): Double =
    cg(a, b, config).residual
