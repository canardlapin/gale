package gale.bench

import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import scala.compiletime.uninitialized

/** Fixed-band storage, varied n and bandwidth. Resetting RHS is included in
  * in-place solve timings; no per-invocation allocation is hidden in setup.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(1)
@State(Scope.Thread)
class BandedCholeskyJmh:
  @Param(Array("300", "1200", "4800"))
  var n: Int = 0
  @Param(Array("4", "16", "64"))
  var bandwidth: Int = 0

  private var bands: DMat = uninitialized
  private var factor: BandedCholesky = uninitialized
  private var rhs: MutableDVec = uninitialized
  private var multiple: DMatBuilder = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    bands = Matrix.tabulate(n, bandwidth + 1): (_, d) =>
      if d == 0 then 2.0 * bandwidth + 1.0 else -0.25
    factor = BandedCholesky.factorLower(bands).orThrow
    rhs = MutableDVec.zeros(n)
    multiple = DMatBuilder.zeros(n, 4)

  @Benchmark
  def factorPacked(): BandedCholesky = BandedCholesky.factorLower(bands).orThrow

  @Benchmark
  def factorAndConditioning(): Double =
    BandedCholesky.factorLower(bands).orThrow.conditioning.conditionNumberUpperBound

  @Benchmark
  def solveInPlace(): Double =
    var i = 0
    while i < n do
      rhs(i) = 1.0
      i += 1
    factor.solveInPlace(rhs).orThrow
    rhs(n / 2)

  @Benchmark
  def solveFourInPlace(): Double =
    multiple.fill(1.0)
    factor.solveInPlace(multiple).orThrow
    multiple(n / 2, 2)
