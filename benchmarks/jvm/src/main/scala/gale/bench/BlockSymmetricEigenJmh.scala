package gale.bench

import gale.linalg.DoubleLinearOperator
import gale.linalg.LinearOperator
import gale.spectral.Eigen
import gale.spectral.EigenOrder
import gale.spectral.EigenSelection
import gale.spectral.SpectralOptions
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Matrix-free block-Krylov baseline with an eight-dimensional repeated top
  * eigenspace target. The diagonal has only twelve distinct roots, so the
  * benchmark measures block expansion, Rayleigh-Ritz solves, and thick restarts
  * without hiding operator cost inside dense multiplication.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 3, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class BlockSymmetricEigenJmh:
  @Param(Array("128", "512", "2048"))
  var n: Int = 0

  private var operator: DoubleLinearOperator = _
  private val selection = EigenSelection.Count(8, EigenOrder.LargestAlgebraic)
  private val options = SpectralOptions(
    tolerance = 1e-8,
    maxIterations = 16,
    subspaceDimension = Some(32)
  )

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    operator = LinearOperator.fromFunction(n, n): (x, into) =>
      var i = 0
      while i < n do
        into(i) = (i % 12).toDouble * x(i)
        i += 1

  @Benchmark
  def repeatedTopEigenpairs(blackhole: Blackhole): Unit =
    blackhole.consume(Eigen.eigSymmetric(operator, n, selection, options))
