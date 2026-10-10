package gale.bench

import scala.compiletime.uninitialized

import gale.linalg.*
import gale.spectral.DenseSpectralKernels
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Kernel-only structure coverage for the D&C dispatch. The forced QL route
  * provides a same-build reference; neither method measures facade diagnostics.
  * Diagonal/identity cases guard against an unnecessary cubic back-transform.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(2)
@State(Scope.Thread)
class SymmetricEigenStructureJmh:
  @Param(Array("128", "512", "1024"))
  var n: Int = 0

  @Param(Array("diagonal", "identity", "dense"))
  var structure: String = ""

  private var a: DMat = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    val rng = new scala.util.Random(600L)
    val raw = Array.fill(n * n)(rng.nextDouble() * 2.0 - 1.0)
    a = Matrix.tabulate(n, n): (r, c) =>
      structure match
        case "diagonal" => if r == c then (n - r).toDouble else 0.0
        case "identity" => if r == c then 1.0 else 0.0
        case "dense" => if r >= c then raw(r * n + c) else raw(c * n + r)
        case other => throw new IllegalArgumentException(s"unknown structure: $other")
    val result = DenseSpectralKernels.symmetricEigen(a, wantVectors = true).toOption.get
    require(result.values.toSeq.forall(_.isFinite))

  @Benchmark def automatic(bh: Blackhole): Unit =
    bh.consume(DenseSpectralKernels.symmetricEigen(a, wantVectors = true).toOption.get)

  @Benchmark def ql(bh: Blackhole): Unit =
    bh.consume(DenseSpectralKernels.symmetricEigen(a, wantVectors = true, divideAndConquer = false).toOption.get)
