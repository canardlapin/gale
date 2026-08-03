package gale.bench

import scala.compiletime.uninitialized

import gale.backend.PureBackend
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Copy-inclusive court for Gale's portable tall-skinny column-pivoted QR.
  *
  * The immutable input is prepared outside the timed boundary. Each invocation
  * includes the mandatory row-major factor copy, exact pivot-norm computation,
  * compact reflector and R construction, rank decision, and owned factor
  * results. Explicit `PureBackend` pins the shared JVM/Scala.js algorithm.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@Threads(1)
@State(Scope.Thread)
class TallPivotedQrJmh:
  @Param(Array("512", "1024", "2048", "4096", "10000"))
  var n: Int = 0

  @Param(Array("3", "5", "6", "8", "16", "24"))
  var p: Int = 0

  private var design: DMat = uninitialized
  private val options = QROptions(pivoting = QRPivoting.Column)

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    require(p <= n, s"tall pivoted-QR court requires n >= p, got n=$n p=$p")
    design = galeMatrix(matrixData(n, p, seed = 0x54535152L + n * 31L + p))
    require(
      design.qr(options)(using PureBackend).diagnostics.rank.contains(p),
      s"fixture is not full rank for n=$n p=$p"
    )

  @Benchmark
  def factorPivotedQr(blackhole: Blackhole): Unit =
    blackhole.consume(design.qr(options)(using PureBackend))
