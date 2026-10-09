package gale.bench

import scala.compiletime.uninitialized

import gale.kernel.DoubleKernels
import gale.linalg.{Axis, DMat, Matrix}
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Gale-only before/after microbench for the max value pass. `ref*` runs a
  * private copy of the `math.max` kernels as they stood before the hill climb
  * (breeze/integration ff03eed); `cur*` runs the live kernel on the same data,
  * so one JMH run gives a paired ratio. `breezeShapeMax` is Breeze 2.1's
  * `max` loop (one `math.max` accumulator) on gale storage.
  *
  * On x86 HotSpot C2 lowers double `Math.max` to a compare-and-branch ladder
  * that does not vectorize, so the live kernel uses an in-order compare scan
  * there (`PlatformMath.scanExtremes`); AArch64 keeps the `fmax` reduction.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 4, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class MaxHillClimbJmh:
  @Param(Array("1024", "65536", "1048576"))
  var n: Int = 0

  private var x: DoubleArray   = uninitialized
  private var side: Int        = 0
  private var a: DMat          = uninitialized
  private var out: DoubleArray = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    x = DoubleArray.adopt(BreezeBenchData.vectorData(n, 5300L))
    side = math.sqrt(n.toDouble).toInt
    a = Matrix.tabulate(side, side)((i, j) => x(i * side + j))
    out = DoubleArray.alloc(side)

  @Benchmark def refMax(): Double         = MaxHillClimbRef.refMax(n, x)
  @Benchmark def breezeShapeMax(): Double = MaxHillClimbRef.breezeShape(n, x)
  @Benchmark def curMax(): Double         = DoubleKernels.dmax(n, x, 0, 1)

  // Column maxima of the side x side row-major matrix: the streamed geometry.
  @Benchmark def refMaxCols(bh: Blackhole): Unit =
    MaxHillClimbRef.refCols(side, x, out)
    bh.consume(out)
  @Benchmark def curMaxCols(bh: Blackhole): Unit = bh.consume(a.max(Axis.Cols))

private object MaxHillClimbRef:
  def refMax(n: Int, x: DoubleArray): Double =
    var m0 = Double.NegativeInfinity
    var m1 = Double.NegativeInfinity
    var m2 = Double.NegativeInfinity
    var m3 = Double.NegativeInfinity
    val limit = n - (n & 3)
    var i = 0
    while i < limit do
      m0 = math.max(m0, x(i))
      m1 = math.max(m1, x(i + 1))
      m2 = math.max(m2, x(i + 2))
      m3 = math.max(m3, x(i + 3))
      i += 4
    var m = math.max(math.max(m0, m1), math.max(m2, m3))
    while i < n do
      m = math.max(m, x(i))
      i += 1
    if m != 0.0 && m == m then m else DoubleKernels.firstZeroOrNaN(x, 0, 1, m)

  def breezeShape(n: Int, x: DoubleArray): Double =
    var m = Double.NegativeInfinity
    var i = 0
    while i < n do
      m = math.max(m, x(i))
      i += 1
    if m != 0.0 && m == m then m else DoubleKernels.firstZeroOrNaN(x, 0, 1, m)

  // The pre-change streamed pass, without the (rarely taken) zero/NaN resolution.
  def refCols(side: Int, x: DoubleArray, out: DoubleArray): Unit =
    var line = 0
    while line < side do
      out(line) = x(line)
      line += 1
    var k = 1
    while k < side do
      val start = k * side
      line = 0
      while line < side do
        out(line) = math.max(out(line), x(start + line))
        line += 1
      k += 1
