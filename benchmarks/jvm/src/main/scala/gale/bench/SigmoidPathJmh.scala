package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseVector as BDV
import gale.bench.BreezeBenchData.*
import gale.kernel.DoubleKernels
import gale.linalg.{DVec, Numerics}
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Attribution bench for the lane-A sigmoid gap (`ReductionBreezeJmh.galeSigmoid` vs `breezeSigmoid`), on the same
  * data (`vectorData(n, 5300)`, uniform on [-1, 1), so no element reaches the `x < -36` branch).
  *
  *   - `breeze`, `gale`: the two lane-A methods.
  *   - `rawHybrid`: the live kernel in a plain loop over a raw array with a fresh output, i.e. Breeze's loop shape;
  *     `gale - rawHybrid` is the facade and loop-shape cost.
  *   - `rawHybrid4`: `rawHybrid` unrolled by four by hand.
  *   - `rawDirect`: Breeze's form `1/(1 + exp(-x))` alone in the same loop. It is not gale's kernel (it overflows
  *     and loses accuracy below -36); `rawDirect - rawHybrid` is the cost of the per-element range test.
  *   - `scanMin`, `scanMin4`, `scanBits`: candidate prescans that would let a bit-identical kernel run the direct
  *     form alone when nothing is below -36, reported per element. They only pay if C2 vectorizes them: one, four
  *     accumulator `Math.min` reductions, and an OR of a sign-and-magnitude bit test (negative with `|x| >= 32`).
  *   - `prescanMin4`, `prescanBits`: those prescans followed by the direct loop (or the hybrid when they fire).
  *     Bit-identical to the live kernel.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SigmoidPathJmh:
  @Param(Array("1024", "65536", "1048576"))
  var n: Int = 0

  private var gx: DVec           = uninitialized
  private var bx: BDV[Double]    = uninitialized
  private var raw: Array[Double] = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    val data = vectorData(n, 5300L)
    gx = galeVector(data)
    bx = breezeVector(data)
    raw = data.clone()

  @Benchmark def breeze(bh: Blackhole): Unit = bh.consume(_root_.breeze.numerics.sigmoid(bx))
  @Benchmark def gale(bh: Blackhole): Unit   = bh.consume(Numerics.sigmoid(gx))

  @Benchmark def rawHybrid(bh: Blackhole): Unit =
    val out = new Array[Double](raw.length)
    SigmoidPathJmh.hybridInto(raw, out, 0, raw.length)
    bh.consume(out)

  @Benchmark def rawDirect(bh: Blackhole): Unit =
    val out = new Array[Double](raw.length)
    SigmoidPathJmh.directInto(raw, out, 0, raw.length)
    bh.consume(out)

  @Benchmark def rawHybrid4(bh: Blackhole): Unit =
    val out = new Array[Double](raw.length)
    SigmoidPathJmh.hybrid4Into(raw, out, 0, raw.length)
    bh.consume(out)

  @Benchmark def scanMin(): Double  = SigmoidPathJmh.min1(raw)
  @Benchmark def scanMin4(): Double = SigmoidPathJmh.min4(raw)
  @Benchmark def scanBits(): Long   = SigmoidPathJmh.bits(raw)

  @Benchmark def prescanMin4(bh: Blackhole): Unit =
    val out = new Array[Double](raw.length)
    if SigmoidPathJmh.min4(raw) >= DoubleKernels.SigmoidDirectFloor then SigmoidPathJmh.directInto(raw, out, 0, raw.length)
    else SigmoidPathJmh.hybridInto(raw, out, 0, raw.length)
    bh.consume(out)

  @Benchmark def prescanBits(bh: Blackhole): Unit =
    val out = new Array[Double](raw.length)
    if SigmoidPathJmh.bits(raw) >= 0L then SigmoidPathJmh.directInto(raw, out, 0, raw.length)
    else SigmoidPathJmh.hybridInto(raw, out, 0, raw.length)
    bh.consume(out)

object SigmoidPathJmh:
  def hybridInto(x: Array[Double], out: Array[Double], from: Int, until: Int): Unit =
    var i = from
    while i < until do
      out(i) = DoubleKernels.sigmoid(x(i))
      i += 1

  def hybrid4Into(x: Array[Double], out: Array[Double], from: Int, until: Int): Unit =
    var i = from
    while i + 3 < until do
      out(i) = DoubleKernels.sigmoid(x(i))
      out(i + 1) = DoubleKernels.sigmoid(x(i + 1))
      out(i + 2) = DoubleKernels.sigmoid(x(i + 2))
      out(i + 3) = DoubleKernels.sigmoid(x(i + 3))
      i += 4
    while i < until do
      out(i) = DoubleKernels.sigmoid(x(i))
      i += 1

  def directInto(x: Array[Double], out: Array[Double], from: Int, until: Int): Unit =
    var i = from
    while i < until do
      out(i) = 1.0 / (1.0 + math.exp(-x(i)))
      i += 1

  def min1(x: Array[Double]): Double =
    var m = Double.PositiveInfinity
    var i = 0
    while i < x.length do
      m = math.min(m, x(i))
      i += 1
    m

  def min4(x: Array[Double]): Double =
    var m0 = Double.PositiveInfinity
    var m1 = m0
    var m2 = m0
    var m3 = m0
    var i = 0
    while i + 3 < x.length do
      m0 = math.min(m0, x(i))
      m1 = math.min(m1, x(i + 1))
      m2 = math.min(m2, x(i + 2))
      m3 = math.min(m3, x(i + 3))
      i += 4
    while i < x.length do
      m0 = math.min(m0, x(i))
      i += 1
    math.min(math.min(m0, m1), math.min(m2, m3))

  /** Negative iff some element is negative with `|x| >= 32` (including -Inf and negative NaNs): adding
    * `2^63 - bits(32.0)` to the magnitude bits carries into bit 63 exactly when `|x| >= 32`, and the AND keeps it
    * only for a set sign bit. A non-negative result proves no element is below -36.
    */
  def bits(x: Array[Double]): Long =
    var acc = 0L
    var i = 0
    while i < x.length do
      val b = java.lang.Double.doubleToRawLongBits(x(i))
      acc |= b & ((b & 0x7fffffffffffffffL) + 0x3fc0000000000000L)
      i += 1
    acc
