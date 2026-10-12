package gale.bench

import scala.compiletime.uninitialized

import gale.kernel.DoubleKernels
import gale.linalg.{Axis, DMat, DVec, Matrix, value}
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
  *
  * x86 candidates, each also run directly so one run compares them on any host:
  * `scanMax` (the scan), `reduceMax` (the `math.max` reduction), `scan4Max` (the
  * scan testing four entries per branch) and `keyMax` (branch-free integer max
  * over order-preserving keys, a candidate for C2 auto-vectorization).
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
  private var v: DVec          = uninitialized
  private var side: Int        = 0
  private var a: DMat          = uninitialized
  private var out: DoubleArray = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    x = DoubleArray.adopt(BreezeBenchData.vectorData(n, 5300L))
    side = math.sqrt(n.toDouble).toInt
    a = Matrix.tabulate(side, side)((i, j) => x(i * side + j))
    out = DoubleArray.alloc(side)
    v = BreezeBenchData.galeVector(BreezeBenchData.vectorData(n, 5300L))
    MaxHillClimbRef.selfCheck(n, x)

  @Benchmark def refMax(): Double         = MaxHillClimbRef.refMax(n, x)
  @Benchmark def breezeShapeMax(): Double = MaxHillClimbRef.breezeShape(n, x)
  @Benchmark def curMax(): Double         = DoubleKernels.dmax(n, x, 0, 1)
  @Benchmark def scanMax(): Double        = DoubleKernels.dmaxScan(n, x, 0, 1)
  @Benchmark def reduceMax(): Double      = DoubleKernels.dmaxReduce(n, x, 0, 1)
  @Benchmark def keyMax(): Double         = MaxHillClimbRef.keyMax(n, x)
  @Benchmark def scan4Max(): Double       = MaxHillClimbRef.scan4(n, x)

  // Facade A/B on the scoreboard's vector: `facadeMax` is `DVec.max`;
  // `generalMax` is its call before the dense entry point (offset and stride
  // read from the vector at run time).
  @Benchmark def facadeMax(): Double = v.max
  @Benchmark def generalMax(): Double =
    DoubleKernels.dmax(v.length, v.data, v.offset.value, v.stride.value)

  // Column maxima of the side x side row-major matrix: the streamed geometry.
  @Benchmark def refMaxCols(bh: Blackhole): Unit =
    MaxHillClimbRef.refCols(side, x, out)
    bh.consume(out)
  @Benchmark def curMaxCols(bh: Blackhole): Unit = bh.consume(a.max(Axis.Cols))

private object MaxHillClimbRef:
  /** Throws unless every max candidate returns the same bits as the scan, on the
    * benchmark data and on lines with NaNs of either sign, signed zeros,
    * infinities and ties planted at the first, middle and last positions.
    */
  def selfCheck(n: Int, data: DoubleArray): Unit =
    def raw(v: Double) = java.lang.Double.doubleToRawLongBits(v)
    def check(len: Int, line: DoubleArray, what: String): Unit =
      val expected = raw(DoubleKernels.dmaxScan(len, line, 0, 1))
      val candidates = Seq(
        "dmax" -> DoubleKernels.dmax(len, line, 0, 1),
        "reduce" -> DoubleKernels.dmaxReduce(len, line, 0, 1),
        "keyMax" -> keyMax(len, line),
        "scan4" -> scan4(len, line)
      )
      for (name, value) <- candidates if raw(value) != expected do
        throw new IllegalStateException(
          s"MaxHillClimbJmh self-check: $name on $what gave ${raw(value).toHexString}, scan ${expected.toHexString}"
        )
    check(n, data, s"bench data n=$n")
    val nanA = java.lang.Double.longBitsToDouble(0x7ff8000000000001L)
    val nanB = java.lang.Double.longBitsToDouble(0xfff8000000000002L)
    val random = new scala.util.Random(20261011L)
    for len <- Seq(1, 3, 4, 5, 8, 63, 64, 65, 1031) do
      val positions = Seq(0, len / 2, len - 1).distinct
      val plants: Seq[(String, Array[Double] => Unit)] = Seq(
        "none" -> (_ => ()),
        "nan" -> (a => a(len / 2) = nanA),
        "neg nan then nan" -> (a => { a(0) = nanB; a(len - 1) = nanA }),
        "+inf" -> (a => a(len - 1) = Double.PositiveInfinity),
        "all -inf" -> (a => java.util.Arrays.fill(a, Double.NegativeInfinity)),
        "-0 then 0" -> (a => { for i <- a.indices do a(i) = -1.0 - i; a(0) = -0.0; a(len - 1) = 0.0 }),
        "0 then -0" -> (a => { for i <- a.indices do a(i) = -1.0 - i; a(0) = 0.0; a(len - 1) = -0.0 }),
        "tie" -> (a => { a(0) = 1e9; a(len - 1) = 1e9 })
      )
      for (name, plant) <- plants; p <- positions do
        val a = Array.fill(len)(random.nextGaussian())
        a(p) = 2.0 * a(p)
        plant(a)
        check(len, DoubleArray.adopt(a), s"$name len=$len")

  // Max over order-preserving integer keys: branch-free integer max reductions
  // (one for the key, one for the magnitude to detect NaN of either sign).
  def keyMax(n: Int, x: DoubleArray): Double =
    var kmax = Long.MinValue
    var amax = 0L
    var i = 0
    while i < n do
      val b = java.lang.Double.doubleToRawLongBits(x(i))
      kmax = math.max(kmax, b ^ ((b >> 63) >>> 1))
      amax = math.max(amax, b & Long.MaxValue)
      i += 1
    if amax > 0x7ff0000000000000L then DoubleKernels.firstZeroOrNaN(x, 0, 1, Double.NaN)
    else
      val m = java.lang.Double.longBitsToDouble(kmax ^ ((kmax >> 63) >>> 1))
      if m != 0.0 then m else DoubleKernels.firstZeroOrNaN(x, 0, 1, m)

  // The scan, four entries per test: one branch when no entry of the four can
  // change the running maximum, the exact in-order update otherwise.
  def scan4(n: Int, x: DoubleArray): Double =
    var m = Double.NegativeInfinity
    var i = 0
    val limit = n - (n & 3)
    while i < limit do
      val v0 = x(i)
      val v1 = x(i + 1)
      val v2 = x(i + 2)
      val v3 = x(i + 3)
      if !(v0 <= m) || !(v1 <= m) || !(v2 <= m) || !(v3 <= m) then
        if !(v0 <= m) then
          if v0 != v0 then return v0
          m = v0
        if !(v1 <= m) then
          if v1 != v1 then return v1
          m = v1
        if !(v2 <= m) then
          if v2 != v2 then return v2
          m = v2
        if !(v3 <= m) then
          if v3 != v3 then return v3
          m = v3
      i += 4
    while i < n do
      val v = x(i)
      if !(v <= m) then
        if v != v then return v
        m = v
      i += 1
    m

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
