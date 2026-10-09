package gale.bench

import scala.compiletime.uninitialized

import gale.kernel.DoubleKernels
import gale.linalg.{Axis, DMat, Matrix}
import gale.syntax.all.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Gale-only before/after microbench for the pure L1, reduction and elementwise
  * kernels. Each `ref*` method runs a private copy of the kernel as it stood
  * before the hill climb (breeze/integration d102321); each `cur*` method runs
  * the live [[DoubleKernels]] kernel on the same data, so one JMH run gives a
  * paired before/after ratio.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 4, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(3)
@State(Scope.Thread)
class L1HillClimbJmh:
  @Param(Array("1024", "65536", "1048576"))
  var n: Int = 0

  private val alpha = 1.5

  private var x: DoubleArray   = uninitialized
  private var y: DoubleArray   = uninitialized
  private var out: DoubleArray = uninitialized
  private var side: Int        = 0
  private var a: DMat          = uninitialized
  private var b: DMat          = uninitialized

  @Setup(Level.Trial)
  def setup(): Unit =
    val rng = new java.util.Random(7L)
    x = DoubleArray.adopt(Array.fill(n)(rng.nextDouble() * 2.0 - 1.0))
    y = DoubleArray.adopt(Array.fill(n)(rng.nextDouble() * 2.0 - 1.0))
    out = DoubleArray.alloc(n)
    side = math.sqrt(n.toDouble).toInt
    a = Matrix.tabulate(side, side)((i, j) => x(i * side + j))
    b = Matrix.tabulate(side, side)((i, j) => y(i * side + j))

  // Each invocation adds alpha*x to y and then subtracts it again, so y stays
  // bounded over the trial and both variants see identical data.
  @Benchmark def refAxpy(): Double =
    L1HillClimbRef.refDaxpy(n, alpha, x, 0, 1, y, 0, 1)
    L1HillClimbRef.refDaxpy(n, -alpha, x, 0, 1, y, 0, 1)
    y(0)
  @Benchmark def curAxpy(): Double =
    DoubleKernels.daxpy(n, alpha, x, 0, 1, y, 0, 1)
    DoubleKernels.daxpy(n, -alpha, x, 0, 1, y, 0, 1)
    y(0)

  @Benchmark def refDot(): Double = L1HillClimbRef.refDdot(n, x, 0, 1, y, 0, 1)
  @Benchmark def curDot(): Double = DoubleKernels.ddot(n, x, 0, 1, y, 0, 1)

  @Benchmark def refNrm2(): Double = L1HillClimbRef.refDnrm2(n, x, 0, 1)
  @Benchmark def curNrm2(): Double = DoubleKernels.dnrm2(n, x, 0, 1)

  @Benchmark def refMax(): Double = x(L1HillClimbRef.refMaxIndex(n, x, 0, 1))
  @Benchmark def curMax(): Double = DoubleKernels.dmax(n, x, 0, 1)

  @Benchmark def refAdd(bh: Blackhole): Unit =
    L1HillClimbRef.refDadd(n, x, 0, 1, y, 0, 1, out, 0, 1)
    bh.consume(out)
  @Benchmark def curAdd(bh: Blackhole): Unit =
    DoubleKernels.dadd(n, x, 0, 1, y, 0, 1, out, 0, 1)
    bh.consume(out)

  // Matrix ops on side x side = sqrt(n) square operands (32, 256, 1024).
  @Benchmark def refHadamard(): DMat = Matrix.tabulate(side, side)((i, j) => a(i, j) * b(i, j))
  @Benchmark def curHadamard(): DMat = a.pointwise * b

  @Benchmark def refMaxCols(bh: Blackhole): Unit =
    val res = DoubleArray.alloc(side)
    L1HillClimbRef.refStreamMax(side, side, x, res)
    bh.consume(res)
  @Benchmark def curMaxCols(bh: Blackhole): Unit = bh.consume(a.max(Axis.Cols))

  @Benchmark def refSigmoid(bh: Blackhole): Unit =
    L1HillClimbRef.refSigmoidInto(n, x, out)
    bh.consume(out)
  @Benchmark def curSigmoid(bh: Blackhole): Unit =
    DoubleKernels.dsigmoidInto(n, x, 0, 1, out, 0, 1)
    bh.consume(out)

  @Benchmark def refExp(bh: Blackhole): Unit =
    L1HillClimbRef.refExpInto(n, x, out)
    bh.consume(out)
  @Benchmark def curExp(bh: Blackhole): Unit =
    DoubleKernels.dexpInto(n, x, 0, 1, out, 0, 1)
    bh.consume(out)

  @Benchmark def refLogSumExp(): Double = L1HillClimbRef.refLogSumExp(n, x)
  @Benchmark def curLogSumExp(): Double = DoubleKernels.dlogSumExp(n, x, 0, 1)
