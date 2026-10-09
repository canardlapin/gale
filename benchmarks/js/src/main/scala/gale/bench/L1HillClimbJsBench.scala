package gale.bench

import gale.kernel.DoubleKernels
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.syntax.all.*
import scala.scalajs.js

/** Scala.js wall-clock pairing of the pre-hill-climb kernels ([[L1HillClimbRef]])
  * with the live [[DoubleKernels]]: alternating rounds, median ns per call.
  * Run with `benchSmokeJSFull`-style full optimisation:
  * `sbt "set benchmarksJS/scalaJSStage := FullOptStage" "benchmarksJS/runMain gale.bench.L1HillClimbJsBench"`.
  */
object L1HillClimbJsBench:
  private val perf = js.Dynamic.global.performance

  def main(args: Array[String]): Unit =
    for n <- Seq(1024, 65536) do
      val rng = new scala.util.Random(7L)
      val x = DoubleArray.alloc(n)
      val y = DoubleArray.alloc(n)
      val out = DoubleArray.alloc(n)
      var i = 0
      while i < n do
        x(i) = rng.nextDouble() * 2.0 - 1.0
        y(i) = rng.nextDouble() * 2.0 - 1.0
        i += 1
      val side = math.sqrt(n.toDouble).toInt
      val a = Matrix.tabulate(side, side)((r, c) => x(r * side + c))
      val b = Matrix.tabulate(side, side)((r, c) => y(r * side + c))
      val maxima = DoubleArray.alloc(side)
      var sink = 0.0
      // Measure in both orders and keep each side's best median, so neither
      // side profits from running first.
      def pair(name: String)(ref: => Double)(cur: => Double): Unit =
        val r1 = measure(ref)
        val c1 = measure(cur)
        val c2 = measure(cur)
        val r2 = measure(ref)
        report(name, n, math.min(r1, r2), math.min(c1, c2))
      pair("axpy") {
        L1HillClimbRef.refDaxpy(n, 1.5, x, 0, 1, y, 0, 1); L1HillClimbRef.refDaxpy(n, -1.5, x, 0, 1, y, 0, 1); y(0)
      } { DoubleKernels.daxpy(n, 1.5, x, 0, 1, y, 0, 1); DoubleKernels.daxpy(n, -1.5, x, 0, 1, y, 0, 1); y(0) }
      pair("dot")(L1HillClimbRef.refDdot(n, x, 0, 1, y, 0, 1))(DoubleKernels.ddot(n, x, 0, 1, y, 0, 1))
      pair("nrm2")(L1HillClimbRef.refDnrm2(n, x, 0, 1))(DoubleKernels.dnrm2(n, x, 0, 1))
      pair("max")(x(L1HillClimbRef.refMaxIndex(n, x, 0, 1)))(DoubleKernels.dmax(n, x, 0, 1))
      pair("add") { L1HillClimbRef.refDadd(n, x, 0, 1, y, 0, 1, out, 0, 1); out(1) } {
        DoubleKernels.dadd(n, x, 0, 1, y, 0, 1, out, 0, 1); out(1)
      }
      pair("hadamard")(Matrix.tabulate(side, side)((r, c) => a(r, c) * b(r, c))(0, 0))((a.pointwise * b)(0, 0))
      pair("maxCols") { L1HillClimbRef.refStreamMax(side, side, x, maxima); maxima(0) }(a.max(Axis.Cols)(0))
      pair("sigmoid") { L1HillClimbRef.refSigmoidInto(n, x, out); out(1) } {
        DoubleKernels.dsigmoidInto(n, x, 0, 1, out, 0, 1); out(1)
      }
      pair("exp") { L1HillClimbRef.refExpInto(n, x, out); out(1) } { DoubleKernels.dexpInto(n, x, 0, 1, out, 0, 1); out(1) }
      pair("logSumExp")(L1HillClimbRef.refLogSumExp(n, x))(DoubleKernels.dlogSumExp(n, x, 0, 1))
      sink += out(0)
      if sink == 42.0 then println(sink)

  /** Median ns per call over 15 rounds of enough calls to fill ~20 ms. */
  private def measure(body: => Double): Double =
    var sink = 0.0
    var calls = 1
    var elapsed = 0.0
    while elapsed < 20.0 do
      calls *= 2
      val start = perf.now().asInstanceOf[Double]
      var k = 0
      while k < calls do
        sink += body
        k += 1
      elapsed = perf.now().asInstanceOf[Double] - start
    val rounds = Array.fill(15) {
      val start = perf.now().asInstanceOf[Double]
      var k = 0
      while k < calls do
        sink += body
        k += 1
      (perf.now().asInstanceOf[Double] - start) * 1e6 / calls
    }
    if sink == 42.0 then println(sink)
    rounds.sorted.apply(7)

  private def report(name: String, n: Int, ref: Double, cur: Double): Unit =
    println(f"$name%-10s n=$n%-7d ref ${ref}%10.0f ns  cur ${cur}%10.0f ns  speedup ${ref / cur}%5.2fx")
