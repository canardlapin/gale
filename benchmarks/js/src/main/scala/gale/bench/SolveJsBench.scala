package gale.bench

import gale.backend.PureBackend
import gale.linalg.*
import scala.scalajs.js

/** Scala.js wall-clock pairing of the pre-`dtrsmLeft` multi-right-hand-side
  * solves ([[SolveRef]]) with the live `LU.solve` / `Cholesky.solve`:
  * alternating rounds, median ns per call. Run with full optimisation:
  * `sbt "set benchmarksJS/scalaJSStage := FullOptStage" "benchmarksJS/runMain gale.bench.SolveJsBench"`.
  */
object SolveJsBench:
  private val perf = js.Dynamic.global.performance

  def main(args: Array[String]): Unit =
    given gale.backend.Backend = PureBackend
    for (n, k) <- Seq((8, 1), (16, 16), (64, 16), (256, 64)) do
      val rng = new scala.util.Random(11L)
      val m = Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)
      val a = Matrix.tabulate(n, n)((i, j) => if i == j then n.toDouble else m(i, j))
      val s = (m * m.t) + Matrix.eye(n) * n.toDouble
      val b = Matrix.tabulate(n, k)((_, _) => rng.nextDouble() * 2.0 - 1.0)
      val lu = a.lu.toOption.get
      val ch = s.cholesky.toOption.get
      def pair(name: String)(ref: => Double)(cur: => Double): Unit =
        val r1 = measure(ref)
        val c1 = measure(cur)
        val c2 = measure(cur)
        val r2 = measure(ref)
        println(
          f"$name%-10s n=$n%-4d k=$k%-4d ref ${math.min(r1, r2)}%11.0f ns  cur ${math.min(c1, c2)}%11.0f ns  speedup ${math.min(r1, r2) / math.min(c1, c2)}%5.2fx"
        )
      pair("luSolve")(SolveRef.luSolve(lu, b)(0, 0))(lu.solve(b).toOption.get(0, 0))
      pair("cholSolve")(SolveRef.choleskySolve(ch, b)(0, 0))(ch.solve(b).toOption.get(0, 0))

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
