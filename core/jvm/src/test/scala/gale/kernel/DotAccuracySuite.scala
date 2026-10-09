package gale.kernel

import gale.numeric.ExactSum
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Accuracy observation for the contiguous dot product, which uses unfused
  * `x*y + acc` accumulators (the JVM-only fused form was ~20% slower). On
  * cancelling, ill-conditioned inputs it compares the unfused kernel and the
  * same four-accumulator loop with `Math.fma` against the exact dot product
  * (`TwoProduct` pieces summed by [[ExactSum]]). The printed medians are
  * recorded in the numerical contract; the only assertion is the classical
  * forward bound `|error| <= gamma_n * sum |x_i y_i|`, which both forms must meet.
  */
class DotAccuracySuite extends munit.FunSuite:
  private val n = 1000
  private val unit = math.ulp(1.0) / 2

  private def fusedDot(x: Array[Double], y: Array[Double]): Double =
    var acc0, acc1, acc2, acc3 = 0.0
    val limit = n - (n & 3)
    var i = 0
    while i < limit do
      acc0 = Math.fma(x(i), y(i), acc0)
      acc1 = Math.fma(x(i + 1), y(i + 1), acc1)
      acc2 = Math.fma(x(i + 2), y(i + 2), acc2)
      acc3 = Math.fma(x(i + 3), y(i + 3), acc3)
      i += 4
    var acc = (acc0 + acc1) + (acc2 + acc3)
    while i < n do
      acc = Math.fma(x(i), y(i), acc)
      i += 1
    acc

  /** Exact dot product, rounded once: each product splits exactly into
    * `p + e` with `e = fma(x, y, -p)`.
    */
  private def exactDot(x: Array[Double], y: Array[Double]): Double =
    val sum = ExactSum.zero()
    var i = 0
    while i < n do
      val p = x(i) * y(i)
      sum.add(p).fold(e => throw new IllegalStateException(e.message), identity)
      sum.add(Math.fma(x(i), y(i), -p)).fold(e => throw new IllegalStateException(e.message), identity)
      i += 1
    sum.value

  /** Random terms spread over `2^±spread` whose last pair cancels all but a
    * fraction `keep` of the running sum (`keep = 0`: cancel to the rounding
    * residual), so the exact result is small next to `sum |x_i y_i|`.
    */
  private def cancelling(seed: Long, spread: Int, keep: Double): (Array[Double], Array[Double]) =
    val random = new scala.util.Random(seed)
    val x = Array.fill(n)(random.nextGaussian() * math.scalb(1.0, random.nextInt(2 * spread + 1) - spread))
    val y = Array.fill(n)(random.nextGaussian())
    x(n - 1) = 0.0
    y(n - 1) = 1.0
    x(n - 1) = -exactDot(x, y) * (1.0 - keep)
    (x, y)

  test("unfused and fused contiguous dot: accuracy on cancelling inputs (observation)") {
    for (spread, keep) <- Seq((0, 1e-3), (20, 1e-3), (20, 1e-8), (0, 0.0), (40, 0.0)) do
      val rows = (0 until 25).map { trial =>
        val (x, y) = cancelling(9000L + 31L * spread + trial, spread, keep)
        val exact = exactDot(x, y)
        val magnitude = x.indices.map(i => math.abs(x(i) * y(i))).sum
        val unfused = DoubleKernels.ddot(n, DoubleArray.adopt(x), 0, 1, DoubleArray.adopt(y), 0, 1)
        val fused = fusedDot(x, y)
        val bound = (n * unit / (1 - n * unit)) * magnitude
        assert(math.abs(unfused - exact) <= bound, s"unfused outside gamma_n bound, spread=$spread trial=$trial")
        assert(math.abs(fused - exact) <= bound, s"fused outside gamma_n bound, spread=$spread trial=$trial")
        // Errors relative to sum |x_i y_i| (condition-free scale).
        (math.abs(unfused - exact) / magnitude, math.abs(fused - exact) / magnitude, 2 * magnitude / math.abs(exact))
      }
      def median(xs: Seq[Double]) = xs.sorted.apply(xs.length / 2)
      println(
        f"[dot-accuracy] spread=2^±$spread%-2d keep=$keep%-6.0e median cond=${median(rows.map(_._3))}%.1e " +
          f"unfused err/sum|xy|=${median(rows.map(_._1))}%.2e fused=${median(rows.map(_._2))}%.2e " +
          f"(unit roundoff ${unit}%.2e)"
      )
  }
