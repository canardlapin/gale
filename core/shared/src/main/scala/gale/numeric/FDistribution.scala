package gale.numeric

/** Rejected distribution parameter or argument. */
enum DistributionError:
  case InvalidParameter(name: String, value: Double)
  case InvalidArgument(name: String, value: Double)
  def message: String = this match
    case InvalidParameter(n, v) => s"$n must be positive (possibly +Infinity), got $v"
    case InvalidArgument(n, v)  => s"$n must not be NaN and must lie in its domain, got $v"

/** Both tails of a regularized incomplete beta integral, each computed without forming `1 - other` when it is the small
  * one: `lower = I_x(a, b)` and `upper = 1 - I_x(a, b) = I_{1-x}(b, a)`.
  */
final case class BetaTails(lower: Double, upper: Double)

/** Regularized incomplete beta function `I_x(a, b)` for real `a, b > 0`. */
object IncompleteBeta:
  import DistributionError.*

  /** Both tails of `I_x(a, b)`. `x` outside `[0, 1]`, NaN `x`, or a non-positive/non-finite shape is rejected. */
  def regularized(a: Double, b: Double, x: Double): Either[DistributionError, BetaTails] =
    if !(a > 0.0) || a.isInfinite then Left(InvalidParameter("a", a))
    else if !(b > 0.0) || b.isInfinite then Left(InvalidParameter("b", b))
    else if !(x >= 0.0 && x <= 1.0) then Left(InvalidArgument("x", x))
    else
      val (lower, upper) = SpecialFunctions.betaTails(a, b, x, 1.0 - x, x - a / (a + b))
      Right(BetaTails(lower, upper))

/** Snedecor F distribution with real degrees of freedom `df1, df2 > 0`; either may be `+Infinity` (the chi-square
  * limits), as in R's `pf`. Non-integer degrees of freedom (Satterthwaite, Kenward-Roger) are supported.
  *
  * Both [[cdf]] and [[survival]] evaluate their own tail through the regularized incomplete beta function, always
  * computing the smaller tail directly, so far upper-tail p-values keep full relative accuracy instead of cancelling in
  * `1 - cdf`. Against R's `pf` over degrees of freedom in `[0.1, 1e6]` the relative error is typically below `1e-12`;
  * the worst observed (about `3e-11`) pairs one degree of freedom above `1e5` with a small other one, where the
  * continued fraction is ill-conditioned. Tails below `Double.MinPositiveValue` underflow to zero, as in R.
  */
final class FDistribution private (val df1: Double, val df2: Double):

  /** `P(F > f)`, the upper-tail p-value. `f <= 0` gives 1, `f = +Infinity` gives 0, NaN propagates. */
  def survival(f: Double): Double = tails(f)._2

  /** `P(F <= f)`. `f <= 0` gives 0, `f = +Infinity` gives 1, NaN propagates. */
  def cdf(f: Double): Double = tails(f)._1

  private def tails(f: Double): (Double, Double) =
    if f.isNaN then (Double.NaN, Double.NaN)
    else if f <= 0.0 then (0.0, 1.0)
    else if f == Double.PositiveInfinity then (1.0, 0.0)
    else if df1.isInfinite && df2.isInfinite then
      if f < 1.0 then (0.0, 1.0) else if f > 1.0 then (1.0, 0.0) else (0.5, 0.5)
    else if df2.isInfinite then
      // F -> chi2(df1) / df1.
      SpecialFunctions.gammaTails(0.5 * df1, 0.5 * df1 * f)
    else if df1.isInfinite then
      // F -> df2 / chi2(df2), so P(F <= f) = P(chi2(df2) >= df2 / f).
      val (lower, upper) = SpecialFunctions.gammaTails(0.5 * df2, 0.5 * df2 / f)
      (upper, lower)
    else
      // P(F <= f) = I_x(df1/2, df2/2) with x = df1 f / (df1 f + df2); x, 1 - x and x - E[x] formed without
      // cancellation. Large df1 * f overflows only in the ratio, which is handled by scaling.
      val (x, y, dx) =
        val p = df1 * f
        if p.isInfinite then
          val s = df2 / df1
          val r = s / f
          (1.0 / (1.0 + r), r / (1.0 + r), 1.0 / (1.0 + r) - 1.0 / (1.0 + s))
        else
          val denom = p + df2
          val n = df1 + df2
          (p / denom, df2 / denom, ((df1 / denom) * (df2 / n)) * (f - 1.0))
      SpecialFunctions.betaTails(0.5 * df1, 0.5 * df2, x, y, dx)

  override def toString: String = s"FDistribution($df1, $df2)"
  override def equals(other: Any): Boolean = other match
    case that: FDistribution => df1 == that.df1 && df2 == that.df2
    case _                   => false
  override def hashCode: Int = (df1, df2).##

object FDistribution:
  import DistributionError.*

  def apply(df1: Double, df2: Double): Either[DistributionError, FDistribution] =
    if !(df1 > 0.0) then Left(InvalidParameter("df1", df1))
    else if !(df2 > 0.0) then Left(InvalidParameter("df2", df2))
    else Right(new FDistribution(df1, df2))

  /** Upper-tail p-value `P(F(df1, df2) > f)`, matching R's `pf(f, df1, df2, lower.tail = FALSE)`. NaN `f` is rejected;
    * for repeated evaluation construct the distribution once and call [[FDistribution.survival]].
    */
  def upperTail(f: Double, df1: Double, df2: Double): Either[DistributionError, Double] =
    if f.isNaN then Left(InvalidArgument("f", f))
    else apply(df1, df2).map(_.survival(f))

/** Internal special-function kernels. Callers validate arguments; these assume finite positive shapes. */
private[numeric] object SpecialFunctions:
  private val MaxIterations = 1_000_000
  private val Epsilon = 1e-16
  private val Tiny = 1e-300
  private val HalfLogTwoPi = 0.9189385332046727417803297

  /** Below this shape `logGamma` is used directly; at or above it the Stirling remainder is. */
  private val StirlingThreshold = 10.0

  // Lanczos approximation (g = 7, n = 9), relative error near 1e-15 for x > 0.
  private val LanczosG = 7.0
  private val Lanczos = Array(
    0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
    12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
  )

  def logGamma(x: Double): Double =
    if x < 0.5 then
      // Reflection keeps tiny positive shapes accurate: lgamma(x) = log(pi / sin(pi x)) - lgamma(1 - x).
      math.log(math.Pi / math.sin(math.Pi * x)) - logGamma(1.0 - x)
    else
      val z = x - 1.0
      var sum = Lanczos(0)
      var i = 1
      while i < Lanczos.length do
        sum += Lanczos(i) / (z + i)
        i += 1
      val t = z + LanczosG + 0.5
      HalfLogTwoPi + (z + 0.5) * math.log(t) - t + math.log(sum)

  /** Stirling remainder `lgamma(z) - [(z - 1/2) log z - z + log(2 pi)/2]` for `z >= 10`. */
  def stirlingRemainder(z: Double): Double =
    val r = 1.0 / z
    val r2 = r * r
    r * (1.0 / 12.0 - r2 * (1.0 / 360.0 - r2 * (1.0 / 1260.0 - r2 * (1.0 / 1680.0 - r2 * (1.0 / 1188.0 -
      r2 * (691.0 / 360360.0 - r2 / 156.0))))))

  /** `log1p(u) - u`, accurate for small `|u|`; `ratio = 1 + u` is supplied exactly by the caller. */
  def log1pMinus(u: Double, ratio: Double): Double =
    if math.abs(u) < 0.25 then
      // -u^2/2 + u^3/3 - ...; |u| < 1/4 needs at most ~25 terms.
      var term = u
      var sum = 0.0
      var k = 2
      var done = false
      while !done do
        term *= -u
        val next = term / k
        sum += next
        done = math.abs(next) <= Epsilon * math.abs(sum)
        k += 1
      sum
    else math.log(ratio) - u

  /** `log(x^a y^b / B(a, b))` with `y = 1 - x` and `dx = x - a / (a + b)` supplied without cancellation. */
  def logBetaPrefactor(a: Double, b: Double, x: Double, y: Double, dx: Double): Double =
    if a >= StirlingThreshold && b >= StirlingThreshold then
      // Centre at the mean x0 = a / (a + b): the large linear terms cancel analytically.
      val n = a + b
      val u = dx * n / a
      val v = -dx * n / b
      a * log1pMinus(u, x * n / a) + b * log1pMinus(v, y * n / b) +
        0.5 * math.log(a / n * b / (2.0 * math.Pi)) - stirlingRemainder(a) - stirlingRemainder(b) +
        stirlingRemainder(n)
    else if b >= StirlingThreshold then smallLarge(a, b, x, y)
    else if a >= StirlingThreshold then smallLarge(b, a, y, x)
    else a * math.log(x) + b * math.log(y) - (logGamma(a) + logGamma(b) - logGamma(a + b))

  /** Prefactor for small `a` and large `b`: lgamma(b) - lgamma(a + b) via Stirling remainders. */
  private def smallLarge(a: Double, b: Double, x: Double, y: Double): Double =
    val n = a + b
    val xn = x * n
    // b log y + (b - 1/2) log1p(a/b) - a, regrouped as b (log1p(w) - w) - x n with 1 + w = y n / b.
    val w = (a - xn) / b
    a * math.log(xn) - xn + b * log1pMinus(w, y * n / b) - 0.5 * math.log1p(a / b) - logGamma(a) -
      stirlingRemainder(b) + stirlingRemainder(n)

  /** Modified Lentz evaluation of the incomplete beta continued fraction for `I_x(a, b)`, with `y = 1 - x`. */
  private def betaContinuedFraction(a: Double, b: Double, x: Double, y: Double): Double =
    val qab = a + b
    val qap = a + 1.0
    val qam = a - 1.0
    var c = 1.0
    // 1 - (a + b) x / (a + 1), rewritten as (1 - b + (a + b) y) / (a + 1) when x > y so the rounded product
    // multiplies the smaller of x and y; otherwise x near 1 with a large shape loses digits to cancellation.
    var d = if x <= y then 1.0 - qab * x / qap else (1.0 - b + qab * y) / qap
    if math.abs(d) < Tiny then d = Tiny
    d = 1.0 / d
    var h = d
    var m = 1
    var done = false
    while !done && m <= MaxIterations do
      val m2 = 2 * m
      var aa = m * (b - m) * x / ((qam + m2) * (a + m2))
      d = 1.0 + aa * d
      if math.abs(d) < Tiny then d = Tiny
      c = 1.0 + aa / c
      if math.abs(c) < Tiny then c = Tiny
      d = 1.0 / d
      h *= d * c
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
      d = 1.0 + aa * d
      if math.abs(d) < Tiny then d = Tiny
      c = 1.0 + aa / c
      if math.abs(c) < Tiny then c = Tiny
      d = 1.0 / d
      val delta = d * c
      h *= delta
      done = math.abs(delta - 1.0) <= Epsilon
      m += 1
    h

  /** `(I_x(a, b), I_y(b, a))` with `y = 1 - x` and `dx = x - a / (a + b)`. The tail on the convergent side of the
    * continued fraction is computed directly; the other is its complement, which is then at least about one half.
    */
  def betaTails(a: Double, b: Double, x: Double, y: Double, dx: Double): (Double, Double) =
    if x <= 0.0 then (0.0, 1.0)
    else if y <= 0.0 then (1.0, 0.0)
    else
      val logPrefactor = logBetaPrefactor(a, b, x, y, dx)
      if x * (a + b + 2.0) < a + 1.0 then
        val lower = math.exp(logPrefactor) * betaContinuedFraction(a, b, x, y) / a
        (lower, 1.0 - lower)
      else
        val upper = math.exp(logPrefactor) * betaContinuedFraction(b, a, y, x) / b
        (1.0 - upper, upper)

  /** `log(z^s e^-z / Gamma(s))`. */
  private def logGammaPrefactor(s: Double, z: Double): Double =
    if s >= StirlingThreshold then
      val u = (z - s) / s
      s * log1pMinus(u, z / s) + 0.5 * math.log(s / (2.0 * math.Pi)) - stirlingRemainder(s)
    else s * math.log(z) - z - logGamma(s)

  /** Regularized incomplete gamma tails `(P(s, z), Q(s, z))` for finite `s > 0`, `z >= 0`. */
  def gammaTails(s: Double, z: Double): (Double, Double) =
    if z <= 0.0 then (0.0, 1.0)
    else if z.isInfinite then (1.0, 0.0)
    else
      val logPrefactor = logGammaPrefactor(s, z)
      if z < s + 1.0 then
        // Series: P = prefactor * sum_k z^k / (s (s+1) ... (s+k)).
        var ap = s
        var term = 1.0 / s
        var sum = term
        var k = 0
        var done = false
        while !done && k < MaxIterations do
          ap += 1.0
          term *= z / ap
          sum += term
          done = math.abs(term) <= Epsilon * math.abs(sum)
          k += 1
        val lower = math.exp(logPrefactor) * sum
        (lower, 1.0 - lower)
      else
        // Modified Lentz continued fraction for Q.
        var b = z + 1.0 - s
        var c = 1.0 / Tiny
        var d = 1.0 / b
        var h = d
        var i = 1
        var done = false
        while !done && i <= MaxIterations do
          val an = -i * (i - s)
          b += 2.0
          d = an * d + b
          if math.abs(d) < Tiny then d = Tiny
          c = b + an / c
          if math.abs(c) < Tiny then c = Tiny
          d = 1.0 / d
          val delta = d * c
          h *= delta
          done = math.abs(delta - 1.0) <= Epsilon
          i += 1
        val upper = math.exp(logPrefactor) * h
        (1.0 - upper, upper)
