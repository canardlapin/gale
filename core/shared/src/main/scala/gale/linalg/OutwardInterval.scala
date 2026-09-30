package gale.linalg

/** A finite, closed real interval whose endpoints are IEEE-754 `Double`s.
  *
  * This is deliberately a small certification primitive. It refuses NaN,
  * infinities, unordered endpoints, a zero-containing divisor, and any
  * intermediate operation which overflows to a non-finite `Double`. Arithmetic
  * expands each finite primitive result by one representable value on the
  * appropriate side with `Math.nextDown`/`Math.nextUp`.
  */
final class OutwardInterval private (val lower: Double, val upper: Double):
  def contains(value: Double): Boolean = value.isFinite && lower <= value && value <= upper

  def +(that: OutwardInterval): Either[OutwardIntervalError, OutwardInterval] =
    OutwardInterval.fromRounded(lower + that.lower, upper + that.upper)

  def -(that: OutwardInterval): Either[OutwardIntervalError, OutwardInterval] =
    OutwardInterval.fromRounded(lower - that.upper, upper - that.lower)

  def *(that: OutwardInterval): Either[OutwardIntervalError, OutwardInterval] =
    OutwardInterval.encloseFour(lower * that.lower, lower * that.upper, upper * that.lower, upper * that.upper)

  def /(that: OutwardInterval): Either[OutwardIntervalError, OutwardInterval] =
    if that.lower <= 0.0 && 0.0 <= that.upper then Left(OutwardIntervalError.ZeroContainingDivisor)
    else OutwardInterval.encloseFour(lower / that.lower, lower / that.upper, upper / that.lower, upper / that.upper)

  /** Enclose the square without using a square root or another transcendental. */
  def square: Either[OutwardIntervalError, OutwardInterval] =
    val crossesZero = lower <= 0.0 && 0.0 <= upper
    val low =
      if crossesZero then 0.0
      else if lower > 0.0 then lower * lower
      else upper * upper
    val high = math.max(lower * lower, upper * upper)
    if !low.isFinite || !high.isFinite then Left(OutwardIntervalError.NonFiniteIntermediate)
    else
      val roundedLower = if crossesZero then 0.0 else java.lang.Math.nextDown(low)
      OutwardInterval.apply(roundedLower, java.lang.Math.nextUp(high))

enum OutwardIntervalError:
  case NonFiniteEndpoint
  case UnorderedEndpoints
  case ZeroContainingDivisor
  case NonFiniteIntermediate
  case LengthMismatch(left: Int, right: Int)

object OutwardInterval:
  def apply(lower: Double, upper: Double): Either[OutwardIntervalError, OutwardInterval] =
    if !lower.isFinite || !upper.isFinite then Left(OutwardIntervalError.NonFiniteEndpoint)
    else if lower > upper then Left(OutwardIntervalError.UnorderedEndpoints)
    else Right(new OutwardInterval(lower, upper))

  def point(value: Double): Either[OutwardIntervalError, OutwardInterval] =
    apply(value, value)

  /** One-pass, constant-extra-storage enclosure of `sum_i left(i) * right(i)`. */
  def streamingDot(left: Array[Double], right: Array[Double]): Either[OutwardIntervalError, OutwardInterval] =
    if left.length != right.length then Left(OutwardIntervalError.LengthMismatch(left.length, right.length))
    else
      var total = point(0.0)
      var index = 0
      while index < left.length && total.isRight do
        total = for
          sum <- total
          x <- point(left(index))
          y <- point(right(index))
          product <- x * y
          next <- sum + product
        yield next
        index += 1
      total

  /** One-pass, constant-extra-storage enclosure of `sum_i abs(values(i))`. */
  def streamingAbsoluteSum(values: Array[Double]): Either[OutwardIntervalError, OutwardInterval] =
    var total = point(0.0)
    var index = 0
    while index < values.length && total.isRight do
      val absolute = if values(index) < 0.0 then -values(index) else values(index)
      total = for
        sum <- total
        term <- point(absolute)
        next <- sum + term
      yield next
      index += 1
    total

  private[linalg] def fromRounded(lower: Double, upper: Double): Either[OutwardIntervalError, OutwardInterval] =
    if !lower.isFinite || !upper.isFinite then Left(OutwardIntervalError.NonFiniteIntermediate)
    else apply(java.lang.Math.nextDown(lower), java.lang.Math.nextUp(upper))

  private[linalg] def encloseFour(a: Double, b: Double, c: Double, d: Double): Either[OutwardIntervalError, OutwardInterval] =
    if !a.isFinite || !b.isFinite || !c.isFinite || !d.isFinite then Left(OutwardIntervalError.NonFiniteIntermediate)
    else fromRounded(math.min(math.min(a, b), math.min(c, d)), math.max(math.max(a, b), math.max(c, d)))
