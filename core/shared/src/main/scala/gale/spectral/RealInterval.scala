package gale.spectral

import gale.linalg.LinAlgError

/** Finite outward-rounded real interval. Construction is checked. */
final class RealInterval private (val lower: Double, val upper: Double):
  def add(that: RealInterval): Either[LinAlgError, RealInterval] = RealInterval.checked(RealInterval.down(lower + that.lower), RealInterval.up(upper + that.upper))
  def subtract(that: RealInterval): Either[LinAlgError, RealInterval] = RealInterval.checked(RealInterval.down(lower - that.upper), RealInterval.up(upper - that.lower))
  def multiply(that: RealInterval): Either[LinAlgError, RealInterval] =
    val xs = Vector(lower * that.lower, lower * that.upper, upper * that.lower, upper * that.upper)
    if xs.exists(x => !x.isFinite) then RealInterval.fail else RealInterval.checked(RealInterval.down(xs.min), RealInterval.up(xs.max))
  def divide(that: RealInterval): Either[LinAlgError, RealInterval] =
    if that.lower <= 0.0 && that.upper >= 0.0 then Left(LinAlgError.InvalidArgument("interval division crosses zero"))
    else
      val xs = Vector(lower / that.lower, lower / that.upper, upper / that.lower, upper / that.upper)
      if xs.exists(x => !x.isFinite) then RealInterval.fail else RealInterval.checked(RealInterval.down(xs.min), RealInterval.up(xs.max))
  def square: Either[LinAlgError, RealInterval] = if lower <= 0.0 && upper >= 0.0 then RealInterval.checked(0.0, RealInterval.up(math.max(lower * lower, upper * upper))) else multiply(this).flatMap(x => RealInterval.checked(math.max(0.0,x.lower),x.upper))
  def sqrt: Either[LinAlgError, RealInterval] =
    if lower < 0.0 then Left(LinAlgError.InvalidArgument("interval square root requires nonnegative lower bound"))
    else for lo <- RealInterval.sqrtDown(lower); hi <- RealInterval.sqrtUp(upper); out <- RealInterval.checked(lo, hi) yield out
object RealInterval:
  def exact(x: Double): Either[LinAlgError, RealInterval] = checked(x, x)
  def checked(lower: Double, upper: Double): Either[LinAlgError, RealInterval] =
    if lower.isFinite && upper.isFinite && lower <= upper then Right(new RealInterval(lower, upper)) else fail
  private[spectral] def fail = Left(LinAlgError.InvalidArgument("interval endpoint is nonfinite or reversed"))
  private[spectral] def down(x: Double): Double = java.lang.Math.nextDown(x)
  private[spectral] def up(x: Double): Double = java.lang.Math.nextUp(x)
  private def sqrtDown(x: Double): Either[LinAlgError, Double] =
    if x == 0.0 then Right(0.0) else
      def seek(y: Double, n: Int): Either[LinAlgError, Double] =
        if !y.isFinite || y < 0.0 then fail
        else if up(y * y) <= x then checked(down(y), down(y)).map(_.lower)
        else if n == 8 then fail else seek(down(y), n + 1)
      seek(math.sqrt(x), 0)
  private def sqrtUp(x: Double): Either[LinAlgError, Double] =
    if x == 0.0 then Right(0.0) else
      def seek(y: Double, n: Int): Either[LinAlgError, Double] =
        if !y.isFinite || y < 0.0 then fail
        else if down(y * y) >= x then checked(up(y), up(y)).map(_.upper)
        else if n == 8 then fail else seek(up(y), n + 1)
      seek(math.sqrt(x), 0)
