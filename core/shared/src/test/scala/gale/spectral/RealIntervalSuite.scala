package gale.spectral
import munit.FunSuite
class RealIntervalSuite extends FunSuite:
  private def rational(x: Double): (BigInt, BigInt) =
    val bits = java.lang.Double.doubleToLongBits(x)
    val exponent = ((bits >>> 52) & 2047L).toInt
    val mantissa = (bits & 0xfffffffffffffL) + (if exponent == 0 then 0L else 1L << 52)
    val n = BigInt(mantissa) * (if bits < 0 then -1 else 1)
    val shift = if exponent == 0 then -1074 else exponent - 1075
    if shift >= 0 then (n << shift, BigInt(1)) else (n, BigInt(1) << -shift)
  private def contains(i: RealInterval, value: (BigInt, BigInt)): Unit =
    val (ln, ld) = rational(i.lower); val (un, ud) = rational(i.upper); val (n, d) = value
    assert(ln * d <= n * ld && n * ud <= un * d, s"exact value not in [${i.lower}, ${i.upper}]")
  private def add(a:(BigInt,BigInt), b:(BigInt,BigInt)) = (a._1*b._2+b._1*a._2,a._2*b._2)
  private def sub(a:(BigInt,BigInt), b:(BigInt,BigInt)) = (a._1*b._2-b._1*a._2,a._2*b._2)
  private def mul(a:(BigInt,BigInt), b:(BigInt,BigInt)) = (a._1*b._1,a._2*b._2)
  private def div(a:(BigInt,BigInt), b:(BigInt,BigInt)) = (a._1*b._2,a._2*b._1)
  test("checked construction and arithmetic enclose exact scalars") {
    val a = RealInterval.exact(2.0).toOption.get; val b = RealInterval.exact(-3.0).toOption.get
    val p = a.multiply(b).toOption.get; assert(p.lower <= -6 && p.upper >= -6)
    assert(RealInterval.checked(Double.NaN, 1).isLeft); assert(a.divide(RealInterval.checked(-1,1).toOption.get).isLeft)
  }
  test("square and sqrt handle zero crossing and subnormal scales") {
    val x = RealInterval.checked(-2, 3).toOption.get.square.toOption.get; assertEquals(x.lower, 0.0); assert(x.upper >= 9)
    // The bounded directed-product verifier deliberately refuses this subnormal
    // court rather than trusting a rounded JS/JVM square.
    assert(RealInterval.exact(java.lang.Double.MIN_VALUE).toOption.get.sqrt.isLeft)
  }
  test("signed subtraction multiplication division and overflow are conservative") {
    val a = RealInterval.checked(-2, 3).toOption.get
    val b = RealInterval.checked(4, 5).toOption.get
    val m = a.multiply(b).toOption.get; assert(m.lower <= -10 && m.upper >= 15)
    val d = a.divide(b).toOption.get; assert(d.lower <= -.5 && d.upper >= .75)
    assert(RealInterval.exact(Double.MaxValue).toOption.get.multiply(RealInterval.exact(2).toOption.get).isLeft)
  }
  test("accepted endpoint courts contain independent exact decimal values") {
    val values = Vector.tabulate(100)(i => ((i - 51).toDouble / 8.0, (i + 3).toDouble / 11.0))
    values.foreach { case (x, y) =>
      val a = RealInterval.exact(x).toOption.get; val b = RealInterval.exact(y).toOption.get
      val rx=rational(x); val ry=rational(y)
      contains(a.add(b).toOption.get, add(rx,ry)); contains(a.subtract(b).toOption.get, sub(rx,ry))
      contains(a.multiply(b).toOption.get, mul(rx,ry)); contains(a.divide(b).toOption.get, div(rx,ry))
    }
  }
