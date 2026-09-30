package gale.linalg

/** The containment oracle is exact dyadic `BigInt` arithmetic, never the
  * production floating-point implementation. */
class OutwardIntervalSuite extends munit.FunSuite:
  private final case class Dyadic(numerator: BigInt, exponent: Int):
    def +(that: Dyadic): Dyadic =
      val exponent0 = math.min(exponent, that.exponent)
      Dyadic((numerator << (exponent - exponent0)) + (that.numerator << (that.exponent - exponent0)), exponent0)
    def unary_- : Dyadic = Dyadic(-numerator, exponent)
    def -(that: Dyadic): Dyadic = this + -that
    def *(that: Dyadic): Dyadic = Dyadic(numerator * that.numerator, exponent + that.exponent)
    def /(that: Dyadic): Rational = Rational(numerator, BigInt(1), exponent) / Rational(that.numerator, BigInt(1), that.exponent)
    def <=(that: Dyadic): Boolean =
      val exponent0 = math.min(exponent, that.exponent)
      (numerator << (exponent - exponent0)) <= (that.numerator << (that.exponent - exponent0))

  private final case class Rational(numerator: BigInt, denominator: BigInt, exponent: Int):
    def /(that: Rational): Rational =
      require(that.numerator != 0)
      val numerator0 = numerator * that.denominator
      val denominator0 = denominator * that.numerator
      if denominator0 < 0 then Rational(-numerator0, -denominator0, exponent - that.exponent)
      else Rational(numerator0, denominator0, exponent - that.exponent)

  private def dyadic(value: Double): Dyadic =
    require(value.isFinite)
    val bits = java.lang.Double.doubleToRawLongBits(value)
    val sign = if (bits < 0L) -BigInt(1) else BigInt(1)
    val exponentBits = ((bits >>> 52) & 0x7ffL).toInt
    val fraction = BigInt(bits & 0x000fffffffffffffL)
    if exponentBits == 0 then Dyadic(sign * fraction, -1074)
    else Dyadic(sign * (fraction + (BigInt(1) << 52)), exponentBits - 1075)

  private def rational(value: Double): Rational =
    val d = dyadic(value)
    Rational(d.numerator, BigInt(1), d.exponent)

  private def lessOrEqual(left: Rational, right: Rational): Boolean =
    val exponent = math.min(left.exponent, right.exponent)
    val leftScaled = (left.numerator * right.denominator) << (left.exponent - exponent)
    val rightScaled = (right.numerator * left.denominator) << (right.exponent - exponent)
    leftScaled <= rightScaled

  private def contains(interval: OutwardInterval, expected: Dyadic): Boolean =
    dyadic(interval.lower) <= expected && expected <= dyadic(interval.upper)

  private def contains(interval: OutwardInterval, expected: Rational): Boolean =
    lessOrEqual(rational(interval.lower), expected) && lessOrEqual(expected, rational(interval.upper))

  private def rectangleProducts(left: OutwardInterval, right: OutwardInterval): Vector[Dyadic] =
    Vector(left.lower, left.upper).flatMap(a => Vector(right.lower, right.upper).map(b => dyadic(a) * dyadic(b)))

  private def rectangleQuotients(left: OutwardInterval, right: OutwardInterval): Vector[Rational] =
    Vector(left.lower, left.upper).flatMap(a => Vector(right.lower, right.upper).map(b => dyadic(a) / dyadic(b)))

  test("construction is finite, ordered, and preserves signed zero") {
    assertEquals(OutwardInterval(Double.NaN, 1.0), Left(OutwardIntervalError.NonFiniteEndpoint))
    assertEquals(OutwardInterval(Double.NegativeInfinity, 1.0), Left(OutwardIntervalError.NonFiniteEndpoint))
    assertEquals(OutwardInterval(1.0, -1.0), Left(OutwardIntervalError.UnorderedEndpoints))
    val zero = OutwardInterval.point(-0.0).toOption.get
    assert(zero.contains(-0.0))
    assert(zero.contains(0.0))
  }

  test("scalar operations enclose exact dyadic values through cancellation and subnormal results") {
    val values = Array(-3.0, -1.5, -java.lang.Double.MIN_VALUE, -0.0, 0.0, java.lang.Double.MIN_VALUE, 1.5, 3.0)
    values.foreach: left =>
      values.foreach: right =>
        val a = OutwardInterval.point(left).toOption.get
        val b = OutwardInterval.point(right).toOption.get
        val sum = (a + b).toOption.get
        assert(contains(sum, dyadic(left) + dyadic(right)), s"sum $left $right")
        val difference = (a - b).toOption.get
        assert(contains(difference, dyadic(left) - dyadic(right)), s"difference $left $right")
        val product = (a * b).toOption.get
        assert(contains(product, dyadic(left) * dyadic(right)), s"product $left $right")
        if right != 0.0 && (left / right).isFinite then
          val quotient = (a / b).toOption.get
          assert(contains(quotient, dyadic(left) / dyadic(right)), s"quotient $left $right")
        else if right != 0.0 then
          assertEquals(a / b, Left(OutwardIntervalError.NonFiniteIntermediate), s"overflowing quotient $left $right")

    val cancellation = (OutwardInterval.point(1.0).toOption.get + OutwardInterval.point(-1.0).toOption.get).toOption.get
    assert(cancellation.lower <= 0.0 && cancellation.upper >= 0.0)
  }

  test("square and streaming enclosures are exact-oracle containing with scalar-only live state by inspection") {
    val crossing = OutwardInterval(-2.0, 3.0).toOption.get.square.toOption.get
    assert(contains(crossing, Dyadic(0, 0)))
    assert(contains(crossing, dyadic(9.0)))

    val left = Array(1.5, -2.0, java.lang.Double.MIN_VALUE, 4.0)
    val right = Array(-2.0, 0.5, 2.0, -0.25)
    val dot = OutwardInterval.streamingDot(left, right).toOption.get
    val exactDot = left.indices.foldLeft(Dyadic(0, 0)) { (sum, index) =>
      sum + dyadic(left(index)) * dyadic(right(index))
    }
    assert(contains(dot, exactDot))

    val absolute = OutwardInterval.streamingAbsoluteSum(Array(-3.0, -0.0, java.lang.Double.MIN_VALUE, 2.5)).toOption.get
    assert(contains(absolute, dyadic(5.5) + dyadic(java.lang.Double.MIN_VALUE)))
  }

  test("fixed-bit randomized finite patterns retain exact containment") {
    var state = 0x9e3779b97f4a7c15L
    def nextFinite(): Double =
      state = state * 6364136223846793005L + 1442695040888963407L
      val sign = state & Long.MinValue
      val exponent = ((state >>> 53) % 80L + 970L) << 52
      val fraction = state & 0x000fffffffffffffL
      java.lang.Double.longBitsToDouble(sign | exponent | fraction)

    var index = 0
    while index < 256 do
      val left = nextFinite()
      val right = nextFinite()
      val a = OutwardInterval.point(left).toOption.get
      val b = OutwardInterval.point(right).toOption.get
      assert(contains((a + b).toOption.get, dyadic(left) + dyadic(right)))
      assert(contains((a * b).toOption.get, dyadic(left) * dyadic(right)))
      if right != 0.0 then assert(contains((a / b).toOption.get, dyadic(left) / dyadic(right)))
      index += 1
  }

  test("nondegenerate multiplication rectangles enclose every exact corner in all sign quadrants") {
    val rectangles = Vector(
      (-4.5, -1.25, -3.5, -0.5),
      (-4.5, -1.25, 0.5, 3.5),
      (1.25, 4.5, -3.5, -0.5),
      (1.25, 4.5, 0.5, 3.5))
    rectangles.foreach: (leftLow, leftHigh, rightLow, rightHigh) =>
      val left = OutwardInterval(leftLow, leftHigh).toOption.get
      val right = OutwardInterval(rightLow, rightHigh).toOption.get
      val product = (left * right).toOption.get
      rectangleProducts(left, right).foreach: exact =>
        assert(contains(product, exact), s"rectangle [$leftLow,$leftHigh] x [$rightLow,$rightHigh] must contain $exact")
  }

  test("nondegenerate positive and negative divisor rectangles enclose every exact rational corner") {
    val cases = Vector(
      (OutwardInterval(-7.5, -2.5).toOption.get, OutwardInterval(1.5, 4.5).toOption.get),
      (OutwardInterval(2.5, 7.5).toOption.get, OutwardInterval(-4.5, -1.5).toOption.get))
    cases.foreach: (left, divisor) =>
      val quotient = (left / divisor).toOption.get
      rectangleQuotients(left, divisor).foreach: exact =>
        assert(contains(quotient, exact), s"quotient [$left] / [$divisor] must contain $exact")
  }

  test("one-sided squares, adjacent maximum expansion, empty streams, late nonfinite values, and signed zero are controlled") {
    Vector(OutwardInterval(2.0, 3.0).toOption.get, OutwardInterval(-3.0, -2.0).toOption.get).foreach: interval =>
      val squared = interval.square.toOption.get
      assert(contains(squared, dyadic(4.0)))
      assert(contains(squared, dyadic(9.0)))
    assertEquals(OutwardInterval.fromRounded(java.lang.Math.nextDown(Double.MaxValue), Double.MaxValue),
      Left(OutwardIntervalError.NonFiniteEndpoint))
    val emptyDot = OutwardInterval.streamingDot(Array.emptyDoubleArray, Array.emptyDoubleArray).toOption.get
    val emptyAbsolute = OutwardInterval.streamingAbsoluteSum(Array.emptyDoubleArray).toOption.get
    assert(contains(emptyDot, dyadic(0.0)))
    assert(contains(emptyAbsolute, dyadic(0.0)))
    assertEquals(OutwardInterval.streamingDot(Array(1.0, Double.NaN), Array(1.0, 1.0)), Left(OutwardIntervalError.NonFiniteEndpoint))
    assertEquals(OutwardInterval.streamingAbsoluteSum(Array(1.0, Double.PositiveInfinity)), Left(OutwardIntervalError.NonFiniteEndpoint))
    val signed = OutwardInterval(-0.0, 0.0).toOption.get
    assertEquals(java.lang.Double.doubleToRawLongBits(signed.lower), java.lang.Double.doubleToRawLongBits(-0.0))
    assertEquals(java.lang.Double.doubleToRawLongBits(signed.upper), java.lang.Double.doubleToRawLongBits(0.0))
  }

  test("overflow, nonfinite input, zero divisors, and incompatible streams refuse") {
    val max = OutwardInterval.point(Double.MaxValue).toOption.get
    assertEquals(max + max, Left(OutwardIntervalError.NonFiniteIntermediate))
    assertEquals(max * max, Left(OutwardIntervalError.NonFiniteIntermediate))
    val crossingZero = OutwardInterval(-1.0, 1.0).toOption.get
    assertEquals(OutwardInterval.point(1.0).toOption.get / crossingZero, Left(OutwardIntervalError.ZeroContainingDivisor))
    assertEquals(OutwardInterval.streamingDot(Array(1.0), Array(1.0, 2.0)), Left(OutwardIntervalError.LengthMismatch(1, 2)))
    assertEquals(OutwardInterval.streamingAbsoluteSum(Array(Double.PositiveInfinity)), Left(OutwardIntervalError.NonFiniteEndpoint))
  }
