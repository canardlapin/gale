package gale.linalg

/** Independent exact rational controls for the streaming certificate. */
class FiniteSymmetricSolveBoundsSuite extends munit.FunSuite:
  private final case class Q(n: BigInt, d: BigInt = 1):
    require(d != 0)
    private val sign = if d < 0 then -1 else 1
    private val divisor = n.gcd(d).abs
    val numerator: BigInt = sign * n / divisor
    val denominator: BigInt = sign * d / divisor
    def +(that: Q): Q = Q(numerator * that.denominator + that.numerator * denominator, denominator * that.denominator)
    def -(that: Q): Q = Q(numerator * that.denominator - that.numerator * denominator, denominator * that.denominator)
    def *(that: Q): Q = Q(numerator * that.numerator, denominator * that.denominator)
    def /(that: Q): Q = Q(numerator * that.denominator, denominator * that.numerator)
    def square: Q = this * this
    def <=(that: Q): Boolean = numerator * that.denominator <= that.numerator * denominator

  private def q(value: Double): Q =
    val bits = java.lang.Double.doubleToRawLongBits(value)
    val sign = if bits < 0L then -BigInt(1) else BigInt(1)
    val exponent = ((bits >>> 52) & 0x7ffL).toInt
    val fraction = BigInt(bits & 0x000fffffffffffffL)
    if exponent == 0 then Q(sign * fraction, BigInt(1) << 1074)
    else
      val significand = sign * (fraction + (BigInt(1) << 52))
      if exponent >= 1075 then Q(significand << (exponent - 1075))
      else Q(significand, BigInt(1) << (1075 - exponent))

  private def interval(lower: Double, upper: Double = Double.NaN): OutwardInterval =
    val high = if upper.isNaN then lower else upper
    OutwardInterval(lower, high).toOption.get

  private def entries(values: (Double, Double, Double)): Iterator[UpperIntervalEntry] =
    Iterator(
      UpperIntervalEntry(0, 0, interval(values._1)),
      UpperIntervalEntry(0, 1, interval(values._2)),
      UpperIntervalEntry(1, 1, interval(values._3))
    )

  private def solve2(a: Q, b: Q, c: Q, d: Q, y0: Q, y1: Q): (Q, Q) =
    val determinant = a * d - b * c
    ((d * y0 - b * y1) / determinant, (a * y1 - c * y0) / determinant)

  test("a known exact solve has an independently checked rational error control") {
    val bound = FiniteSymmetricSolveBounds
      .certify(2, entries((2.0, 1.0, 2.0)), Vector(interval(1.0), interval(0.0)), Vector(0.5, 0.0))
      .toOption
      .get
    val exact = solve2(Q(2), Q(1), Q(1), Q(2), Q(1), Q(0))
    val errorSquared = (exact._1 - Q(1, 2)).square + exact._2.square
    // Case-class equality compares the raw constructor fields (45/324 here).
    // Compare canonical exact values so equivalent fractions remain equal.
    assertEquals((errorSquared.numerator, errorSquared.denominator), (BigInt(5), BigInt(36)))
    assert(bound.gershgorinMarginLower > 0.0 && bound.gershgorinMarginLower <= 1.0)
    assert(bound.residualL1Upper >= 0.5)
    assert(bound.solutionError2Upper < 0.51)
    assert(errorSquared <= q(bound.solutionError2Upper).square)
    assertEquals(bound.upperTriangleEntryCount, 3L)
  }

  test("nondegenerate dyadic matrix and rhs corners are bounded by independent exact solves") {
    val stream = Iterator(
      UpperIntervalEntry(0, 0, interval(3.75, 4.25)),
      UpperIntervalEntry(0, 1, interval(-0.25, 0.25)),
      UpperIntervalEntry(1, 1, interval(3.75, 4.25))
    )
    val bound = FiniteSymmetricSolveBounds
      .certify(2, stream, Vector(interval(-1.0, 1.0), interval(-0.5, 0.5)), Vector(0.0, 0.0))
      .toOption
      .get
    val matrixCorners =
      Vector(3.75, 4.25).flatMap(a => Vector(-0.25, 0.25).flatMap(b => Vector(3.75, 4.25).map(d => (a, b, d))))
    val rhsCorners = Vector(-1.0, 1.0).flatMap(y0 => Vector(-0.5, 0.5).map(y1 => (y0, y1)))
    matrixCorners.foreach: (a, b, d) =>
      rhsCorners.foreach: (y0, y1) =>
        val solved = solve2(q(a), q(b), q(b), q(d), q(y0), q(y1))
        val squaredNorm = solved._1.square + solved._2.square
        assert(squaredNorm <= q(bound.solutionError2Upper).square, s"corner ($a, $b, $d; $y0, $y1)")
  }

  test("non-dominant, singular, indefinite, and interval-singular enclosures refuse honestly") {
    assert(
      FiniteSymmetricSolveBounds
        .certify(2, entries((1.0, 2.0, 5.0)), Vector(interval(0.0), interval(0.0)), Vector(0.0, 0.0))
        .isLeft
    )
    assert(
      FiniteSymmetricSolveBounds
        .certify(2, entries((1.0, 1.0, 1.0)), Vector(interval(0.0), interval(0.0)), Vector(0.0, 0.0))
        .isLeft
    )
    assert(
      FiniteSymmetricSolveBounds
        .certify(2, entries((-1.0, 0.0, -1.0)), Vector(interval(0.0), interval(0.0)), Vector(0.0, 0.0))
        .isLeft
    )
    assert(
      FiniteSymmetricSolveBounds
        .certify(1, Iterator(UpperIntervalEntry(0, 0, interval(0.0, 1.0))), Vector(interval(0.0)), Vector(0.0))
        .isLeft
    )
  }

  test("canonical completeness, shapes, and finite candidate validation fail before a false certificate") {
    val rhs = Vector(interval(0.0), interval(0.0))
    assertEquals(
      FiniteSymmetricSolveBounds.certify(0, Iterator.empty, Vector.empty, Vector.empty),
      Left(FiniteSolveBoundError.NonPositiveSize(0))
    )
    assertEquals(
      FiniteSymmetricSolveBounds.certify(2, Iterator.empty, Vector(interval(0.0)), Vector(0.0, 0.0)),
      Left(FiniteSolveBoundError.RhsLengthMismatch(2, 1))
    )
    assertEquals(
      FiniteSymmetricSolveBounds.certify(2, Iterator.empty, rhs, Vector(0.0)),
      Left(FiniteSolveBoundError.CandidateLengthMismatch(2, 1))
    )
    assertEquals(
      FiniteSymmetricSolveBounds.certify(2, Iterator.empty, rhs, Vector(Double.NaN, 0.0)),
      Left(FiniteSolveBoundError.NonFiniteCandidate(0))
    )
    assert(
      FiniteSymmetricSolveBounds
        .certify(2, Iterator(UpperIntervalEntry(0, 0, interval(2.0))), rhs, Vector(0.0, 0.0))
        .isLeft
    )
    assert(
      FiniteSymmetricSolveBounds
        .certify(2, Iterator(UpperIntervalEntry(0, 1, interval(0.0))), rhs, Vector(0.0, 0.0))
        .isLeft
    )
    assert(
      FiniteSymmetricSolveBounds
        .certify(
          1,
          Iterator(UpperIntervalEntry(0, 0, interval(1.0)), UpperIntervalEntry(0, 0, interval(0.0))),
          Vector(interval(0.0)),
          Vector(0.0)
        )
        .isLeft
    )
  }

  test("subnormal, signed cancellation, and nonfinite intermediate controls are fail closed") {
    val tiny = java.lang.Math.scalb(1.0, -30)
    val tinyBound = FiniteSymmetricSolveBounds
      .certify(2, entries((1.0, 0.0, tiny)), Vector(interval(0.0), interval(0.0)), Vector(0.0, 0.0))
      .toOption
      .get
    assert(tinyBound.gershgorinMarginLower > 0.0)
    assert(tinyBound.gershgorinMarginLower <= tiny)

    // Outward subtraction loses a strictly positive margin at the smallest
    // subnormal. Refusal is honest; it must not produce a zero-margin bound.
    assert(
      FiniteSymmetricSolveBounds
        .certify(
          1,
          Iterator(UpperIntervalEntry(0, 0, interval(java.lang.Double.MIN_VALUE))),
          Vector(interval(0.0)),
          Vector(0.0)
        )
        .isLeft
    )

    val cancellation = FiniteSymmetricSolveBounds
      .certify(2, entries((2.0, -1.0, 2.0)), Vector(interval(1.0), interval(1.0)), Vector(1.0, 1.0))
      .toOption
      .get
    assert(cancellation.residualL1Upper >= 0.0)
    assert(
      FiniteSymmetricSolveBounds
        .certify(1, Iterator(UpperIntervalEntry(0, 0, interval(Double.MaxValue))), Vector(interval(0.0)), Vector(2.0))
        .isLeft
    )
  }
