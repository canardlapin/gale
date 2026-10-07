package gale.numeric

class ExactSumSuite extends munit.FunSuite:

  // ExactSum promises bit-exact, correctly rounded totals, so these assertions compare doubles exactly.

  test("cancellation is exact where naive summation loses the small terms") {
    assertEquals(sumOf(1e16, 1.0, -1e16), 1.0)
    assertEquals(sumOf(1e308, 1e308, -1e308, -1e308, 3.5), 3.5)
    assertEquals(sumOf(Double.MinPositiveValue, -Double.MinPositiveValue), 0.0)
    assertEquals(sumOf(), 0.0)
  }

  test("totals round once to nearest, ties to even") {
    val ulpAtOne = math.ulp(1.0)
    assertEquals(sumOf(1.0, ulpAtOne / 2.0), 1.0)
    assertEquals(sumOf(1.0, ulpAtOne / 2.0, Double.MinPositiveValue), 1.0 + ulpAtOne)
    assertEquals(sumOf(1.0 + ulpAtOne, ulpAtOne / 2.0), 1.0 + 2.0 * ulpAtOne)
    assertEquals(sumOf(-1.0, -ulpAtOne / 2.0, -Double.MinPositiveValue), -(1.0 + ulpAtOne))
    assertEquals(sumOf(Double.MinPositiveValue, Double.MinPositiveValue), 2.0 * Double.MinPositiveValue)
    assertEquals(sumOf(Double.MaxValue, Double.MaxValue), Double.PositiveInfinity)
  }

  test("totals are independent of order and of how partial sums are merged") {
    val values = Vector.tabulate(257) { index =>
      val raw = math.sin(index.toDouble * 12.9898) * 43758.5453
      (raw - math.floor(raw) - 0.5) * math.pow(10.0, (index % 9 - 4).toDouble)
    }
    val reference = sumOf(values*)
    assertEquals(sumOf(values.reverse*), reference)
    (1 to 17).foreach { width =>
      val parts = values
        .grouped(width)
        .map { group =>
          val part = ExactSum.zero()
          group.foreach(value => accepted(part.add(value)))
          part
        }
        .toVector
      val merged = ExactSum.zero()
      parts.reverse.foreach(part => accepted(merged.addAll(part)))
      assertEquals(merged.value, reference, s"width=$width")
    }
  }

  test("non-finite inputs follow IEEE addition") {
    assertEquals(sumOf(1.0, Double.PositiveInfinity), Double.PositiveInfinity)
    assert(sumOf(Double.PositiveInfinity, Double.NegativeInfinity).isNaN)
    assert(sumOf(1.0, Double.NaN).isNaN)
  }

  test("boundary totals: subnormals, large ties, window edges, overflow threshold, and signed zero") {
    val minNormal = java.lang.Double.MIN_NORMAL
    val tiny = Double.MinPositiveValue
    // Subnormal totals are exact; the first rounding happens at 2^-1021, where one ulp is two units.
    assertEquals(sumOf(minNormal, -tiny), minNormal - tiny)
    assertEquals(sumOf(tiny, tiny, tiny), 3.0 * tiny)
    assertEquals(sumOf(2.0 * minNormal, tiny), 2.0 * minNormal)
    assertEquals(sumOf(2.0 * minNormal, 3.0 * tiny), 2.0 * minNormal + 4.0 * tiny)
    // Ties at large magnitudes.
    assertEquals(sumOf(math.pow(2.0, 1000), math.pow(2.0, 947)), math.pow(2.0, 1000))
    assertEquals(
      sumOf(math.pow(2.0, 1000) + math.pow(2.0, 948), math.pow(2.0, 947)),
      math.pow(2.0, 1000) + math.pow(2.0, 949)
    )
    // Overflow happens exactly at (2^1024 - 2^970), the midpoint above MaxValue.
    val halfTop = math.ulp(Double.MaxValue) / 2.0
    assertEquals(sumOf(Double.MaxValue, halfTop), Double.PositiveInfinity)
    assertEquals(sumOf(-Double.MaxValue, -halfTop), Double.NegativeInfinity)
    assertEquals(sumOf(Double.MaxValue, halfTop, -tiny), Double.MaxValue)
    assertEquals(sumOf(-Double.MaxValue, -halfTop, tiny), -Double.MaxValue)
    // Exact zero is always +0.0, including a sum of negative zeros.
    Vector(Vector(-0.0), Vector(-0.0, -0.0), Vector(5.0, -5.0), Vector(-tiny, tiny)).foreach { values =>
      assertEquals(java.lang.Double.doubleToRawLongBits(sumOf(values*)), 0L, values.toString)
    }
    // The integer fast path ends at 2^62 units; probe ties and near-ties on both sides of each window edge.
    Vector(53, 61, 62, 63, 64).foreach { exponent =>
      val edge = powerOfTwo(exponent - 1074)
      val ulp = math.ulp(edge)
      Vector(
        Vector(edge, ulp / 2.0),
        Vector(edge, ulp / 2.0, tiny),
        Vector(edge, ulp / 2.0, -tiny),
        Vector(edge, ulp, ulp / 2.0),
        Vector(edge, -ulp / 4.0),
        Vector(edge, -ulp / 4.0, -tiny),
        Vector(edge, -edge / 2.0, ulp / 8.0),
        Vector(-edge, -ulp / 2.0, -tiny)
      ).foreach(values => assertMatchesReference(values, s"edge 2^$exponent units"))
    }
  }

  test("randomized totals match an independent BigInteger reference rounded to nearest, ties to even") {
    val random = new scala.util.Random(0x5ca1af1L)
    val tiny = Double.MinPositiveValue
    def mantissa(): Long = (random.nextLong() >>> 11) | (1L << 52)
    def clustered(base: Int): Double =
      val exponent = math.max(-1074, math.min(971, base - random.nextInt(70)))
      val magnitude = mantissa().toDouble * powerOfTwo(exponent)
      if random.nextBoolean() then magnitude else -magnitude
    def anyFinite(): Double =
      var candidate = java.lang.Double.longBitsToDouble(random.nextLong())
      while !candidate.isFinite do candidate = java.lang.Double.longBitsToDouble(random.nextLong())
      candidate
    (0 until 600).foreach { trial =>
      val values: Vector[Double] = trial % 6 match
        case 0 => Vector.fill(1 + random.nextInt(12))(anyFinite())
        case 1 =>
          val base = random.nextInt(2040) - 1070
          Vector.fill(2 + random.nextInt(20))(clustered(base))
        case 2 =>
          val x = clustered(random.nextInt(1900) - 1000)
          val half = math.ulp(x) / 2.0
          Vector(x, if random.nextBoolean() then half else -half) ++
            Vector.fill(random.nextInt(3))(if random.nextBoolean() then tiny else -tiny)
        case 3 => Vector.fill(1 + random.nextInt(30))(clustered(-1000 - random.nextInt(74)))
        case 4 =>
          val top = Double.MaxValue - math.ulp(Double.MaxValue) * random.nextInt(4).toDouble
          Vector(top, math.ulp(Double.MaxValue) / 2.0, clustered(970 - random.nextInt(60)))
        case _ =>
          val x = clustered(random.nextInt(1900) - 1000)
          Vector(x, -x, clustered(random.nextInt(1900) - 1000), clustered(random.nextInt(1900) - 1000))
      assertMatchesReference(values, s"trial $trial")
    }
  }

  test("copy storage and merge sources remain independent after subsequent mutation") {
    val original = accumulator(Vector(1e16, 1.0, -1e16))
    val copied = original.copy()
    val receiver = accumulator(Vector(2.0))
    accepted(receiver.addAll(original))
    assertEquals(receiver.value, 3.0)
    assertEquals(receiver.finiteTerms, 4L)
    assertEquals(original.value, 1.0)
    assertEquals(original.finiteTerms, 3L)
    accepted(original.add(4.0))
    accepted(copied.add(8.0))
    assertEquals(original.value, 5.0)
    assertEquals(copied.value, 9.0)
    assertEquals(receiver.value, 3.0)
    assertEquals(original.finiteTerms, 4L)
    assertEquals(copied.finiteTerms, 4L)
  }

  test("reading a rounded value preserves exact state for continued cancellation and merging") {
    val source = accumulator(Vector(1e16, 1.0))
    assertEquals(source.value, 1e16)
    assertEquals(source.value, 1e16)
    val copy = source.copy()
    accepted(source.add(-1e16))
    assertEquals(source.value, 1.0)
    val receiver = accumulator(Vector(-1e16))
    accepted(receiver.addAll(copy))
    assertEquals(receiver.value, 1.0)
    assertEquals(copy.value, 1e16)
  }

  test("self-merge doubles exact state, finite counts, and IEEE special channels") {
    val sum = accumulator(Vector(1e16, 1.0, -1e16))
    accepted(sum.addAll(sum))
    assertEquals(sum.value, 2.0)
    assertEquals(sum.finiteTerms, 6L)
    val special = accumulator(Vector(2.0, Double.PositiveInfinity))
    accepted(special.addAll(special))
    assertEquals(special.value, Double.PositiveInfinity)
    assertEquals(special.finiteTerms, 2L)
  }

  test("mixed-exponent partial sums and balanced merge trees retain independent oracle totals") {
    val random = new scala.util.Random(0x6a1eL)
    (0 until 80).foreach { trial =>
      val values = Vector.fill(31) {
        var candidate = java.lang.Double.longBitsToDouble(random.nextLong())
        while !candidate.isFinite do candidate = java.lang.Double.longBitsToDouble(random.nextLong())
        candidate
      }
      val expected = sumOf(values*)
      assertMatchesReference(values, s"merge trial $trial")
      var layer = random.shuffle(values).grouped(3).map(accumulator).toVector
      while layer.length > 1 do
        layer = random
          .shuffle(layer)
          .grouped(2)
          .map { group =>
            val merged = group.head.copy()
            group.drop(1).foreach(part => accepted(merged.addAll(part)))
            merged
          }
          .toVector
      assertEquals(layer.head.value, expected, s"merge trial $trial")
      assertEquals(layer.head.finiteTerms, values.count(_ != 0.0).toLong)
    }
  }

  test("finite term capacity is admitted exactly and refused additions and merges leave state unchanged") {
    val limit = accumulator(Vector(1.0))
    var doublings = 0
    while doublings < 60 do
      accepted(limit.addAll(limit))
      doublings += 1
    assertEquals(limit.finiteTerms, ExactSum.MaxTerms)
    assertEquals(limit.value, powerOfTwo(60))
    val before = limit.copy()
    assertEquals(limit.add(-1.0), Left(ExactSumError.CapacityExceeded(ExactSum.MaxTerms, 1L)))
    assertEquals(limit.addAll(limit), Left(ExactSumError.CapacityExceeded(ExactSum.MaxTerms, ExactSum.MaxTerms)))
    assertEquals(limit.value, before.value)
    assertEquals(limit.finiteTerms, before.finiteTerms)
    val receiver = accumulator(Vector(-1.0))
    assertEquals(receiver.addAll(limit), Left(ExactSumError.CapacityExceeded(1L, ExactSum.MaxTerms)))
    assertEquals(receiver.value, -1.0)
    assertEquals(receiver.finiteTerms, 1L)
    assertEquals(limit.value, before.value)
    val previous = accumulator(Vector(1.0))
    var bit = 1
    while bit < 60 do
      accepted(previous.addAll(previous))
      bit += 1
    accepted(previous.addAll(previous.copy()))
    assertEquals(previous.finiteTerms, ExactSum.MaxTerms)
  }

  test("zeros and nonfinite channels do not consume finite capacity") {
    val sum = accumulator(Vector(-0.0, 0.0, Double.PositiveInfinity, Double.NaN))
    assertEquals(sum.finiteTerms, 0L)
    assert(sum.value.isNaN)
    val receiver = accumulator(Vector(1.0))
    accepted(receiver.addAll(sum))
    assertEquals(receiver.finiteTerms, 1L)
    assert(receiver.value.isNaN)
    val full = accumulator(Vector(1.0))
    var bit = 0
    while bit < 60 do
      accepted(full.addAll(full))
      bit += 1
    accepted(full.add(-0.0))
    accepted(full.add(Double.PositiveInfinity))
    accepted(full.addAll(sum))
    assertEquals(full.finiteTerms, ExactSum.MaxTerms)
    assert(full.value.isNaN)
  }

  test("finite cancellation retains term counts and finite overflow can cancel back to a finite total") {
    val cancellation = accumulator(Vector(1.0, -1.0))
    assertEquals(cancellation.value, 0.0)
    assertEquals(cancellation.finiteTerms, 2L)
    val overflow = accumulator(Vector(Double.MaxValue, Double.MaxValue))
    assertEquals(overflow.value, Double.PositiveInfinity)
    accepted(overflow.add(-Double.MaxValue))
    assertEquals(overflow.value, Double.MaxValue)
    assertEquals(overflow.finiteTerms, 3L)
  }

  private def accepted(result: Either[ExactSumError, Unit]): Unit =
    assertEquals(result, Right(()))

  test("ratios normalize overflowing totals and weights before rounding") {
    val max = Double.MaxValue
    assertEquals(accumulator(Vector(max, max)).ratio(accumulator(Vector(1.0, 1.0))), Right(max))
    assertEquals(accumulator(Vector(-max, -max)).ratio(accumulator(Vector(1.0, 1.0))), Right(-max))
    assertEquals(accumulator(Vector(max, max)).ratio(accumulator(Vector(max, max))), Right(1.0))
    assertEquals(accumulator(Vector(max, max, -max, -max, 3.0)).ratio(accumulator(Vector(2.0))), Right(1.5))
    // Unequal positive weights: 2 * 3 + 8 * 1, divided by 3 + 1.
    assertEquals(accumulator(Vector(6.0, 8.0)).ratio(accumulator(Vector(3.0, 1.0))), Right(3.5))
    assertEquals(accumulator(Vector(1.0)).ratio(accumulator(Vector(3.0))), Right(1.0 / 3.0))
    assertEquals(accumulator(Vector(1.0)).ratio(accumulator(Vector(-2.0))), Right(-0.5))
    val full = accumulator(Vector(max))
    (0 until 60).foreach(_ => accepted(full.addAll(full)))
    assertEquals(full.ratio(full), Right(1.0))
    assertEquals(full.finiteTerms, ExactSum.MaxTerms)
  }

  test("ratios round ties, subnormals, signed underflow and overflow exactly") {
    val tiny = Double.MinPositiveValue
    val cases = Vector(
      (Vector(tiny), Vector(2.0), 0.0),
      (Vector(3.0 * tiny), Vector(2.0), 2.0 * tiny),
      (Vector(5.0 * tiny), Vector(2.0), 2.0 * tiny),
      (Vector(java.lang.Double.MIN_NORMAL, -tiny), Vector(1.0), java.lang.Double.MIN_NORMAL - tiny),
      (Vector(2.0, math.ulp(1.0)), Vector(2.0), 1.0),
      (Vector(2.0, math.ulp(1.0), tiny), Vector(2.0), 1.0 + math.ulp(1.0)),
      (Vector(2.0, 3.0 * math.ulp(1.0)), Vector(2.0), 1.0 + 2.0 * math.ulp(1.0)),
      (Vector(Double.MaxValue, math.ulp(Double.MaxValue) / 2.0), Vector(1.0), Double.PositiveInfinity),
      (Vector(Double.MaxValue, math.ulp(Double.MaxValue) / 2.0, -tiny), Vector(1.0), Double.MaxValue),
      (Vector(Double.MaxValue), Vector(tiny), Double.PositiveInfinity)
    )
    cases.foreach { (n, d, expected) =>
      assertEquals(accumulator(n).ratio(accumulator(d)), Right(expected))
      val negative = accumulator(n.map(-_)).ratio(accumulator(d)).toOption.get
      assertEquals(java.lang.Double.doubleToRawLongBits(negative), java.lang.Double.doubleToRawLongBits(-expected))
    }
    assertEquals(accumulator(Vector(1.0, -1.0)).ratio(accumulator(Vector(-1.0))), Right(0.0))
  }

  test("ratio refusals and successful reads preserve both builders") {
    val n = accumulator(Vector(Double.MaxValue, Double.MaxValue))
    val d = accumulator(Vector(1.0, 1.0))
    assertEquals(n.ratio(d), Right(Double.MaxValue))
    assertEquals(n.ratio(d), Right(Double.MaxValue))
    assertEquals(n.value, Double.PositiveInfinity)
    assertEquals(d.value, 2.0)
    assertEquals(n.finiteTerms, 2L)
    assertEquals(d.finiteTerms, 2L)
    accepted(n.add(-Double.MaxValue))
    accepted(d.add(2.0))
    assertEquals(n.ratio(d), Right(Double.MaxValue / 4.0))
    Vector(Vector.empty[Double], Vector(-0.0), Vector(1.0, -1.0)).foreach { values =>
      assertEquals(n.ratio(accumulator(values)), Left(ExactSumError.ZeroNormalizer))
    }
    Vector(Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity).foreach { value =>
      val special = accumulator(Vector(value))
      assertEquals(n.ratio(special), Left(ExactSumError.NonFiniteRatio))
      assertEquals(special.ratio(d), Left(ExactSumError.NonFiniteRatio))
      assertEquals(special.ratio(ExactSum.zero()), Left(ExactSumError.NonFiniteRatio))
      assertEquals(special.finiteTerms, 0L)
    }
  }

  test("ratios match independent rational distances across permutations and merge trees") {
    val random = new scala.util.Random(0x19ea1L)
    def finite(): Double =
      var x = java.lang.Double.longBitsToDouble(random.nextLong())
      while !x.isFinite do x = java.lang.Double.longBitsToDouble(random.nextLong())
      x
    def merged(values: Vector[Double]): ExactSum =
      var layer = random.shuffle(values).grouped(2).map(accumulator).toVector
      while layer.length > 1 do
        layer = random
          .shuffle(layer)
          .grouped(2)
          .map { group =>
            val total = group.head.copy()
            group.drop(1).foreach(part => accepted(total.addAll(part)))
            total
          }
          .toVector
      layer.head
    (0 until 400).foreach { trial =>
      val ns = Vector.fill(1 + random.nextInt(8))(finite())
      val ds = Vector.fill(1 + random.nextInt(8))(math.abs(finite())) :+ Double.MinPositiveValue
      val n = ns.foldLeft(java.math.BigInteger.ZERO)((s, x) => s.add(units(x)))
      val d = ds.foldLeft(java.math.BigInteger.ZERO)((s, x) => s.add(units(x)))
      val actual = accumulator(ns).ratio(accumulator(ds)).toOption.get
      assertEquals(merged(ns).ratio(merged(ds)), Right(actual), s"merge trial $trial")
      assertEquals(accumulator(ns.reverse).ratio(accumulator(ds.reverse)), Right(actual), s"order trial $trial")
      // Compare exact rational distances to adjacent doubles; this oracle does no binary quotient rounding.
      val scaled = n.shiftLeft(1074)
      if scaled.abs().compareTo(d.multiply(OverflowThreshold)) >= 0 then
        assertEquals(actual, if n.signum() > 0 then Double.PositiveInfinity else Double.NegativeInfinity)
      else
        assert(actual.isFinite, s"trial $trial")
        val distance = scaled.subtract(d.multiply(units(actual))).abs()
        Vector(nextUp(actual), nextDown(actual)).filter(_.isFinite).foreach { neighbour =>
          val other = scaled.subtract(d.multiply(units(neighbour))).abs()
          val order = distance.compareTo(other)
          assert(order <= 0, s"trial $trial: $neighbour is nearer than $actual")
          if order == 0 then assertEquals(java.lang.Double.doubleToRawLongBits(actual) & 1L, 0L)
        }
    }
  }

  private def accumulator(values: Seq[Double]): ExactSum =
    val sum = ExactSum.zero()
    values.foreach(value => accepted(sum.add(value)))
    sum

  private def sumOf(values: Double*): Double =
    val sum = ExactSum.zero()
    values.foreach(value => accepted(sum.add(value)))
    sum.value

  // Independent reference: the exact total in units of 2^-1074 as a BigInteger, built from the IEEE bit fields.
  // java.math.BigDecimal is deliberately avoided; its double conversions are not exact on Scala.js.
  private val OverflowThreshold =
    java.math.BigInteger.ONE.shiftLeft(2098).subtract(java.math.BigInteger.ONE.shiftLeft(2044))

  private def units(value: Double): java.math.BigInteger =
    val bits = java.lang.Double.doubleToRawLongBits(value)
    val exponent = ((bits >>> 52) & 0x7ffL).toInt
    val fraction = bits & 0xfffffffffffffL
    val magnitude =
      if exponent == 0 then java.math.BigInteger.valueOf(fraction)
      else java.math.BigInteger.valueOf(fraction | (1L << 52)).shiftLeft(exponent - 1)
    if bits < 0L then magnitude.negate() else magnitude

  /** Exact 2^exponent for exponent in [-1074, 1023], built from bits so it is portable to Scala.js. */
  private def powerOfTwo(exponent: Int): Double =
    if exponent >= -1022 then java.lang.Double.longBitsToDouble((exponent + 1023).toLong << 52)
    else java.lang.Double.longBitsToDouble(1L << (exponent + 1074))

  private def nextUp(value: Double): Double =
    if value == 0.0 then Double.MinPositiveValue
    else
      val bits = java.lang.Double.doubleToRawLongBits(value)
      java.lang.Double.longBitsToDouble(if value > 0.0 then bits + 1L else bits - 1L)

  private def nextDown(value: Double): Double = -nextUp(-value)

  private def assertMatchesReference(values: Vector[Double], clue: String): Unit =
    val exact = values.foldLeft(java.math.BigInteger.ZERO)((total, value) => total.add(units(value)))
    val actual = sumOf(values*)
    val context = s"$clue values=$values actual=$actual"
    if exact.abs().compareTo(OverflowThreshold) >= 0 then
      assertEquals(actual, if exact.signum() > 0 then Double.PositiveInfinity else Double.NegativeInfinity, context)
    else if exact.signum() == 0 then assertEquals(java.lang.Double.doubleToRawLongBits(actual), 0L, context)
    else
      assert(actual.isFinite, context)
      assertEquals(math.signum(actual), exact.signum().toDouble, context)
      val distance = exact.subtract(units(actual)).abs()
      Vector(nextUp(actual), nextDown(actual)).filter(_.isFinite).foreach { neighbour =>
        val other = exact.subtract(units(neighbour)).abs()
        val order = distance.compareTo(other)
        assert(order <= 0, s"$context: $neighbour is nearer")
        if order == 0 then
          assertEquals(java.lang.Double.doubleToRawLongBits(actual) & 1L, 0L, s"$context: tie not rounded to even")
      }
