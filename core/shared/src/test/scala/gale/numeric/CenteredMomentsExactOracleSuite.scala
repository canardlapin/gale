package gale.numeric

class CenteredMomentsExactOracleSuite extends munit.FunSuite:
  private def ok[E, A](value: Either[E, A]): A = value.fold(e => fail(e.toString), identity)
  // Independent full fixed-unit accumulation; it has no compact-shift state.
  private def units(value: Double): BigInt =
    val bits = java.lang.Double.doubleToRawLongBits(value)
    val exponent = ((bits >>> 52) & 2047L).toInt
    val fraction = bits & 0xfffffffffffffL
    val magnitude = if exponent == 0 then BigInt(fraction)
      else BigInt(fraction | (1L << 52)) << (exponent - 1)
    if bits < 0L then -magnitude else magnitude
  private def bits(value: Double): Long = java.lang.Double.doubleToRawLongBits(value)
  private def check(rows: Vector[Array[Double]], weights: Vector[Long], degrees: Int): Unit =
    val dimension = rows.head.length
    val mass = weights.sum
    val sums = Vector.tabulate(dimension)(r => rows.indices.foldLeft(BigInt(0))((s, i) => s + weights(i) * units(rows(i)(r))))
    val means = sums.map(s => ExactDyadic.roundRatio(s, BigInt(mass) << 1074))
    val covariance = Vector.tabulate(dimension * dimension): cell =>
      val r = cell / dimension; val c = cell % dimension
      val products = rows.indices.foldLeft(BigInt(0))((s, i) => s + weights(i) * units(rows(i)(r)) * units(rows(i)(c)))
      ExactDyadic.roundRatio(mass * products - sums(r) * sums(c), (BigInt(mass) * (mass - degrees)) << 2148)
    val state = ok(CenteredMoments.make(dimension))
    rows.indices.foreach(i => ok(state.add(rows(i), weights(i))))
    val result = state.finish(degrees)
    if !(means ++ covariance).forall(_.isFinite) then assert(result.isLeft)
    else
      val summary = ok(result)
      assertEquals(summary.mean.toSeq.map(bits).toVector, means.map(bits))
      assertEquals(summary.covariance.valuesRowMajor.map(bits).toVector, covariance.map(bits))
      val left = ok(CenteredMoments.make(dimension)); val right = ok(CenteredMoments.make(dimension))
      rows.indices.reverse.foreach(i => ok((if i % 2 == 0 then left else right).add(rows(i), weights(i))))
      val untouched = right.copy()
      ok(left.merge(right))
      val merged = ok(left.finish(degrees))
      assertEquals(merged.mean.toSeq.map(bits).toVector, means.map(bits))
      assertEquals(merged.covariance.valuesRowMajor.map(bits).toVector, covariance.map(bits))
      assertEquals(right.frequency, untouched.frequency)
      if right.frequency > degrees then
        assertEquals(ok(right.finish(degrees)).mean.toSeq.map(bits).toVector, ok(untouched.finish(degrees)).mean.toSeq.map(bits).toVector)

  test("independent fixed-unit oracle across subnormal, normal, zero and widely separated exponents"):
    val random = new scala.util.Random(417L)
    val exponents = Vector(0, 1, 512, 1000, 1023, 1040, 1533)
    for trial <- 0 until 40 do
      val rows = Vector.tabulate(9): i =>
        Array.tabulate(4): r =>
          if (i + r + trial) % 11 == 0 then (if i % 2 == 0 then 0.0 else -0.0)
          else
            val exponent = exponents((i + r + trial) % exponents.size)
            val sign = if random.nextBoolean() then Long.MinValue else 0L
            java.lang.Double.longBitsToDouble(sign | (exponent.toLong << 52) | (random.nextLong() & 0xfffffffffffffL))
      val weights = Vector.tabulate(9)(i => if i == 3 then 0L else (i + 1L) * (trial + 1L))
      check(rows, weights, 0); check(rows, weights, 1)

  test("full-range constants and frequency weights preserve exact cancellation at final normalization"):
    for value <- Vector(Double.MaxValue, -Double.MaxValue, java.lang.Double.longBitsToDouble(1L), 0.0) do
      check(Vector.fill(3)(Array(value, -value)), Vector(1L << 58, 1L << 58, 1L << 58), 1)
    val large = java.lang.Double.longBitsToDouble((1023L + 500) << 52)
    val tiny = java.lang.Double.longBitsToDouble(1L)
    check(Vector(Array(large, tiny), Array(-large, -tiny), Array(0.0, tiny)), Vector(2L, 2L, 1L), 1)
