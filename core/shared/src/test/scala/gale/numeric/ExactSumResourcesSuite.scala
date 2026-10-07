package gale.numeric

class ExactSumResourcesSuite extends munit.FunSuite:
  private def ok[E, A](value: Either[E, A]): A = value.fold(e => fail(e.toString), identity)

  test("resource inspection shares immutable facts without rounding or normalization") {
    val sum = ExactSum.zero()
    ok(sum.add(1e16)); ok(sum.add(1.0)); ok(sum.add(-1e16))
    val before = sum.finiteTerms
    val facts = sum.resources
    assert(facts eq ExactSum.resources)
    assert(facts eq sum.copy().resources)
    assertEquals(sum.finiteTerms, before)
    assertEquals(sum.value, 1.0)
    assert(facts.stateBytes > 0)
    assertEquals(facts.copyAdditionalBytes, facts.stateBytes)
    assertEquals(facts.valueScratchBytes, facts.stateBytes)
    assert(facts.ratioScratchBytes > facts.valueScratchBytes * 2)
  }

  test("capacity-sized maximum and tiny states retain the same estimates and exact ratios") {
    val numerator = ExactSum.zero()
    val normalizer = ExactSum.zero()
    ok(numerator.add(Double.MaxValue)); ok(normalizer.add(1.0))
    val doubling = 63 - java.lang.Long.numberOfLeadingZeros(ExactSum.MaxTerms)
    (0 until doubling).foreach { _ => ok(numerator.addAll(numerator)); ok(normalizer.addAll(normalizer)) }
    assertEquals(numerator.finiteTerms, ExactSum.MaxTerms)
    assertEquals(ok(numerator.ratio(normalizer)), Double.MaxValue)
    assertEquals(numerator.resources, ExactSum.zero().resources)
    assertEquals(numerator.resources.maximumExactIntegerBits, 2098 + doubling)
    assert(numerator.resources.maximumRatioOperandBits >= numerator.resources.maximumExactIntegerBits + 1074)
    assert(numerator.add(1.0).isLeft)
    assertEquals(ok(numerator.ratio(normalizer)), Double.MaxValue)
  }

  test("signed cancellation, subnormal, overflow and special states remain unmodified") {
    val numerator = ExactSum.zero(); val denominator = ExactSum.zero()
    ok(numerator.add(Double.MinPositiveValue)); ok(denominator.add(2.0))
    val before = (numerator.value, denominator.value, numerator.finiteTerms, denominator.finiteTerms)
    assertEquals(ok(numerator.ratio(denominator)), 0.0)
    val resources = numerator.resources
    assertEquals((numerator.value, denominator.value, numerator.finiteTerms, denominator.finiteTerms), before)
    ok(numerator.add(-Double.MinPositiveValue))
    assertEquals(ok(numerator.ratio(denominator)), 0.0)
    ok(numerator.add(Double.MaxValue)); ok(numerator.add(Double.MaxValue))
    assertEquals(ok(numerator.ratio(denominator)), Double.MaxValue)
    assertEquals(numerator.resources, resources)
    ok(numerator.add(Double.PositiveInfinity))
    assert(numerator.ratio(denominator).isLeft)
    assertEquals(numerator.resources, resources)
  }
