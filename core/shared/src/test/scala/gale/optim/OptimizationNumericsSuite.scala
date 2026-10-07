package gale.optim

final class OptimizationNumericsSuite extends munit.FunSuite:
  test("scaled norms preserve zeros, mixed magnitudes, subnormals and overflow"):
    val fixtures = Vector(
      Array(0.0, -0.0),
      Array(3.0, 4.0),
      Array(1e155, -1e155),
      Array(1e-170, -1e-170),
      Array(1e-200, 1e200, -2e199),
      Array(Double.MinPositiveValue, Double.MinPositiveValue),
      Array(Double.MaxValue, Double.MaxValue)
    )
    fixtures.foreach: values =>
      val expected = values.foldLeft(0.0)(Math.hypot)
      val actual = OptimNumerics.norm(values)
      if expected.isInfinite then assert(actual.isInfinite)
      else assert(Math.abs(actual - expected) <= 4.0 * Math.ulp(expected), s"$actual != $expected")
    assert(OptimNumerics.norm(Array(Double.NaN, 1.0)).isNaN)
