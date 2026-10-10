package gale.parity

import gale.parity.GenerateBreezeGoldens.{Arr, Case, compare}

class GoldenSpectralFreshnessSuite extends munit.FunSuite:
  private def fixture(sigma: Array[Double]): Case =
    Case("svd", "svd/synthetic", Vector(
      Arr("A", 3, 2, Array(sigma(0), 0.0, 0.0, sigma(1), 0.0, 0.0), true),
      Arr("sigma", 2, 1, sigma, false),
      Arr("U", 3, 2, Array(1.0, 0.0, 0.0, 1.0, 0.0, 0.0), false),
      Arr("VT", 2, 2, Array(1.0, 0.0, 0.0, 1.0), false)
    ))

  private def replace(c: Case, key: String, values: Array[Double]): Case =
    c.copy(arrays = c.arrays.map(a => if a.key == key then a.copy(values = values) else a))

  test("freshness accepts paired signs but rejects a one-sided sign change") {
    val old = fixture(Array(3.0, 1.0))
    val left = replace(old, "U", Array(-1.0, 0.0, 0.0, 1.0, 0.0, 0.0))
    assert(compare(Vector(old), Vector(left)).exists(_.contains("reconstruction")))
    val paired = replace(left, "VT", Array(-1.0, 0.0, 0.0, 1.0))
    assertEquals(compare(Vector(old), Vector(paired)), Vector.empty)
  }

  test("freshness accepts common repeated-cluster rotations and rejects corrupt ambient rows") {
    val old = fixture(Array(3.0, 3.0))
    val left = replace(old, "U", Array(0.6, -0.8, 0.8, 0.6, 0.0, 0.0))
    val rotated = replace(left, "VT", Array(0.6, 0.8, -0.8, 0.6))
    assertEquals(compare(Vector(old), Vector(rotated)), Vector.empty)
    val corrupt = replace(old, "U", Array(1.0, 0.0, 0.0, 1.0, 1e-3, 0.0))
    assert(compare(Vector(old), Vector(corrupt)).nonEmpty)
  }

  test("freshness permits a different valid economy null-space completion") {
    val old = fixture(Array(3.0, 0.0))
    val other = replace(old, "U", Array(1.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    assertEquals(compare(Vector(old), Vector(other)), Vector.empty)
  }
