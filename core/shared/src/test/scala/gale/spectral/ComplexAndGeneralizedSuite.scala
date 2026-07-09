package gale.spectral

class ComplexAndGeneralizedSuite extends munit.FunSuite:
  test("Complex exposes magnitude, conjugate, and scalar division") {
    val z = Complex(3.0, 4.0)
    assertEqualsDouble(z.magnitude, 5.0, 1e-12)
    assertEquals(z.conjugate, Complex(3.0, -4.0))
    assertEquals(z / 2.0, Complex(1.5, 2.0))
    assertEquals(Complex.real(2.0), Complex(2.0, 0.0))
    assertEquals(Complex.Zero, Complex(0.0, 0.0))
  }

  test("GeneralizedEigenvalue: a finite pair reports value = alpha / beta") {
    val ge = GeneralizedEigenvalue(Complex(6.0, -2.0), 2.0)
    assert(ge.isFinitePair)
    assert(!ge.isInfinite)
    assertEquals(ge.value, Complex(3.0, -1.0))
  }

  test("GeneralizedEigenvalue: beta == 0 marks an infinite eigenvalue") {
    val inf = GeneralizedEigenvalue(Complex(1.0, 0.0), 0.0)
    assert(inf.isInfinite)
    assert(!inf.isFinitePair)
    assert(!inf.value.magnitude.isFinite)
  }

  test("GeneralizedSingularValue classifies finite / infinite / zero") {
    assertEqualsDouble(GeneralizedSingularValue.Finite(2.5).value, 2.5, 0.0)
    assert(GeneralizedSingularValue.Infinite.value.isPosInfinity)
    assertEqualsDouble(GeneralizedSingularValue.Zero.value, 0.0, 0.0)
  }
