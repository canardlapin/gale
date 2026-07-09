package gale.spectral

import gale.linalg.*

class NonsymmetricEigenPackingSuite extends munit.FunSuite:
  // Spectrum: 2 (real), 1 + 3i / 1 − 3i (conjugate pair), −5 (real). The pair's
  // positive-imaginary member comes first, per the real-Schur packing invariant.
  private val re = Vec(2.0, 1.0, 1.0, -5.0)
  private val im = Vec(0.0, 3.0, -3.0, 0.0)

  // Right eigenvector packing:
  //   col0 = real eigenvector of λ = 2
  //   col1 = v_re, col2 = v_im of the conjugate pair
  //   col3 = real eigenvector of λ = −5
  private val rightPacked = Matrix.dense(4, 4)(
    1.0, 0.5, 0.1, 0.0,
    0.0, 0.6, 0.2, 0.0,
    0.0, 0.0, 0.7, 0.0,
    0.0, 0.0, 0.8, 1.0
  )

  private def decomp(left: Option[DMat]): NonsymmetricEigenDecomposition =
    new NonsymmetricEigenDecomposition(
      re,
      im,
      rightPacked,
      left,
      SpectralDiagnostics(4, 4, Vec(0.0, 0.0, 0.0, 0.0), 0.0, 0)
    )

  private def assertVec(actual: DVec, expected: Seq[Double])(using munit.Location): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < expected.length do
      assertEqualsDouble(actual(i), expected(i), 1e-12)
      i += 1

  /** Build a decomposition from raw spectrum and packing, with all-converged
    * diagnostics sized to the spectrum.
    */
  private def build(reV: Seq[Double], imV: Seq[Double], packed: DMat): NonsymmetricEigenDecomposition =
    new NonsymmetricEigenDecomposition(
      Vec(reV*),
      Vec(imV*),
      packed,
      None,
      SpectralDiagnostics(reV.length, reV.length, DVec.zeros(reV.length), 0.0, 0)
    )

  test("eigenvalue reads the structure-of-arrays real/imaginary parts") {
    val d = decomp(None)
    assertEquals(d.size, 4)
    assertEquals(d.eigenvalue(0), Complex(2.0, 0.0))
    assertEquals(d.eigenvalue(1), Complex(1.0, 3.0))
    assertEquals(d.eigenvalue(2), Complex(1.0, -3.0))
    assertEquals(d.eigenvalue(3), Complex(-5.0, 0.0))
  }

  test("isRealPair distinguishes real eigenvalues from conjugate-pair members") {
    val d = decomp(None)
    assert(d.isRealPair(0))
    assert(!d.isRealPair(1))
    assert(!d.isRealPair(2))
    assert(d.isRealPair(3))
  }

  test("eigenvector decodes real-Schur packing, including conjugate signs") {
    val d = decomp(None)

    // Real eigenvalues: single column, zero imaginary part.
    val (re0, im0) = d.eigenvector(0)
    assertVec(re0, Seq(1.0, 0.0, 0.0, 0.0))
    assertVec(im0, Seq(0.0, 0.0, 0.0, 0.0))
    val (re3, im3) = d.eigenvector(3)
    assertVec(re3, Seq(0.0, 0.0, 0.0, 1.0))
    assertVec(im3, Seq(0.0, 0.0, 0.0, 0.0))

    // Positive-imaginary member: v_re + i·v_im over the two adjacent columns.
    val (re1, im1) = d.eigenvector(1)
    assertVec(re1, Seq(0.5, 0.6, 0.0, 0.0))
    assertVec(im1, Seq(0.1, 0.2, 0.7, 0.8))

    // Negative-imaginary member (conjugate): same real part, negated imaginary part.
    val (re2, im2) = d.eigenvector(2)
    assertVec(re2, Seq(0.5, 0.6, 0.0, 0.0))
    assertVec(im2, Seq(-0.1, -0.2, -0.7, -0.8))
  }

  test("left eigenvectors decode when present and are unsupported when absent") {
    val withLeft = decomp(Some(rightPacked))
    val (lre, lim) = withLeft.leftEigenvector(1)
    assertVec(lre, Seq(0.5, 0.6, 0.0, 0.0))
    assertVec(lim, Seq(0.1, 0.2, 0.7, 0.8))

    intercept[LinAlgError.UnsupportedOperation] {
      decomp(None).leftEigenvector(0)
    }
  }

  test("accessors bounds-check and report vectors-not-computed") {
    val d = decomp(None)
    intercept[LinAlgError.IndexOutOfBounds](d.eigenvalue(4))
    intercept[LinAlgError.IndexOutOfBounds](d.eigenvector(-1))

    val valuesOnly = new NonsymmetricEigenDecomposition(
      re,
      im,
      DMat.zeros(4, 0),
      None,
      SpectralDiagnostics(4, 4, Vec(0.0, 0.0, 0.0, 0.0), 0.0, 0)
    )
    intercept[LinAlgError.UnsupportedOperation](valuesOnly.eigenvector(0))
  }

  test("re/im length mismatch is rejected at construction") {
    intercept[IllegalArgumentException] {
      new NonsymmetricEigenDecomposition(
        Vec(1.0, 2.0),
        Vec(0.0),
        DMat.zeros(2, 0),
        None,
        SpectralDiagnostics(2, 2, Vec(0.0, 0.0), 0.0, 0)
      )
    }
  }

  test("trailing conjugate pair decodes col(i+1) at the last column") {
    // Spectrum: 3, −2 (real), then 1 ± 5i occupying the final two columns.
    val d = build(
      Seq(3.0, -2.0, 1.0, 1.0),
      Seq(0.0, 0.0, 5.0, -5.0),
      Matrix.dense(4, 4)(
        1.0, 0.0, 0.2, 0.6,
        0.0, 1.0, 0.3, 0.7,
        0.0, 0.0, 0.4, 0.8,
        0.0, 0.0, 0.5, 0.9
      )
    )
    assertEquals(d.eigenvalue(2), Complex(1.0, 5.0))
    assertEquals(d.eigenvalue(3), Complex(1.0, -5.0))
    val (re2, im2) = d.eigenvector(2)
    assertVec(re2, Seq(0.2, 0.3, 0.4, 0.5))
    assertVec(im2, Seq(0.6, 0.7, 0.8, 0.9))
    val (re3, im3) = d.eigenvector(3)
    assertVec(re3, Seq(0.2, 0.3, 0.4, 0.5))
    assertVec(im3, Seq(-0.6, -0.7, -0.8, -0.9))
  }

  test("leading conjugate pair decodes col(i-1) at index 0/1") {
    // Spectrum: 4 ± 6i occupying the first two columns, then −7 (real).
    val d = build(
      Seq(4.0, 4.0, -7.0),
      Seq(6.0, -6.0, 0.0),
      Matrix.dense(3, 3)(
        0.1, 0.4, 0.0,
        0.2, 0.5, 0.0,
        0.3, 0.6, 1.0
      )
    )
    assertEquals(d.eigenvalue(0), Complex(4.0, 6.0))
    assertEquals(d.eigenvalue(1), Complex(4.0, -6.0))
    val (re0, im0) = d.eigenvector(0)
    assertVec(re0, Seq(0.1, 0.2, 0.3))
    assertVec(im0, Seq(0.4, 0.5, 0.6))
    val (re1, im1) = d.eigenvector(1)
    assertVec(re1, Seq(0.1, 0.2, 0.3))
    assertVec(im1, Seq(-0.4, -0.5, -0.6))
  }

  test("two adjacent conjugate pairs decode without crossing pair boundaries") {
    // Spectrum: 2 ± 1i (cols 0,1) then −3 ± 8i (cols 2,3); distinct blocks so a
    // boundary-crossing decode would be caught.
    val d = build(
      Seq(2.0, 2.0, -3.0, -3.0),
      Seq(1.0, -1.0, 8.0, -8.0),
      Matrix.dense(4, 4)(
        1.0, 1.2, 0.0, 0.0,
        1.1, 1.3, 0.0, 0.0,
        0.0, 0.0, 2.0, 2.2,
        0.0, 0.0, 2.1, 2.3
      )
    )
    val (re0, im0) = d.eigenvector(0)
    assertVec(re0, Seq(1.0, 1.1, 0.0, 0.0))
    assertVec(im0, Seq(1.2, 1.3, 0.0, 0.0))
    val (re1, im1) = d.eigenvector(1)
    assertVec(re1, Seq(1.0, 1.1, 0.0, 0.0))
    assertVec(im1, Seq(-1.2, -1.3, 0.0, 0.0))
    val (re2, im2) = d.eigenvector(2)
    assertVec(re2, Seq(0.0, 0.0, 2.0, 2.1))
    assertVec(im2, Seq(0.0, 0.0, 2.2, 2.3))
    val (re3, im3) = d.eigenvector(3)
    assertVec(re3, Seq(0.0, 0.0, 2.0, 2.1))
    assertVec(im3, Seq(0.0, 0.0, -2.2, -2.3))
  }

  test("malformed packing is rejected at construction") {
    // Split pair: the positive-imaginary member is not followed by its conjugate.
    intercept[IllegalArgumentException] {
      build(Seq(1.0, 1.0, 1.0), Seq(3.0, 0.0, -3.0), DMat.zeros(3, 0))
    }
    // Negative-imaginary member first.
    intercept[IllegalArgumentException] {
      build(Seq(1.0, 1.0), Seq(-2.0, 2.0), DMat.zeros(2, 0))
    }
    // Adjacent but mismatched magnitudes (not actually conjugates).
    intercept[IllegalArgumentException] {
      build(Seq(1.0, 1.0), Seq(3.0, -4.0), DMat.zeros(2, 0))
    }
    // Packed columns neither 0 nor size.
    intercept[IllegalArgumentException] {
      build(Seq(1.0, -1.0), Seq(0.0, 0.0), DMat.zeros(2, 1))
    }
  }
