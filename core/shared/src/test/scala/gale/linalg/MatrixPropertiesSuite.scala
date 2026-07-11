package gale.linalg

class MatrixPropertiesSuite extends munit.FunSuite:
  private val symmetric = Matrix.dense(3, 3)(
    4.0, 1.0, -2.0,
    1.0, 5.0, 0.5,
    -2.0, 0.5, 6.0
  )

  test("assume wrappers are zero-cost subtypes and compose") {
    val sym: Symmetric[DMat] = symmetric.assumeSymmetric
    val spd: PositiveDefinite[DMat] = sym.assumePositiveDefinite
    val lower: LowerTriangular[DMat] = symmetric.assumeLowerTriangular
    val upper: UpperTriangular[DMat] = symmetric.assumeUpperTriangular

    summon[Symmetric[DMat] <:< DMat]
    summon[PositiveDefinite[DMat] <:< Symmetric[DMat]]
    summon[LowerTriangular[DMat] <:< DMat]
    summon[UpperTriangular[DMat] <:< DMat]
    assert(sym.value eq symmetric)
    assert(spd.value eq symmetric)
    assert(lower.value eq symmetric)
    assert(upper.value eq symmetric)
    assertEquals((spd * Vec(1.0, 0.0, 0.0)).length, 3)
  }

  test("verifySymmetric uses a scale-aware tolerance") {
    val nearly = Matrix.dense(2, 2)(
      1.0e12, 3.0e8,
      3.0e8 + 1.0e-4, 2.0e12
    )

    assert(nearly.verifySymmetric(tolerance = 1.0e-12).isRight)
    assert(nearly.verifySymmetric(tolerance = 1.0e-14).isLeft)
    assert(symmetric.verifySymmetric(tolerance = 0.0).isRight)
  }

  test("verifySymmetric rejects structural, tolerance, and non-finite violations") {
    val rectangular = Matrix.zeros(2, 3)
    val asymmetric = symmetric.updated(0, 1, 9.0)
    val nonFinite = symmetric.updated(0, 1, Double.NaN).updated(1, 0, Double.NaN)

    assert(rectangular.verifySymmetric().left.exists(_.isInstanceOf[LinAlgError.NonSquareMatrix]))
    assert(asymmetric.verifySymmetric().left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(nonFinite.verifySymmetric().left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(symmetric.verifySymmetric(tolerance = -1.0).left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(symmetric.verifySymmetric(tolerance = Double.PositiveInfinity).left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
  }

  test("verifyPositiveDefinite proves symmetry as well as Cholesky success") {
    val spd = symmetric.verifyPositiveDefinite.orThrow
    summon[PositiveDefinite[DMat] <:< Symmetric[DMat]]
    assert(spd.value eq symmetric)
    assert(spd.cholesky.isRight)

    val indefinite = Matrix.dense(2, 2)(
      1.0, 2.0,
      2.0, 1.0
    )
    assertEquals(indefinite.verifyPositiveDefinite, Left(LinAlgError.NotPositiveDefinite(1)))

    // Cholesky itself reads only the lower triangle; property verification must
    // not certify this matrix merely because that triangle is SPD.
    val nonsymmetric = Matrix.dense(2, 2)(
      2.0, 99.0,
      0.0, 2.0
    )
    assert(nonsymmetric.cholesky.isRight)
    assert(nonsymmetric.verifyPositiveDefinite.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
  }

  test("positive-definite verification is invariant to positive common scaling") {
    val tiny = Matrix.tabulate(symmetric.rows, symmetric.cols)((i, j) => symmetric(i, j) * 1.0e-200)
    val huge = Matrix.tabulate(symmetric.rows, symmetric.cols)((i, j) => symmetric(i, j) * 1.0e200)

    assert(tiny.verifyPositiveDefinite.isRight)
    assert(huge.verifyPositiveDefinite.isRight)
  }

  test("triangular verification checks the forbidden triangle with tolerance") {
    val lower = Matrix.dense(3, 3)(
      1.0, 0.0, 0.0,
      2.0, 3.0, 0.0,
      4.0, 5.0, 6.0
    )
    val upper = lower.t
    assert(lower.verifyLowerTriangular.isRight)
    assert(upper.verifyUpperTriangular.isRight)
    assert(lower.verifyUpperTriangular.isLeft)
    assert(upper.verifyLowerTriangular.isLeft)

    val almostLower = lower.updated(0, 2, 1.0e-13)
    assert(almostLower.verifyLowerTriangular.isLeft)
    assert(almostLower.verifyLowerTriangular(tolerance = 1.0e-12).isRight)
    assert(lower.updated(0, 2, Double.NaN).verifyLowerTriangular.isLeft)
    assert(Matrix.zeros(2, 3).verifyLowerTriangular.left.exists(_.isInstanceOf[LinAlgError.NonSquareMatrix]))
  }
