package gale.linalg

class CholeskySuite extends munit.FunSuite:
  test("Cholesky reconstructs a symmetric positive-definite matrix") {
    val A = Matrix.dense(2, 2)(
      4.0, 1.0,
      1.0, 3.0
    )

    val ch = A.cholesky.orThrow
    val reconstructed = ch.lower * ch.lower.t

    assertMatrixClose(reconstructed, A, 1e-12)
  }

  test("Cholesky solve has small residual") {
    val A = Matrix.dense(2, 2)(
      4.0, 1.0,
      1.0, 3.0
    )
    val b = Vec(1.0, 2.0)

    val x = A.cholesky.orThrow.solve(b).orThrow
    val r = A * x - b

    assert(norm(r) < 1e-12)
    assert(math.abs(x(0) - 1.0 / 11.0) < 1e-12)
    assert(math.abs(x(1) - 7.0 / 11.0) < 1e-12)
  }

  test("Cholesky reports non-positive-definite and non-square inputs") {
    val indefinite = Matrix.dense(2, 2)(
      1.0, 2.0,
      2.0, 1.0
    )
    val rectangular = Matrix.zeros(2, 3)

    assertEquals(indefinite.cholesky, Left(LinAlgError.NotPositiveDefinite(1)))
    assert(rectangular.cholesky.left.exists(_.isInstanceOf[LinAlgError.NonSquareMatrix]))
  }

  test("Cholesky solves matrix right-hand sides in one result") {
    val A = Matrix.dense(3, 3)(
      6.0, 2.0, 1.0,
      2.0, 5.0, 2.0,
      1.0, 2.0, 4.0
    )
    val expected = Matrix.dense(3, 2)(
      1.0, -2.0,
      3.0, 0.5,
      -1.0, 4.0
    )

    val actual = A.cholesky.orThrow.solve(A * expected).orThrow

    assertMatrixClose(actual, expected, 1e-11)
  }

  test("explicit Cholesky pivot tolerance rejects numerically tiny pivots") {
    val A = Matrix.dense(2, 2)(
      1.0, 0.0,
      0.0, 1.0e-14
    )

    assert(A.cholesky.isRight)
    assertEquals(
      A.cholesky(CholeskyOptions(pivotTolerance = 1.0e-12)),
      Left(LinAlgError.NotPositiveDefinite(1))
    )
  }

  private def assertMatrixClose(actual: DMat, expected: DMat, tolerance: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        assert(math.abs(actual(i, j) - expected(i, j)) < tolerance)
        j += 1
      i += 1
