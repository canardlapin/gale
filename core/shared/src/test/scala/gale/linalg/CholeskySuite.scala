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

  test("non-finite vector and strided matrix RHS are rejected without mutating inputs or factor") {
    val factor = Matrix.eye(2).cholesky.orThrow
    val lowerBefore = Vector.tabulate(4)(i => factor.lower(i / 2, i % 2))
    for bad <- Vector(Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity) do
      val vector = Vec(1.0, bad)
      assertEquals(factor.solve(vector), Left(LinAlgError.InvalidArgument("non-finite Cholesky right-hand side")))
      assertEquals(vector(0), 1.0)
      assertEquals(java.lang.Double.doubleToRawLongBits(vector(1)), java.lang.Double.doubleToRawLongBits(bad))
      val matrix = Matrix(2, 2)(1.0, bad, 2.0, 3.0).t
      assertEquals(factor.solve(matrix), Left(LinAlgError.InvalidArgument("non-finite Cholesky right-hand side")))
      assertEquals(java.lang.Double.doubleToRawLongBits(matrix(1, 0)), java.lang.Double.doubleToRawLongBits(bad))
    assertEquals(Vector.tabulate(4)(i => factor.lower(i / 2, i % 2)), lowerBefore)
    assertEquals(factor.solve(Vec(2.0, 3.0)).orThrow.toSeq, Seq(2.0, 3.0))
  }

  test("finite-input overflow cannot publish a successful vector or matrix solution") {
    val factor = Matrix(1, 1)(1e-200).cholesky.orThrow
    val vector = Vec(Double.MaxValue)
    val matrix = Matrix(1, 1)(Double.MaxValue)
    assertEquals(factor.solve(vector), Left(LinAlgError.InvalidArgument("non-finite Cholesky solution")))
    assertEquals(factor.solve(matrix), Left(LinAlgError.InvalidArgument("non-finite Cholesky solution")))
    assertEquals(vector(0), Double.MaxValue)
    assertEquals(matrix(0, 0), Double.MaxValue)
    assertEqualsDouble(factor.solve(Vec(1e-200)).orThrow(0), 1.0, 1e-14)
    assertEqualsDouble(factor.solve(Matrix(1, 1)(1e-200)).orThrow(0, 0), 1.0, 1e-14)
  }

  test("shape failures retain precedence over RHS numerical validation") {
    val factor = Matrix.eye(2).cholesky.orThrow
    assert(factor.solve(Vec(Double.NaN)).left.exists(_.isInstanceOf[LinAlgError.DimensionMismatch]))
    assert(factor.solve(Matrix(1, 1)(Double.NaN)).left.exists(_.isInstanceOf[LinAlgError.DimensionMismatch]))
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
