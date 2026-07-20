package gale.linalg

import gale.backend.DenseDoubleFactorizations
import scala.compiletime.testing.typeCheckErrors

class FactorizationCapabilitySuite extends munit.FunSuite:
  private val square = Matrix.dense(3, 3)(
    6.0, 2.0, 1.0,
    2.0, 5.0, 2.0,
    1.0, 2.0, 4.0
  )

  private val tall = Matrix.dense(5, 3)(
    1.0, 0.0, 2.0,
    1.0, 1.0, -1.0,
    1.0, 2.0, 0.5,
    1.0, 3.0, 1.5,
    1.0, 4.0, -0.5
  )

  test("LU and Cholesky satisfy the exact-solve vector and matrix laws") {
    val expectedVector = Vec(1.0, -2.0, 3.0)
    val expectedMatrix = Matrix.dense(3, 2)(
      1.0, -2.0,
      3.0, 0.5,
      -1.0, 4.0
    )
    val vectorRhs = square * expectedVector
    // Form the same product transposed, then take its O(1) transpose view. The
    // logical RHS is 3x2 but its storage is non-contiguous row-major.
    val matrixRhs = (expectedMatrix.t * square.t).t
    val factors: Seq[ExactSolveFactor] = Seq(square.lu.orThrow, square.cholesky.orThrow)

    factors.foreach: factor =>
      assertEquals(factor.size, 3)
      assert(factor.diagnostics.isSuccess)
      assertVectorClose(factor.solve(vectorRhs).orThrow, expectedVector, 1e-11)
      assertMatrixClose(factor.solve(matrixRhs).orThrow, expectedMatrix, 1e-11)
  }

  test("LU matrix right-hand sides preserve pivoting and typed dimension errors") {
    val pivoted = Matrix.dense(3, 3)(
      0.0, 2.0, 1.0,
      1.0, -2.0, 0.0,
      3.0, 1.0, 4.0
    )
    val expected = Matrix.dense(3, 2)(
      1.0, -1.0,
      2.0, 0.5,
      -3.0, 4.0
    )
    val factor: ExactSolveFactor = pivoted.lu.orThrow

    assertMatrixClose(factor.solve(pivoted * expected).orThrow, expected, 1e-11)
    assert(factor.solve(Vec(1.0, 2.0)).left.exists(_.isInstanceOf[LinAlgError.DimensionMismatch]))
    assert(factor.solve(Matrix.zeros(2, 4)).left.exists(_.isInstanceOf[LinAlgError.DimensionMismatch]))
  }

  test("QR satisfies least-squares vector and matrix laws") {
    val expectedVector = Vec(1.0, -2.0, 3.0)
    val expectedMatrix = Matrix.dense(3, 2)(
      1.0, -2.0,
      3.0, 0.5,
      -1.0, 4.0
    )
    val factor: LeastSquaresFactor = tall.qr

    assertEquals(factor.observationCount, 5)
    assertEquals(factor.coefficientCount, 3)
    assert(factor.diagnostics.isSuccess)
    assertVectorClose(factor.solveLeastSquares(tall * expectedVector).orThrow, expectedVector, 1e-10)
    assertMatrixClose(factor.solveLeastSquares(tall * expectedMatrix).orThrow, expectedMatrix, 1e-10)
    assert(
      factor
        .solveLeastSquares(Vec(1.0, 2.0))
        .left
        .exists(_.isInstanceOf[LinAlgError.DimensionMismatch])
    )
  }

  test("rank deficiency remains a typed least-squares failure") {
    val deficient = Matrix.dense(3, 2)(
      1.0, 2.0,
      2.0, 4.0,
      3.0, 6.0
    )
    val factor: LeastSquaresFactor = deficient.qr
    assertEquals(
      factor.solveLeastSquares(Vec(1.0, 2.0, 3.0)),
      Left(LinAlgError.RankDeficient(1, 2))
    )
  }

  test("backend-produced dense factors retain the same capabilities") {
    val provider = new DenseDoubleFactorizations:
      def lu(a: DMat): Either[LinAlgError, LU] = DenseDecompositions.lu(a)
      def cholesky(a: DMat): Either[LinAlgError, Cholesky] = DenseDecompositions.cholesky(a)
      def qr(a: DMat): Either[LinAlgError, QR] = Right(DenseDecompositions.qr(a))

    val exact: ExactSolveFactor = provider.lu(square).orThrow
    val exactSpd: ExactSolveFactor = provider.cholesky(square).orThrow
    val leastSquares: LeastSquaresFactor = provider.qr(tall).orThrow
    assertVectorClose(exact.solve(square * Vec(2.0, -1.0, 0.5)).orThrow, Vec(2.0, -1.0, 0.5), 1e-11)
    assertVectorClose(exactSpd.solve(square * Vec(2.0, -1.0, 0.5)).orThrow, Vec(2.0, -1.0, 0.5), 1e-11)
    assertVectorClose(
      leastSquares.solveLeastSquares(tall * Vec(2.0, -1.0, 0.5)).orThrow,
      Vec(2.0, -1.0, 0.5),
      1e-10
    )
  }

  test("compile-time API keeps exact, least-squares, and concrete extras distinct") {
    val qrIsNotExact = typeCheckErrors("""
      import gale.linalg.*
      def requireExact(factor: ExactSolveFactor): Unit = ()
      def invalid(qr: QR): Unit = requireExact(qr)
    """)
    val luIsNotLeastSquares = typeCheckErrors("""
      import gale.linalg.*
      def requireLeastSquares(factor: LeastSquaresFactor): Unit = ()
      def invalid(lu: LU): Unit = requireLeastSquares(lu)
    """)
    val determinantIsConcrete = typeCheckErrors("""
      import gale.linalg.*
      def invalid(factor: ExactSolveFactor) = factor.det
    """)
    val orthogonalFactorIsConcrete = typeCheckErrors("""
      import gale.linalg.*
      def invalid(factor: LeastSquaresFactor) = factor.q
    """)

    assert(qrIsNotExact.nonEmpty)
    assert(luIsNotLeastSquares.nonEmpty)
    assert(determinantIsConcrete.nonEmpty)
    assert(orthogonalFactorIsConcrete.nonEmpty)
  }

  private def assertVectorClose(actual: DVec, expected: DVec, tolerance: Double): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assert(math.abs(actual(i) - expected(i)) <= tolerance, s"index $i")
      i += 1

  private def assertMatrixClose(actual: DMat, expected: DMat, tolerance: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assert(math.abs(actual(row, col) - expected(row, col)) <= tolerance, s"($row,$col)")
        col += 1
      row += 1
