package gale.linalg

import gale.TestAccess

class TriangularSolveSuite extends munit.FunSuite:
  private def assertClose(actual: DVec, expected: Seq[Double], tol: Double = 1e-12): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < expected.length do
      assert(math.abs(actual(i) - expected(i)) <= tol, s"at $i: ${actual(i)} != ${expected(i)}")
      i += 1

  test("lower solves a hand-computed 3x3 forward substitution") {
    // L y = b, L lower-triangular (non-unit diagonal).
    val l = Matrix.dense(3, 3)(
      2.0, 0.0, 0.0,
      6.0, 3.0, 0.0,
      4.0, 2.0, 5.0
    )
    val b = Vec(4.0, 24.0, 31.0)
    assertClose(TriangularSolve.lower(l, b).orThrow, Seq(2.0, 4.0, 3.0))
  }

  test("upper solves a hand-computed 3x3 back substitution") {
    val u = Matrix.dense(3, 3)(
      5.0, 2.0, 4.0,
      0.0, 3.0, 6.0,
      0.0, 0.0, 2.0
    )
    val b = Vec(21.0, 24.0, 6.0)
    assertClose(TriangularSolve.upper(u, b).orThrow, Seq(1.0, 2.0, 3.0))
  }

  test("a zero diagonal is a Left(SingularMatrix) at that pivot") {
    val l = Matrix.dense(3, 3)(
      2.0, 0.0, 0.0,
      6.0, 0.0, 0.0,
      4.0, 2.0, 5.0
    )
    val b = Vec(4.0, 24.0, 31.0)
    TriangularSolve.lower(l, b) match
      case Left(LinAlgError.SingularMatrix(index)) => assertEquals(index, 1)
      case other                                   => fail(s"expected Left(SingularMatrix(1)), got $other")

    val u = Matrix.dense(3, 3)(
      5.0, 2.0, 4.0,
      0.0, 3.0, 6.0,
      0.0, 0.0, 0.0
    )
    TriangularSolve.upper(u, Vec(21.0, 24.0, 6.0)) match
      case Left(LinAlgError.SingularMatrix(index)) => assertEquals(index, 2)
      case other                                   => fail(s"expected Left(SingularMatrix(2)), got $other")
  }

  test("a strided-view matrix and vector give the same solution as contiguous storage") {
    // L realised as a strided submatrix view: values interleaved with sentinels.
    val lBacking = TestAccess.doubleArray(
      2.0, -1.0, 0.0, -1.0, 0.0, -1.0,
      6.0, -1.0, 3.0, -1.0, 0.0, -1.0,
      4.0, -1.0, 2.0, -1.0, 5.0, -1.0
    )
    val lStrided = TestAccess.mat(lBacking, offset = 0, rows = 3, cols = 3, rowStride = 6, colStride = 2)
    val bStrided = TestAccess.stridedCopy(Vec(4.0, 24.0, 31.0), stride = 3)

    val contiguous =
      TriangularSolve
        .lower(Matrix.dense(3, 3)(2.0, 0.0, 0.0, 6.0, 3.0, 0.0, 4.0, 2.0, 5.0), Vec(4.0, 24.0, 31.0))
        .orThrow
    val strided = TriangularSolve.lower(lStrided, bStrided).orThrow
    assertClose(strided, contiguous.toSeq)
  }

  test("dimension guards reject non-square and length-mismatched inputs") {
    TriangularSolve.lower(Matrix.zeros(2, 3), Vec(1.0, 2.0)) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected NonSquareMatrix, got $other")
    TriangularSolve.upper(Matrix.eye(3), Vec(1.0, 2.0)) match
      case Left(_: LinAlgError.DimensionMismatch) => ()
      case other                                  => fail(s"expected DimensionMismatch, got $other")
  }
