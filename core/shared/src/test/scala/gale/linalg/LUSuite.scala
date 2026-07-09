package gale.linalg

class LUSuite extends munit.FunSuite:
  test("LU solves a dense square system with small residual") {
    val A = Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0,
      7.0, 8.0, 10.0
    )
    val b = Vec(3.0, 3.0, 4.0)

    val x = A.solve(b).orThrow
    val r = A * x - b

    assert(norm(r) < 1e-10)
  }

  test("LU uses partial pivoting") {
    val A = Matrix.dense(2, 2)(
      0.0, 1.0,
      2.0, 3.0
    )
    val b = Vec(1.0, 5.0)

    val lu = A.lu.orThrow
    val x = lu.solve(b).orThrow

    assertEquals(lu.pivots.toArray.toSeq, Seq(1, 0))
    assert(math.abs(x(0) - 1.0) < 1e-12)
    assert(math.abs(x(1) - 1.0) < 1e-12)
  }

  test("LU determinant comes from packed factors and pivot parity") {
    val A = Matrix.dense(2, 2)(
      1.0, 2.0,
      3.0, 4.0
    )

    assert(math.abs(A.lu.orThrow.det.orThrow + 2.0) < 1e-12)
  }

  test("singular and non-square matrices return diagnostics") {
    val singular = Matrix.dense(2, 2)(
      1.0, 2.0,
      2.0, 4.0
    )
    val rectangular = Matrix.zeros(2, 3)

    assertEquals(singular.lu, Left(LinAlgError.SingularMatrix(1)))
    assert(rectangular.lu.left.exists(_.isInstanceOf[LinAlgError.NonSquareMatrix]))
  }
