package gale.linalg

class LinearOperatorSuite extends munit.FunSuite:
  test("DMat acts as a DoubleLinearOperator") {
    val A = Matrix.dense(2, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0
    )
    val y = MutableVec.zeros(2)

    A.applyTo(Vec(1.0, 2.0, 1.0), y)
    assertEquals(y.asVec.toSeq, Seq(8.0, 20.0))

    val z = MutableVec.zeros(3)
    A.transposeApplyTo(Vec(1.0, -1.0), z)
    assertEquals(z.asVec.toSeq, Seq(-3.0, -3.0, -3.0))
  }

  test("matrix-free operators validate shape and write into destination") {
    val op =
      LinearOperator.fromFunction(2, 2): (x, into) =>
        into(0) = 2.0 * x(0)
        into(1) = 3.0 * x(1)

    assertEquals((op * Vec(4.0, 5.0)).toSeq, Seq(8.0, 15.0))

    intercept[LinAlgError.VectorLengthMismatch] {
      op * Vec(1.0, 2.0, 3.0)
    }
  }

  test("operators apply to matrix right-hand sides in one batch result") {
    val operator = Matrix.dense(2, 3)(
      1.0, 2.0, 0.0,
      -1.0, 0.0, 3.0
    )
    val input = Matrix.dense(3, 2)(
      1.0, 4.0,
      2.0, 5.0,
      3.0, 6.0
    )

    val actual = operator.applyTo(input).orThrow

    assertMatrixClose(actual, operator * input, 1e-12)
    assert(operator.applyTo(Matrix.zeros(4, 1)).left.exists(_.isInstanceOf[LinAlgError.DimensionMismatch]))
  }

  test("adjoint, composition, scaling, and restrictions obey operator laws") {
    val a = Matrix.dense(3, 2)(
      1.0, 2.0,
      0.0, -1.0,
      3.0, 1.0
    )
    val b = Matrix.dense(2, 3)(
      2.0, 0.0, 1.0,
      -1.0, 4.0, 0.0
    )
    val x = Vec(0.5, -2.0, 3.0)
    val y = Vec(1.0, -1.0, 2.0)

    val composed = a.compose(b).orThrow
    assertVectorClose(composed(x), a * (b * x), 1e-12)
    assert(math.abs(composed(x).dot(y) - x.dot(composed.adjoint(y))) < 1e-12)
    assertVectorClose(composed.scaled(-2.0)(x), composed(x) * -2.0, 1e-12)

    val rowRestricted = a.restrictRows(Vector(2, 0)).orThrow
    assertVectorClose(rowRestricted(Vec(2.0, -1.0)), Vec(5.0, 0.0), 1e-12)
    val colRestricted = a.restrictColumns(Vector(1)).orThrow
    assertVectorClose(colRestricted(Vec(2.0)), Vec(4.0, -2.0, 2.0), 1e-12)
    assert(a.restrictRows(Vector(0, 0)).isLeft)
  }

  test("block diagonal and rectangular block operators match dense assembly") {
    val a = Matrix.dense(2, 2)(
      1.0, 2.0,
      3.0, 4.0
    )
    val b = Matrix.dense(1, 1)(5.0)
    val diagonal = LinearOperator.blockDiagonal(Vector(a, b)).orThrow
    assertVectorClose(diagonal(Vec(1.0, 2.0, 3.0)), Vec(5.0, 11.0, 15.0), 1e-12)
    assertVectorClose(diagonal.adjoint(Vec(1.0, 2.0, 3.0)), Vec(7.0, 10.0, 15.0), 1e-12)

    val topLeft = Matrix.dense(1, 2)(1.0, 2.0)
    val topRight = Matrix.dense(1, 1)(3.0)
    val bottomLeft = Matrix.dense(2, 2)(
      4.0, 5.0,
      6.0, 7.0
    )
    val bottomRight = Matrix.dense(2, 1)(8.0, 9.0)
    val blocked = LinearOperator.block(
      Vector(
        Vector(topLeft, topRight),
        Vector(bottomLeft, bottomRight)
      )
    ).orThrow
    val dense = Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 8.0,
      6.0, 7.0, 9.0
    )
    val input = Vec(1.0, -2.0, 0.5)

    assertVectorClose(blocked(input), dense * input, 1e-12)
    assertVectorClose(blocked.adjoint(input), dense.t * input, 1e-12)
  }

  private def assertVectorClose(actual: DVec, expected: DVec, tolerance: Double): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assert(math.abs(actual(i) - expected(i)) <= tolerance, s"index $i: ${actual(i)} != ${expected(i)}")
      i += 1

  private def assertMatrixClose(actual: DMat, expected: DMat, tolerance: Double): Unit =
    assertEquals((actual.rows, actual.cols), (expected.rows, expected.cols))
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assert(math.abs(actual(row, col) - expected(row, col)) <= tolerance)
        col += 1
      row += 1
