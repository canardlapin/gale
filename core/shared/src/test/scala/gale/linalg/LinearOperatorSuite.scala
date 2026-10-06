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

  test("allocating operator application snapshots a retained destination") {
    var retained: MutableDVec = null
    val operator =
      LinearOperator.fromFunction(2, 2): (input, output) =>
        retained = output
        output(0) = input(0) + 1.0
        output(1) = input(1) + 2.0

    val result = operator(Vec(3.0, 4.0))
    retained(0) = 99.0

    assertEquals(result.toSeq, Seq(4.0, 6.0))
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

  test("Kronecker operators match dense products without materialization") {
    val left = Matrix.dense(2, 3)(
      1.0, 2.0, -1.0,
      0.5, 0.0, 3.0
    )
    val right = Matrix.dense(3, 2)(
      2.0, 0.0,
      -1.0, 4.0,
      0.5, 1.0
    )
    val operator = LinearOperator.kronecker(left, right).orThrow
    val dense = left.kron(right)
    val input = Vec(1.0, -2.0, 0.5, 3.0, -1.0, 2.0)
    val dualInput = Vec(1.0, 0.0, -1.0, 2.0, 0.5, -0.5)

    assertEquals((operator.rows, operator.cols), (dense.rows, dense.cols))
    assert(operator.adjoint.isInstanceOf[KroneckerLinearOperator])
    assertVectorClose(operator(input), dense * input, 1e-12)
    assertVectorClose(operator.adjoint(dualInput), dense.t * dualInput, 1e-12)
  }

  test("Kronecker operator rejects an unstorable vector shape") {
    val huge = LinearOperator.fromFunctions(Int.MaxValue, 1)((_, _) => (), (_, _) => ())
    val small = LinearOperator.fromFunctions(2, 1)((_, _) => (), (_, _) => ())

    assert(LinearOperator.kronecker(huge, small).left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
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

  test("empty row restrictions are zero-sized forward maps with zero adjoints") {
    var calls = 0
    val source = LinearOperator.fromFunctions(2, 3)(
      (_, _) => { calls += 1; throw new AssertionError("empty restriction evaluated forward") },
      (_, _) => { calls += 1; throw new AssertionError("empty restriction evaluated adjoint") }
    )
    val restricted = source.restrictRows(Vector.empty).orThrow
    assertEquals((restricted.rows, restricted.cols), (0, 3))
    assertEquals(restricted(Vec(1.0, 2.0, 3.0)).length, 0)
    val output = MutableVec.zeros(3)
    output(0) = Double.NaN
    output(1) = Double.PositiveInfinity
    output(2) = -7.0
    restricted.transposeApplyTo(Vec.zeros(0), output)
    assertEquals(output.asVec.toSeq, Seq(0.0, 0.0, 0.0))
    assertEquals(calls, 0)
    intercept[LinAlgError.VectorLengthMismatch] { restricted(Vec(1.0)) }
    assertEquals(calls, 0)
  }

  test("empty column restrictions are zero maps with zero-sized adjoints") {
    var calls = 0
    val source = LinearOperator.fromFunctions(2, 3)(
      (_, _) => { calls += 1; throw new AssertionError("empty restriction evaluated forward") },
      (_, _) => { calls += 1; throw new AssertionError("empty restriction evaluated adjoint") }
    )
    val restricted = source.restrictColumns(Vector.empty).orThrow
    assertEquals((restricted.rows, restricted.cols), (2, 0))
    val output = MutableVec.zeros(2)
    output(0) = Double.NaN
    output(1) = Double.PositiveInfinity
    restricted.applyTo(Vec.zeros(0), output)
    assertEquals(output.asVec.toSeq, Seq(0.0, 0.0))
    assertEquals(restricted.adjoint(Vec(1.0, 2.0)).length, 0)
    assertEquals(calls, 0)
    intercept[LinAlgError.VectorLengthMismatch] { restricted.applyTo(Vec.zeros(0), MutableVec.zeros(1)) }
    assertEquals(calls, 0)
  }

  test("empty restrictions compose and batch consistently with empty dense matrices") {
    val source = Matrix(2, 3)(1, 2, 3, 4, 5, 6)
    val rows = source.restrictRows(Vector.empty).orThrow
    val columns = source.restrictColumns(Vector.empty).orThrow
    val forwardRows = rows.applyTo(Matrix.eye(3)).orThrow
    assertEquals((forwardRows.rows, forwardRows.cols), (0, 3))
    val forwardColumns = columns.applyTo(Matrix.zeros(0, 4)).orThrow
    assertEquals((forwardColumns.rows, forwardColumns.cols), (2, 4))
    assertEquals(forwardColumns.valuesRowMajor, Vector.fill(8)(0.0))
    val composed = rows.andThen(columns).orThrow
    assertEquals(composed(Vec(1, 2, 3)).toSeq, Seq(0.0, 0.0))
    assertEquals(source.restrictRows(Vector(-1)), Left(LinAlgError.IndexOutOfBounds(-1, 2)))
    assert(source.restrictColumns(Vector(1, 1)).isLeft)
  }
