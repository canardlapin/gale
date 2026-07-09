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
