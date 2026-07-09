package gale.linalg

import gale.TestAccess

class DenseCoreSuite extends munit.FunSuite:
  test("Vec.apply returns DVec with dot, norm, arithmetic, and slices") {
    val x = Vec(1.0, 2.0, 3.0)
    val y = Vec(4.0, 5.0, 6.0)

    assert(x.isInstanceOf[DVec])
    assertEquals(x.length, 3)
    assertEquals(x.dot(y), 32.0)
    assert(math.abs(x.norm2 - math.sqrt(14.0)) < 1e-12)
    assertEquals((x + y).toSeq, Seq(5.0, 7.0, 9.0))
    assertEquals((y - x).toSeq, Seq(3.0, 3.0, 3.0))
    assertEquals((2.0 * x).toSeq, Seq(2.0, 4.0, 6.0))
    assertEquals(x.slice(1, 3).toSeq, Seq(2.0, 3.0))
  }

  test("Matrix.dense returns row-major DMat with row, column, transpose, and updates") {
    val A = Matrix.dense(2, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0
    )

    assert(A.isInstanceOf[DMat])
    assertEquals(A.rows, 2)
    assertEquals(A.cols, 3)
    assertEquals(A(1, 2), 6.0)
    assertEquals(A.row(1).toSeq, Seq(4.0, 5.0, 6.0))
    assertEquals(A.col(1).toSeq, Seq(2.0, 5.0))
    assertEquals(A.t.rows, 3)
    assertEquals(A.t.cols, 2)
    assertEquals(A.t(0, 1), 4.0)

    val B = A.updated(0, 1, 9.0)
    assertEquals(A(0, 1), 2.0)
    assertEquals(B(0, 1), 9.0)
  }

  test("matrix-vector and matrix-matrix products use dense kernels") {
    val A = Matrix.dense(2, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0
    )
    val x = Vec(1.0, 2.0, 1.0)

    assertEquals((A * x).toSeq, Seq(8.0, 20.0))

    val B = Matrix.dense(3, 2)(
      1.0, 2.0,
      0.0, 1.0,
      1.0, 0.0
    )
    val C = A * B
    assertEquals(C.rows, 2)
    assertEquals(C.cols, 2)
    assertEquals(C.valuesRowMajor, Seq(4.0, 4.0, 10.0, 13.0))
  }

  test("mulInto handles row-major and transpose-view layouts") {
    val A = Matrix.dense(2, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0
    )

    val rowMajorBuffer = MutableVec.from(Vec.fill(4)(-7.0))
    val rowMajorOut = TestAccess.mutableViewInto(rowMajorBuffer, offset = 1, length = 2, stride = 1)
    A.mulInto(Vec(1.0, 2.0, 1.0), rowMajorOut)
    assertEquals(rowMajorBuffer.asVec.toSeq, Seq(-7.0, 8.0, 20.0, -7.0))

    val transposeBuffer = MutableVec.from(Vec.fill(5)(-7.0))
    val transposeOut = TestAccess.mutableViewInto(transposeBuffer, offset = 1, length = 3, stride = 1)
    A.t.mulInto(Vec(2.0, -1.0), transposeOut)
    assertEquals(transposeBuffer.asVec.toSeq, Seq(-7.0, -2.0, -1.0, 0.0, -7.0))
  }

  test("shape errors are reported before kernel entry") {
    val A = Matrix.eye(3)
    val x = Vec(1.0, 2.0)

    intercept[LinAlgError.VectorLengthMismatch] {
      A * x
    }

    intercept[LinAlgError.VectorLengthMismatch] {
      Vec(1.0, 2.0, 3.0).dot(x)
    }
  }
