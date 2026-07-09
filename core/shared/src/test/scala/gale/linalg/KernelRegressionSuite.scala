package gale.linalg

import gale.TestAccess
import gale.sparse.Sparse

class KernelRegressionSuite extends munit.FunSuite:
  // Item 1: with beta == 0 the destination must be assigned, not read-and-scaled.
  // A NaN already sitting in the output buffer would poison the result via
  // 0.0 * NaN == NaN if the kernel reads y before overwriting it.
  test("mulInto with beta==0 ignores NaN already in the destination (all layouts)") {
    // A == [[1, 2], [3, 4]] realised three different storage ways; A * [1, 1] == [3, 7].
    val expected = Seq(3.0, 7.0)
    val x = Vec(1.0, 1.0)

    // (a) contiguous row-major -> dgemvRowMajor
    val rowMajor = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)

    // (b) column-major view -> dgemvColMajor (A.t of a row-major matrix)
    val colMajor = Matrix.dense(2, 2)(1.0, 3.0, 2.0, 4.0).t

    // (c) strided submatrix view -> dgemv (rowStride != 1 and colStride != 1)
    val backing = TestAccess.doubleArray(1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0)
    val strided = TestAccess.mat(backing, offset = 0, rows = 2, cols = 2, rowStride = 4, colStride = 2)

    for (label, a) <- Seq("row-major" -> rowMajor, "col-major" -> colMajor, "strided" -> strided) do
      val y = MutableDVec.zeros(2)
      y(0) = Double.NaN
      y(1) = Double.NaN
      a.mulInto(x, y)
      assert(y(0).isFinite && y(1).isFinite, s"$label produced non-finite output: ${y.asVec.toSeq}")
      assertEquals(y.asVec.toSeq, expected, s"$label mismatch")
  }

  // Item 1 continued: dgemm (matrix-matrix) has the same beta==0 read-not-assign hazard.
  test("dgemm with beta==0 ignores NaN already in the destination") {
    val a = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
    val b = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val c = DMat.zeros(2, 2)
    TestAccess.poisonWithNaN(c)
    TestAccess.gemm(a, b, c)
    assertEquals(c.valuesRowMajor, Seq(1.0, 2.0, 3.0, 4.0))
  }

  // Item 3: writing into a destination that shares storage with the source
  // vector corrupts the computation mid-flight. The kernels must reject it
  // rather than silently produce garbage.
  test("DMat.mulInto rejects a destination aliasing the source vector") {
    val a = Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0,
      7.0, 8.0, 9.0
    )
    val y = MutableDVec.zeros(3)
    y(0) = 1.0
    y(1) = 2.0
    y(2) = 3.0
    val aliased = y.asVec // shares y's backing storage
    intercept[LinAlgError.UnsupportedOperation] {
      a.mulInto(aliased, y)
    }
  }

  test("CSR.mulInto and tMulInto reject a destination aliasing the source vector") {
    val a =
      Sparse
        .coo(3, 3)
        .add(0, 0, 1.0)
        .add(1, 1, 2.0)
        .add(2, 0, 3.0)
        .add(2, 2, 4.0)
        .toCSR()
    val y = MutableDVec.zeros(3)
    y(0) = 1.0
    y(1) = 2.0
    y(2) = 3.0
    intercept[LinAlgError.UnsupportedOperation] {
      a.mulInto(y.asVec, y)
    }
    val z = MutableDVec.zeros(3)
    z(0) = 1.0
    z(1) = 2.0
    z(2) = 3.0
    intercept[LinAlgError.UnsupportedOperation] {
      a.tMulInto(z.asVec, z)
    }
  }
