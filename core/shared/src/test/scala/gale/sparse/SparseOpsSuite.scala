package gale.sparse

import gale.linalg.*

class SparseOpsSuite extends munit.FunSuite:
  test("CSR supports add, subtract, scale, mapValues, and zipValues") {
    val A =
      Sparse
        .coo(2, 3)
        .add(0, 0, 1.0)
        .add(0, 2, 2.0)
        .add(1, 1, 3.0)
        .toCSR()
    val B =
      Sparse
        .coo(2, 3)
        .add(0, 2, 5.0)
        .add(1, 0, 7.0)
        .toCSR()

    assertEquals((A + B).toDense().valuesRowMajor, Seq(1.0, 0.0, 7.0, 7.0, 3.0, 0.0))
    assertEquals((A - B).toDense().valuesRowMajor, Seq(1.0, 0.0, -3.0, -7.0, 3.0, 0.0))
    assertEquals((A * 2.0).toDense().valuesRowMajor, Seq(2.0, 0.0, 4.0, 0.0, 6.0, 0.0))
    assertEquals(A.mapValues(_ + 1.0).toDense().valuesRowMajor, Seq(2.0, 0.0, 3.0, 0.0, 4.0, 0.0))
    assertEquals(A.zipValues(B)(_ * _).toDense().valuesRowMajor, Seq(0.0, 0.0, 10.0, 0.0, 0.0, 0.0))
  }

  test("CSR supports transpose multiply, sparse-dense multiply, diagonal, and dense guard") {
    val A =
      Sparse
        .coo(2, 3)
        .add(0, 0, 1.0)
        .add(0, 2, 2.0)
        .add(1, 1, 3.0)
        .toCSR()
    val y = MutableVec.zeros(3)
    A.tMulInto(Vec(10.0, 20.0), y)

    val B = Matrix.dense(3, 2)(
      1.0, 2.0,
      3.0, 4.0,
      5.0, 6.0
    )

    assertEquals(y.asVec.toSeq, Seq(10.0, 60.0, 20.0))
    assertEquals((A * B).valuesRowMajor, Seq(11.0, 14.0, 9.0, 12.0))
    assertEquals(A.diagonal.toSeq, Seq(1.0, 3.0))
    intercept[LinAlgError.UnsupportedOperation] {
      A.toDense(maxEntries = 3)
    }
  }
