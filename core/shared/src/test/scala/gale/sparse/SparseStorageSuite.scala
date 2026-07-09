package gale.sparse

import gale.linalg.*

class SparseStorageSuite extends munit.FunSuite:
  test("CSR and CSC expose element, row, column, transpose, trace, and density") {
    val csr =
      Sparse
        .coo(3, 3)
        .add(0, 0, 1.0)
        .add(1, 1, 2.0)
        .add(2, 0, 3.0)
        .add(2, 2, 4.0)
        .toCSR()
    val csc = csr.toCSC

    assertEquals(csr.nnz, 4)
    assertEquals(csr.density, 4.0 / 9.0)
    assertEquals(csr(2, 0), 3.0)
    assertEquals(csr.row(2).toSeq, Seq(3.0, 0.0, 4.0))
    assertEquals(csr.col(0).toSeq, Seq(1.0, 0.0, 3.0))
    assertEquals(csr.trace, 7.0)
    assertEquals(csr.t.toCSR.toDense().valuesRowMajor, csr.toDense().t.valuesRowMajor)
    assertEquals((csc * Vec(1.0, 2.0, 3.0)).toSeq, (csr * Vec(1.0, 2.0, 3.0)).toSeq)
  }

  test("structured sparse matrices multiply directly") {
    val diagonal = Sparse.diagonal(2.0, 0.0, -1.0)
    val identity = Sparse.identity(3)
    val zero = Sparse.zero(2, 3)
    val permutation = Sparse.permutation(2, 0, 1)

    assertEquals((diagonal * Vec(3.0, 4.0, 5.0)).toSeq, Seq(6.0, 0.0, -5.0))
    assertEquals((identity * Vec(3.0, 4.0, 5.0)).toSeq, Seq(3.0, 4.0, 5.0))
    assertEquals((zero * Vec(3.0, 4.0, 5.0)).toSeq, Seq(0.0, 0.0))
    assertEquals((permutation * Vec(10.0, 20.0, 30.0)).toSeq, Seq(30.0, 10.0, 20.0))
    assertEquals((permutation.t * Vec(10.0, 20.0, 30.0)).toSeq, Seq(20.0, 30.0, 10.0))
    assertEquals(diagonal.nnz, 2)
    assertEquals(zero.nnz, 0)
  }
