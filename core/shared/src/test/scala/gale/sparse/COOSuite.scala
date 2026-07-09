package gale.sparse

import gale.linalg.*

class COOSuite extends munit.FunSuite:
  test("COO builder canonicalizes sorted duplicates") {
    val builder =
      Sparse
        .coo(rows = 3, cols = 3)
        .add(2, 1, 4.0)
        .add(0, 2, 1.0)
        .add(0, 2, 3.0)
        .add(1, 1, 0.0)

    assert(!builder.sortedIndices)
    assert(!builder.hasCanonicalFormat)

    builder.canonicalize(DuplicatePolicy.Sum).pruneZeros
    val coo = builder.toCOO()

    assert(builder.sortedIndices)
    assert(builder.hasCanonicalFormat)
    assertEquals(coo.entries, Seq(COOEntry(0, 2, 4.0), COOEntry(2, 1, 4.0)))
  }

  test("COO duplicate policies support last and error") {
    val last =
      Sparse
        .coo(2, 2)
        .add(0, 0, 1.0)
        .add(0, 0, 7.0)
        .toCOO(DuplicatePolicy.Last)

    assertEquals(last(0, 0), 7.0)

    intercept[LinAlgError.InvalidArgument] {
      Sparse.coo(2, 2).add(0, 0, 1.0).add(0, 0, 2.0).toCOO(DuplicatePolicy.Error)
    }
  }

  test("COO converts to CSR and CSC") {
    val coo =
      Sparse
        .coo(2, 3)
        .add(0, 1, 2.0)
        .add(1, 0, 3.0)
        .add(1, 2, 4.0)
        .toCOO()
    val csr = coo.toCSR
    val csc = coo.toCSC

    assertEquals((csr * Vec(1.0, 2.0, 3.0)).toSeq, Seq(4.0, 15.0))
    assertEquals((csc * Vec(1.0, 2.0, 3.0)).toSeq, Seq(4.0, 15.0))
    assertEquals(csr.toDense().valuesRowMajor, Seq(0.0, 2.0, 0.0, 3.0, 0.0, 4.0))
  }
