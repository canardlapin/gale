package gale.laws

import gale.linalg.*
import gale.sparse.*

/** Canonicalization laws over the public sparse API. Publicly-constructed
  * matrices never carry duplicate coordinates, so the reachable non-canonical
  * feature is an explicitly-stored zero; these laws exercise that and the
  * already-canonical case.
  */
class CanonicalizationLawSuite extends munit.FunSuite:
  // An explicit zero at (0, 1) survives COO -> CSR construction.
  private def csrWithExplicitZero: CSR =
    Sparse.coo(2, 3).add(0, 0, 1.0).add(0, 1, 0.0).add(1, 2, 3.0).toCSR()

  private def cscWithExplicitZero: CSC =
    Sparse.coo(2, 3).add(0, 0, 1.0).add(0, 1, 0.0).add(1, 2, 3.0).toCSC()

  test("CSR: canonicalize prunes explicit zeros and reports canonical format") {
    val a = csrWithExplicitZero
    assert(!a.hasCanonicalFormat, "an explicit zero is not canonical")
    val c = a.canonicalize
    assert(c.hasCanonicalFormat, "hasCanonicalFormat must hold after canonicalize")
    // Dense image is preserved (the zero contributes nothing).
    SparseLaws.densifiesTo(c.toDense(), a.toDense())
    // Idempotence at the dense level.
    SparseLaws.canonicalizeIsIdempotentDense(c.toDense(), c.canonicalize.toDense())
  }

  test("CSC: canonicalize prunes explicit zeros and reports canonical format") {
    val a = cscWithExplicitZero
    assert(!a.hasCanonicalFormat, "an explicit zero is not canonical")
    val c = a.canonicalize
    assert(c.hasCanonicalFormat, "hasCanonicalFormat must hold after canonicalize")
    SparseLaws.densifiesTo(c.toDense(), a.toDense())
    SparseLaws.canonicalizeIsIdempotentDense(c.toDense(), c.canonicalize.toDense())
  }

  test("an already-canonical CSR reports canonical format and canonicalize is a no-op") {
    val a = Sparse.coo(2, 2).add(0, 0, 1.0).add(1, 1, 2.0).toCSR()
    assert(a.hasCanonicalFormat)
    SparseLaws.densifiesTo(a.canonicalize.toDense(), a.toDense())
  }

  test("COO.canonicalize sums duplicates, prunes cancellations, and is idempotent at the dense level") {
    // Builder path is the only public source of duplicates, and it sums them at
    // toCOO time; COO.canonicalize additionally prunes zeros (incl. cancellations).
    val a =
      Sparse
        .coo(2, 2)
        .add(0, 0, 3.0)
        .add(0, 0, -3.0) // cancels to zero
        .add(1, 1, 4.0)
        .toCOO()
    val c = a.canonicalize
    // (0,0) cancelled to zero must be gone; only (1,1)=4 remains.
    assertEquals(c.nnz, 1)
    assertEquals(c.toCSR.toDense().valuesRowMajor, Seq(0.0, 0.0, 0.0, 4.0))
    assertEquals(c.canonicalize.toCSR.toDense().valuesRowMajor, c.toCSR.toDense().valuesRowMajor)
  }
