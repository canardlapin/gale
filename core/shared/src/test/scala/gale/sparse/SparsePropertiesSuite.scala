package gale.sparse

import gale.TestAccess
import gale.linalg.*

class SparsePropertiesSuite extends munit.FunSuite:
  private def messyCsr: CSR =
    TestAccess.csr(
      rows = 2,
      cols = 4,
      rowPtr = Seq(0, 4, 6),
      colIdx = Seq(3, 0, 3, 1, 2, 0),
      values = Seq(5.0, 2.0, -3.0, 0.0, 4.0, 1.0)
    )

  test("assumeCanonicalSparse is zero-cost and preserves the concrete type") {
    val csr = messyCsr
    val assumed: CanonicalSparse[CSR] = csr.assumeCanonicalSparse
    summon[CanonicalSparse[CSR] <:< CSR]
    assert(assumed.value eq csr)
    assertEquals(assumed.rows, csr.rows)
  }

  test("verifyCanonicalSparse accepts canonical CSR, CSC, and COO") {
    val csr = messyCsr.canonicalize
    val csc = csr.toCSC
    val coo = Sparse.coo(2, 3).add(0, 0, 1.0).add(1, 2, 2.0).toCOO()

    assert(csr.verifyCanonicalSparse.isRight)
    assert(csc.verifyCanonicalSparse.isRight)
    assert(coo.verifyCanonicalSparse.isRight)
  }

  test("verifyCanonicalSparse rejects noncanonical forms") {
    val csr = messyCsr
    val csc = csr.t
    val canonicalCoo = Sparse.coo(2, 3).add(0, 2, 1.0).add(1, 0, 2.0).toCOO()
    val noncanonicalCoo = canonicalCoo.t

    assert(csr.verifyCanonicalSparse.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(csc.verifyCanonicalSparse.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(noncanonicalCoo.verifyCanonicalSparse.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
  }

  test("canonicalization enables verification without changing the represented operator") {
    val original = messyCsr
    val canonical = original.canonicalize.verifyCanonicalSparse.orThrow
    val x = Vec(1.0, 2.0, 3.0, 4.0)

    assertEquals((canonical * x).toSeq, (original * x).toSeq)
    assert(canonical.value.hasCanonicalFormat)
  }
