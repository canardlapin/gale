package gale.parity

import gale.linalg.*
import gale.parity.NumpyScipyFixtures.*
import gale.parity.NumpyScipySupport.*
import gale.sparse.direct.*

/** Sparse-direct factorization: lock the empty-provider contract and the
  * SciPy SuperLU mathematical target.
  *
  * Gale ships no sparse LU / Cholesky / QR. Breeze's SuiteSparse path is not
  * an honest public reference for this tree. SciPy `splu` (and R
  * `Matrix::lu`) is the intended differential target when a JVM provider
  * lands. Until then:
  *
  *   - `SparseDirectProvider.none` advertises no capability and cannot create
  *     a workspace;
  *   - Gale's dense LU on the same numbers must already match SciPy SuperLU,
  *     so a future provider has a pinned `Ax = b` answer.
  */
class SparseDirectParitySuite extends munit.FunSuite:

  test("default sparse-direct provider has no factorization capability") {
    assertEquals(SparseDirectProvider.none.name, "none")
    assertEquals(SparseDirectProvider.none.capabilities, Set.empty)
    assertEquals(SparseDirect.capabilities, Set.empty)
    SparseDirect.newWorkspace() match
      case Left(LinAlgError.UnsupportedOperation(op)) =>
        assert(op.contains("sparse direct"), s"unexpected message: $op")
      case other =>
        fail(s"expected UnsupportedOperation, got $other")
  }

  test("SciPy SuperLU solution matches Gale dense LU on the same sparse pattern") {
    for ref <- sparseLu do
      val gx = galeMatrix(ref.a).solve(galeVector(ref.b)).orThrow
      assertVecClose(gx, ref.x, 1e-10, s"${ref.name} dense LU vs scipy.sparse.linalg.splu")
  }
