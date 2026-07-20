package gale.sparse

import gale.TestAccess
import gale.linalg.*

/** Equivalence tests for the P5 sparse rewrites: single-pass CSR column
  * extraction (K2b), the two-pointer `zipValues` / in-place `mapValues` (K3),
  * direct CSC kernels, zero-copy transpose, and COO scatter (K4).
  */
class SparsePerfSuite extends munit.FunSuite:
  private def csrOf(rows: Int, cols: Int, entries: (Int, Int, Double)*): CSR =
    val builder = Sparse.coo(rows, cols)
    entries.foreach { case (r, c, v) => builder.add(r, c, v) }
    builder.toCSR()

  // dense [[1,0,2,0],[0,3,0,0],[4,0,5,6],[0,0,0,7]]
  private def sample: CSR =
    csrOf(4, 4, (0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0), (2, 0, 4.0), (2, 2, 5.0), (2, 3, 6.0), (3, 3, 7.0))

  private def assertVecApprox(actual: Seq[Double], expected: Seq[Double], tol: Double = 1e-12): Unit =
    assertEquals(actual.length, expected.length)
    actual.zip(expected).zipWithIndex.foreach { case ((a, e), i) =>
      assert(math.abs(a - e) <= tol, s"index $i: $a != $e")
    }

  // ---- K2b: CSR.col single pass -------------------------------------------

  test("CSR.col matches the dense column for every index") {
    val a = sample
    val dense = a.toDense()
    var j = 0
    while j < a.cols do
      assertEquals(a.col(j).toSeq, dense.col(j).toSeq, s"column $j")
      j += 1
  }

  // ---- K3: zipValues / mapValues ------------------------------------------

  private def denseZip(a: CSR, b: CSR, f: (Double, Double) => Double): Seq[Double] =
    val da = a.toDense()
    val db = b.toDense()
    Matrix.tabulate(a.rows, a.cols)((i, j) => f(da(i, j), db(i, j))).valuesRowMajor

  test("CSR.zipValues (two-pointer merge) matches dense elementwise for +, -, *") {
    val a = csrOf(3, 4, (0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0), (2, 3, 4.0))
    val b = csrOf(3, 4, (0, 2, 5.0), (0, 3, 6.0), (1, 1, -3.0), (2, 0, 7.0))

    assertEquals(a.zipValues(b)(_ + _).toDense().valuesRowMajor, denseZip(a, b, _ + _))
    assertEquals(a.zipValues(b)(_ - _).toDense().valuesRowMajor, denseZip(a, b, _ - _))
    assertEquals(a.zipValues(b)(_ * _).toDense().valuesRowMajor, denseZip(a, b, _ * _))

    // (1,1): 3 + (-3) == 0 must be pruned, not stored.
    val sum = a.zipValues(b)(_ + _)
    assert(!TestAccess.colIdx(sum).isEmpty)
    assertEquals(sum(1, 1), 0.0)
    // Result stays canonical: recovering it via toCSC/toCSR round-trip agrees.
    assertEquals(sum.toDense().valuesRowMajor, sum.toCSC.toCSR.toDense().valuesRowMajor)
  }

  test("CSR.mapValues preserves structure and explicit zeros without re-sorting") {
    val a = sample
    val dense = a.toDense()

    // f(0) == 0: mapping equals mapping the dense matrix.
    assertEquals(
      a.mapValues(_ * 2.0).toDense().valuesRowMajor,
      Matrix.tabulate(a.rows, a.cols)((i, j) => dense(i, j) * 2.0).valuesRowMajor
    )

    // f(0) != 0: only stored entries are mapped; structural zeros stay zero and
    // the nonzero count is unchanged.
    val shifted = a.mapValues(_ + 10.0)
    assertEquals(shifted.nnz, a.nnz)
    assertEquals(
      shifted.toDense().valuesRowMajor,
      Matrix.tabulate(a.rows, a.cols)((i, j) => if dense(i, j) != 0.0 then dense(i, j) + 10.0 else 0.0).valuesRowMajor
    )

    // Numeric updates retain the pattern even when every stored value becomes
    // zero; only the explicitly structure-changing prune removes those entries.
    val zeroed = a.mapValues(_ * 0.0)
    assertEquals(zeroed.nnz, a.nnz)
    assertEquals(zeroed.pattern, a.pattern)
    assertEquals(zeroed.pruneZeros.nnz, 0)
  }

  // ---- K4: direct CSC kernels ---------------------------------------------

  test("CSC.mulInto and tMulInto match the dense products") {
    val csr = sample
    val csc = csr.toCSC
    val dense = csr.toDense()
    val x = Vec(1.0, 2.0, 3.0, 4.0)

    val y = MutableVec.zeros(csc.rows)
    csc.mulInto(x, y)
    assertVecApprox(y.asVec.toSeq, (dense * x).toSeq)
    assertVecApprox((csc * x).toSeq, (dense * x).toSeq)

    val yt = MutableVec.zeros(csc.cols)
    csc.tMulInto(x, yt)
    assertVecApprox(yt.asVec.toSeq, (dense.t * x).toSeq)
  }

  test("CSC.mulInto and tMulInto reject aliased and mismatched arguments") {
    val csc = sample.toCSC
    val y = MutableVec.zeros(4)
    intercept[LinAlgError.UnsupportedOperation] {
      csc.mulInto(y.asVec, y)
    }
    intercept[LinAlgError.DimensionMismatch] {
      csc.mulInto(Vec(1.0, 2.0), MutableVec.zeros(4))
    }
    intercept[LinAlgError.DimensionMismatch] {
      csc.tMulInto(Vec(1.0, 2.0), MutableVec.zeros(4))
    }
  }

  // ---- K4: zero-copy transpose --------------------------------------------

  test("CSR.t is a zero-copy CSC that transposes correctly") {
    val csr = csrOf(2, 3, (0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0))
    val t = csr.t
    assert(t.isInstanceOf[CSC])
    assertEquals(t.toCSR.toDense().valuesRowMajor, csr.toDense().t.valuesRowMajor)
    // Values array is shared, not copied.
    assert(
      TestAccess.sameStorage(TestAccess.csrValues(csr), TestAccess.cscValues(t)),
      "CSR.t copied the value array instead of sharing it"
    )
  }

  test("CSC.t is a zero-copy CSR that transposes correctly") {
    val csc = csrOf(3, 2, (0, 0, 1.0), (1, 1, 2.0), (2, 0, 3.0)).toCSC
    val t = csc.t
    assert(t.isInstanceOf[CSR])
    assertEquals(t.toDense().valuesRowMajor, csc.toCSR.toDense().t.valuesRowMajor)
    assert(
      TestAccess.sameStorage(TestAccess.cscValues(csc), TestAccess.csrValues(t)),
      "CSC.t copied the value array instead of sharing it"
    )
  }

  // ---- K4: COO direct scatter ---------------------------------------------

  test("COO.applyTo and transposeApplyTo scatter directly, summing duplicates") {
    // Includes a duplicate (0,0) entry to prove the scatter sums, matching dense.
    val coo =
      Sparse
        .coo(3, 4)
        .add(0, 0, 1.0)
        .add(0, 0, 2.0) // duplicate -> sums to 3 at (0,0)
        .add(0, 2, -1.0)
        .add(1, 3, 4.0)
        .add(2, 1, 5.0)
        .toCOO()
    val dense = coo.toCSR.toDense()

    val x = Vec(1.0, 2.0, 3.0, 4.0)
    val y = MutableVec.zeros(3)
    coo.applyTo(x, y)
    assertVecApprox(y.asVec.toSeq, (dense * x).toSeq)

    val xt = Vec(1.0, 2.0, 3.0)
    val yt = MutableVec.zeros(4)
    coo.transposeApplyTo(xt, yt)
    assertVecApprox(yt.asVec.toSeq, (dense.t * xt).toSeq)
  }
