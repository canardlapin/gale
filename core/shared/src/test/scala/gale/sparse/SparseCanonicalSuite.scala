package gale.sparse

import gale.TestAccess
import gale.linalg.*

/** Full canonicalization path (sort + sum duplicates + prune zeros) exercised on
  * deliberately non-canonical CSR/CSC input built through [[TestAccess]], which
  * the public API never produces. Correctness is checked through the
  * sum-semantics matvec and an explicit expected dense, because a
  * duplicate-bearing matrix's `toDense` is last-write-wins, not the sum.
  */
class SparseCanonicalSuite extends munit.FunSuite:
  // Row 0 (cols 4): unsorted, a duplicate column (3), and an explicit zero (col 1).
  // Row 1: two entries in order.
  private def messyCsr: CSR =
    TestAccess.csr(
      rows = 2,
      cols = 4,
      rowPtr = Seq(0, 4, 6),
      colIdx = Seq(3, 0, 3, 1, 2, 0),
      values = Seq(5.0, 2.0, -3.0, 0.0, 4.0, 1.0)
    )

  // Sum-semantics dense: row0 col0=2, col3=5+(-3)=2; row1 col0=1, col2=4.
  private val expectedDense = Seq(2.0, 0.0, 0.0, 2.0, 1.0, 0.0, 4.0, 0.0)

  private def assertDense(actual: DMat, expected: Seq[Double]): Unit =
    assertEquals(actual.valuesRowMajor, expected)

  test("a non-canonical CSR reports hasCanonicalFormat == false") {
    assert(!messyCsr.hasCanonicalFormat)
  }

  test("CSR.canonicalize sorts, sums duplicates, and prunes zeros") {
    val c = messyCsr.canonicalize
    assert(c.hasCanonicalFormat, "canonicalize must yield canonical format")
    assertDense(c.toDense(), expectedDense)
    // Structure: row 0 -> {0: 2, 3: 2}; row 1 -> {0: 1, 2: 4}. No zero at col 1.
    assertEquals(TestAccess.rowPtr(c), Seq(0, 2, 4))
    assertEquals(TestAccess.colIdx(c), Seq(0, 3, 0, 2))
  }

  test("CSR.canonicalize is idempotent (dense and structure)") {
    val once = messyCsr.canonicalize
    val twice = once.canonicalize
    assertEquals(twice.toDense().valuesRowMajor, once.toDense().valuesRowMajor)
    assertEquals(TestAccess.rowPtr(twice), TestAccess.rowPtr(once))
    assertEquals(TestAccess.colIdx(twice), TestAccess.colIdx(once))
  }

  test("canonicalize preserves the sum-semantics matvec") {
    val x = Vec(1.0, 2.0, 3.0, 4.0)
    val before = messyCsr * x
    val after = messyCsr.canonicalize * x
    assertEquals(after.toSeq, before.toSeq)
  }

  test("CSR.sortedIndices orders columns within each row without summing or dropping") {
    val s = messyCsr.sortedIndices
    assertEquals(TestAccess.rowPtr(s), Seq(0, 4, 6))
    assertEquals(TestAccess.colIdx(s), Seq(0, 1, 3, 3, 0, 2))
    // Values follow their columns; the duplicate 3s and the zero survive.
    assertEquals(TestAccess.csrValues(s).toArray.toSeq, Seq(2.0, 0.0, 5.0, -3.0, 1.0, 4.0))
  }

  test("CSR.pruneZeros drops only exact zeros; prune(absBelow) drops small magnitudes") {
    val pruned = messyCsr.pruneZeros
    // The explicit zero at (0, 1) is gone; everything else (incl. dup 3s) stays.
    assertEquals(TestAccess.colIdx(pruned), Seq(3, 0, 3, 2, 0))
    val small = TestAccess.csr(1, 3, Seq(0, 3), Seq(0, 1, 2), Seq(1.0, 0.05, -2.0))
    assertEquals(TestAccess.colIdx(small.prune(absBelow = 0.1)), Seq(0, 2))
  }

  test("CSC canonicalization mirrors CSR via the zero-copy transpose") {
    // Column 0 (rows): unsorted, duplicate row 2, explicit zero row 1.
    val messyCsc =
      TestAccess.csc(
        rows = 4,
        cols = 2,
        colPtr = Seq(0, 4, 6),
        rowIdx = Seq(3, 0, 3, 1, 2, 0),
        values = Seq(5.0, 2.0, -3.0, 0.0, 4.0, 1.0)
      )
    assert(!messyCsc.hasCanonicalFormat)
    val c = messyCsc.canonicalize
    assert(c.hasCanonicalFormat)
    // Same numbers as the CSR case, transposed: expectedDense as a 4x2 column layout.
    val x = Vec(1.0, 2.0)
    assertEquals((c * x).toSeq, (messyCsc * x).toSeq)
  }
