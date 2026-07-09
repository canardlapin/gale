package gale.sparse

import gale.TestAccess
import gale.linalg.*

class CanonicalizationSuite extends munit.FunSuite:
  // Item 6: a canonical COO stores entries in row-major sorted order. Transposing
  // swaps (row, col), which turns that into column-major order, so the transpose
  // is generally NOT canonical. The flag must not lie about it.
  private val a =
    Sparse
      .coo(3, 3)
      .add(0, 0, 1.0)
      .add(0, 2, 2.0)
      .add(1, 1, 3.0)
      .add(2, 0, 4.0)
      .add(2, 2, 5.0)
      .toCOO() // canonical, row-major sorted

  test("a builder-produced COO is canonical") {
    assert(a.hasCanonicalFormat)
  }

  test("COO.t does not falsely advertise canonical format") {
    assert(!a.t.hasCanonicalFormat, "transpose falsely claims canonical format")
  }

  test("CSR built from the transpose still has sorted columns within each row") {
    val csr = a.t.toCSR
    val rowPtr = TestAccess.rowPtr(csr)
    val colIdx = TestAccess.colIdx(csr)
    var row = 0
    while row < csr.rows do
      var p = rowPtr(row)
      val end = rowPtr(row + 1)
      var prev = -1
      while p < end do
        val c = colIdx(p)
        assert(c > prev, s"row $row column indices not strictly increasing at p=$p")
        prev = c
        p += 1
      row += 1
  }

  test("CSR built from the transpose matches the dense transpose") {
    val actual = a.t.toCSR.toDense().valuesRowMajor
    val expected = a.toCSR.toDense().t.valuesRowMajor
    assertEquals(actual, expected)
    assertEquals(actual, Seq(1.0, 0.0, 4.0, 0.0, 3.0, 0.0, 2.0, 0.0, 5.0))
  }
