package gale.sparse

import gale.TestAccess
import gale.linalg.*

/** CSC arithmetic/reduction parity with CSR (the two agree because CSC reuses the
  * CSR code through the zero-copy transpose), plus the CSR/COO apply-agreement
  * contract.
  */
class CSCParitySuite extends munit.FunSuite:
  private val entriesA = Seq((0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0), (2, 0, 4.0), (2, 2, 5.0))
  private val entriesB = Seq((0, 0, 10.0), (0, 1, -1.0), (1, 1, 7.0), (2, 2, 2.0))

  private def cooOf(entries: Seq[(Int, Int, Double)]): COO =
    Sparse.fromCOO(3, 3, entries)

  private val csrA = cooOf(entriesA).toCSR
  private val csrB = cooOf(entriesB).toCSR
  private val cscA = cooOf(entriesA).toCSC
  private val cscB = cooOf(entriesB).toCSC

  private def assertSameDense(csc: CSC, csr: CSR, clue: String): Unit =
    assertEquals(csc.toDense().valuesRowMajor, csr.toDense().valuesRowMajor, clue)

  test("CSC + / - / scale match the CSR result") {
    assertSameDense(cscA + cscB, csrA + csrB, "add")
    assertSameDense(cscA - cscB, csrA - csrB, "sub")
    assertSameDense(cscA * 2.5, csrA * 2.5, "scale")
  }

  test("CSC mapValues and zipValues match the CSR result") {
    assertSameDense(cscA.mapValues(_ * -1.0), csrA.mapValues(_ * -1.0), "mapValues")
    assertSameDense(cscA.zipValues(cscB)(_ + _), csrA.zipValues(csrB)(_ + _), "zipValues")
  }

  test("CSC toDense is contiguous row-major and equals the CSR densification") {
    val dense = cscA.toDense()
    assert(dense.isContiguousRowMajor, "CSC.toDense must be contiguous row-major")
    assertEquals(dense.valuesRowMajor, csrA.toDense().valuesRowMajor)
  }

  test("CSC diagonal and trace match the CSR result") {
    assertEquals(cscA.diagonal.toSeq, csrA.diagonal.toSeq)
    assertEquals(cscA.trace, csrA.trace)
  }

  test("CSC row and col accessors match the CSR result") {
    var i = 0
    while i < 3 do
      assertEquals(cscA.row(i).toSeq, csrA.row(i).toSeq, s"row $i")
      assertEquals(cscA.col(i).toSeq, csrA.col(i).toSeq, s"col $i")
      i += 1
  }

  test("CSR.apply and COO.apply agree on a canonical matrix") {
    val coo = cooOf(entriesA)
    assert(coo.hasCanonicalFormat)
    var i = 0
    while i < 3 do
      var j = 0
      while j < 3 do
        assertEquals(csrA(i, j), coo(i, j), s"apply($i,$j)")
        j += 1
      i += 1
  }

  test("COO.apply sums duplicate coordinates in non-canonical form") {
    // Two entries at (0, 0): a non-canonical COO the public builder never yields.
    val coo = TestAccess.coo(2, 2, Seq((0, 0, 1.0), (0, 0, 2.0), (1, 1, 4.0)), canonical = false)
    assertEquals(coo(0, 0), 3.0)
    assertEquals(coo(1, 1), 4.0)
    // The canonical CSR of the same matrix reports the summed value once.
    assertEquals(coo.canonicalize.toCSR(0, 0), 3.0)
  }
