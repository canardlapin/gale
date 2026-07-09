package gale.sparse

import gale.TestAccess
import gale.linalg.*

class SparseViewSuite extends munit.FunSuite:
  // Item 2: sparse apply/mul kernels must honour the destination's offset and
  // stride, and must never touch cells outside the strided view. Each
  // destination below is a stride-2 slice starting at offset 1 of a larger
  // buffer whose other cells are sentinels that must survive untouched.
  private val Sentinel = -999.0

  private def stridedDest(outLen: Int): MutableDVec =
    val backing = TestAccess.filled(2 * outLen + 1, Sentinel)
    TestAccess.mutableVec(backing, offset = 1, length = outLen, stride = 2)

  private def assertStrided(dest: MutableDVec, expected: Seq[Double]): Unit =
    assertEquals(dest.asVec.toSeq, expected)
    val backing = TestAccess.backingSnapshot(dest)
    var i = 0
    while i < backing.length do
      if i % 2 == 0 then
        assertEquals(backing(i), Sentinel, s"neighbour cell $i was overwritten")
      i += 1

  private val csr =
    Sparse
      .coo(2, 3)
      .add(0, 0, 1.0)
      .add(0, 2, 2.0)
      .add(1, 1, 3.0)
      .toCSR() // dense [[1,0,2],[0,3,0]]

  test("CSR.mulInto writes into a strided destination, neighbours untouched") {
    val x = Vec(1.0, 2.0, 3.0)
    val dest = stridedDest(2)
    csr.mulInto(x, dest)
    val denseExpected = (csr.toDense() * x).toSeq
    assertEquals(denseExpected, Seq(7.0, 6.0))
    assertStrided(dest,denseExpected)
  }

  test("CSR.tMulInto writes into a strided destination, neighbours untouched") {
    val x = Vec(10.0, 20.0)
    val dest = stridedDest(3)
    csr.tMulInto(x, dest)
    val denseExpected = (csr.toDense().t * x).toSeq
    assertEquals(denseExpected, Seq(10.0, 60.0, 20.0))
    assertStrided(dest,denseExpected)
  }

  test("Diagonal.applyTo writes into a strided destination, neighbours untouched") {
    val d = Sparse.diagonal(2.0, 3.0, 4.0)
    val x = Vec(1.0, 1.0, 1.0)
    val dest = stridedDest(3)
    d.applyTo(x, dest)
    assertStrided(dest,Seq(2.0, 3.0, 4.0))
  }

  test("Identity.applyTo writes into a strided destination, neighbours untouched") {
    val id = Sparse.identity(3)
    val x = Vec(5.0, 6.0, 7.0)
    val dest = stridedDest(3)
    id.applyTo(x, dest)
    assertStrided(dest,Seq(5.0, 6.0, 7.0))
  }

  test("Permutation.applyTo writes into a strided destination, neighbours untouched") {
    val p = Sparse.permutation(2, 0, 1)
    val x = Vec(10.0, 20.0, 30.0)
    val dest = stridedDest(3)
    p.applyTo(x, dest)
    // into(row) = x(columnsByRow(row)); columnsByRow = [2, 0, 1]
    assertStrided(dest,Seq(30.0, 10.0, 20.0))
  }
