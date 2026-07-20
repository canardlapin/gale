package gale.sparse

import gale.TestAccess
import gale.linalg.*

class SparseWorkspaceSuite extends munit.FunSuite:
  private def messy: CSR =
    TestAccess.csr(
      rows = 3,
      cols = 5,
      rowPtr = Seq(0, 5, 8, 10),
      colIdx = Seq(3, 1, 3, 2, 2, 4, 0, 4, 2, 2),
      values = Seq(2.0, 1.0, -2.0, 5.0, -5.0, 7.0, 3.0, 1.0, 4.0, -1.0)
    )

  test("CSR reports exact Double and index canonicalization scratch") {
    val requirement = messy.canonicalizeScratchRequirement.toOption.get
    assertEquals(requirement.doubleElements, messy.nnz)
    assertEquals(requirement.indexElements, messy.nnz)

    val empty = TestAccess.csr(2, 4, Seq(0, 0, 0), Seq.empty, Seq.empty)
    assertEquals(empty.canonicalizeScratchRequirement, Right(ScratchRequirement.empty))
  }

  test("canonicalizeWith matches the allocating facade, reuses both regions, and owns results") {
    val input = messy
    val expected = input.canonicalize
    val workspace = DenseWorkspace.empty

    val first = input.canonicalizeWith(workspace)
    assertEquals(first.toDense().valuesRowMajor, expected.toDense().valuesRowMajor)
    assertEquals(TestAccess.rowPtr(first), TestAccess.rowPtr(expected))
    assertEquals(TestAccess.colIdx(first), TestAccess.colIdx(expected))
    assert(first.hasCanonicalFormat)
    assertEquals(workspace.doubleCapacity, input.nnz)
    assertEquals(workspace.indexCapacity, input.nnz)
    val values = first.toDense().valuesRowMajor
    val doubleBacking = TestAccess.workBacking(workspace)
    val indexBacking = TestAccess.indexBacking(workspace)

    val second = input.canonicalizeWith(workspace)
    assert(TestAccess.sameStorage(doubleBacking, TestAccess.workBacking(workspace)))
    assert(TestAccess.sameIndexStorage(indexBacking, TestAccess.indexBacking(workspace)))
    assertEquals(second.toDense().valuesRowMajor, values)
    assertEquals(first.toDense().valuesRowMajor, values)
  }

  test("zero-nnz canonicalization does not grow an empty workspace") {
    val input = TestAccess.csr(3, 2, Seq(0, 0, 0, 0), Seq.empty, Seq.empty)
    val workspace = DenseWorkspace.empty
    val result = input.canonicalizeWith(workspace)
    assertEquals(result.nnz, 0)
    assert(result.hasCanonicalFormat)
    assertEquals(workspace.doubleCapacity, 0)
    assertEquals(workspace.indexCapacity, 0)
  }
