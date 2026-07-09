package gale.linalg

import gale.TestAccess

class WorkspaceSuite extends munit.FunSuite:
  test("empty workspace starts with no scratch and grows on demand") {
    val workspace = DenseWorkspace.empty
    assertEquals(workspace.workCapacity, 0)

    workspace.work(5)
    assertEquals(workspace.workCapacity, 5)

    // A smaller request keeps the larger buffer.
    workspace.work(3)
    assertEquals(workspace.workCapacity, 5)
  }

  test("qrWith uses the length-m reflector scratch and reuses it without aliasing results") {
    val A = Matrix.dense(3, 3)(
      12.0, -51.0, 4.0,
      6.0, 167.0, -68.0,
      -4.0, 24.0, -41.0
    )
    val workspace = DenseWorkspace.empty

    val first = A.qrWith(workspace)
    // The only genuine QR scratch is the length-rows Householder reflector, so
    // the previously-empty workspace grew to exactly 3 (not m*n or a tau buffer).
    assertEquals(workspace.workCapacity, 3)
    val backingAfterFirst = TestAccess.workBacking(workspace)
    val q1 = first.q.valuesRowMajor
    val r1 = first.r.valuesRowMajor

    val second = A.qrWith(workspace)
    val backingAfterSecond = TestAccess.workBacking(workspace)

    // Genuine reuse: the second call reused the same backing scratch array
    // rather than allocating a fresh one.
    assert(
      TestAccess.sameStorage(backingAfterFirst, backingAfterSecond),
      "qrWith reallocated scratch instead of reusing the workspace buffer"
    )
    // No scratch was aliased into the results: the first factorization is
    // byte-for-byte unchanged after the second call ran through the same buffer.
    assertEquals(first.q.valuesRowMajor, q1)
    assertEquals(first.r.valuesRowMajor, r1)
    assertEquals(second.diagnostics.rank, Some(3))
  }
