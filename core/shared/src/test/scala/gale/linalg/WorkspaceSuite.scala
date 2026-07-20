package gale.linalg

import gale.TestAccess

class WorkspaceSuite extends munit.FunSuite:
  private def requirement(doubles: Long, indices: Long): ScratchRequirement =
    ScratchRequirement.checked(doubles, indices).toOption.get

  test("scratch requirements are checked at zero, negative, and address-space boundaries") {
    assertEquals(ScratchRequirement.checked(0L, 0L), Right(ScratchRequirement.empty))
    assert(ScratchRequirement.checked(-1L, 0L).isLeft)
    assert(ScratchRequirement.checked(0L, -1L).isLeft)
    assert(ScratchRequirement.checked(Int.MaxValue.toLong + 1L, 0L).isLeft)
    assert(ScratchRequirement.checked(0L, Int.MaxValue.toLong + 1L).isLeft)

    val maximum = requirement(Int.MaxValue.toLong, Int.MaxValue.toLong)
    assert(maximum.simultaneous(requirement(1L, 0L)).isLeft)
    assert(maximum.simultaneous(requirement(0L, 1L)).isLeft)
  }

  test("simultaneous sums, alternative takes maxima, and both compositions obey their laws") {
    val a = requirement(2L, 7L)
    val b = requirement(5L, 3L)
    val c = requirement(11L, 13L)

    assertEquals(a.simultaneous(b), Right(requirement(7L, 10L)))
    assertEquals(a.alternative(b), requirement(5L, 7L))
    assertEquals(a.simultaneous(b), b.simultaneous(a))
    assertEquals(a.alternative(b), b.alternative(a))
    assertEquals(
      a.simultaneous(b).flatMap(_.simultaneous(c)),
      b.simultaneous(c).flatMap(a.simultaneous)
    )
    assertEquals(a.alternative(b).alternative(c), a.alternative(b.alternative(c)))
    assertEquals(a.simultaneous(ScratchRequirement.empty), Right(a))
    assertEquals(a.alternative(ScratchRequirement.empty), a)
  }

  test("empty workspace starts with no scratch and grows on demand") {
    val workspace = DenseWorkspace.empty
    assertEquals(workspace.workCapacity, 0)
    assertEquals(workspace.indexCapacity, 0)

    workspace.work(5)
    assertEquals(workspace.workCapacity, 5)

    // A smaller request keeps the larger buffer.
    workspace.work(3)
    assertEquals(workspace.workCapacity, 5)
  }

  test("reserve grows primitive regions independently and retains larger backing identities") {
    val workspace = DenseWorkspace.empty
    workspace.reserve(requirement(7L, 5L))
    assertEquals(workspace.doubleCapacity, 7)
    assertEquals(workspace.indexCapacity, 5)
    val doubles0 = TestAccess.workBacking(workspace)
    val indices0 = TestAccess.indexBacking(workspace)

    workspace.reserve(requirement(3L, 4L))
    assert(TestAccess.sameStorage(doubles0, TestAccess.workBacking(workspace)))
    assert(TestAccess.sameIndexStorage(indices0, TestAccess.indexBacking(workspace)))

    workspace.reserve(requirement(9L, 4L))
    assertEquals(workspace.doubleCapacity, 9)
    assertEquals(workspace.indexCapacity, 5)
    assert(!TestAccess.sameStorage(doubles0, TestAccess.workBacking(workspace)))
    assert(TestAccess.sameIndexStorage(indices0, TestAccess.indexBacking(workspace)))
    val doubles1 = TestAccess.workBacking(workspace)

    workspace.reserve(requirement(8L, 12L))
    assertEquals(workspace.doubleCapacity, 9)
    assertEquals(workspace.indexCapacity, 12)
    assert(TestAccess.sameStorage(doubles1, TestAccess.workBacking(workspace)))
    assert(!TestAccess.sameIndexStorage(indices0, TestAccess.indexBacking(workspace)))
  }

  test("QR reports its exact branch-specific scratch requirement") {
    assertEquals(DenseWorkspace.qrRequirement(3, 2), Right(requirement(3L, 0L)))
    assertEquals(DenseWorkspace.qrRequirement(129, 97), Right(requirement(8256L, 0L)))
    assert(DenseWorkspace.qrRequirement(-1, 2).isLeft)
    assert(DenseWorkspace.qrRequirement(2, -1).isLeft)
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

  test("blocked qrWith reuses panel, WY, and trailing-update workspace") {
    val rng = new scala.util.Random(2026071003L)
    val rows = 129
    val cols = 97
    val A = Matrix.dense(rows, cols, Seq.fill(rows * cols)(rng.nextDouble() * 2.0 - 1.0))
    val workspace = DenseWorkspace.empty

    val first = A.qrWith(workspace)
    val capacity = workspace.workCapacity
    val backing = TestAccess.workBacking(workspace)
    assert(capacity > rows)

    val second = A.qrWith(workspace)
    assertEquals(workspace.workCapacity, capacity)
    assert(TestAccess.sameStorage(backing, TestAccess.workBacking(workspace)))
    assertEquals(first.diagnostics.rank, Some(cols))
    assertEquals(second.diagnostics.rank, Some(cols))
  }
