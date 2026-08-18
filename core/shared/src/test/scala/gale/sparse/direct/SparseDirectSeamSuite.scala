package gale.sparse.direct

import gale.linalg.LinAlgError

/** Empty-provider lock for the shared sparse-direct seam. Runs on JVM, Scala.js,
  * and the experimental Wasm lane. A capable provider is still a later import;
  * this suite only pins that the default `given` advertises nothing.
  */
class SparseDirectSeamSuite extends munit.FunSuite:

  test("default provider advertises no capability and cannot create a workspace") {
    assertEquals(SparseDirectProvider.none.name, "none")
    assertEquals(SparseDirectProvider.none.capabilities, Set.empty)
    assertEquals(SparseDirect.capabilities, Set.empty)
    assertEquals(SparseDirectProvider.validationErrors(SparseDirectProvider.none), Nil)

    val report = SparseDirectProvider.current
    assertEquals(report.name, "none")
    assertEquals(report.capabilities, Set.empty)
    assert(!report.available)
    assertEquals(report.config, gale.backend.BackendConfig.singleThreaded)

    SparseDirect.newWorkspace() match
      case Left(LinAlgError.UnsupportedOperation(op)) =>
        assert(op.contains("sparse direct"), s"unexpected message: $op")
        assert(!op.contains("JVM"), s"error text must be platform-neutral: $op")
        assert(op.contains("no provider is installed"), s"unexpected message: $op")
      case other =>
        fail(s"expected UnsupportedOperation, got $other")
  }

  test("no factorization family is advertised until a capable provider is imported") {
    for family <- SparseDirectFactorization.values do
      assert(
        !SparseDirect.capabilities.contains(family.requiredCapability),
        s"$family leaked onto the default provider"
      )
  }
