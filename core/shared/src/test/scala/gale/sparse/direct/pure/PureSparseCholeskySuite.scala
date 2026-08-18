package gale.sparse.direct.pure

import gale.linalg.*
import gale.sparse.*
import gale.sparse.direct.*

class PureSparseCholeskySuite extends munit.FunSuite:
  import gale.sparse.direct.pure.given

  private def csr(rows: Int, cols: Int, entries: (Int, Int, Double)*): CSR =
    val builder = Sparse.coo(rows, cols)
    entries.foreach { case (row, col, value) => builder.add(row, col, value) }
    builder.toCSR()

  private def csrFromDense(a: DMat): CSR =
    val builder = Sparse.coo(a.rows, a.cols)
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        val v = a(i, j)
        if v != 0.0 then builder.add(i, j, v)
        j += 1
      i += 1
    builder.toCSR()

  private def residual(a: CSR, x: DVec, b: DVec): Double =
    val ax = MutableDVec.zeros(a.rows)
    a.mulInto(x, ax)
    (ax.toVec - b).norm2

  private def assertVecClose(actual: DVec, expected: DVec, tolerance: Double, clue: String): Unit =
    assertEquals(actual.length, expected.length, clue)
    var i = 0
    while i < actual.length do
      assert(
        math.abs(actual(i) - expected(i)) <= tolerance,
        s"$clue: index $i: ${actual(i)} != ${expected(i)}"
      )
      i += 1

  test("pure given advertises Cholesky and not LU or QR") {
    assertEquals(SparseDirectProvider.pure.name, "pure")
    assert(SparseDirect.capabilities.contains(SparseDirectCapability.Cholesky))
    assert(SparseDirect.capabilities.contains(SparseDirectCapability.UserOrdering))
    assert(SparseDirect.capabilities.contains(SparseDirectCapability.TransposeSolve))
    assert(SparseDirect.capabilities.contains(SparseDirectCapability.MultipleRhs))
    assert(!SparseDirect.capabilities.contains(SparseDirectCapability.LU))
    assert(!SparseDirect.capabilities.contains(SparseDirectCapability.QR))
    assertEquals(SparseDirectProvider.validationErrors(SparseDirectProvider.pure), Nil)
  }

  test("default given stays none when the pure given is not imported at the call") {
    assertEquals(SparseDirect.capabilities(using SparseDirectProvider.none), Set.empty)
  }

  test("Natural Cholesky solve matches dense Cholesky on a small SPD matrix") {
    val dense = Matrix.dense(3, 3)(6.0, 2.0, 1.0, 2.0, 5.0, 2.0, 1.0, 2.0, 4.0)
    val matrix = csrFromDense(dense)
    val workspace = SparseDirect.newWorkspace().orThrow
    val analysis =
      SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.Cholesky, workspace, SparseDirectOrdering.Natural).orThrow
    val factor = SparseDirect.factor(analysis, matrix, workspace).orThrow
    val b = Vec(1.0, 0.0, -1.0)
    val x = SparseDirect.solve(factor, b, workspace).orThrow.solution
    val expected = dense.cholesky.orThrow.solve(b).orThrow
    assertVecClose(x, expected, 1e-12, "Natural vs dense")
    assert(residual(matrix, x, b) < 1e-12, s"CSR residual ${residual(matrix, x, b)}")
    workspace.close()
    factor.close()
    analysis.close()
  }

  test("ProviderDefault minimum degree and User ordering still solve with a small residual") {
    val dense = Matrix.dense(4, 4)(
      10.0, 1.0, 0.0, 2.0,
      1.0, 8.0, 1.0, 0.0,
      0.0, 1.0, 7.0, 1.0,
      2.0, 0.0, 1.0, 9.0
    )
    val matrix = csrFromDense(dense)
    val b = Vec(1.0, -2.0, 0.5, 3.0)
    val expected = dense.cholesky.orThrow.solve(b).orThrow

    val defaultWs = SparseDirect.newWorkspace().orThrow
    val defaultAnalysis =
      SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.Cholesky, defaultWs, SparseDirectOrdering.ProviderDefault).orThrow
    val defaultFactor = SparseDirect.factor(defaultAnalysis, matrix, defaultWs).orThrow
    val xDefault = SparseDirect.solve(defaultFactor, b, defaultWs).orThrow.solution
    assertVecClose(xDefault, expected, 1e-11, "minimum degree")
    assertEquals(defaultAnalysis.columnPermutation.rows, 4)
    assert(defaultAnalysis.diagnostics.deterministic)
    assertEquals(defaultAnalysis.diagnostics.ordering, SparseDirectOrdering.ProviderDefault)

    val userPerm = Sparse.permutation(3, 1, 0, 2)
    val userWs = SparseDirect.newWorkspace().orThrow
    val userAnalysis =
      SparseDirect
        .analyze(matrix.pattern, SparseDirectFactorization.Cholesky, userWs, SparseDirectOrdering.User(userPerm))
        .orThrow
    val userFactor = SparseDirect.factor(userAnalysis, matrix, userWs).orThrow
    val xUser = SparseDirect.solve(userFactor, b, userWs).orThrow.solution
    assertVecClose(xUser, expected, 1e-11, "user ordering")
    assertEquals(userAnalysis.columnPermutation.toIndexSeq, userPerm.toIndexSeq)
  }

  test("symbolic analysis is reusable across rebind and independent of analysis close") {
    val dense0 = Matrix.dense(2, 2)(4.0, 1.0, 1.0, 3.0)
    val matrix0 = csrFromDense(dense0)
    val workspace = SparseDirect.newWorkspace().orThrow
    val analysis =
      SparseDirect.analyze(matrix0.pattern, SparseDirectFactorization.Cholesky, workspace, SparseDirectOrdering.Natural).orThrow
    val factor0 = SparseDirect.factor(analysis, matrix0, workspace).orThrow
    val rebound = analysis.pattern.bind(Array(5.0, 2.0, 2.0, 4.0)).orThrow
    val factor1 = SparseDirect.factor(analysis, rebound, workspace).orThrow
    val b = Vec(1.0, 2.0)
    val x0 = SparseDirect.solve(factor0, b, workspace).orThrow.solution
    val x1 = SparseDirect.solve(factor1, b, workspace).orThrow.solution
    assertVecClose(x0, dense0.cholesky.orThrow.solve(b).orThrow, 1e-12, "first values")
    val dense1 = Matrix.dense(2, 2)(5.0, 2.0, 2.0, 4.0)
    assertVecClose(x1, dense1.cholesky.orThrow.solve(b).orThrow, 1e-12, "rebound values")
    analysis.close()
    assert(SparseDirect.factor(analysis, matrix0, workspace).isLeft)
    val still = SparseDirect.solve(factor0, b, workspace).orThrow.solution
    assertVecClose(still, x0, 1e-12, "factor outlives analysis")
  }

  test("transpose and multiple-RHS solves honor the SPD contract") {
    val dense = Matrix.dense(2, 2)(4.0, 1.0, 1.0, 3.0)
    val matrix = csrFromDense(dense)
    val workspace = SparseDirect.newWorkspace().orThrow
    val analysis = SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.Cholesky, workspace).orThrow
    val factor = SparseDirect.factor(analysis, matrix, workspace).orThrow
    val b = Vec(1.0, 2.0)
    val normal = SparseDirect.solve(factor, b, workspace).orThrow
    val transposed = SparseDirect.solve(factor, b, workspace, SparseSolveOperation.Transpose).orThrow
    assertVecClose(transposed.solution, normal.solution, 1e-12, "A^T = A")
    assertEquals(transposed.diagnostics.operation, SparseSolveOperation.Transpose)

    val rhs = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val solved = SparseDirect.solve(factor, rhs, workspace).orThrow
    val expected = dense.cholesky.orThrow.solve(rhs).orThrow
    var i = 0
    while i < 2 do
      var j = 0
      while j < 2 do
        assert(math.abs(solved.solution(i, j) - expected(i, j)) < 1e-12, s"($i,$j)")
        j += 1
      i += 1
    assertEquals(solved.diagnostics.rightHandSides, 2)
  }

  test("indefinite, empty, and LU requests are typed failures") {
    val indefinite = csr(2, 2, (0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 1.0))
    val workspace = SparseDirect.newWorkspace().orThrow
    val analysis = SparseDirect.analyze(indefinite.pattern, SparseDirectFactorization.Cholesky, workspace).orThrow
    analysis.factorNumeric(indefinite, workspace) match
      case Left(_: LinAlgError.NotPositiveDefinite) => ()
      case other                                    => fail(s"expected NotPositiveDefinite, got $other")

    val empty = Sparse.coo(0, 0).toCSR()
    val emptyWs = SparseDirect.newWorkspace().orThrow
    val emptyAnalysis = SparseDirect.analyze(empty.pattern, SparseDirectFactorization.Cholesky, emptyWs).orThrow
    val emptyFactor = SparseDirect.factor(emptyAnalysis, empty, emptyWs).orThrow
    val emptySolved = SparseDirect.solve(emptyFactor, Vec(), emptyWs).orThrow
    assertEquals(emptySolved.solution.length, 0)

    SparseDirect.analyze(indefinite.pattern, SparseDirectFactorization.LU, workspace) match
      case Left(LinAlgError.UnsupportedOperation(op)) =>
        assert(op.contains("LU"), s"unexpected message: $op")
      case other =>
        fail(s"expected UnsupportedOperation for LU, got $other")
  }

  test("closed workspace and closed factor are rejected before a solve") {
    val matrix = csr(1, 1, (0, 0, 2.0))
    val workspace = SparseDirect.newWorkspace().orThrow
    val analysis = SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.Cholesky, workspace).orThrow
    val factor = SparseDirect.factor(analysis, matrix, workspace).orThrow
    factor.close()
    assert(SparseDirect.solve(factor, Vec(1.0), workspace).isLeft)
    val openFactor = SparseDirect.factor(analysis, matrix, workspace).orThrow
    workspace.close()
    assert(SparseDirect.solve(openFactor, Vec(1.0), workspace).isLeft)
    workspace.close()
  }

  test("Natural arrowhead fill that exceeds the documented guard is InvalidArgument") {
    val n = 4100
    val builder = Sparse.coo(n, n)
    var i = 0
    while i < n do
      builder.add(i, i, 1.0)
      if i > 0 then builder.add(i, 0, 1.0)
      i += 1
    val matrix = builder.toCSR()
    val workspace = SparseDirect.newWorkspace().orThrow
    SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.Cholesky, workspace, SparseDirectOrdering.Natural) match
      case Left(LinAlgError.InvalidArgument(message)) =>
        assert(message.contains("fill estimate"), s"unexpected message: $message")
        assert(message.contains("exceeds guard"), s"unexpected message: $message")
      case other =>
        fail(s"expected InvalidArgument fill guard, got $other")
  }

  test("fill-guard constants match the documented rule") {
    assert(SparseCholeskyKernels.exceedsFill(100, SparseCholeskyKernels.AbsoluteNnzCap + 1))
    assert(!SparseCholeskyKernels.exceedsFill(100, 100L * SparseCholeskyKernels.FillFactor))
    assertEquals(SparseCholeskyKernels.fillCap(100), SparseCholeskyKernels.AbsoluteNnzCap)
    assertEquals(SparseCholeskyKernels.fillCap(200_000), 200_000L * SparseCholeskyKernels.FillFactor)
  }
