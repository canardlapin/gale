package gale.sparse.direct

import gale.backend.BackendConfig
import gale.linalg.*
import gale.sparse.*
import java.util.concurrent.atomic.AtomicInteger
import scala.concurrent.Await
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.concurrent.duration.*

class SparseDirectProviderSuite extends munit.FunSuite:
  private def csr(rows: Int, cols: Int, entries: (Int, Int, Double)*): CSR =
    val builder = Sparse.coo(rows, cols)
    entries.foreach { case (row, col, value) => builder.add(row, col, value) }
    builder.toCSR()

  private def identityPermutation(size: Int): Permutation =
    Sparse.permutation((0 until size)*)

  private def assertVecApprox(actual: DVec, expected: DVec, tolerance: Double = 1e-12): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assert(math.abs(actual(i) - expected(i)) <= tolerance, s"index $i: ${actual(i)} != ${expected(i)}")
      i += 1

  private def assertMatApprox(actual: DMat, expected: DMat, tolerance: Double = 1e-12): Unit =
    assertEquals(actual.shape, expected.shape)
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assert(
          math.abs(actual(row, col) - expected(row, col)) <= tolerance,
          s"($row,$col): ${actual(row, col)} != ${expected(row, col)}"
        )
        col += 1
      row += 1

  private final class TestWorkspace(val provider: SparseDirectProvider) extends SparseDirectWorkspace:
    private var closed = false
    def isClosed: Boolean = closed
    def close(): Unit = closed = true

  private final class TestAnalysis(
      val provider: TestProvider,
      val pattern: CSRPattern,
      val factorization: SparseDirectFactorization,
      val ordering: SparseDirectOrdering
  ) extends SparseDirectSymbolicAnalysis:
    private var closed = false
    val columnPermutation: Permutation = identityPermutation(pattern.cols)
    val diagnostics: SparseSymbolicDiagnostics =
      SparseSymbolicDiagnostics(
        provider.name,
        factorization,
        pattern.nnz,
        predictedFactorNnz = Some(pattern.nnz.toLong),
        ordering,
        deterministic = true
      )

    def isClosed: Boolean = closed
    def close(): Unit = closed = true

    def factorNumeric(
        matrix: CSR,
        workspace: SparseDirectWorkspace
    ): Either[LinAlgError, SparseDirectNumericFactor] =
      provider.factorCalls.incrementAndGet()
      val dense = matrix.toDense()
      DenseDecompositions.lu(dense).flatMap: normal =>
        DenseDecompositions.lu(dense.t).map: transpose =>
          new TestFactor(provider, normal, transpose, dense.rows)

  private final class TestFactor(
      val provider: TestProvider,
      normal: LU,
      transpose: LU,
      val inputRows: Int
  ) extends SparseDirectNumericFactor:
    val inputCols: Int = inputRows
    val factorization: SparseDirectFactorization = SparseDirectFactorization.LU
    private var closed = false
    val solveCalls = new AtomicInteger(0)
    val rowPermutation: Permutation = Sparse.permutation(normal.pivots.toIndexSeq*)
    val columnPermutation: Permutation = identityPermutation(inputCols)
    val diagnostics: SparseNumericDiagnostics =
      SparseNumericDiagnostics(
        provider.name,
        factorization,
        factorNnz = None,
        rank = Some(inputRows),
        pivotCount = Some(normal.pivots.toIndexSeq.zipWithIndex.count { case (pivot, index) => pivot != index }),
        reciprocalConditionEstimate = None
      )

    def isClosed: Boolean = closed
    def close(): Unit = closed = true
    def rhsRows(operation: SparseSolveOperation): Int = inputRows
    def solutionRows(operation: SparseSolveOperation): Int = inputRows

    private def selected(operation: SparseSolveOperation): LU =
      if operation == SparseSolveOperation.Normal then normal else transpose

    def solveVectorInto(
        rhs: DVec,
        destination: MutableDVec,
        operation: SparseSolveOperation,
        workspace: SparseDirectWorkspace
    ): Either[LinAlgError, SparseSolveDiagnostics] =
      solveCalls.incrementAndGet()
      selected(operation).solve(rhs).map: solution =>
        var i = 0
        while i < solution.length do
          destination(i) = solution(i)
          i += 1
        SparseSolveDiagnostics(provider.name, operation, 1, residualNorm = None, refinementSteps = 0)

    def solveMatrixInto(
        rhs: DMat,
        destination: DMatBuilder,
        operation: SparseSolveOperation,
        workspace: SparseDirectWorkspace
    ): Either[LinAlgError, SparseSolveDiagnostics] =
      solveCalls.incrementAndGet()
      selected(operation).solve(rhs).map: solution =>
        var row = 0
        while row < solution.rows do
          var col = 0
          while col < solution.cols do
            destination(row, col) = solution(row, col)
            col += 1
          row += 1
        SparseSolveDiagnostics(provider.name, operation, rhs.cols, residualNorm = None, refinementSteps = 0)

  private final class TestProvider extends SparseDirectProvider:
    val name: String = "deterministic-test-lu"
    val capabilities: Set[SparseDirectCapability] = Set(
      SparseDirectCapability.LU,
      SparseDirectCapability.TransposeSolve,
      SparseDirectCapability.MultipleRhs
    )
    val config: BackendConfig = BackendConfig.singleThreaded
    val analyzeCalls = new AtomicInteger(0)
    val factorCalls = new AtomicInteger(0)

    def createWorkspace(): Either[LinAlgError, SparseDirectWorkspace] =
      Right(new TestWorkspace(this))

    def analyze(
        pattern: CSRPattern,
        factorization: SparseDirectFactorization,
        ordering: SparseDirectOrdering,
        workspace: SparseDirectWorkspace
    ): Either[LinAlgError, SparseDirectSymbolicAnalysis] =
      analyzeCalls.incrementAndGet()
      Right(new TestAnalysis(this, pattern, factorization, ordering))

  test("capability-less default reports unavailable and never implies implementation") {
    val report = SparseDirectProvider.current(using SparseDirectProvider.none)
    assertEquals(report.name, "none")
    assertEquals(report.capabilities, Set.empty)
    assert(!report.available)
    assert(SparseDirect.newWorkspace()(using SparseDirectProvider.none).isLeft)
  }

  test("symbolic analysis is reusable across changing numeric values and exposes owned diagnostics") {
    val provider = new TestProvider
    val matrix = csr(2, 2, (0, 0, 4.0), (0, 1, 1.0), (1, 0, 2.0), (1, 1, 3.0))
    val workspace = SparseDirect.newWorkspace()(using provider).toOption.get
    val analysis =
      SparseDirect
        .analyze(matrix.pattern, SparseDirectFactorization.LU, workspace, SparseDirectOrdering.Natural)(using provider)
        .toOption
        .get
    assertEquals(analysis.diagnostics.inputNnz, 4)
    assertEquals(analysis.diagnostics.ordering, SparseDirectOrdering.Natural)
    assert(analysis.diagnostics.deterministic)
    val permutationCopy = analysis.columnPermutation.toArray
    permutationCopy(0) = 1
    assertEquals(analysis.columnPermutation.toIndexSeq, IndexedSeq(0, 1))

    val factor0 = SparseDirect.factor(analysis, matrix, workspace).toOption.get
    val rebound = analysis.pattern.bind(Array(5.0, 2.0, 1.0, 4.0)).toOption.get
    val factor1 = SparseDirect.factor(analysis, rebound, workspace).toOption.get
    assertEquals(provider.factorCalls.get(), 2)

    val rhs = Vec(1.0, 2.0)
    val expected0 = DenseDecompositions.lu(matrix.toDense()).toOption.get.solve(rhs).toOption.get
    val expected1 = DenseDecompositions.lu(rebound.toDense()).toOption.get.solve(rhs).toOption.get
    assertVecApprox(SparseDirect.solve(factor0, rhs, workspace).toOption.get.solution, expected0)
    assertVecApprox(SparseDirect.solve(factor1, rhs, workspace).toOption.get.solution, expected1)
  }

  test("normal, transpose, destination, and multiple-RHS solves honor the typed contract") {
    val provider = new TestProvider
    val matrix = csr(2, 2, (0, 0, 4.0), (0, 1, 1.0), (1, 0, 2.0), (1, 1, 3.0))
    val workspace = SparseDirect.newWorkspace()(using provider).toOption.get
    val analysis = SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.LU, workspace)(using provider).toOption.get
    val factor = SparseDirect.factor(analysis, matrix, workspace).toOption.get.asInstanceOf[TestFactor]
    val rhs = Vec(1.0, 2.0)

    val expectedTranspose = DenseDecompositions.lu(matrix.toDense().t).toOption.get.solve(rhs).toOption.get
    val transposed =
      SparseDirect.solve(factor, rhs, workspace, SparseSolveOperation.Transpose).toOption.get
    assertVecApprox(transposed.solution, expectedTranspose)
    assertEquals(transposed.diagnostics.operation, SparseSolveOperation.Transpose)

    val destination = MutableVec.zeros(2)
    assert(SparseDirect.solveInto(factor, rhs, destination, workspace).isRight)
    val expectedNormal = DenseDecompositions.lu(matrix.toDense()).toOption.get.solve(rhs).toOption.get
    assertVecApprox(destination.asVec, expectedNormal)

    val rhsMatrix = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val expectedMatrix = DenseDecompositions.lu(matrix.toDense()).toOption.get.solve(rhsMatrix).toOption.get
    val solvedMatrix = SparseDirect.solve(factor, rhsMatrix, workspace).toOption.get
    assertMatApprox(solvedMatrix.solution, expectedMatrix)
    assertEquals(solvedMatrix.diagnostics.rightHandSides, 2)

    val calls = factor.solveCalls.get()
    val wrongDestination = MutableVec.zeros(1)
    assert(SparseDirect.solveInto(factor, rhs, wrongDestination, workspace).isLeft)
    assertEquals(factor.solveCalls.get(), calls, "provider solve ran after facade dimension rejection")
  }

  test("changed patterns, wrong workspaces, unsupported capabilities, and closed resources are typed failures") {
    val provider = new TestProvider
    val matrix = csr(2, 2, (0, 0, 4.0), (0, 1, 1.0), (1, 0, 2.0), (1, 1, 3.0))
    val workspace = SparseDirect.newWorkspace()(using provider).toOption.get
    val analysis = SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.LU, workspace)(using provider).toOption.get

    val changed = csr(2, 2, (0, 0, 4.0), (1, 1, 3.0))
    val factorCalls = provider.factorCalls.get()
    assert(SparseDirect.factor(analysis, changed, workspace).isLeft)
    assertEquals(provider.factorCalls.get(), factorCalls)

    assert(SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.Cholesky, workspace)(using provider).isLeft)
    assert(SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.QR, workspace)(using provider).isLeft)
    assert(
      SparseDirect
        .analyze(
          matrix.pattern,
          SparseDirectFactorization.LU,
          workspace,
          SparseDirectOrdering.User(identityPermutation(2))
        )(using provider)
        .isLeft
    )
    val noncanonical = CSRPattern.checked(2, 2, Array(0, 2, 2), Array(1, 0)).toOption.get
    assert(SparseDirect.analyze(noncanonical, SparseDirectFactorization.LU, workspace)(using provider).isLeft)
    val rectangular = csr(2, 3, (0, 0, 1.0)).pattern
    assert(SparseDirect.analyze(rectangular, SparseDirectFactorization.LU, workspace)(using provider).isLeft)

    val anotherProvider = new TestProvider
    val wrongWorkspace = SparseDirect.newWorkspace()(using anotherProvider).toOption.get
    assert(SparseDirect.factor(analysis, matrix, wrongWorkspace).isLeft)

    val factor = SparseDirect.factor(analysis, matrix, workspace).toOption.get
    analysis.close()
    assert(SparseDirect.factor(analysis, matrix, workspace).isLeft)
    // A completed numeric factor owns an independent lifetime.
    assert(SparseDirect.solve(factor, Vec(1.0, 2.0), workspace).isRight)
    factor.close()
    assert(SparseDirect.solve(factor, Vec(1.0, 2.0), workspace).isLeft)
    workspace.close()
    assert(SparseDirect.solve(factor, Vec(1.0, 2.0), workspace).isLeft)
    workspace.close() // close is idempotent
  }

  test("one immutable factor supports concurrent solves with distinct workspaces") {
    val provider = new TestProvider
    val matrix = csr(2, 2, (0, 0, 4.0), (0, 1, 1.0), (1, 0, 2.0), (1, 1, 3.0))
    val analysisWorkspace = SparseDirect.newWorkspace()(using provider).toOption.get
    val analysis =
      SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.LU, analysisWorkspace)(using provider).toOption.get
    val factor = SparseDirect.factor(analysis, matrix, analysisWorkspace).toOption.get
    val expected = DenseDecompositions.lu(matrix.toDense()).toOption.get.solve(Vec(1.0, 2.0)).toOption.get
    val workspaces = Vector.fill(8)(SparseDirect.newWorkspace()(using provider).toOption.get)
    val futures = workspaces.map: workspace =>
      Future(SparseDirect.solve(factor, Vec(1.0, 2.0), workspace).toOption.get.solution)
    val results = Await.result(Future.sequence(futures), 10.seconds)
    results.foreach(result => assertVecApprox(result, expected))
  }

  test("facade rejects a provider that returns a handle owned by another provider") {
    val delegate = new TestProvider
    val malformed = new SparseDirectProvider:
      val name: String = "malformed-owner"
      val capabilities: Set[SparseDirectCapability] = Set(SparseDirectCapability.LU)
      val config: BackendConfig = BackendConfig.singleThreaded
      def createWorkspace(): Either[LinAlgError, SparseDirectWorkspace] = Right(new TestWorkspace(this))
      def analyze(
          pattern: CSRPattern,
          factorization: SparseDirectFactorization,
          ordering: SparseDirectOrdering,
          workspace: SparseDirectWorkspace
      ): Either[LinAlgError, SparseDirectSymbolicAnalysis] =
        Right(new TestAnalysis(delegate, pattern, factorization, ordering))

    val matrix = csr(2, 2, (0, 0, 1.0), (1, 1, 1.0))
    val workspace = SparseDirect.newWorkspace()(using malformed).toOption.get
    assert(SparseDirect.analyze(matrix.pattern, SparseDirectFactorization.LU, workspace)(using malformed).isLeft)
  }
