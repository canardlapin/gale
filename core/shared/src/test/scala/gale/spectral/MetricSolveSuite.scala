package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.linalg.Vec
import gale.linalg.assumePositiveDefinite
import gale.linalg.orThrow
import gale.solvers.SolverConfig
import gale.solvers.ToleranceMode

class MetricSolveSuite extends munit.FunSuite:

  private val diagonal: DMat =
    DMat.tabulate(4, 4)((i, j) => if i == j then (i + 1).toDouble else 0.0)

  private def assertVec(actual: DVec, expected: IndexedSeq[Double], tolerance: Double = 1e-12): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEqualsDouble(actual(i), expected(i), tolerance)
      i += 1

  test("direct adapter reuses a caller-created factor and owns its solution") {
    val factor = diagonal.cholesky.orThrow
    val solve = LinearSolveOperator.direct(factor)
    val rhs = Vec(2.0, 6.0, 12.0, 20.0)
    val result = solve.solve(rhs).orThrow

    assertEquals(solve.kind, LinearSolveKind.DirectFactor)
    assertVec(result.solution, IndexedSeq(2.0, 3.0, 4.0, 5.0))
    assertEquals(
      result.diagnostics,
      LinearSolveDiagnostics(
        converged = true,
        iterations = 0,
        residualNorm = None,
        operatorApplications = 0L
      )
    )
  }

  test("CG adapter reports convergence, residual, and exact operator applications") {
    val solve =
      LinearSolveOperator
        .conjugateGradient(
          diagonal.assumePositiveDefinite,
          SolverConfig(tolerance = 1e-12, maxIterations = 20),
          toleranceMode = ToleranceMode.RelativeToRhs
        )
        .orThrow
    val result = solve.solve(Vec(2.0, 6.0, 12.0, 20.0)).orThrow

    assert(result.diagnostics.converged)
    var i = 0
    while i < result.solution.length do
      assertEqualsDouble(result.solution(i), (i + 2).toDouble, 1e-10)
      i += 1
    assert(result.diagnostics.residualNorm.exists(_ <= 1e-10))
    assertEquals(
      result.diagnostics.operatorApplications,
      result.diagnostics.iterations.toLong + 1L
    )
  }

  test("iterative non-convergence is a result, not a structural Left") {
    val solve =
      LinearSolveOperator
        .conjugateGradient(
          diagonal.assumePositiveDefinite,
          SolverConfig(tolerance = 0.0, maxIterations = 0)
        )
        .orThrow
    val result = solve.solve(Vec(1.0, 1.0, 1.0, 1.0)).orThrow

    assert(!result.diagnostics.converged)
    assertEquals(result.diagnostics.iterations, 0)
    assertEquals(result.diagnostics.operatorApplications, 1L)
    assert(result.diagnostics.residualNorm.exists(_ > 0.0))
  }

  test("solve boundary validates dimensions, backend diagnostics, and returned ownership") {
    var retained: Option[DVec] = None
    val malformed =
      LinearSolveOperator
        .backendProvided(2): rhs =>
          retained = Some(rhs)
          Right(
            LinearSolveResult(
              Vec(1.0),
              LinearSolveDiagnostics(true, 0, None, 0L)
            )
          )
        .orThrow

    assert(malformed.solve(Vec(1.0)).isLeft)
    malformed.solve(Vec(1.0, 2.0)) match
      case Left(_: LinAlgError.VectorLengthMismatch) => ()
      case other                                     => fail(s"expected malformed result rejection, got $other")
    assert(retained.nonEmpty)
  }

  test("backend solve validation rejects malformed diagnostics and non-finite solutions") {
    val valid = LinearSolveDiagnostics(converged = true, iterations = 0, residualNorm = None, operatorApplications = 0L)

    def solveWith(solution: DVec, diagnostics: LinearSolveDiagnostics): Either[LinAlgError, LinearSolveResult] =
      LinearSolveOperator
        .backendProvided(2)(_ => Right(LinearSolveResult(solution, diagnostics)))
        .orThrow
        .solve(Vec(1.0, 2.0))

    val malformed = Seq(
      solveWith(Vec(1.0, 2.0), valid.copy(iterations = -1)),
      solveWith(Vec(1.0, 2.0), valid.copy(operatorApplications = -1L)),
      solveWith(Vec(1.0, 2.0), valid.copy(residualNorm = Some(Double.NaN))),
      solveWith(Vec(1.0, 2.0), valid.copy(residualNorm = Some(-1.0))),
      solveWith(Vec(1.0, 2.0), valid.copy(converged = false, residualNorm = None)),
      solveWith(Vec(Double.PositiveInfinity, 2.0), valid)
    )

    malformed.foreach:
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected malformed backend solve rejection, got $other")

    assert(LinearSolveOperator.backendProvided(-1)(_ => fail("negative-size backend must not run")).isLeft)

    val throwing =
      LinearSolveOperator
        .backendProvided(2)(_ => throw LinAlgError.InvalidArgument("provider failure"))
        .orThrow
    throwing.solve(Vec(1.0, 2.0)) match
      case Left(LinAlgError.InvalidArgument("provider failure")) => ()
      case other => fail(s"expected retained typed provider failure, got $other")
  }

  test("CG construction rejects invalid operator and iteration policy before solving") {
    assert(
      LinearSolveOperator
        .conjugateGradient(Matrix.zeros(2, 3).assumePositiveDefinite)
        .isLeft
    )
    assert(
      LinearSolveOperator
        .conjugateGradient(diagonal.assumePositiveDefinite, SolverConfig(tolerance = Double.NaN))
        .isLeft
    )
    assert(
      LinearSolveOperator
        .conjugateGradient(diagonal.assumePositiveDefinite, SolverConfig(tolerance = -1.0))
        .isLeft
    )
    assert(
      LinearSolveOperator
        .conjugateGradient(diagonal.assumePositiveDefinite, SolverConfig(maxIterations = -1))
        .isLeft
    )
  }

  test("metric binding rejects mismatched solves and represents B-inverse action") {
    val metric = diagonal.assumePositiveDefinite
    val solve = LinearSolveOperator.direct(diagonal.cholesky.orThrow)
    val bound = MetricSolveOperator.bind(metric, solve).orThrow

    assertEquals(bound.size, 4)
    assertVec(
      bound.solve(Vec(2.0, 6.0, 12.0, 20.0)).orThrow.solution,
      IndexedSeq(2.0, 3.0, 4.0, 5.0)
    )

    val wrong =
      LinearSolveOperator
        .backendProvided(3): rhs =>
          Right(LinearSolveResult(rhs, LinearSolveDiagnostics(true, 0, None, 0L)))
        .orThrow
    MetricSolveOperator.bind(metric, wrong) match
      case Left(_: LinAlgError.VectorLengthMismatch) => ()
      case other                                     => fail(s"expected metric/solve mismatch, got $other")
  }

  test("B-inverse A is self-adjoint in the B inner product") {
    val b = DMat.tabulate(3, 3)((i, j) => if i == j then (i + 2).toDouble else 0.0)
    val entries = Array(
      Array(4.0, 1.0, -0.5),
      Array(1.0, 5.0, 0.25),
      Array(-0.5, 0.25, 6.0)
    )
    val a = Matrix.tabulate(3, 3)((i, j) => entries(i)(j))
    val metric = MetricSolveOperator
      .bind(b.assumePositiveDefinite, LinearSolveOperator.direct(b.cholesky.orThrow))
      .orThrow
    val x = Vec(0.5, -1.0, 2.0)
    val y = Vec(1.5, 0.25, -0.75)
    val tx = metric.solve(a * x).orThrow.solution
    val ty = metric.solve(a * y).orThrow.solution

    val left = x.dot(b * ty)
    val right = tx.dot(b * y)
    assertEqualsDouble(left, right, 1e-12)
  }

  test("LinearSolvePlan resolves an explicit solve or an advertised backend") {
    val factorPlan = LinearSolvePlan.direct(diagonal.cholesky.orThrow)
    val explicit =
      LinearSolvePlan
        .resolve(factorPlan, diagonal, None, 0.5)(using SpectralBackend.none)
        .orThrow
    assertEquals(explicit.kind, LinearSolveKind.DirectFactor)

    val backendSolve =
      LinearSolveOperator
        .backendProvided(4): rhs =>
          Right(LinearSolveResult(rhs, LinearSolveDiagnostics(true, 0, None, 0L)))
        .orThrow
    given SpectralBackend = new SpectralBackend:
      def name: String = "solve-test"
      def capabilities: Set[SpectralCapability] = Set(SpectralCapability.ShiftInvertSolve)
      override def shiftInvertSolve(
          a: DMat,
          b: Option[DMat],
          sigma: Double
      ): Either[LinAlgError, LinearSolveOperator] =
        Right(backendSolve)

    val resolved = LinearSolvePlan.resolve(LinearSolvePlan.Backend, diagonal, None, 0.5).orThrow
    assertEquals(resolved.kind, LinearSolveKind.BackendProvided)
  }

  test("LinearSolvePlan rejects invalid shapes, shifts, capabilities, and returned sizes") {
    val solve3 =
      LinearSolveOperator
        .backendProvided(3): rhs =>
          Right(LinearSolveResult(rhs, LinearSolveDiagnostics(true, 0, None, 0L)))
        .orThrow

    assert(LinearSolvePlan.resolve(LinearSolvePlan.Use(solve3), diagonal, None, 0.0)(using SpectralBackend.none).isLeft)
    assert(
      LinearSolvePlan
        .resolve(LinearSolvePlan.Use(solve3), Matrix.zeros(2, 3), None, 0.0)(using SpectralBackend.none)
        .isLeft
    )
    assert(
      LinearSolvePlan
        .resolve(LinearSolvePlan.Use(solve3), diagonal, Some(Matrix.eye(3)), 0.0)(using SpectralBackend.none)
        .isLeft
    )
    assert(
      LinearSolvePlan
        .resolve(LinearSolvePlan.Use(solve3), diagonal, None, Double.NaN)(using SpectralBackend.none)
        .isLeft
    )
    assert(LinearSolvePlan.resolve(LinearSolvePlan.Backend, diagonal, None, 0.0)(using SpectralBackend.none).isLeft)

    val wrongSizeBackend = new SpectralBackend:
      def name: String = "wrong-size-solve"
      def capabilities: Set[SpectralCapability] = Set(SpectralCapability.ShiftInvertSolve)
      override def shiftInvertSolve(
          a: DMat,
          b: Option[DMat],
          sigma: Double
      ): Either[LinAlgError, LinearSolveOperator] =
        Right(solve3)

    assert(LinearSolvePlan.resolve(LinearSolvePlan.Backend, diagonal, None, 0.0)(using wrongSizeBackend).isLeft)
  }
