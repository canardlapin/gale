package gale.spectral

import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.ExactSolveFactor
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import gale.linalg.PositiveDefinite
import gale.linalg.Shape
import gale.linalg.Rows
import gale.linalg.Cols
import gale.solvers.IterativeSolvers
import gale.solvers.Preconditioner
import gale.solvers.SolverConfig
import gale.solvers.ToleranceMode

/** Provenance of an executable linear solve.
  *
  * This is an ADT rather than a string tag so routing and diagnostics cannot
  * silently invent new, unhandled solve modes.
  */
enum LinearSolveKind:
  case DirectFactor, ConjugateGradient, BackendProvided

/** Diagnostics for one application of a [[LinearSolveOperator]].
  *
  * `operatorApplications` counts applications of the system operator made by
  * this solve. A reusable direct factor reports zero because triangular solves
  * do not reapply the original operator. `residualNorm = None` means the solve
  * provider did not measure a residual; it never means a measured zero.
  */
final case class LinearSolveDiagnostics(
    converged: Boolean,
    iterations: Int,
    residualNorm: Option[Double],
    operatorApplications: Long
)

/** One independently owned solution and its per-call diagnostics. */
final case class LinearSolveResult(
    solution: DVec,
    diagnostics: LinearSolveDiagnostics
)

/** Aggregated inner-solve work attached to a spectral result. */
final case class LinearSolveSummary(
    solves: Int,
    converged: Int,
    iterations: Long,
    operatorApplications: Long,
    worstResidualNorm: Option[Double]
):
  private[spectral] def append(diagnostics: LinearSolveDiagnostics): LinearSolveSummary =
    val worst =
      (worstResidualNorm, diagnostics.residualNorm) match
        case (None, next)             => next
        case (current, None)          => current
        case (Some(a), Some(b))       => Some(math.max(a, b))
    LinearSolveSummary(
      solves = solves + 1,
      converged = converged + (if diagnostics.converged then 1 else 0),
      iterations = iterations + diagnostics.iterations.toLong,
      operatorApplications = operatorApplications + diagnostics.operatorApplications,
      worstResidualNorm = worst
    )

object LinearSolveSummary:
  val Empty: LinearSolveSummary =
    LinearSolveSummary(0, 0, 0L, 0L, None)

/** An executable solve capability for one fixed square linear system.
  *
  * This deliberately does not extend `DoubleLinearOperator`: applying a matrix
  * and solving a system are different operations with different failure and
  * work-accounting contracts. Implementations receive a snapshot of the right
  * hand side. [[solve]] validates and snapshots the returned solution, so a
  * later factor/workspace reuse cannot mutate a published result.
  *
  * Structural failures are `Left(LinAlgError)`. Iterative non-convergence is a
  * successful call carrying the best iterate with
  * `diagnostics.converged == false`; an outer algorithm decides whether that
  * condition is fatal.
  */
trait LinearSolveOperator:
  def size: Int
  def kind: LinearSolveKind

  protected def solveImpl(rhs: DVec): Either[LinAlgError, LinearSolveResult]

  final def solve(rhs: DVec): Either[LinAlgError, LinearSolveResult] =
    if size < 0 then
      Left(LinAlgError.InvalidArgument(s"linear solve size must be non-negative, got $size"))
    else if rhs.length != size then
      Left(LinAlgError.VectorLengthMismatch(size, rhs.length))
    else
      try
        solveImpl(rhs.copy).flatMap(validateResult)
      catch case error: LinAlgError => Left(error)

  private def validateResult(result: LinearSolveResult): Either[LinAlgError, LinearSolveResult] =
    val diagnostics = result.diagnostics
    if result.solution.length != size then
      Left(LinAlgError.VectorLengthMismatch(size, result.solution.length))
    else if diagnostics.iterations < 0 then
      Left(LinAlgError.InvalidArgument(s"linear solve reported negative iterations: ${diagnostics.iterations}"))
    else if diagnostics.operatorApplications < 0L then
      Left(
        LinAlgError.InvalidArgument(
          s"linear solve reported negative operator applications: ${diagnostics.operatorApplications}"
        )
      )
    else if diagnostics.residualNorm.exists(value => !value.isFinite || value < 0.0) then
      Left(
        LinAlgError.InvalidArgument(
          s"linear solve residual must be finite and non-negative when measured, got ${diagnostics.residualNorm}"
        )
      )
    else if !diagnostics.converged && diagnostics.residualNorm.isEmpty then
      Left(LinAlgError.InvalidArgument("a non-converged linear solve must report its residual norm"))
    else
      var i = 0
      while i < result.solution.length do
        if !result.solution(i).isFinite then
          return Left(LinAlgError.InvalidArgument(s"linear solve returned a non-finite solution entry at $i"))
        i += 1
      Right(result.copy(solution = result.solution.copy))

object LinearSolveOperator:

  /** Adapt a caller-created reusable factor. No matrix is factorized here. */
  def direct(factor: ExactSolveFactor): LinearSolveOperator =
    new LinearSolveOperator:
      val size: Int = factor.size
      val kind: LinearSolveKind = LinearSolveKind.DirectFactor

      protected def solveImpl(rhs: DVec): Either[LinAlgError, LinearSolveResult] =
        factor.solve(rhs).map: solution =>
          LinearSolveResult(
            solution,
            LinearSolveDiagnostics(
              converged = true,
              iterations = 0,
              residualNorm = None,
              operatorApplications = 0L
            )
          )

  /** Build a repeated conjugate-gradient solve for a caller-asserted SPD system.
    *
    * Every call owns its CG workspace. The application count includes CG's
    * initial residual evaluation and every Krylov step. A custom preconditioner
    * remains responsible for its own concurrency contract.
    */
  def conjugateGradient[A <: DoubleLinearOperator](
      operator: PositiveDefinite[A],
      config: SolverConfig = SolverConfig(),
      preconditioner: Preconditioner = Preconditioner.Identity,
      toleranceMode: ToleranceMode = ToleranceMode.RelativeToRhs
  ): Either[LinAlgError, LinearSolveOperator] =
    if operator.rows != operator.cols then
      Left(LinAlgError.NonSquareMatrix(Shape(Rows(operator.rows), Cols(operator.cols))))
    else if !config.tolerance.isFinite || config.tolerance < 0.0 then
      Left(
        LinAlgError.InvalidArgument(
          s"linear solve tolerance must be finite and non-negative, got ${config.tolerance}"
        )
      )
    else if config.maxIterations < 0 then
      Left(
        LinAlgError.InvalidArgument(
          s"linear solve maxIterations must be non-negative, got ${config.maxIterations}"
        )
      )
    else
      Right(
        new LinearSolveOperator:
          val size: Int = operator.rows
          val kind: LinearSolveKind = LinearSolveKind.ConjugateGradient

          protected def solveImpl(rhs: DVec): Either[LinAlgError, LinearSolveResult] =
            var applications = 0L
            val counted = new DoubleLinearOperator:
              val rows: Int = operator.rows
              val cols: Int = operator.cols

              def applyTo(x: DVec, into: MutableDVec): Unit =
                applications += 1L
                operator.applyTo(x, into)

            val result =
              IterativeSolvers.cg(
                counted,
                rhs,
                config,
                preconditioner,
                initial = None,
                toleranceMode = toleranceMode
              )
            Right(
              LinearSolveResult(
                result.x,
                LinearSolveDiagnostics(
                  converged = result.converged,
                  iterations = result.iterations,
                  residualNorm = Some(result.residual),
                  operatorApplications = applications
                )
              )
            )
      )

  /** Build an optional-backend solve while retaining gale's validation and
    * ownership boundary around every returned result.
    */
  def backendProvided(
      systemSize: Int
  )(
      solveOne: DVec => Either[LinAlgError, LinearSolveResult]
  ): Either[LinAlgError, LinearSolveOperator] =
    if systemSize < 0 then
      Left(LinAlgError.InvalidArgument(s"linear solve size must be non-negative, got $systemSize"))
    else
      Right(
        new LinearSolveOperator:
          val size: Int = systemSize
          val kind: LinearSolveKind = LinearSolveKind.BackendProvided

          protected def solveImpl(rhs: DVec): Either[LinAlgError, LinearSolveResult] =
            solveOne(rhs)
      )

/** A solve explicitly bound to the positive-definite metric it inverts.
  *
  * If `A` is symmetric and every solve is the action of `B^-1`, then
  * `T = B^-1 A` is self-adjoint in the `B` inner product:
  * `<x, T y>_B = x^T A y = (A x)^T y = <T x, y>_B`.
  * The wrapper records the identity of that metric contract and rejects a
  * dimensionally incompatible solver; it never forms `B^-1`.
  */
final class MetricSolveOperator[B <: DoubleLinearOperator] private (
    val metric: PositiveDefinite[B],
    val solver: LinearSolveOperator
):
  def size: Int = solver.size

  def solve(rhs: DVec): Either[LinAlgError, LinearSolveResult] =
    solver.solve(rhs)

object MetricSolveOperator:
  def bind[B <: DoubleLinearOperator](
      metric: PositiveDefinite[B],
      solver: LinearSolveOperator
  ): Either[LinAlgError, MetricSolveOperator[B]] =
    if metric.rows != metric.cols then
      Left(LinAlgError.NonSquareMatrix(Shape(Rows(metric.rows), Cols(metric.cols))))
    else if solver.size != metric.rows then
      Left(LinAlgError.VectorLengthMismatch(metric.rows, solver.size))
    else Right(new MetricSolveOperator(metric, solver))
