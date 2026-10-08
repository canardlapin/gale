package gale.optim

import gale.linalg.{DMat, DVec}

/** Constraint values and analytic Jacobian at the same point. Equality rows come first; all remaining rows mean c <= 0.
  * The Jacobian has one row per constraint and one column per parameter.
  */
final case class ConstraintEvaluation(values: DVec, jacobian: DMat)

/** Smooth constraints on a single-column parameter matrix. Callbacks must be defined outside the nonlinear feasible
  * set; only optional box bounds are enforced at every evaluation.
  */
trait NonlinearConstraints:
  def variableRows: Int
  def equalities: Int
  def inequalities: Int
  def evaluate(at: DMat): Either[FirstOrderError, ConstraintEvaluation]

object NonlinearConstraints:
  def apply(rows: Int, equalityCount: Int, inequalityCount: Int)(
      callback: DMat => Either[FirstOrderError, ConstraintEvaluation]
  ): Either[FirstOrderError, NonlinearConstraints] =
    if rows <= 0 || equalityCount < 0 || inequalityCount < 0 ||
      equalityCount.toLong + inequalityCount.toLong > Int.MaxValue
    then
      Left(FirstOrderError.InvalidConfiguration("constraint dimensions must be non-negative with positive parameters"))
    else
      Right(new NonlinearConstraints:
        val variableRows = rows
        val equalities = equalityCount
        val inequalities = inequalityCount
        def evaluate(at: DMat): Either[FirstOrderError, ConstraintEvaluation] = callback(at))

enum AugmentedLagrangianStatus:
  case Converged, IterationLimit, InnerIterationLimit, LineSearchFailed, EvaluationLimit, Cancelled, PenaltyLimit,
    NumericalStagnation

/** Absolute tolerances. Constraint values/Jacobian rows are divided by constraintScales (default 1); returned values
  * and multipliers retain original units. Stationarity and complementarity use the original objective and multiplier
  * units. roundoffAwareLineSearch permits ALM-only approximate Wolfe acceptance when observed and predicted value
  * changes are within eight ulps; strong curvature and final original-problem tolerances remain required.
  */
final case class AugmentedLagrangianConfig(
    maxIterations: Int = 100,
    maxInnerIterations: Int = 500,
    feasibilityTolerance: Double = 1e-8,
    scaledFeasibilityTolerance: Double = 1e-8,
    stationarityTolerance: Double = 1e-7,
    complementarityTolerance: Double = 1e-8,
    initialPenalty: Double = 10.0,
    penaltyGrowth: Double = 10.0,
    feasibilityContraction: Double = 0.5,
    maximumPenalty: Double = 1e12,
    multiplierBound: Double = 1e12,
    initialInnerTolerance: Double = 1e-2,
    innerToleranceContraction: Double = 0.1,
    historySize: Int = 10,
    maxLineSearch: Int = 40,
    constraintScales: Option[DVec] = None,
    control: SolverControl = SolverControl(),
    adaptiveInnerTolerance: Boolean = true,
    reuseCurvatureScale: Boolean = true,
    roundoffAwareLineSearch: Boolean = true
)

/** First-order diagnostics, not a certificate of a minimum or of infeasibility. Stationarity is the infinity norm of
  * the unit box-projected Lagrangian gradient. Complementarity is max |mu_i * g_i|; feasibility also includes box
  * violations.
  */
final case class ConstrainedDiagnostics(
    feasibility: Double,
    scaledFeasibility: Double,
    stationarity: Double,
    complementarity: Double,
    dualFeasibility: Double
):
  def satisfies(config: AugmentedLagrangianConfig): Boolean =
    feasibility <= config.feasibilityTolerance && scaledFeasibility <= config.scaledFeasibilityTolerance && stationarity <= config.stationarityTolerance &&
      complementarity <= config.complementarityTolerance && dualFeasibility == 0.0

/** One attempted inner solve. User callback counts exclude the initial outer sample and include any endpoint recheck.
  * Augmented evaluations also include cache hits. Retention uses control.traceCapacity; interrupted attempts are kept.
  * Scale fields describe the scalar carried between subproblems (always 1 when reuse is disabled), not inner history.
  * approximateSteps counts accepted inner steps that used the roundoff-bounded approximate-decrease criterion.
  */
final case class AugmentedLagrangianIteration(
    outerIteration: Int,
    penalty: Double,
    tolerance: Double,
    initialCurvatureScale: Double,
    finalCurvatureScale: Double,
    iterations: Int,
    augmentedEvaluations: Long,
    objectiveEvaluations: Long,
    constraintEvaluations: Long,
    status: Option[FirstOrderStoppingStatus],
    endpointChecked: Boolean,
    approximateSteps: Int = 0
)

final case class AugmentedLagrangianResult(
    primal: DMat,
    objective: Double,
    constraints: DVec,
    multipliers: DVec,
    diagnostics: ConstrainedDiagnostics,
    status: AugmentedLagrangianStatus,
    iterations: Int,
    innerIterations: Long,
    penalty: Double,
    evaluations: EvaluationCounts,
    trace: Vector[OptimizationProgress],
    settings: AugmentedLagrangianConfig,
    lastInnerStatus: Option[FirstOrderStoppingStatus],
    innerTrace: Vector[AugmentedLagrangianIteration] = Vector.empty
)
