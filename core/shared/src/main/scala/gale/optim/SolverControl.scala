package gale.optim

import gale.linalg.{DMat, DoubleLinearOperator}
import scala.util.control.NonFatal

/** Counts attempted callbacks, including rejected trials and final checks. A fused call contributes one callback, one
  * value, and one gradient.
  */
final case class EvaluationCounts(
    callbacks: Int = 0,
    values: Int = 0,
    gradients: Int = 0,
    proximals: Int = 0,
    projections: Int = 0,
    forwards: Int = 0,
    adjoints: Int = 0,
    residuals: Int = 0,
    jacobians: Int = 0
)

final case class OptimizationProgress(
    iteration: Int,
    objective: Double,
    residual: Double,
    evaluations: EvaluationCounts
)

final case class SolverControl(
    maxEvaluations: Int = Int.MaxValue,
    cancelled: () => Boolean = () => false,
    progress: Option[OptimizationProgress => Boolean] = None,
    traceCapacity: Int = 0
):
  private[optim] def validate: Either[FirstOrderError, Unit] =
    if maxEvaluations <= 0 || traceCapacity < 0 then
      Left(FirstOrderError.InvalidConfiguration("evaluation budget must be positive and trace capacity non-negative"))
    else Right(())

enum AlgorithmSettings:
  case FixedStep
  case Accelerated(initialLipschitz: Double, backtrackingFactor: Double, maxBacktracks: Int, restart: Boolean)
  case LimitedMemory(historySize: Int, c1: Double, c2: Double, maxLineSearch: Int, maximumStep: Double)
  case BoundedLimitedMemory(historySize: Int, c1: Double, c2: Double, maxLineSearch: Int)
  case DampedLeastSquares(initialDamping: Double, acceptanceRatio: Double)

private[optim] enum EvaluationKind:
  case Value, Gradient, Fused, Proximal, Projection, Forward, Adjoint, Residual, Jacobian

private[optim] final class OptimizationExecution(val control: SolverControl):
  private var work = EvaluationCounts()
  private var history = Vector.empty[OptimizationProgress]
  private var stopped = false
  private var hookError: Option[FirstOrderError] = None
  var last: Option[FirstOrderSolution] = None
  def counts: EvaluationCounts = work
  def trace: Vector[OptimizationProgress] = history
  def continue: Boolean = !stopped

  def invoke[A](kind: EvaluationKind)(body: => Either[FirstOrderError, A]): Either[FirstOrderError, A] =
    try
      if stopped || control.cancelled() then Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled))
      else if work.callbacks >= control.maxEvaluations then
        Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.EvaluationLimit))
      else
        work = kind match
          case EvaluationKind.Value    => work.copy(callbacks = work.callbacks + 1, values = work.values + 1)
          case EvaluationKind.Gradient => work.copy(callbacks = work.callbacks + 1, gradients = work.gradients + 1)
          case EvaluationKind.Fused    =>
            work.copy(callbacks = work.callbacks + 1, values = work.values + 1, gradients = work.gradients + 1)
          case EvaluationKind.Proximal   => work.copy(callbacks = work.callbacks + 1, proximals = work.proximals + 1)
          case EvaluationKind.Projection =>
            work.copy(callbacks = work.callbacks + 1, projections = work.projections + 1)
          case EvaluationKind.Forward  => work.copy(callbacks = work.callbacks + 1, forwards = work.forwards + 1)
          case EvaluationKind.Adjoint  => work.copy(callbacks = work.callbacks + 1, adjoints = work.adjoints + 1)
          case EvaluationKind.Residual => work.copy(callbacks = work.callbacks + 1, residuals = work.residuals + 1)
          case EvaluationKind.Jacobian => work.copy(callbacks = work.callbacks + 1, jacobians = work.jacobians + 1)
        body
    catch
      case NonFatal(error) =>
        Left(FirstOrderError.OracleFailure("callback", Option(error.getMessage).getOrElse(error.getClass.getName)))

  def evaluate(objective: DifferentiableObjective, at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
    invoke(EvaluationKind.Fused)(objective.evaluate(at)).flatMap: result =>
      for
        _ <- OptimNumerics.finiteScalar(result.value, "objective")
        _ <- OptimNumerics.like(result.gradient, at)
      yield result

  def notify(iteration: Int, objective: Double, residual: Double): Boolean =
    try
      if control.progress.nonEmpty || control.traceCapacity > 0 then
        val event = OptimizationProgress(iteration, objective, residual, work)
        if control.traceCapacity > 0 then history = (history :+ event).takeRight(control.traceCapacity)
        if control.progress.exists(callback => !callback(event)) then stopped = true
      if control.cancelled() then stopped = true
      !stopped
    catch
      case NonFatal(error) =>
        hookError = Some(
          FirstOrderError.OracleFailure(
            "progress/cancellation callback",
            Option(error.getMessage).getOrElse(error.getClass.getName)
          )
        )
        stopped = true
        false

  def record(result: FirstOrderSolution): Unit =
    last = Some(result)
    notify(result.certificate.iterations, result.objective, result.certificate.primalResidual)
    ()

  def finish(result: Either[FirstOrderError, FirstOrderSolution]): Either[FirstOrderError, FirstOrderSolution] =
    val selected = result match
      case Left(FirstOrderError.ExecutionStopped(status)) =>
        last.map(_.copy(status = status)).toRight(FirstOrderError.ExecutionStopped(status))
      case other => other
    hookError match
      case Some(error) => Left(error)
      case None        => selected.map(solution => solution.copy(evaluations = work, trace = history))

  def smooth(source: SmoothObjective): SmoothObjective = new SmoothObjective:
    def variableRows: Int = source.variableRows
    def lipschitz: Double = source.lipschitz
    def value(at: DMat): Either[FirstOrderError, Double] = invoke(EvaluationKind.Value)(source.value(at))
    def gradient(at: DMat): Either[FirstOrderError, DMat] = invoke(EvaluationKind.Gradient)(source.gradient(at))

  def term(source: ProximalTerm): ProximalTerm = new ProximalTerm:
    def variableRows: Int = source.variableRows
    def value(at: DMat): Either[FirstOrderError, Double] = invoke(EvaluationKind.Value)(source.value(at))
    def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] =
      invoke(EvaluationKind.Proximal)(source.proximal(at, step))

  def projection(source: ProjectionSet): ProjectionSet = new ProjectionSet:
    def variableRows: Int = source.variableRows
    def project(at: DMat): Either[FirstOrderError, DMat] = invoke(EvaluationKind.Projection)(source.project(at))

  def primal(source: ProximalObjective): ProximalObjective = new ProximalObjective:
    def variableRows: Int = source.variableRows
    def value(at: DMat): Either[FirstOrderError, Double] = invoke(EvaluationKind.Value)(source.value(at))
    def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] =
      invoke(EvaluationKind.Proximal)(source.proximal(at, step))

  def functional(source: LinearCompositeFunctional): LinearCompositeFunctional = new LinearCompositeFunctional:
    def targetRows: Int = source.targetRows
    def value(at: DMat): Either[FirstOrderError, Double] = invoke(EvaluationKind.Value)(source.value(at))
    def proximalConjugate(at: DMat, step: Double): Either[FirstOrderError, DMat] =
      invoke(EvaluationKind.Proximal)(source.proximalConjugate(at, step))

  def forward(source: DoubleLinearOperator, at: DMat, adjoint: Boolean): Either[FirstOrderError, DMat] =
    invoke(if adjoint then EvaluationKind.Adjoint else EvaluationKind.Forward):
      (if adjoint then source.transposeApplyTo(at) else source.applyTo(at)).left
        .map(FirstOrderError.OperatorFailure("optimization operator", _))
