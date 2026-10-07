package gale.optim

import gale.linalg.{DMat, DVec}

/** Value and analytic gradient evaluated together at the same matrix. */
final case class ObjectiveEvaluation(value: Double, gradient: DMat)

/** Smooth objective that does not require a global Lipschitz bound. */
trait DifferentiableObjective:
  def variableRows: Int
  def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation]

object DifferentiableObjective:
  def apply(
      rows: Int
  )(callback: DMat => Either[FirstOrderError, ObjectiveEvaluation]): Either[FirstOrderError, DifferentiableObjective] =
    if rows <= 0 then Left(FirstOrderError.InvalidConfiguration("objective rows must be positive"))
    else
      Right(new DifferentiableObjective:
        val variableRows = rows
        def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] = callback(at))

  /** Adapt an existing bound-based objective without changing its callbacks. */
  def fromSmooth(objective: SmoothObjective): DifferentiableObjective =
    new DifferentiableObjective:
      def variableRows: Int = objective.variableRows
      def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
        for
          value <- objective.value(at)
          gradient <- objective.gradient(at)
        yield ObjectiveEvaluation(value, gradient)

  /** Single-column adapter. No mutable storage is shared with a caller. */
  def vector(
      rows: Int
  )(callback: DVec => Either[FirstOrderError, (Double, DVec)]): Either[FirstOrderError, DifferentiableObjective] =
    apply(rows): at =>
      if at.cols != 1 then Left(FirstOrderError.InvalidConfiguration("vector objectives require one column"))
      else
        callback(at.col(0)).map: (value, gradient) =>
          ObjectiveEvaluation(value, DMat.tabulate(gradient.length, 1)((r, _) => gradient(r)))

object SmoothObjectives:
  def apply(rows: Int, bound: Double)(
      objective: DMat => Either[FirstOrderError, Double],
      derivative: DMat => Either[FirstOrderError, DMat]
  ): Either[FirstOrderError, SmoothObjective] =
    if rows <= 0 || !bound.isFinite || bound <= 0.0 then
      Left(
        FirstOrderError.InvalidConfiguration(
          "smooth objective requires positive rows and a finite positive Lipschitz bound"
        )
      )
    else
      Right(new SmoothObjective:
        val variableRows = rows
        val lipschitz = bound
        def value(at: DMat): Either[FirstOrderError, Double] = objective(at)
        def gradient(at: DMat): Either[FirstOrderError, DMat] = derivative(at))

final case class DirectionalGradientCheck(
    analytic: Double,
    finiteDifference: Double,
    absoluteError: Double,
    threshold: Double
):
  def passed: Boolean = absoluteError <= threshold

object GradientCheck:
  /** Central differences are a debugging diagnostic, not a proof or a solver gradient substitute. The direction is used
    * exactly as supplied.
    */
  def directional(
      objective: DifferentiableObjective,
      at: DMat,
      direction: DMat,
      step: Double = 1e-5,
      absoluteTolerance: Double = 1e-6,
      relativeTolerance: Double = 1e-5
  ): Either[FirstOrderError, DirectionalGradientCheck] =
    import OptimNumerics.*
    if !step.isFinite || step <= 0.0 || !absoluteTolerance.isFinite || absoluteTolerance < 0.0 ||
      !relativeTolerance.isFinite || relativeTolerance < 0.0
    then Left(FirstOrderError.InvalidConfiguration("gradient-check step/tolerances are invalid"))
    else
      val execution = new OptimizationExecution(SolverControl())
      for
        _ <- validate(at, objective.variableRows)
        _ <- like(direction, at)
        center <- execution.evaluate(objective, at)
        plus <- affine(at, direction, step)
        minus <- affine(at, direction, -step)
        _ <-
          if identical(plus, at) || identical(minus, at) then
            Left(FirstOrderError.NumericalFailure("gradient-check step is not representable"))
          else Right(())
        high <- execution.evaluate(objective, plus)
        low <- execution.evaluate(objective, minus)
        analytic = dot(center.gradient, direction)
        measured = (high.value - low.value) / (2.0 * step)
        threshold = absoluteTolerance + relativeTolerance * Math.max(Math.abs(analytic), Math.abs(measured))
        _ <- finiteScalar(analytic, "analytic directional derivative")
        _ <- finiteScalar(measured, "finite-difference derivative")
        _ <- finiteScalar(threshold, "gradient-check threshold")
      yield DirectionalGradientCheck(analytic, measured, Math.abs(analytic - measured), threshold)
