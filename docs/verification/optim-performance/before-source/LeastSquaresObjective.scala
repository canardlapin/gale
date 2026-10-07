package gale.optim

import gale.linalg.{DMat, DVec}

/** Residual model for `0.5 * ||r(x)||²`. Parameters form one vector; the Jacobian is residuals by parameters.
  * Construction never evaluates user callbacks. Finite differences are explicitly selected and every residual call made
  * by the solver is included in the evaluation budget.
  */
final class LeastSquaresObjective private (
    val parameterCount: Int,
    val residualCount: Int,
    private[optim] val residualFunction: DVec => Either[FirstOrderError, DVec],
    private[optim] val jacobianFunction: Option[DVec => Either[FirstOrderError, DMat]],
    private[optim] val differenceStep: Double
):
  /** Fixed nonnegative observation weights: transforms residual/Jacobian rows by sqrt(weight). */
  def weighted(weights: DVec): Either[FirstOrderError, LeastSquaresObjective] =
    if weights.length != residualCount then
      Left(FirstOrderError.InvalidConfiguration("weight count must equal residual count"))
    else if (0 until weights.length).exists(i => !weights(i).isFinite || weights(i) < 0.0) then
      Left(FirstOrderError.InvalidConfiguration("least-squares weights must be finite and nonnegative"))
    else
      val roots = DVec.tabulate(weights.length)(i => Math.sqrt(weights(i)))
      Right(
        new LeastSquaresObjective(
          parameterCount,
          residualCount,
          x =>
            residualFunction(x).flatMap: r =>
              LeastSquaresObjective
                .checkResidual(r, residualCount)
                .map(_ => DVec.tabulate(residualCount)(i => roots(i) * r(i))),
          jacobianFunction.map(callback =>
            (x: DVec) =>
              callback(x).flatMap: j =>
                LeastSquaresObjective
                  .checkJacobian(j, residualCount, parameterCount)
                  .map(_ => DMat.tabulate(residualCount, parameterCount)((r, c) => roots(r) * j(r, c)))
          ),
          differenceStep
        )
      )

object LeastSquaresObjective:
  def apply(parameters: Int, residuals: Int)(
      residual: DVec => Either[FirstOrderError, DVec],
      jacobian: DVec => Either[FirstOrderError, DMat]
  ): Either[FirstOrderError, LeastSquaresObjective] =
    dimensions(parameters, residuals).map(_ =>
      new LeastSquaresObjective(parameters, residuals, residual, Some(jacobian), 0.0)
    )

  def finiteDifferences(parameters: Int, residuals: Int, relativeStep: Double = 6e-6)(
      residual: DVec => Either[FirstOrderError, DVec]
  ): Either[FirstOrderError, LeastSquaresObjective] =
    for
      _ <- dimensions(parameters, residuals)
      _ <-
        if relativeStep.isFinite && relativeStep > 0.0 then Right(())
        else Left(FirstOrderError.InvalidConfiguration("finite-difference step must be finite and positive"))
    yield new LeastSquaresObjective(parameters, residuals, residual, None, relativeStep)

  private def dimensions(parameters: Int, residuals: Int): Either[FirstOrderError, Unit] =
    if parameters <= 0 || residuals <= 0 then
      Left(FirstOrderError.InvalidConfiguration("least squares requires positive parameter and residual counts"))
    else Right(())

  private[optim] def checkResidual(value: DVec, size: Int): Either[FirstOrderError, Unit] =
    if value.length != size then Left(FirstOrderError.ShapeMismatch("residual", size, value.length))
    else
      var i = 0
      while i < size do
        if !value(i).isFinite then return Left(FirstOrderError.NonFiniteValue("residual", i, value(i)))
        i += 1
      Right(())

  private[optim] def checkJacobian(value: DMat, rows: Int, cols: Int): Either[FirstOrderError, Unit] =
    if value.rows != rows || value.cols != cols then
      Left(FirstOrderError.InvalidConfiguration(s"Jacobian must have shape ${rows}x${cols}"))
    else OptimNumerics.finiteMatrix(value)

  private[optim] def residual(
      model: LeastSquaresObjective,
      at: DVec,
      execution: OptimizationExecution
  ): Either[FirstOrderError, DVec] =
    execution
      .invoke(EvaluationKind.Residual)(model.residualFunction(at.copy))
      .flatMap: value =>
        checkResidual(value, model.residualCount).map(_ => value.copy)

  private[optim] def jacobian(
      model: LeastSquaresObjective,
      at: DVec,
      execution: OptimizationExecution
  ): Either[FirstOrderError, DMat] =
    model.jacobianFunction match
      case Some(callback) =>
        execution
          .invoke(EvaluationKind.Jacobian)(callback(at.copy))
          .flatMap: value =>
            checkJacobian(value, model.residualCount, model.parameterCount)
              .map(_ => DMat.tabulate(value.rows, value.cols)((r, c) => value(r, c)))
      case None =>
        val entries = new Array[Double](model.residualCount * model.parameterCount)
        var column = 0
        while column < model.parameterCount do
          val h = model.differenceStep * Math.max(1.0, Math.abs(at(column)))
          val high = at(column) + h
          val low = at(column) - h
          val width = high - low
          if !high.isFinite || !low.isFinite || !width.isFinite || width <= 0.0 || high == at(column) || low == at(
              column
            )
          then return Left(FirstOrderError.NumericalFailure("finite-difference perturbation is not representable"))
          val plus =
            residual(model, DVec.tabulate(at.length)(i => if i == column then high else at(i)), execution) match
              case Left(error)  => return Left(error)
              case Right(value) => value
          val minus =
            residual(model, DVec.tabulate(at.length)(i => if i == column then low else at(i)), execution) match
              case Left(error)  => return Left(error)
              case Right(value) => value
          var row = 0
          while row < model.residualCount do
            entries(row * model.parameterCount + column) = (plus(row) - minus(row)) / width
            row += 1
          column += 1
        val result = DMat.fromArrayRowMajor(model.residualCount, model.parameterCount, entries)
        OptimNumerics.finiteMatrix(result).map(_ => result)
