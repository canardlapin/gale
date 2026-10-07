package gale.optim

import gale.backend.Backend
import gale.linalg.{DMat, DVec}

/** Reusable dense smooth objectives with a single-column parameter contract. Inputs are copied on construction so later
  * mutation through unsafe aliases cannot change the objective seen by a solver.
  */
object DenseObjectives:
  def leastSquares(design: DMat, response: DVec)(using
      backend: Backend
  ): Either[FirstOrderError, DifferentiableObjective] =
    for
      _ <- validateDesign(design, response)
      copiedDesign = copy(design)
      copiedResponse = response.copy
    yield new DifferentiableObjective:
      val variableRows = copiedDesign.cols
      def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
        for
          weight <- singleColumn(at, variableRows)
          fitted = copiedDesign * weight
          residual = DVec.tabulate(copiedDesign.rows)(i => fitted(i) - copiedResponse(i))
          _ <- finite(residual, "least-squares residual")
          squared = dot(residual, residual)
          _ <- scalar(squared, "least-squares squared residual")
          rawGradient = copiedDesign.t * residual
          _ <- finite(rawGradient, "least-squares gradient")
          scale = 1.0 / copiedDesign.rows
          gradient = DMat.tabulate(variableRows, 1)((row, _) => rawGradient(row) * scale)
          _ <- OptimNumerics.finiteMatrix(gradient)
          value = 0.5 * squared * scale
          _ <- scalar(value, "least-squares objective")
        yield ObjectiveEvaluation(value, gradient)

  def logistic(design: DMat, labels: DVec, l2: Double = 0.0)(using
      backend: Backend
  ): Either[FirstOrderError, DifferentiableObjective] =
    for
      _ <- validateDesign(design, labels)
      _ <-
        if l2.isFinite && l2 >= 0.0 then Right(())
        else Left(FirstOrderError.InvalidConfiguration("logistic L2 weight must be finite and non-negative"))
      _ <- validateLabels(labels)
      copiedDesign = copy(design)
      copiedLabels = labels.copy
    yield new DifferentiableObjective:
      val variableRows = copiedDesign.cols
      def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
        for
          weight <- singleColumn(at, variableRows)
          fitted = copiedDesign * weight
          lossAndResidual <- logisticTerms(fitted, copiedLabels)
          rawGradient = copiedDesign.t * lossAndResidual._2
          _ <- finite(rawGradient, "logistic gradient")
          weightNorm = dot(weight, weight)
          _ <- scalar(weightNorm, "logistic weight norm")
          scale = 1.0 / copiedDesign.rows
          gradient = DMat.tabulate(variableRows, 1): (row, _) =>
            rawGradient(row) * scale + l2 * weight(row)
          _ <- OptimNumerics.finiteMatrix(gradient)
          value = lossAndResidual._1 * scale + 0.5 * l2 * weightNorm
          _ <- scalar(value, "logistic objective")
        yield ObjectiveEvaluation(value, gradient)

  private def validateDesign(design: DMat, response: DVec): Either[FirstOrderError, Unit] =
    if design.rows <= 0 || design.cols <= 0 then
      Left(FirstOrderError.InvalidConfiguration("dense design must have positive rows and columns"))
    else if response.length != design.rows then
      Left(FirstOrderError.InvalidConfiguration("dense response length must equal design rows"))
    else for _ <- OptimNumerics.finiteMatrix(design); _ <- finite(response, "dense response") yield ()

  private def validateLabels(labels: DVec): Either[FirstOrderError, Unit] =
    var i = 0
    while i < labels.length do
      if labels(i) != -1.0 && labels(i) != 1.0 then
        return Left(FirstOrderError.InvalidConfiguration("logistic labels must be exactly -1 or +1"))
      i += 1
    Right(())

  private def singleColumn(at: DMat, rows: Int): Either[FirstOrderError, DVec] =
    if at.cols != 1 then
      Left(FirstOrderError.InvalidConfiguration("dense objectives require a single parameter column"))
    else OptimNumerics.validate(at, rows).map(_ => at.col(0))

  private def logisticTerms(fitted: DVec, labels: DVec): Either[FirstOrderError, (Double, DVec)] =
    val residual = DVec.newBuilder(fitted.length)
    var loss = 0.0
    var i = 0
    while i < fitted.length do
      val margin = labels(i) * fitted(i)
      if !margin.isFinite then return Left(FirstOrderError.NonFiniteValue("logistic margin", i, margin))
      val tail = Math.exp(-Math.abs(margin))
      if margin >= 0.0 then
        loss += Math.log1p(tail)
        residual(i) = -labels(i) * tail / (1.0 + tail)
      else
        loss += -margin + Math.log1p(tail)
        residual(i) = -labels(i) / (1.0 + tail)
      if !loss.isFinite || !residual(i).isFinite then
        return Left(FirstOrderError.NumericalFailure("logistic loss arithmetic is non-finite"))
      i += 1
    Right(loss -> residual.result())

  private def copy(matrix: DMat): DMat =
    DMat.tabulate(matrix.rows, matrix.cols)((row, column) => matrix(row, column))

  private def finite(vector: DVec, context: String): Either[FirstOrderError, Unit] =
    var i = 0
    while i < vector.length do
      if !vector(i).isFinite then return Left(FirstOrderError.NonFiniteValue(context, i, vector(i)))
      i += 1
    Right(())

  private def scalar(value: Double, context: String): Either[FirstOrderError, Unit] =
    if value.isFinite then Right(()) else Left(FirstOrderError.NumericalFailure(s"$context is non-finite"))

  private def dot(left: DVec, right: DVec): Double =
    var value = 0.0
    var i = 0
    while i < left.length do
      value += left(i) * right(i)
      i += 1
    value
