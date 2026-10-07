package gale.optim

import gale.linalg.DMat

/** Bounded strong-Wolfe search used internally by smooth solvers.
  *
  * Bracketing doubles an effective step; zoom uses bisection. Bisection is deliberately conservative: it avoids
  * interpolation arithmetic becoming a second source of non-finite trial points on portable runtimes.
  */
private[optim] object StrongWolfe:
  final case class Accepted(point: DMat, evaluation: ObjectiveEvaluation, step: Double, directionalDerivative: Double)

  def search(
      objective: DifferentiableObjective,
      origin: DMat,
      originEvaluation: ObjectiveEvaluation,
      direction: DMat,
      execution: OptimizationExecution,
      c1: Double,
      c2: Double,
      maxTrials: Int,
      maximumStep: Double,
      feasibleBoundary: Boolean = false,
      bounds: Option[MatrixBoxBounds] = None
  ): Either[FirstOrderError, Accepted] =
    import OptimNumerics.*
    val initialDerivative = dot(originEvaluation.gradient, direction)
    if !initialDerivative.isFinite || initialDerivative >= 0.0 then
      Left(FirstOrderError.NumericalFailure("strong-Wolfe search requires a finite descent direction"))
    else
      def trial(step: Double): Either[FirstOrderError, Option[Accepted]] =
        if !step.isFinite || step <= 0.0 then Right(None)
        else
          val proposed = affine(origin, direction, step).map: raw =>
            bounds.fold(raw)(box =>
              DMat.tabulate(raw.rows, raw.cols)((r, c) =>
                Math.max(box.lower(r, c), Math.min(box.upper(r, c), raw(r, c)))
              )
            )
          proposed match
            case Left(_)                                  => Right(None)
            case Right(point) if identical(point, origin) => Right(None)
            case Right(point)                             =>
              execution.evaluate(objective, point) match
                case Right(evaluation) =>
                  val derivative = dot(evaluation.gradient, direction)
                  if evaluation.value.isFinite && derivative.isFinite then
                    Right(Some(Accepted(point, evaluation, step, derivative)))
                  else Right(None)
                case Left(error @ FirstOrderError.OracleFailure(_, _))                             => Left(error)
                case Left(error @ FirstOrderError.ExecutionStopped(_))                             => Left(error)
                case Left(_: FirstOrderError.NonFiniteValue | _: FirstOrderError.NumericalFailure) => Right(None)
                case Left(error)                                                                   => Left(error)

      def sufficient(candidate: Accepted): Boolean =
        candidate.evaluation.value <= originEvaluation.value + c1 * candidate.step * initialDerivative

      def curvature(candidate: Accepted): Boolean =
        Math.abs(candidate.directionalDerivative) <= -c2 * initialDerivative

      def zoom(low: Accepted, initialHigh: Double, used: Int): Either[FirstOrderError, Accepted] =
        var lower = low
        var upper = initialHigh
        var attempts = used
        while attempts < maxTrials do
          val width = upper - lower.step
          val step = lower.step + 0.5 * width
          if !width.isFinite || Math.abs(width) <= Math.ulp(
              Math.max(1.0, Math.abs(lower.step))
            ) || step == lower.step || step == upper
          then return Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))
          trial(step) match
            case Left(error)            => return Left(error)
            case Right(None)            => upper = step
            case Right(Some(candidate)) =>
              if sufficient(candidate) && curvature(candidate) then return Right(candidate)
              else if !sufficient(candidate) || candidate.evaluation.value > lower.evaluation.value then upper = step
              else
                if candidate.directionalDerivative * (upper - lower.step) >= 0.0 then upper = lower.step
                lower = candidate
          attempts += 1
        Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))

      var previous = Accepted(origin, originEvaluation, 0.0, initialDerivative)
      var step = Math.min(1.0, maximumStep)
      var attempt = 0
      while attempt < maxTrials do
        trial(step) match
          case Left(error)            => return Left(error)
          case Right(None)            => return zoom(previous, step, attempt + 1)
          case Right(Some(candidate)) =>
            if sufficient(candidate) && (curvature(candidate) || (feasibleBoundary && step == maximumStep)) then
              return Right(candidate)
            else if !sufficient(candidate) || (attempt > 0 && candidate.evaluation.value > previous.evaluation.value)
            then return zoom(previous, step, attempt + 1)
            else if candidate.directionalDerivative >= 0.0 then return zoom(candidate, previous.step, attempt + 1)
            else
              previous = candidate
              val next = Math.min(maximumStep, step * 2.0)
              if next <= step || !next.isFinite then
                return Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))
              step = next
        attempt += 1
      Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))
