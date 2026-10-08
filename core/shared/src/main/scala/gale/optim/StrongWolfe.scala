package gale.optim

import gale.linalg.DMat

/** Bounded strong-Wolfe search used internally by smooth solvers.
  *
  * Bracketing doubles an effective step. Zoom uses bisection by default. ALM inner solves opt into cached
  * derivative/value interpolation with interior guards and a forced midpoint after poor contraction; nonfinite
  * endpoints use bisection. ALM may also opt into a roundoff-bounded approximate-decrease test while retaining strong
  * curvature. Public smooth solvers use strict sufficient decrease.
  */
private[optim] object StrongWolfe:
  val ignoreApproximateAcceptance: () => Unit = () => ()

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
      bounds: Option[MatrixBoxBounds] = None,
      interpolate: Boolean = false,
      approximateWolfe: Boolean = false,
      onApproximateAcceptance: () => Unit = ignoreApproximateAcceptance
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

      // Hager-Zhang approximate Wolfe near value resolution, with the existing strong-curvature test retained.
      // Restrict both observed and predicted changes to a few ulps; no unit-sized floor or relative 1e-6 allowance.
      val valueResolution = 8.0 * Math.ulp(originEvaluation.value)
      def acceptApproximate(candidate: Accepted): Boolean =
        val acceptable = approximateWolfe && c1 < 0.5 &&
          Math.abs(candidate.evaluation.value - originEvaluation.value) <= valueResolution &&
          Math.abs(candidate.step * initialDerivative) <= valueResolution &&
          candidate.directionalDerivative >= c2 * initialDerivative &&
          candidate.directionalDerivative <= (2.0 * c1 - 1.0) * initialDerivative && curvature(candidate)
        if acceptable then onApproximateAcceptance()
        acceptable

      def zoom(
          low: Accepted,
          initialHigh: Double,
          highValue: Option[Accepted],
          used: Int
      ): Either[FirstOrderError, Accepted] =
        var lower = low
        var upper = initialHigh
        var upperValue = highValue
        var attempts = used
        var bisect = false
        while attempts < maxTrials do
          val width = upper - lower.step
          val midpoint = lower.step + 0.5 * width
          val interpolated =
            if !interpolate || bisect then midpoint
            else
              upperValue.fold(midpoint): high =>
                // Normalize derivatives so a finite secant does not overflow their difference.
                val scale = Math.max(Math.abs(lower.directionalDerivative), Math.abs(high.directionalDerivative))
                val a = lower.directionalDerivative / scale
                val b = high.directionalDerivative / scale
                val secant = -a / (b - a)
                val fraction =
                  if secant.isFinite && secant > 0.0 && secant < 1.0 then secant
                  else
                    val slope = (high.evaluation.value - lower.evaluation.value) / width
                    -lower.directionalDerivative / (2.0 * (slope - lower.directionalDerivative))
                val candidate = lower.step + fraction * width
                // A proposal near the high endpoint can waste the budget when constraint activity changes.
                // Keep tiny proposals near low (needed by stiff quadratics), but bisect the high half.
                if fraction.isFinite && fraction > 0.0 && fraction <= 0.5 && candidate.isFinite &&
                  candidate > Math.min(lower.step, upper) && candidate < Math.max(lower.step, upper)
                then candidate
                else midpoint
          val step = interpolated
          if !width.isFinite || Math.abs(width) <= Math.ulp(
              Math.max(1.0, Math.abs(lower.step))
            ) || step == lower.step || step == upper
          then return Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))
          trial(step) match
            case Left(error) => return Left(error)
            case Right(None) =>
              upper = step
              upperValue = None
            case Right(Some(candidate)) =>
              if (sufficient(candidate) && curvature(candidate)) || acceptApproximate(candidate) then
                return Right(candidate)
              else if !sufficient(candidate) || candidate.evaluation.value > lower.evaluation.value then
                upper = step
                upperValue = Some(candidate)
              else
                if candidate.directionalDerivative * (upper - lower.step) >= 0.0 then
                  upper = lower.step
                  upperValue = Some(lower)
                lower = candidate
          // Permit very small useful interpolants, but never let repeated endpoint-hugging trials stall zoom.
          bisect = Math.abs(upper - lower.step) > 0.66 * Math.abs(width)
          attempts += 1
        Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))

      var previous = Accepted(origin, originEvaluation, 0.0, initialDerivative)
      var step = Math.min(1.0, maximumStep)
      var attempt = 0
      while attempt < maxTrials do
        trial(step) match
          case Left(error)            => return Left(error)
          case Right(None)            => return zoom(previous, step, None, attempt + 1)
          case Right(Some(candidate)) =>
            if (sufficient(candidate) && (curvature(
                candidate
              ) || (feasibleBoundary && step == maximumStep))) || acceptApproximate(candidate)
            then return Right(candidate)
            else if !sufficient(candidate) || (attempt > 0 && candidate.evaluation.value > previous.evaluation.value)
            then return zoom(previous, step, Some(candidate), attempt + 1)
            else if candidate.directionalDerivative >= 0.0 then
              return zoom(candidate, previous.step, Some(previous), attempt + 1)
            else
              previous = candidate
              val next = Math.min(maximumStep, step * 2.0)
              if next <= step || !next.isFinite then
                return Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))
              step = next
        attempt += 1
      Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))
