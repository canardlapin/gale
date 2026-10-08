package gale.optim

import gale.linalg.{DMat, DVec}
import gale.numeric.ExactSum
import scala.util.control.NonFatal

/** Safeguarded Powell-Hestenes-Rockafellar augmented Lagrangian with L-BFGS/B inner solves. This is a local smooth
  * constrained solver. Convergence requires original-problem first-order diagnostics, never just inner convergence.
  */
object AugmentedLagrangian:
  import OptimNumerics.*

  private final case class Sample(point: DMat, objective: ObjectiveEvaluation, constraints: ConstraintEvaluation)

  /** Recover cancellation in the sum of rounded products before declaring stationarity. */
  private def lagrangianGradient(at: Sample, weights: Array[Double]): Array[Double] =
    tabulate(at.point.rows): j =>
      var sum = at.objective.gradient(j, 0)
      var correction = 0.0
      var magnitude = Math.abs(sum)
      var finiteProducts = true
      var i = 0
      while i < weights.length do
        val product = weights(i) * at.constraints.jacobian(i, j)
        finiteProducts &&= product.isFinite
        magnitude += Math.abs(product)
        val adjusted = product - correction
        val next = sum + adjusted
        correction = (next - sum) - adjusted
        sum = next
        i += 1
      if !finiteProducts then Double.NaN
      else if magnitude == 0.0 then 0.0
      else if Math.abs(sum) <= 32.0 * Math.ulp(1.0) * magnitude || !sum.isFinite then
        val exact = ExactSum.zero()
        val _ = exact.add(at.objective.gradient(j, 0))
        i = 0
        while i < weights.length do
          val _ = exact.add(weights(i) * at.constraints.jacobian(i, j))
          i += 1
        exact.value
      else sum

  private def finiteConstraints(values: DVec): Either[FirstOrderError, Unit] =
    var i = 0
    while i < values.length do
      if !values(i).isFinite then return Left(FirstOrderError.NumericalFailure("non-finite constraint value"))
      i += 1
    Right(())

  def minimize(
      objective: DifferentiableObjective,
      constraints: NonlinearConstraints,
      initial: DMat,
      bounds: Option[MatrixBoxBounds] = None,
      config: AugmentedLagrangianConfig = AugmentedLagrangianConfig(),
      initialMultipliers: Option[DVec] = None
  ): Either[FirstOrderError, AugmentedLagrangianResult] =
    val count = constraints.equalities.toLong + constraints.inequalities.toLong
    val positives = Seq(
      config.feasibilityTolerance,
      config.scaledFeasibilityTolerance,
      config.stationarityTolerance,
      config.complementarityTolerance,
      config.initialPenalty,
      config.maximumPenalty,
      config.multiplierBound,
      config.initialInnerTolerance
    )
    if constraints.variableRows != objective.variableRows || initial.cols != 1 || constraints.equalities < 0 ||
      constraints.inequalities < 0 || count > Int.MaxValue || config.maxIterations < 0 ||
      config.maxInnerIterations <= 0 || config.historySize <= 0 || config.maxLineSearch <= 0 ||
      positives.exists(v => !v.isFinite || v <= 0.0) || config.maximumPenalty < config.initialPenalty ||
      !config.penaltyGrowth.isFinite || config.penaltyGrowth <= 1.0 ||
      !config.feasibilityContraction.isFinite || config.feasibilityContraction <= 0.0 || config.feasibilityContraction >= 1.0 ||
      !config.innerToleranceContraction.isFinite || config.innerToleranceContraction <= 0.0 || config.innerToleranceContraction >= 1.0
    then Left(FirstOrderError.InvalidConfiguration("invalid augmented Lagrangian dimensions, limits or tolerances"))
    else
      val m = count.toInt
      val scales = config.constraintScales.fold(Array.fill(m)(1.0))(_.toArray)
      val multipliers = initialMultipliers.fold(new Array[Double](m))(_.toArray)
      if scales.length != m || scales.exists(v => !v.isFinite || v <= 0.0) || multipliers.length != m ||
        multipliers.indices.exists(i =>
          !multipliers(i).isFinite || (i >= constraints.equalities && multipliers(i) < 0.0)
        )
      then Left(FirstOrderError.InvalidConfiguration("invalid constraint scales or initial multipliers"))
      else
        for
          _ <- validate(initial, objective.variableRows)
          _ <- config.control.validate
          _ <-
            if multipliers.indices.forall(i => (multipliers(i) * scales(i)).isFinite) then Right(())
            else Left(FirstOrderError.InvalidConfiguration("scaled initial multipliers overflow"))
          _ <-
            if bounds.forall(_.contains(initial)) then Right(())
            else Left(FirstOrderError.InvalidConfiguration("initial point must satisfy the box bounds"))
          result <- solve(objective, constraints, initial, bounds, config, scales, multipliers)
        yield result

  private def solve(
      objective: DifferentiableObjective,
      constraints: NonlinearConstraints,
      initial: DMat,
      bounds: Option[MatrixBoxBounds],
      config: AugmentedLagrangianConfig,
      scales: Array[Double],
      originalMultipliers: Array[Double]
  ): Either[FirstOrderError, AugmentedLagrangianResult] =
    val execution = new OptimizationExecution(config.control)
    val n = initial.rows
    val m = scales.length
    val equalities = constraints.equalities
    val lower = bounds.fold(Array.fill(n)(Double.NegativeInfinity))(b => array(b.lower))
    val upper = bounds.fold(Array.fill(n)(Double.PositiveInfinity))(b => array(b.upper))
    var cache: Option[Sample] = None
    var iterations = 0
    var innerIterations = 0L
    var lastInnerStatus: Option[FirstOrderStoppingStatus] = None
    var penalty = config.initialPenalty
    var trace = Vector.empty[OptimizationProgress]
    var innerTrace = Vector.empty[AugmentedLagrangianIteration]
    val curvatureScale = if config.reuseCurvatureScale then Some(new CurvatureScale) else None

    def sample(at: DMat): Either[FirstOrderError, Sample] =
      cache match
        case Some(value) if identical(value.point, at) => Right(value)
        case _                                         =>
          for
            f <- execution.evaluate(objective, at)
            c <-
              if m == 0 then Right(ConstraintEvaluation(DVec.zeros(0), DMat.zeros(0, n)))
              else execution.invoke(EvaluationKind.Jacobian)(constraints.evaluate(at))
            _ <-
              if c.values.length == m && c.jacobian.rows == m && c.jacobian.cols == n then Right(())
              else Left(FirstOrderError.InvalidConfiguration("constraint values/Jacobian have incorrect dimensions"))
            _ <- finiteMatrix(c.jacobian)
            _ <- finiteConstraints(c.values)
          yield
            val value = Sample(at, f, c)
            cache = Some(value)
            value

    def diagnostics(at: Sample, multipliers: Array[Double]): Either[FirstOrderError, ConstrainedDiagnostics] =
      val gradient = lagrangianGradient(at, multipliers)
      var feasibility = 0.0
      var scaledFeasibility = 0.0
      var complementarity = 0.0
      var dualFeasibility = 0.0
      var i = 0
      while i < m do
        val c = at.constraints.values(i)
        val violation = if i < equalities then Math.abs(c) else Math.max(0.0, c)
        feasibility = Math.max(feasibility, violation)
        scaledFeasibility = Math.max(scaledFeasibility, violation / scales(i))
        if i >= equalities then
          complementarity = Math.max(complementarity, Math.abs(multipliers(i) * c))
          dualFeasibility = Math.max(dualFeasibility, -multipliers(i))
        i += 1
      val x = array(at.point)
      i = 0
      while i < n do
        val violation = Math.max(0.0, Math.max(lower(i) - x(i), x(i) - upper(i)))
        feasibility = Math.max(feasibility, violation)
        scaledFeasibility = Math.max(scaledFeasibility, violation)
        i += 1
      val pg = LBFGSB.projected(x, gradient, lower, upper)
      val stationarity = pg.foldLeft(0.0)((acc, v) => Math.max(acc, Math.abs(v)))
      if !allFinite(gradient) || !allFinite(
          Array(feasibility, scaledFeasibility, stationarity, complementarity, dualFeasibility)
        )
      then Left(FirstOrderError.NumericalFailure("constrained diagnostic overflow"))
      else Right(ConstrainedDiagnostics(feasibility, scaledFeasibility, stationarity, complementarity, dualFeasibility))

    var current = sample(initial) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    var multipliers = originalMultipliers.clone()
    var residuals = diagnostics(current, multipliers) match
      case Left(error)  => return Left(error)
      case Right(value) => value

    def result(status: AugmentedLagrangianStatus): AugmentedLagrangianResult =
      AugmentedLagrangianResult(
        current.point,
        current.objective.value,
        current.constraints.values,
        DVec.fromArray(multipliers),
        residuals,
        status,
        iterations,
        innerIterations,
        penalty,
        execution.counts,
        trace,
        config,
        lastInnerStatus,
        innerTrace
      )

    def stopped(error: FirstOrderError): Either[FirstOrderError, AugmentedLagrangianResult] = error match
      case FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.EvaluationLimit) =>
        Right(result(AugmentedLagrangianStatus.EvaluationLimit))
      case FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled) =>
        Right(result(AugmentedLagrangianStatus.Cancelled))
      case other => Left(other)

    def notify(): Either[FirstOrderError, Boolean] =
      try
        val event = OptimizationProgress(
          iterations,
          current.objective.value,
          Math.max(residuals.stationarity, Math.max(residuals.scaledFeasibility, residuals.complementarity)),
          execution.counts
        )
        if config.control.traceCapacity > 0 then trace = (trace :+ event).takeRight(config.control.traceCapacity)
        val proceed = config.control.progress.forall(_(event))
        Right(proceed && !config.control.cancelled())
      catch
        case NonFatal(error) =>
          Left(
            FirstOrderError.OracleFailure(
              "progress/cancellation callback",
              Option(error.getMessage).getOrElse(error.getClass.getName)
            )
          )

    var previousViolation = Double.PositiveInfinity
    var innerTolerance = Math.max(config.stationarityTolerance, config.initialInnerTolerance)
    while true do
      notify() match
        case Left(error)  => return Left(error)
        case Right(false) => return Right(result(AugmentedLagrangianStatus.Cancelled))
        case Right(true)  => ()
      if residuals.satisfies(config) then return Right(result(AugmentedLagrangianStatus.Converged))
      if iterations >= config.maxIterations then return Right(result(AugmentedLagrangianStatus.IterationLimit))

      // Safeguard only the estimates fed into the NEXT subproblem. Diagnostics use the actual update, not clipped values.
      if multipliers.indices.exists(i => !(multipliers(i) * scales(i)).isFinite) then
        return Left(FirstOrderError.NumericalFailure("scaled multipliers overflow"))
      val estimates = tabulate(m): i =>
        val value = multipliers(i) * scales(i)
        Math.max(if i < equalities then -config.multiplierBound else 0.0, Math.min(config.multiplierBound, value))
      val rho = penalty
      val beforeCounts = execution.counts
      val initialCurvature = curvatureScale.fold(1.0)(_.value)
      var augmentedEvaluations = 0L
      var approximateSteps = 0
      val recordApproximate: () => Unit = () => approximateSteps += 1
      def recordInner(solved: Option[FirstOrderSolution], checked: Boolean): Unit =
        if config.control.traceCapacity > 0 then
          val after = execution.counts
          innerTrace = (innerTrace :+ AugmentedLagrangianIteration(
            iterations + (if checked then 0 else 1),
            rho,
            innerTolerance,
            initialCurvature,
            curvatureScale.fold(1.0)(_.value),
            solved.fold(0)(_.certificate.iterations),
            augmentedEvaluations,
            after.values - beforeCounts.values,
            after.jacobians - beforeCounts.jacobians,
            solved.map(_.status),
            checked,
            approximateSteps
          )).takeRight(config.control.traceCapacity)
      val augmented = new DifferentiableObjective:
        val variableRows = n
        def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
          augmentedEvaluations += 1
          sample(at).flatMap: base =>
            val weights = new Array[Double](m)
            var value = base.objective.value
            var i = 0
            while i < m do
              val c = base.constraints.values(i) / scales(i)
              val shifted = estimates(i) + rho * c
              if !c.isFinite || !shifted.isFinite then value = Double.NaN
              val active = i < equalities || shifted > 0.0
              // Subtract the multiplier-only constant algebraically, avoiding the difference of two large squares.
              if active then value += c * (estimates(i) + 0.5 * rho * c)
              else value -= 0.5 * estimates(i) * (estimates(i) / rho)
              if active then weights(i) = shifted / scales(i)
              i += 1
            val gradient = lagrangianGradient(base, weights)
            if !value.isFinite || !allFinite(gradient) then
              Left(FirstOrderError.NumericalFailure("augmented objective overflow"))
            else Right(ObjectiveEvaluation(value, matrix(gradient, n, 1)))

      val tolerance = FirstOrderTolerance.from(innerTolerance, 0.0).toOption.get
      val innerControl = SolverControl(cancelled = config.control.cancelled)
      val inner = bounds match
        case Some(box) =>
          LBFGSB.minimizeWithScale(
            augmented,
            box,
            current.point,
            LBFGSBConfig(
              maxIterations = config.maxInnerIterations,
              tolerance = tolerance,
              historySize = config.historySize,
              maxLineSearch = config.maxLineSearch,
              control = innerControl
            ),
            curvatureScale,
            interpolate = true,
            approximateWolfe = config.roundoffAwareLineSearch,
            onApproximateAcceptance = recordApproximate
          )
        case None =>
          LBFGS.minimizeWithScale(
            augmented,
            current.point,
            LBFGSConfig(
              maxIterations = config.maxInnerIterations,
              tolerance = tolerance,
              historySize = config.historySize,
              maxLineSearch = config.maxLineSearch,
              control = innerControl
            ),
            curvatureScale,
            interpolate = true,
            approximateWolfe = config.roundoffAwareLineSearch,
            onApproximateAcceptance = recordApproximate
          )
      val solved = inner match
        case Left(error) =>
          recordInner(None, false)
          return stopped(error)
        case Right(value) => value
      innerIterations += solved.certificate.iterations.toLong
      lastInnerStatus = Some(solved.status)
      // A stopped inner solve may have discarded its last rejected trial. Preserve the previous checked outer point if
      // its returned endpoint is not cached and the global budget cannot evaluate it.
      val next = sample(solved.primal) match
        case Left(error) =>
          recordInner(Some(solved), false)
          return stopped(error)
        case Right(value) => value
      var violation = 0.0
      val updated = new Array[Double](m)
      var i = 0
      while i < m do
        val c = next.constraints.values(i) / scales(i)
        val shifted = estimates(i) + rho * c
        updated(i) = (if i < equalities then shifted else Math.max(0.0, shifted)) / scales(i)
        val v = if i < equalities then c else Math.max(c, -estimates(i) / rho)
        violation = Math.max(violation, Math.abs(v))
        i += 1
      if !allFinite(updated) || !violation.isFinite then
        return Left(FirstOrderError.NumericalFailure("multiplier update overflow"))
      val checked = diagnostics(next, updated) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      current = next
      multipliers = updated
      residuals = checked
      iterations += 1
      recordInner(Some(solved), true)

      solved.status match
        case FirstOrderStoppingStatus.Cancelled       => return Right(result(AugmentedLagrangianStatus.Cancelled))
        case FirstOrderStoppingStatus.EvaluationLimit => return Right(result(AugmentedLagrangianStatus.EvaluationLimit))
        case FirstOrderStoppingStatus.LineSearchFailed if !residuals.satisfies(config) =>
          return Right(result(AugmentedLagrangianStatus.LineSearchFailed))
        case FirstOrderStoppingStatus.IterationLimit if !residuals.satisfies(config) =>
          return Right(result(AugmentedLagrangianStatus.InnerIterationLimit))
        case FirstOrderStoppingStatus.NumericalStagnation if !residuals.satisfies(config) =>
          return Right(result(AugmentedLagrangianStatus.NumericalStagnation))
        case _ => ()
      if !residuals.satisfies(config) && violation > config.feasibilityContraction * previousViolation then
        if penalty >= config.maximumPenalty then return Right(result(AugmentedLagrangianStatus.PenaltyLimit))
        penalty =
          if penalty > config.maximumPenalty / config.penaltyGrowth then config.maximumPenalty
          else penalty * config.penaltyGrowth
      if penalty != rho then curvatureScale.foreach(_.value = 1.0)
      previousViolation = violation
      val nextTolerance =
        if config.adaptiveInnerTolerance && solved.certificate.iterations > 0 &&
          (residuals.feasibility > config.feasibilityTolerance ||
            residuals.scaledFeasibility > config.scaledFeasibilityTolerance ||
            residuals.complementarity > config.complementarityTolerance)
        then Math.min(innerTolerance, config.innerToleranceContraction * violation)
        else innerTolerance * config.innerToleranceContraction
      innerTolerance = Math.max(config.stationarityTolerance, nextTolerance)
    Left(FirstOrderError.NumericalFailure("unreachable augmented Lagrangian state"))
