package gale.optim

import gale.linalg.DMat

/** Backtracking accelerated proximal gradient for a convex smooth objective and convex term. */
final case class AcceleratedConfig(
    maxIterations: Int = 10000,
    tolerance: FirstOrderTolerance = FirstOrderTolerance.strict,
    initialLipschitz: Double = 1.0,
    backtrackingFactor: Double = 2.0,
    maxBacktracks: Int = 40,
    restart: Boolean = true,
    control: SolverControl = SolverControl()
)

object AcceleratedProximal:
  def minimize(
      objective: DifferentiableObjective,
      term: ProximalTerm,
      initial: DMat,
      config: AcceleratedConfig = AcceleratedConfig()
  ): Either[FirstOrderError, FirstOrderSolution] =
    import OptimNumerics.*
    if config.maxIterations < 0 || config.maxBacktracks <= 0 ||
      !config.initialLipschitz.isFinite || config.initialLipschitz <= 0.0 ||
      !config.backtrackingFactor.isFinite || config.backtrackingFactor <= 1.0
    then
      Left(
        FirstOrderError.InvalidConfiguration(
          "acceleration needs nonnegative iterations, positive finite Lipschitz estimate, backtracking factor > 1 and positive backtracking budget"
        )
      )
    else if objective.variableRows != term.variableRows then
      Left(FirstOrderError.ShapeMismatch("proximal term", objective.variableRows, term.variableRows))
    else
      for
        _ <- validate(initial, objective.variableRows)
        _ <- config.control.validate
        result <-
          val execution = new OptimizationExecution(config.control)
          execution.finish(solve(objective, execution.term(term), initial, config, execution))
      yield result

  private def solve(
      objective: DifferentiableObjective,
      term: ProximalTerm,
      initial: DMat,
      config: AcceleratedConfig,
      execution: OptimizationExecution
  ): Either[FirstOrderError, FirstOrderSolution] =
    import OptimNumerics.*
    var point = initial
    var evaluation = execution.evaluate(objective, point) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    var termValue = term.value(point).flatMap(value => finiteScalar(value, "proximal objective").map(_ => value)) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    var total = evaluation.value + termValue
    if !total.isFinite then return Left(FirstOrderError.NumericalFailure("combined objective overflow"))
    var extrapolated = point
    var momentum = 1.0
    var lipschitz = config.initialLipschitz
    var iteration = 0
    var reference = 1.0
    var change = 0.0
    var stoppedMoving = false
    val settings = FirstOrderSettings(
      FirstOrderMethod.AcceleratedProximalGradient,
      config.maxIterations,
      config.tolerance,
      1.0,
      1.0,
      algorithm = AlgorithmSettings
        .Accelerated(config.initialLipschitz, config.backtrackingFactor, config.maxBacktracks, config.restart),
      maxEvaluations = config.control.maxEvaluations
    )

    def proximal(at: DMat, gradient: DMat, step: Double): Either[FirstOrderError, DMat] =
      for
        trial <- affine(at, gradient, -step)
        next <- term.proximal(trial, step)
        _ <- like(next, at)
      yield next

    // Cached accepted evaluations are reused for the residual and next iteration.
    while true do
      val step = 1.0 / lipschitz
      if !step.isFinite || step <= 0.0 then
        return Left(FirstOrderError.NumericalFailure("accelerated step is not finite and positive"))
      // Fixed-point evidence needs a positive step, not another objective evaluation.
      // Cap it independently of the update to prevent huge steps diluting the residual.
      val verificationStep = Math.min(1.0, step)
      val verificationFixed = proximal(point, evaluation.gradient, verificationStep) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      val (residual, resolution) = mapping(point, verificationFixed, verificationStep) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      if iteration == 0 then reference = Math.max(1.0, residual)
      val threshold = config.tolerance.threshold(reference)
      if !threshold.isFinite then return Left(FirstOrderError.NumericalFailure("residual threshold overflow"))
      val status =
        if residual <= threshold && resolution <= threshold then FirstOrderStoppingStatus.Converged
        else if stoppedMoving || (step == verificationStep && identical(point, verificationFixed)) then
          FirstOrderStoppingStatus.NumericalStagnation
        else FirstOrderStoppingStatus.IterationLimit
      val latest = result(
        point,
        total,
        residual,
        iteration,
        status,
        settings.copy(primalResidualScale = reference),
        execution,
        step,
        resolution,
        change
      )
      execution.record(latest)
      if !execution.continue then return Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled))
      if status != FirstOrderStoppingStatus.IterationLimit || iteration >= config.maxIterations then
        return Right(latest)

      var base = extrapolated
      var baseEvaluation =
        if identical(base, point) then evaluation
        else
          execution.evaluate(objective, base) match
            case Left(error)  => return Left(error)
            case Right(value) => value
      var restarted = false
      var accepted = false
      var candidate = point
      var candidateEvaluation = evaluation
      var candidateTerm = termValue
      var trialCount = 0
      while !accepted && trialCount < config.maxBacktracks do
        val trialStep = 1.0 / lipschitz
        if !trialStep.isFinite || trialStep <= 0.0 then
          return Left(FirstOrderError.NumericalFailure("backtracking step underflow"))
        candidate = (if identical(base, point) && trialStep == verificationStep then Right(verificationFixed)
                     else proximal(base, baseEvaluation.gradient, trialStep)) match
          case Left(error)  => return Left(error)
          case Right(value) => value
        candidateEvaluation = execution.evaluate(objective, candidate) match
          case Left(error)  => return Left(error)
          case Right(value) => value
        val displacement = difference(candidate, base)
        val linear = dot(array(baseEvaluation.gradient), displacement)
        val quadratic = .5 * lipschitz * dot(displacement, displacement)
        val bound = baseEvaluation.value + linear + quadratic
        val slack = 8.0 * Math.ulp(Math.max(Math.abs(baseEvaluation.value), Math.abs(candidateEvaluation.value)))
        if !bound.isFinite || !slack.isFinite then
          return Left(FirstOrderError.NumericalFailure("majorization arithmetic overflow"))
        if candidateEvaluation.value <= bound + slack then
          candidateTerm =
            term.value(candidate).flatMap(value => finiteScalar(value, "proximal objective").map(_ => value)) match
              case Left(error)  => return Left(error)
              case Right(value) => value
          val candidateTotal = candidateEvaluation.value + candidateTerm
          if !candidateTotal.isFinite then return Left(FirstOrderError.NumericalFailure("combined objective overflow"))
          val totalSlack = 8.0 * Math.ulp(Math.max(Math.abs(total), Math.abs(candidateTotal)))
          if config.restart && !restarted && !identical(base, point) && candidateTotal > total + totalSlack then
            base = point
            baseEvaluation = evaluation
            momentum = 1.0
            restarted = true
          else accepted = true
        else
          lipschitz *= config.backtrackingFactor
          if !lipschitz.isFinite then return Left(FirstOrderError.NumericalFailure("Lipschitz estimate overflow"))
        trialCount += 1
      if !accepted then return Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed))
      // O'Donoghue-Candes gradient restart uses the accepted update, without a retry.
      // Normalize both differences before their dot product to keep this sign test finite.
      val update = difference(candidate, point)
      val descent = difference(base, candidate)
      var updateScale = 0.0
      var descentScale = 0.0
      var k = 0
      while k < update.length do
        updateScale = Math.max(updateScale, Math.abs(update(k)))
        descentScale = Math.max(descentScale, Math.abs(descent(k)))
        k += 1
      var restartDot = 0.0
      if config.restart && updateScale > 0.0 && descentScale > 0.0 then
        if !updateScale.isFinite || !descentScale.isFinite then
          return Left(FirstOrderError.NumericalFailure("restart displacement overflow"))
        k = 0
        while k < update.length do
          restartDot += (update(k) / updateScale) * (descent(k) / descentScale)
          k += 1
      val nextMomentum =
        if config.restart && restartDot > 0.0 then 1.0
        else .5 * (1.0 + Math.sqrt(1.0 + 4.0 * momentum * momentum))
      val weight = if nextMomentum == 1.0 then 0.0 else (momentum - 1.0) / nextMomentum
      extrapolated =
        if weight == 0.0 then candidate
        else DMat.tabulate(point.rows, point.cols)((r, c) => candidate(r, c) + weight * (candidate(r, c) - point(r, c)))
      stoppedMoving = identical(base, point) && identical(candidate, point)
      finiteMatrix(extrapolated) match
        case Left(error) => return Left(error)
        case _           => ()
      change = Math.abs(candidateEvaluation.value + candidateTerm - total)
      if !change.isFinite then return Left(FirstOrderError.NumericalFailure("objective change overflow"))
      point = candidate
      evaluation = candidateEvaluation
      termValue = candidateTerm
      total = evaluation.value + termValue
      momentum = nextMomentum
      iteration += 1
    Left(FirstOrderError.NumericalFailure("unreachable accelerated state"))
