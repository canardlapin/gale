package gale.optim

import gale.linalg.{DMat, DVec}

/** Configuration for the portable, unconstrained limited-memory BFGS solver. */
final case class LBFGSConfig(
    maxIterations: Int = 1000,
    tolerance: FirstOrderTolerance = FirstOrderTolerance.strict,
    historySize: Int = 10,
    c1: Double = 1e-4,
    c2: Double = .9,
    maxLineSearch: Int = 40,
    maximumStep: Double = 1e20,
    control: SolverControl = SolverControl()
)

object LBFGS:
  private final case class Pair(s: Array[Double], y: Array[Double], rho: Double)

  def minimize(
      objective: DifferentiableObjective,
      initial: DMat,
      config: LBFGSConfig = LBFGSConfig()
  ): Either[FirstOrderError, FirstOrderSolution] =
    import OptimNumerics.*
    validateConfig(config).flatMap: _ =>
      for
        _ <- validate(initial, objective.variableRows)
        _ <- config.control.validate
        result <- solve(objective, initial, config)
      yield result

  def minimize(
      objective: DifferentiableObjective,
      initial: DVec,
      config: LBFGSConfig
  ): Either[FirstOrderError, FirstOrderSolution] =
    minimize(objective, DMat.tabulate(initial.length, 1)((row, _) => initial(row)), config)

  private def validateConfig(config: LBFGSConfig): Either[FirstOrderError, Unit] =
    if config.maxIterations < 0 then
      Left(FirstOrderError.InvalidConfiguration("L-BFGS iterations must be non-negative"))
    else if config.historySize <= 0 then
      Left(FirstOrderError.InvalidConfiguration("L-BFGS history size must be positive"))
    else if !config.c1.isFinite || !config.c2.isFinite || config.c1 <= 0.0 || config.c1 >= config.c2 || config.c2 >= 1.0
    then Left(FirstOrderError.InvalidConfiguration("L-BFGS requires 0 < c1 < c2 < 1"))
    else if config.maxLineSearch <= 0 then
      Left(FirstOrderError.InvalidConfiguration("L-BFGS line-search limit must be positive"))
    else if !config.maximumStep.isFinite || config.maximumStep <= 0.0 then
      Left(FirstOrderError.InvalidConfiguration("L-BFGS maximum step must be finite and positive"))
    else Right(())

  private def solve(
      objective: DifferentiableObjective,
      initial: DMat,
      config: LBFGSConfig
  ): Either[FirstOrderError, FirstOrderSolution] =
    import OptimNumerics.*
    val execution = new OptimizationExecution(config.control)
    val settings = FirstOrderSettings(
      FirstOrderMethod.LBFGS,
      config.maxIterations,
      config.tolerance,
      1.0,
      0.0,
      algorithm = AlgorithmSettings
        .LimitedMemory(config.historySize, config.c1, config.c2, config.maxLineSearch, config.maximumStep),
      maxEvaluations = config.control.maxEvaluations
    )
    def solution(
        point: DMat,
        evaluation: ObjectiveEvaluation,
        iterations: Int,
        status: FirstOrderStoppingStatus,
        scale: Double,
        step: Double
    ): FirstOrderSolution =
      result(
        point,
        evaluation.value,
        normInf(evaluation.gradient),
        iterations,
        status,
        settings.copy(primalResidualScale = scale),
        execution,
        step
      )

    val initialEvaluation = execution.evaluate(objective, initial) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    val initialScale = Math.max(1.0, normInf(initialEvaluation.gradient))
    if !config.tolerance.threshold(initialScale).isFinite then
      return Left(FirstOrderError.NumericalFailure("gradient tolerance overflow"))
    var current = initial
    var currentEvaluation = initialEvaluation
    var acceptedStep = 0.0
    var iteration = 0
    var history = Vector.empty[Pair]
    var latest = solution(
      current,
      currentEvaluation,
      iteration,
      FirstOrderStoppingStatus.IterationLimit,
      initialScale,
      acceptedStep
    )
    execution.record(latest)
    while iteration <= config.maxIterations do
      if !execution.continue then
        return execution.finish(Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled)))
      val residual = normInf(currentEvaluation.gradient)
      if residual <= config.tolerance.threshold(initialScale) then
        val converged = solution(
          current,
          currentEvaluation,
          iteration,
          FirstOrderStoppingStatus.Converged,
          initialScale,
          acceptedStep
        )
        return execution.finish(Right(converged))
      if iteration >= config.maxIterations then return execution.finish(Right(latest))
      val directionArray = twoLoop(array(currentEvaluation.gradient), history)
      var direction = matrix(directionArray, current.rows, current.cols)
      val directionalDerivative = dot(currentEvaluation.gradient, direction)
      if !finiteMatrix(direction).isRight || !directionalDerivative.isFinite || directionalDerivative >= 0.0 then
        direction = matrix(array(currentEvaluation.gradient).map(-_), current.rows, current.cols)
        history = Vector.empty
      StrongWolfe.search(
        objective,
        current,
        currentEvaluation,
        direction,
        execution,
        config.c1,
        config.c2,
        config.maxLineSearch,
        config.maximumStep
      ) match
        case Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed)) =>
          val failed = latest.copy(status = FirstOrderStoppingStatus.LineSearchFailed)
          return execution.finish(Right(failed))
        case Left(error @ FirstOrderError.ExecutionStopped(_)) => return execution.finish(Left(error))
        case Left(error)                                       => return execution.finish(Left(error))
        case Right(found)                                      =>
          val s = difference(found.point, current)
          val y = difference(found.evaluation.gradient, currentEvaluation.gradient)
          val curvature = dot(s, y)
          val scale = Math.sqrt(Math.max(0.0, dot(s, s))) * Math.sqrt(Math.max(0.0, dot(y, y)))
          if curvature.isFinite && (1.0 / curvature).isFinite && scale.isFinite && curvature > 1e-12 * scale then
            history = (history :+ Pair(s, y, 1.0 / curvature)).takeRight(config.historySize)
          current = found.point
          currentEvaluation = found.evaluation
          acceptedStep = found.step
          iteration += 1
          latest = solution(
            current,
            currentEvaluation,
            iteration,
            FirstOrderStoppingStatus.IterationLimit,
            initialScale,
            acceptedStep
          )
          execution.record(latest)
    execution.finish(Right(latest))

  /** Standard O(mn) two-loop recursion for an inverse-Hessian times gradient. */
  private def twoLoop(gradient: Array[Double], history: Vector[Pair]): Array[Double] =
    val q = gradient.clone()
    val alpha = new Array[Double](history.size)
    var i = history.size - 1
    while i >= 0 do
      val pair = history(i)
      alpha(i) = pair.rho * OptimNumerics.dot(pair.s, q)
      var j = 0
      while j < q.length do
        q(j) -= alpha(i) * pair.y(j)
        j += 1
      i -= 1
    var gamma = 1.0
    history.lastOption.foreach: pair =>
      val yy = OptimNumerics.dot(pair.y, pair.y)
      val sy = OptimNumerics.dot(pair.s, pair.y)
      if yy.isFinite && yy > 0.0 && sy.isFinite && sy > 0.0 then gamma = sy / yy
    var j = 0
    while j < q.length do
      q(j) *= gamma
      j += 1
    i = 0
    while i < history.size do
      val pair = history(i)
      val beta = pair.rho * OptimNumerics.dot(pair.y, q)
      j = 0
      while j < q.length do
        q(j) += pair.s(j) * (alpha(i) - beta)
        j += 1
      i += 1
    j = 0
    while j < q.length do
      q(j) = -q(j)
      j += 1
    q
