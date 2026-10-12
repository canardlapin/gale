package gale.optim

import gale.linalg.{DMat, DVec}
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

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

/** Positive scalar inverse curvature shared only between related solves; correction pairs are never shared. */
private[optim] final class CurvatureScale:
  var value: Double = 1.0
  def update(sy: Double, yy: Double): Unit =
    val next = sy / yy
    if next.isFinite && next > 0.0 && (1.0 / next).isFinite then value = next

object LBFGS:
  /** One correction pair with its curvature `sy = sᵀy` (so `rho = 1 / sy`) and `yy = yᵀy`. */
  private final class Pair(val s: Array[Double], val y: Array[Double], val sy: Double, val yy: Double):
    val rho: Double = 1.0 / sy

  /** The newest `capacity` correction pairs, oldest first, in a ring. */
  private final class History(capacity: Int):
    private val pairs = new Array[Pair](capacity)
    private var start = 0
    private var count = 0
    def size: Int = count
    def apply(index: Int): Pair = pairs((start + index) % capacity)
    def clear(): Unit = count = 0
    def append(pair: Pair): Unit =
      if count < capacity then
        pairs((start + count) % capacity) = pair
        count += 1
      else
        pairs(start) = pair
        start = (start + 1) % capacity

  def minimize(
      objective: DifferentiableObjective,
      initial: DMat,
      config: LBFGSConfig = LBFGSConfig()
  ): Either[FirstOrderError, FirstOrderSolution] =
    minimizeWithScale(
      objective,
      initial,
      config,
      None,
      scaleSteepestStep = true,
      interpolate = false,
      approximateWolfe = false,
      onApproximateAcceptance = StrongWolfe.ignoreApproximateAcceptance
    )

  private[optim] def minimizeWithScale(
      objective: DifferentiableObjective,
      initial: DMat,
      config: LBFGSConfig,
      curvatureScale: Option[CurvatureScale],
      scaleSteepestStep: Boolean,
      interpolate: Boolean,
      approximateWolfe: Boolean,
      onApproximateAcceptance: () => Unit
  ): Either[FirstOrderError, FirstOrderSolution] =
    import OptimNumerics.*
    validateConfig(config).flatMap: _ =>
      for
        _ <- validate(initial, objective.variableRows)
        _ <- config.control.validate
        result <- solve(
          objective,
          initial,
          config,
          curvatureScale,
          scaleSteepestStep,
          interpolate,
          approximateWolfe,
          onApproximateAcceptance
        )
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
      config: LBFGSConfig,
      curvatureScale: Option[CurvatureScale],
      scaleSteepestStep: Boolean,
      interpolate: Boolean,
      approximateWolfe: Boolean,
      onApproximateAcceptance: () => Unit
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
        residual: Double,
        iterations: Int,
        status: FirstOrderStoppingStatus,
        scale: Double,
        step: Double
    ): FirstOrderSolution =
      result(
        point,
        evaluation.value,
        residual,
        iterations,
        status,
        settings.copy(primalResidualScale = scale),
        execution,
        step
      )

    val initialEvaluation = execution.evaluate(objective, initial) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    val initialResidual = normInf(initialEvaluation.gradient)
    val initialScale = Math.max(1.0, initialResidual)
    if !config.tolerance.threshold(initialScale).isFinite then
      return Left(FirstOrderError.NumericalFailure("gradient tolerance overflow"))
    var current = initial
    var currentEvaluation = initialEvaluation
    var residual = initialResidual
    var acceptedStep = 0.0
    var iteration = 0
    val history = new History(config.historySize)
    val alpha = new Array[Double](config.historySize)
    var latest = solution(
      current,
      currentEvaluation,
      residual,
      iteration,
      FirstOrderStoppingStatus.IterationLimit,
      initialScale,
      acceptedStep
    )
    execution.record(latest)
    while iteration <= config.maxIterations do
      if !execution.continue then
        return execution.finish(Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled)))
      if residual <= config.tolerance.threshold(initialScale) then
        val converged = solution(
          current,
          currentEvaluation,
          residual,
          iteration,
          FirstOrderStoppingStatus.Converged,
          initialScale,
          acceptedStep
        )
        return execution.finish(Right(converged))
      if iteration >= config.maxIterations then return execution.finish(Right(latest))
      val directionData = twoLoop(currentEvaluation.gradient, history, alpha, curvatureScale.fold(1.0)(_.value))
      var direction = DMat.fromDoubleArrayOwned(current.rows, current.cols, directionData)
      var directionalDerivative = dot(currentEvaluation.gradient, direction)
      if !finiteData(directionData) || !directionalDerivative.isFinite || directionalDerivative >= 0.0 then
        direction = matrix(array(currentEvaluation.gradient).map(-_), current.rows, current.cols)
        directionalDerivative = dot(currentEvaluation.gradient, direction)
        history.clear()
        curvatureScale.foreach(_.value = 1.0)
      // With no correction pairs and no carried curvature, the direction is the raw `-g`, so a unit step has the
      // gradient's units. Try `1 / ‖d‖₂` (Breeze `LBFGS.determineStepSize`, Liu-Nocedal `STP1`): a unit-length
      // displacement. Once a pair exists, `gamma = sᵀy / yᵀy` scales `d` and the unit step is the natural trial.
      val initialStep =
        if scaleSteepestStep && history.size == 0 && curvatureScale.isEmpty then
          val inverse = 1.0 / norm(array(direction))
          if inverse.isFinite && inverse > 0.0 then inverse else 1.0
        else 1.0
      StrongWolfe.searchFrom(
        objective,
        current,
        currentEvaluation,
        direction,
        directionalDerivative,
        execution,
        config.c1,
        config.c2,
        config.maxLineSearch,
        config.maximumStep,
        interpolate = interpolate,
        approximateWolfe = approximateWolfe,
        onApproximateAcceptance = onApproximateAcceptance,
        initialStep = initialStep
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
          val yy = dot(y, y)
          val scale = Math.sqrt(Math.max(0.0, dot(s, s))) * Math.sqrt(Math.max(0.0, yy))
          if curvature.isFinite && (1.0 / curvature).isFinite && scale.isFinite && curvature > 1e-12 * scale then
            history.append(Pair(s, y, curvature, yy))
            curvatureScale.foreach(_.update(curvature, yy))
          current = found.point
          currentEvaluation = found.evaluation
          residual = normInf(currentEvaluation.gradient)
          acceptedStep = found.step
          iteration += 1
          latest = solution(
            current,
            currentEvaluation,
            residual,
            iteration,
            FirstOrderStoppingStatus.IterationLimit,
            initialScale,
            acceptedStep
          )
          execution.record(latest)
    execution.finish(Right(latest))

  private def finiteData(values: DoubleArray): Boolean =
    var i = 0
    while i < values.length do
      if !values(i).isFinite then return false
      i += 1
    true

  /** Standard O(mn) two-loop recursion for the descent direction `-H g`, returned as owned row-major storage.
    *
    * Each pass over the iterate fuses one pair's update with the next pair's dot product (and the initial scaling and
    * final negation into the adjacent passes). Every element still receives the same rounded operations in the same
    * order as the unfused recursion `alpha_i = rho_i s_iᵀq; q -= alpha_i y_i` (newest first), `q *= gamma`,
    * `beta_i = rho_i y_iᵀq; q += (alpha_i - beta_i) s_i` (oldest first), `-q`, so the result is bit-identical to it.
    */
  private def twoLoop(gradient: DMat, history: History, alpha: Array[Double], initialCurvature: Double): DoubleArray =
    val q = OptimNumerics.array(gradient)
    val n = q.length
    val out = DoubleArray.alloc(n)
    val k = history.size
    var gamma = initialCurvature
    if k > 0 then
      val newest = history(k - 1)
      if newest.yy.isFinite && newest.yy > 0.0 && newest.sy.isFinite && newest.sy > 0.0 then
        gamma = newest.sy / newest.yy
    if k == 0 then
      var j = 0
      while j < n do
        out(j) = -(q(j) * gamma)
        j += 1
      return out
    // Newest to oldest; `projection` is s_iᵀq for the pair about to be applied.
    var projection = OptimNumerics.dot(history(k - 1).s, q)
    var i = k - 1
    while i >= 0 do
      val pair = history(i)
      val a = pair.rho * projection
      alpha(i) = a
      val y = pair.y
      // After the oldest pair, scale by gamma and start y_0ᵀq for the forward loop.
      val next = if i > 0 then history(i - 1).s else history(0).y
      val scale = i == 0
      var sum = 0.0
      var j = 0
      while j < n do
        var v = q(j) - a * y(j)
        if scale then v *= gamma
        q(j) = v
        sum += next(j) * v
        j += 1
      projection = sum
      i -= 1
    // Oldest to newest; `projection` is y_iᵀq for the pair about to be applied.
    i = 0
    while i < k do
      val pair = history(i)
      val c = alpha(i) - pair.rho * projection
      val s = pair.s
      if i < k - 1 then
        val next = history(i + 1).y
        var sum = 0.0
        var j = 0
        while j < n do
          val v = q(j) + s(j) * c
          q(j) = v
          sum += next(j) * v
          j += 1
        projection = sum
      else
        var j = 0
        while j < n do
          out(j) = -(q(j) + s(j) * c)
          j += 1
      i += 1
    out
