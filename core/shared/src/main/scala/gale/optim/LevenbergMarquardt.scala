package gale.optim

import gale.backend.Backend
import gale.linalg.{DMat, DVec, LinAlgError, QROptions, QRPivoting}
import gale.numeric.ExactSum

final case class LevenbergMarquardtConfig(
    maxIterations: Int = 200,
    tolerance: FirstOrderTolerance = FirstOrderTolerance.strict,
    initialDamping: Double = 1e-3,
    acceptanceRatio: Double = 1e-4,
    parameterScale: Option[DVec] = None,
    control: SolverControl = SolverControl()
)

/** All residual/Jacobian values include any fixed weights supplied by the objective. Stationarity is raw ||Jᵀr||∞; a
  * converged nonzero residual does not establish a root, identifiability, or a global minimum.
  */
final case class LeastSquaresSolution(
    solution: FirstOrderSolution,
    residuals: DVec,
    jacobian: DMat,
    damping: Double,
    rejectedSteps: Int,
    scaledGradientNorm: Double
):
  def parameters: DVec = solution.primal.col(0)
  def status: FirstOrderStoppingStatus = solution.status
  def objective: Double = solution.objective

/** Dense, unconstrained LM. Each damped step is solved by pivoted QR of the augmented scaled Jacobian; JᵀJ is never
  * formed. Iterations count proposed steps, including rejected trials. Only stationarity establishes convergence.
  */
object LevenbergMarquardt:
  private final case class Evaluation(residual: DVec, jacobian: DMat, value: Double, gradient: DMat)

  def minimize(
      objective: LeastSquaresObjective,
      initial: DVec,
      config: LevenbergMarquardtConfig = LevenbergMarquardtConfig()
  )(using Backend): Either[FirstOrderError, LeastSquaresSolution] =
    val scaleValid = config.parameterScale.forall(s =>
      s.length == objective.parameterCount && (0 until s.length).forall(i =>
        s(i).isFinite && s(i) > 0.0 && (1.0 / s(i)).isFinite
      )
    )
    if config.maxIterations < 0 || !config.initialDamping.isFinite || config.initialDamping <= 0.0 ||
      !config.acceptanceRatio.isFinite || config.acceptanceRatio < 0.0 || config.acceptanceRatio >= 1.0 || !scaleValid
    then Left(FirstOrderError.InvalidConfiguration("invalid LM iteration, damping, acceptance, or parameter scale"))
    else
      for
        _ <- LeastSquaresObjective.checkResidual(initial, objective.parameterCount)
        _ <- config.control.validate
        result <- solve(objective, initial.copy, config)
      yield result

  private def cost(residual: DVec): Either[FirstOrderError, Double] =
    val norm = residual.norm2
    val value = (0.5 * norm) * norm
    OptimNumerics.finiteScalar(value, "least-squares cost").map(_ => value)

  private def evaluate(residual: DVec, jacobian: DMat, value: Double): Either[FirstOrderError, Evaluation] =
    val gradient = DMat.tabulate(jacobian.cols, 1): (column, _) =>
      var sum = 0.0
      var correction = 0.0
      var magnitude = 0.0
      var row = 0
      while row < jacobian.rows do
        val product = jacobian(row, column) * residual(row)
        magnitude += Math.abs(product)
        val adjusted = product - correction
        val next = sum + adjusted
        correction = (next - sum) - adjusted
        sum = next
        row += 1
      // Resolve severe cancellation before certifying stationarity. Products themselves are rounded doubles.
      if Math.abs(sum) <= 32.0 * Math.ulp(1.0) * magnitude || !sum.isFinite then
        val exact = ExactSum.zero()
        row = 0
        while row < jacobian.rows do
          val _ = exact.add(jacobian(row, column) * residual(row))
          row += 1
        exact.value
      else sum
    OptimNumerics.finiteMatrix(gradient).map(_ => Evaluation(residual, jacobian, value, gradient))

  private def solve(model: LeastSquaresObjective, start: DVec, config: LevenbergMarquardtConfig)(using
      Backend
  ): Either[FirstOrderError, LeastSquaresSolution] =
    import OptimNumerics.*
    val execution = new OptimizationExecution(config.control)
    val first = for
      r <- LeastSquaresObjective.residual(model, start, execution)
      j <- LeastSquaresObjective.jacobian(model, start, execution)
      value <- cost(r).flatMap(value => evaluate(r, j, value))
    yield value
    var evaluation = first match
      case Left(error)  => return Left(error)
      case Right(value) => value
    var point = start
    val scale = new Array[Double](model.parameterCount)
    def updateScale(): Unit =
      var c = 0
      while c < scale.length do
        config.parameterScale match
          case Some(s) => scale(c) = 1.0 / s(c)
          case None    =>
            val norm = evaluation.jacobian.col(c).norm2
            scale(c) = Math.max(scale(c), if norm > 0.0 then norm else 1.0)
        c += 1
    updateScale()
    val reference = Math.max(1.0, normInf(evaluation.gradient))
    val threshold = config.tolerance.threshold(reference)
    if !threshold.isFinite then return Left(FirstOrderError.NumericalFailure("LM gradient tolerance overflow"))
    val settings = FirstOrderSettings(
      FirstOrderMethod.LevenbergMarquardt,
      config.maxIterations,
      config.tolerance,
      1.0,
      0.0,
      primalResidualScale = reference,
      algorithm = AlgorithmSettings.DampedLeastSquares(config.initialDamping, config.acceptanceRatio),
      maxEvaluations = config.control.maxEvaluations
    )
    var damping = config.initialDamping
    var growth = 2.0
    var iterations = 0
    var rejected = 0
    var stepNorm = 0.0
    var change = 0.0
    def current(status: FirstOrderStoppingStatus): FirstOrderSolution =
      result(
        DMat.tabulate(point.length, 1)((r, _) => point(r)),
        evaluation.value,
        normInf(evaluation.gradient),
        iterations,
        status,
        settings,
        execution,
        stepNorm,
        change = change
      )
    def finish(answer: Either[FirstOrderError, FirstOrderSolution]): Either[FirstOrderError, LeastSquaresSolution] =
      val coherent = answer match
        case Left(FirstOrderError.ExecutionStopped(status)) => Right(current(status))
        case other                                          => other
      execution
        .finish(coherent)
        .flatMap: solved =>
          var scaled = 0.0
          var i = 0
          while i < scale.length do
            scaled = Math.max(scaled, Math.abs(evaluation.gradient(i, 0) / scale(i)))
            i += 1
          finiteScalar(scaled, "scaled LM gradient").map(_ =>
            LeastSquaresSolution(solved, evaluation.residual, evaluation.jacobian, damping, rejected, scaled)
          )
    def reject(): Boolean =
      rejected += 1
      damping *= growth
      growth *= 2.0
      damping.isFinite && growth.isFinite
    while true do
      val status =
        if normInf(evaluation.gradient) <= threshold then FirstOrderStoppingStatus.Converged
        else FirstOrderStoppingStatus.IterationLimit
      val latest = current(status)
      execution.record(latest)
      if !execution.continue then
        return finish(Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled)))
      if status == FirstOrderStoppingStatus.Converged || iterations >= config.maxIterations then
        return finish(Right(latest))
      val m = model.residualCount
      val n = model.parameterCount
      val augmented = DMat.newBuilder(m + n, n)
      var row = 0
      while row < m do
        var column = 0
        while column < n do
          val entry = evaluation.jacobian(row, column) / scale(column)
          if !entry.isFinite then
            return finish(Left(FirstOrderError.NonFiniteValue("numerical value", row * n + column, entry)))
          augmented(row, column) = entry
          column += 1
        row += 1
      val rootDamping = Math.sqrt(damping)
      var column = 0
      while column < n do
        augmented(m + column, column) = rootDamping
        column += 1
      val rhs = DVec.tabulate(m + n)(r => if r < m then -evaluation.residual(r) else 0.0)
      val solved = augmented.consumeQR(QROptions(pivoting = QRPivoting.Column)).solveLeastSquares(rhs)
      iterations += 1
      solved match
        case Left(_: LinAlgError.RankDeficient) =>
          if !reject() then return finish(Left(FirstOrderError.NumericalFailure("LM damping overflow")))
        case Left(error) => return finish(Left(FirstOrderError.OperatorFailure("LM augmented QR", error)))
        case Right(q)    =>
          val trial = DVec.tabulate(n)(i => point(i) + q(i) / scale(i))
          LeastSquaresObjective.checkResidual(trial, n) match
            case Left(_) =>
              if !reject() then return finish(Left(FirstOrderError.NumericalFailure("LM trial and damping overflow")))
            case Right(_) =>
              val displacement = Array.tabulate(n)(i => trial(i) - point(i))
              stepNorm = norm(displacement)
              if !stepNorm.isFinite then
                return finish(Left(FirstOrderError.NumericalFailure("LM displacement overflow")))
              if stepNorm == 0.0 then return finish(Right(current(FirstOrderStoppingStatus.NumericalStagnation)))
              // Prediction uses the actual rounded displacement, not the unrounded QR solution.
              var linear = 0.0
              val image = new ScaledNorm
              var c = 0
              while c < n do
                linear += evaluation.gradient(c, 0) * displacement(c)
                c += 1
              var r = 0
              while r < m do
                var value = 0.0
                c = 0
                while c < n do
                  value += evaluation.jacobian(r, c) * displacement(c)
                  c += 1
                image.add(value)
                r += 1
              val imageNorm = image.result
              val predicted = -linear - (0.5 * imageNorm) * imageNorm
              if !predicted.isFinite || predicted <= 0.0 then
                if !reject() then
                  return finish(Left(FirstOrderError.NumericalFailure("LM model reduction is unresolved")))
              else
                val candidateResidual =
                  LeastSquaresObjective.residual(model, trial, execution).flatMap(r => cost(r).map(r -> _))
                candidateResidual match
                  case Left(_: FirstOrderError.NonFiniteValue) =>
                    if !reject() then return finish(Left(FirstOrderError.NumericalFailure("LM damping overflow")))
                  case Left(error)              => return finish(Left(error))
                  case Right((residual, value)) =>
                    val actual = evaluation.value - value
                    val ratio = actual / predicted
                    if actual > 0.0 && ratio >= config.acceptanceRatio then
                      val candidate = for
                        j <- LeastSquaresObjective.jacobian(model, trial, execution)
                        e <- evaluate(residual, j, value)
                      yield e
                      candidate match
                        case Left(error) => return finish(Left(error))
                        case Right(next) =>
                          point = trial
                          evaluation = next
                          change = actual
                          updateScale()
                          val factor = Math.max(1.0 / 3.0, 1.0 - Math.pow(2.0 * Math.min(ratio, 1.0) - 1.0, 3))
                          damping = Math.max(java.lang.Double.MIN_NORMAL, damping * factor)
                          growth = 2.0
                    else if !reject() then return finish(Left(FirstOrderError.NumericalFailure("LM damping overflow")))
    finish(Left(FirstOrderError.NumericalFailure("unreachable LM state")))
