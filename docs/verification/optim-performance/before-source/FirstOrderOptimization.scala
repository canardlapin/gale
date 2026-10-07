package gale.optim

import gale.linalg.DMat
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError

/** Portable optimization methods selected from the structure exposed by a downstream objective.
  */
enum FirstOrderMethod:
  case ProximalGradient
  case ProjectedGradient
  case SmoothCompositePrimalDual
  case LinearCompositePrimalDual
  case ExactLinearReduction
  case AcceleratedProximalGradient
  case LBFGS
  case LBFGSB
  case LevenbergMarquardt

enum SolverMethodRequest:
  case Automatic
  case Require(method: FirstOrderMethod)

enum FirstOrderError:
  case InvalidConfiguration(detail: String)
  case ShapeMismatch(context: String, expectedRows: Int, actualRows: Int)
  case NonFiniteValue(context: String, index: Int, value: Double)
  case OracleFailure(context: String, detail: String)
  case OperatorFailure(context: String, cause: LinAlgError)
  case MissingCapability(
      requested: Option[FirstOrderMethod],
      compatible: Vector[FirstOrderMethod],
      available: Set[FirstOrderMethod]
  )
  case InvalidReduction(detail: String)
  case NumericalFailure(detail: String)
  case ExecutionStopped(status: FirstOrderStoppingStatus)

  def message: String =
    this match
      case InvalidConfiguration(detail)                     => detail
      case ShapeMismatch(context, expectedRows, actualRows) =>
        s"$context expected $expectedRows rows, got $actualRows"
      case NonFiniteValue(context, index, value) =>
        s"$context value at linear index $index is not finite: $value"
      case OracleFailure(context, detail)                      => s"$context failed: $detail"
      case OperatorFailure(context, cause)                     => s"$context failed: ${cause.getMessage}"
      case MissingCapability(requested, compatible, available) =>
        val requestedText = requested.fold("automatic selection")(_.toString)
        s"$requestedText is unavailable; compatible methods are ${compatible.mkString(", ")}; " +
          s"available methods are ${available.toVector.sortBy(_.ordinal).mkString(", ")}"
      case InvalidReduction(detail) => detail
      case NumericalFailure(detail) => detail
      case ExecutionStopped(status) => s"execution stopped: $status"

final case class FirstOrderCapabilities private (methods: Set[FirstOrderMethod]):
  def supports(method: FirstOrderMethod): Boolean =
    methods.contains(method)

  def select(
      compatible: Vector[FirstOrderMethod],
      request: SolverMethodRequest
  ): Either[FirstOrderError, FirstOrderMethod] =
    val distinct = compatible.distinct
    request match
      case SolverMethodRequest.Require(method) =>
        if distinct.contains(method) && supports(method) then Right(method)
        else Left(FirstOrderError.MissingCapability(Some(method), distinct, methods))
      case SolverMethodRequest.Automatic =>
        FirstOrderCapabilities.preference
          .find(method => distinct.contains(method) && supports(method))
          .toRight(FirstOrderError.MissingCapability(None, distinct, methods))

object FirstOrderCapabilities:
  private val preference = Vector(
    FirstOrderMethod.ExactLinearReduction,
    FirstOrderMethod.LevenbergMarquardt,
    FirstOrderMethod.LBFGSB,
    FirstOrderMethod.LBFGS,
    FirstOrderMethod.AcceleratedProximalGradient,
    FirstOrderMethod.ProximalGradient,
    FirstOrderMethod.ProjectedGradient,
    FirstOrderMethod.SmoothCompositePrimalDual,
    FirstOrderMethod.LinearCompositePrimalDual
  )

  def from(methods: Set[FirstOrderMethod]): Either[FirstOrderError, FirstOrderCapabilities] =
    if methods.nonEmpty then Right(FirstOrderCapabilities(methods))
    else Left(FirstOrderError.InvalidConfiguration("solver capabilities must not be empty"))

  val portable: FirstOrderCapabilities =
    FirstOrderCapabilities(FirstOrderMethod.values.toSet)

final case class FirstOrderTolerance private (
    absolute: Double,
    relative: Double
):
  def threshold(scale: Double): Double =
    absolute + relative * Math.abs(scale)

object FirstOrderTolerance:
  def from(absolute: Double, relative: Double): Either[FirstOrderError, FirstOrderTolerance] =
    if !absolute.isFinite || absolute < 0.0 then
      Left(FirstOrderError.InvalidConfiguration(s"absolute tolerance must be finite and non-negative, got $absolute"))
    else if !relative.isFinite || relative < 0.0 then
      Left(FirstOrderError.InvalidConfiguration(s"relative tolerance must be finite and non-negative, got $relative"))
    else if absolute == 0.0 && relative == 0.0 then
      Left(FirstOrderError.InvalidConfiguration("at least one stopping tolerance must be positive"))
    else Right(FirstOrderTolerance(absolute, relative))

  val strict: FirstOrderTolerance =
    FirstOrderTolerance(1e-8, 1e-8)

final case class FirstOrderConfig private (
    maxIterations: Int,
    tolerance: FirstOrderTolerance,
    stepSafety: Double,
    extrapolation: Double,
    control: SolverControl = SolverControl()
)

object FirstOrderConfig:
  def from(
      maxIterations: Int,
      tolerance: FirstOrderTolerance,
      stepSafety: Double = 0.99,
      extrapolation: Double = 1.0,
      control: SolverControl = SolverControl()
  ): Either[FirstOrderError, FirstOrderConfig] =
    if maxIterations <= 0 then
      Left(FirstOrderError.InvalidConfiguration(s"iteration budget must be positive, got $maxIterations"))
    else if !stepSafety.isFinite || stepSafety <= 0.0 || stepSafety >= 1.0 then
      Left(FirstOrderError.InvalidConfiguration(s"step safety must be finite and in (0, 1), got $stepSafety"))
    else if !extrapolation.isFinite || extrapolation < 0.0 || extrapolation > 1.0 then
      Left(FirstOrderError.InvalidConfiguration(s"extrapolation must be finite and in [0, 1], got $extrapolation"))
    else control.validate.map(_ => FirstOrderConfig(maxIterations, tolerance, stepSafety, extrapolation, control))

  val portable: FirstOrderConfig =
    FirstOrderConfig(10000, FirstOrderTolerance.strict, 0.99, 1.0)

final case class ValueSummary(
    rows: Int,
    columns: Int,
    sum: Double,
    squaredNorm: Double,
    maxAbs: Double
)

object ValueSummary:
  def from(value: DMat): ValueSummary =
    var sum = 0.0
    var squaredNorm = 0.0
    var maxAbs = 0.0
    var row = 0
    while row < value.rows do
      var column = 0
      while column < value.cols do
        val current = value(row, column)
        sum += current
        squaredNorm += current * current
        maxAbs = Math.max(maxAbs, Math.abs(current))
        column += 1
      row += 1
    ValueSummary(value.rows, value.cols, sum, squaredNorm, maxAbs)

final case class FirstOrderSettings(
    method: FirstOrderMethod,
    maxIterations: Int,
    tolerance: FirstOrderTolerance,
    stepSafety: Double,
    extrapolation: Double,
    primalResidualScale: Double = 1.0,
    dualResidualScale: Double = 1.0,
    algorithm: AlgorithmSettings = AlgorithmSettings.FixedStep,
    maxEvaluations: Int = Int.MaxValue
)

/** A stopping certificate bound to the exact returned primal and dual values.
  *
  * Residuals describe the implemented first-order conditions under `settings`; they do not assert application-specific
  * statistical validity.
  */
final case class FirstOrderCertificate(
    primal: ValueSummary,
    dual: Option[ValueSummary],
    objective: Double,
    primalResidual: Double,
    dualResidual: Double,
    objectiveChange: Double,
    iterations: Int,
    settings: FirstOrderSettings,
    primalResolution: Double = 0.0,
    dualResolution: Double = 0.0,
    private[optim] val boundPrimal: Option[DMat] = None,
    private[optim] val boundDual: Option[DMat] = None
):
  require(objective.isFinite, "certified objective must be finite")
  require(primalResidual.isFinite && primalResidual >= 0.0, "primal residual must be finite and non-negative")
  require(dualResidual.isFinite && dualResidual >= 0.0, "dual residual must be finite and non-negative")
  require(objectiveChange.isFinite && objectiveChange >= 0.0, "objective change must be finite and non-negative")
  require(iterations >= 0, "iteration count must be non-negative")

  require(primalResolution.isFinite && primalResolution >= 0.0, "primal resolution must be finite and non-negative")
  require(dualResolution.isFinite && dualResolution >= 0.0, "dual resolution must be finite and non-negative")

  /** Exact content binding; summaries alone never establish identity. */
  def binds(primalValue: DMat, dualValue: Option[DMat]): Boolean =
    boundPrimal.exists(sameValues(_, primalValue)) && ((boundDual, dualValue) match
      case (None, None)              => true
      case (Some(left), Some(right)) => sameValues(left, right)
      case _                         => false)

enum FirstOrderStoppingStatus:
  case Converged
  case IterationLimit
  case NumericalStagnation
  case EvaluationLimit
  case Cancelled
  case LineSearchFailed

final case class FirstOrderSolution(
    primal: DMat,
    dual: Option[DMat],
    objective: Double,
    status: FirstOrderStoppingStatus,
    certificate: FirstOrderCertificate,
    evaluations: EvaluationCounts = EvaluationCounts(),
    primalStep: Double = 0.0,
    dualStep: Double = 0.0,
    trace: Vector[OptimizationProgress] = Vector.empty
):
  /** Verification mappings cap steps at one; iteration steps remain separately visible. */
  def primalResidualStep: Double =
    certificate.settings.method match
      case FirstOrderMethod.LBFGS | FirstOrderMethod.LevenbergMarquardt => 0.0
      case FirstOrderMethod.LBFGSB                                      => 1.0
      case _                                                            => Math.min(1.0, primalStep)
  def dualResidualStep: Double = Math.min(1.0, dualStep)

  require(certificate.binds(primal, dual), "solver certificate must bind the returned values")
  require(objective == certificate.objective, "solution and certificate objectives must agree")
  require(
    status != FirstOrderStoppingStatus.Converged ||
      (certificate.primalResidual <= certificate.settings.tolerance
        .threshold(certificate.settings.primalResidualScale) &&
        certificate.dualResidual <= certificate.settings.tolerance.threshold(certificate.settings.dualResidualScale) &&
        certificate.primalResolution <= certificate.settings.tolerance
          .threshold(certificate.settings.primalResidualScale) &&
        certificate.dualResolution <= certificate.settings.tolerance.threshold(certificate.settings.dualResidualScale)),
    "convergence requires resolved final-point residuals"
  )

/** Differentiable objective with a caller-supplied gradient Lipschitz upper bound. The entire matrix is one variable;
  * columns may be coupled. Values must be finite at every evaluated point. Convexity is a premise of convex convergence
  * guarantees, not a property verified by this interface.
  */
trait SmoothObjective:
  def variableRows: Int
  def lipschitz: Double
  def value(at: DMat): Either[FirstOrderError, Double]
  def gradient(at: DMat): Either[FirstOrderError, DMat]

/** Directly proximable term over a matrix of one or more parameter columns. */
trait ProximalTerm:
  def variableRows: Int
  def value(at: DMat): Either[FirstOrderError, Double]
  def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat]

/** Feasible set with an exact projection supplied by the caller. Convex sets admit convex convergence guarantees.
  * Nonconvex projections remain usable, but the result reports only local computed fixed-point evidence.
  */
trait ProjectionSet:
  def variableRows: Int
  def project(at: DMat): Either[FirstOrderError, DMat]

/** Complete proximable primal objective for linear-composite splitting. */
trait ProximalObjective:
  def variableRows: Int
  def value(at: DMat): Either[FirstOrderError, Double]
  def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat]

/** Functional `g` applied after a bounded linear operator.
  *
  * Primal-dual iteration uses the proximal map of the convex conjugate `g*`.
  */
trait LinearCompositeFunctional:
  def targetRows: Int
  def value(at: DMat): Either[FirstOrderError, Double]
  def proximalConjugate(at: DMat, step: Double): Either[FirstOrderError, DMat]

/** A linear operator paired with a caller-certified induced-norm upper bound.
  *
  * Gale validates the shape of the claim, but does not estimate or prove the bound. An underestimate invalidates the
  * primal-dual step-size premise.
  */
final case class BoundedLinearOperator private (
    linearOperator: DoubleLinearOperator,
    normUpperBound: Double
)

object BoundedLinearOperator:
  def from(
      linearOperator: DoubleLinearOperator,
      normUpperBound: Double
  ): Either[FirstOrderError, BoundedLinearOperator] =
    if linearOperator.rows < 0 || linearOperator.cols < 0 then
      Left(
        FirstOrderError.InvalidConfiguration(
          s"operator dimensions must be non-negative, got ${linearOperator.rows}x${linearOperator.cols}"
        )
      )
    else if !normUpperBound.isFinite || normUpperBound < 0.0 then
      Left(
        FirstOrderError.InvalidConfiguration(
          s"operator norm bound must be finite and non-negative, got $normUpperBound"
        )
      )
    else Right(BoundedLinearOperator(linearOperator, normUpperBound))

/** Portable proximal, projected, and primal-dual first-order solvers.
  *
  * Iteration exhaustion returns a successful [[FirstOrderSolution]] with [[FirstOrderStoppingStatus.IterationLimit]];
  * it is never relabeled as convergence.
  */
object FirstOrderSolvers:

  def proximalGradient(
      smooth: SmoothObjective,
      term: ProximalTerm,
      initial: DMat,
      config: FirstOrderConfig = FirstOrderConfig.portable
  ): Either[FirstOrderError, FirstOrderSolution] =
    controlled(config): execution =>
      proximalGradientImpl(execution.smooth(smooth), execution.term(term), initial, config, execution)

  def projectedGradient(
      smooth: SmoothObjective,
      feasible: ProjectionSet,
      initial: DMat,
      config: FirstOrderConfig = FirstOrderConfig.portable
  ): Either[FirstOrderError, FirstOrderSolution] =
    controlled(config): execution =>
      projectedGradientImpl(execution.smooth(smooth), execution.projection(feasible), initial, config, execution)

  def linearCompositePrimalDual(
      primalObjective: ProximalObjective,
      functional: LinearCompositeFunctional,
      operator: BoundedLinearOperator,
      initial: DMat,
      config: FirstOrderConfig = FirstOrderConfig.portable,
      initialDual: Option[DMat] = None
  ): Either[FirstOrderError, FirstOrderSolution] =
    controlled(config): execution =>
      linearCompositePrimalDualImpl(
        execution.primal(primalObjective),
        execution.functional(functional),
        operator,
        initial,
        config,
        execution,
        initialDual
      )

  def smoothCompositePrimalDual(
      smooth: SmoothObjective,
      direct: ProximalTerm,
      functional: LinearCompositeFunctional,
      operator: BoundedLinearOperator,
      initial: DMat,
      config: FirstOrderConfig = FirstOrderConfig.portable,
      initialDual: Option[DMat] = None
  ): Either[FirstOrderError, FirstOrderSolution] =
    controlled(config): execution =>
      smoothCompositePrimalDualImpl(
        execution.smooth(smooth),
        execution.term(direct),
        execution.functional(functional),
        operator,
        initial,
        config,
        execution,
        initialDual
      )

  private def controlled(config: FirstOrderConfig)(
      run: OptimizationExecution => Either[FirstOrderError, FirstOrderSolution]
  ): Either[FirstOrderError, FirstOrderSolution] =
    config.control.validate.flatMap: _ =>
      val execution = new OptimizationExecution(config.control)
      execution.finish(run(execution))

  private def proximalGradientImpl(
      smooth: SmoothObjective,
      term: ProximalTerm,
      initial: DMat,
      config: FirstOrderConfig,
      execution: OptimizationExecution
  ): Either[FirstOrderError, FirstOrderSolution] =
    if smooth.variableRows != term.variableRows then
      Left(FirstOrderError.ShapeMismatch("proximal term", smooth.variableRows, term.variableRows))
    else if !smooth.lipschitz.isFinite || smooth.lipschitz <= 0.0 then
      Left(
        FirstOrderError.InvalidConfiguration(
          s"smooth Lipschitz bound must be finite and positive, got ${smooth.lipschitz}"
        )
      )
    else
      for
        _ <- validateInitial(initial, smooth.variableRows)
        _ <- validateStep(config.stepSafety / smooth.lipschitz)
        initialObjective <- combinedValue(smooth, term, initial)
        result <- iterateSingle(
          FirstOrderMethod.ProximalGradient,
          initial,
          initialObjective,
          config,
          execution,
          config.stepSafety / smooth.lipschitz,
          (value, updateStep) =>
            for
              gradient <- smooth.gradient(value)
              _ <- validateLike("smooth gradient", gradient, value)
              trial <- affineInput(value, gradient, -updateStep)
              next <- term.proximal(trial, updateStep)
              _ <- validateLike("proximal result", next, value)
              objective <- combinedValue(smooth, term, next)
            yield next -> objective
        )
      yield result

  private def projectedGradientImpl(
      smooth: SmoothObjective,
      feasible: ProjectionSet,
      initial: DMat,
      config: FirstOrderConfig,
      execution: OptimizationExecution
  ): Either[FirstOrderError, FirstOrderSolution] =
    if smooth.variableRows != feasible.variableRows then
      Left(FirstOrderError.ShapeMismatch("projection set", smooth.variableRows, feasible.variableRows))
    else if !smooth.lipschitz.isFinite || smooth.lipschitz <= 0.0 then
      Left(
        FirstOrderError.InvalidConfiguration(
          s"smooth Lipschitz bound must be finite and positive, got ${smooth.lipschitz}"
        )
      )
    else
      val step = config.stepSafety / smooth.lipschitz
      for
        _ <- validateInitial(initial, smooth.variableRows)
        _ <- validateStep(config.stepSafety / smooth.lipschitz)
        initialObjective <- smooth.value(initial).flatMap(validateScalar("smooth objective", _))
        result <- iterateSingle(
          FirstOrderMethod.ProjectedGradient,
          initial,
          initialObjective,
          config,
          execution,
          step,
          (value, updateStep) =>
            for
              gradient <- smooth.gradient(value)
              _ <- validateLike("smooth gradient", gradient, value)
              trial <- affineInput(value, gradient, -updateStep)
              next <- feasible.project(trial)
              _ <- validateLike("projection result", next, value)
              objective <- smooth.value(next).flatMap(validateScalar("smooth objective", _))
            yield next -> objective
        )
      yield result

  private def linearCompositePrimalDualImpl(
      primalObjective: ProximalObjective,
      functional: LinearCompositeFunctional,
      operator: BoundedLinearOperator,
      initial: DMat,
      config: FirstOrderConfig,
      execution: OptimizationExecution,
      initialDual: Option[DMat]
  ): Either[FirstOrderError, FirstOrderSolution] =
    if config.extrapolation != 1.0 then
      Left(FirstOrderError.InvalidConfiguration("primal-dual convergence requires extrapolation = 1"))
    else if operator.linearOperator.cols != primalObjective.variableRows then
      Left(
        FirstOrderError
          .ShapeMismatch("linear-composite source", primalObjective.variableRows, operator.linearOperator.cols)
      )
    else if operator.linearOperator.rows != functional.targetRows then
      Left(
        FirstOrderError.ShapeMismatch("linear-composite target", functional.targetRows, operator.linearOperator.rows)
      )
    else
      for
        _ <- validateInitial(initial, primalObjective.variableRows)
        _ <- initialDual.fold[Either[FirstOrderError, Unit]](Right(()))(value =>
          if value.cols != initial.cols then Left(FirstOrderError.InvalidConfiguration("initial dual column mismatch"))
          else validateMatrix("initial dual", value, operator.linearOperator.rows)
        )
        mappedInitial <- forward(operator.linearOperator, initial, "linear-composite forward", execution, false)
        initialObjective <- compositeValue(primalObjective, functional, initial, mappedInitial)
        result <- iteratePrimalDual(
          primalObjective,
          functional,
          operator,
          initial,
          initialObjective,
          config,
          execution,
          initialDual
        )
      yield result

  /** Condat--Vu splitting for `f(x) + h(x) + g(Kx)`.
    *
    * `f` is differentiable with a certified Lipschitz bound, `h` has an exact proximal map, and `g` has an exact
    * conjugate proximal map. The chosen steps satisfy the coupled condition `1/tau - sigma * ||K||^2 > L/2` from the
    * supplied upper bound, including the zero-operator case. Only extrapolation 1 is admitted. Values must be finite at
    * all evaluated points; general extended-valued composite indicators are not supported.
    */
  private def smoothCompositePrimalDualImpl(
      smooth: SmoothObjective,
      direct: ProximalTerm,
      functional: LinearCompositeFunctional,
      operator: BoundedLinearOperator,
      initial: DMat,
      config: FirstOrderConfig,
      execution: OptimizationExecution,
      initialDual: Option[DMat]
  ): Either[FirstOrderError, FirstOrderSolution] =
    if config.extrapolation != 1.0 then
      Left(FirstOrderError.InvalidConfiguration("primal-dual convergence requires extrapolation = 1"))
    else if smooth.variableRows != direct.variableRows then
      Left(FirstOrderError.ShapeMismatch("direct proximal term", smooth.variableRows, direct.variableRows))
    else if operator.linearOperator.cols != smooth.variableRows then
      Left(FirstOrderError.ShapeMismatch("smooth-composite source", smooth.variableRows, operator.linearOperator.cols))
    else if operator.linearOperator.rows != functional.targetRows then
      Left(
        FirstOrderError.ShapeMismatch("smooth-composite target", functional.targetRows, operator.linearOperator.rows)
      )
    else if !smooth.lipschitz.isFinite || smooth.lipschitz <= 0.0 then
      Left(
        FirstOrderError.InvalidConfiguration(
          s"smooth Lipschitz bound must be finite and positive, got ${smooth.lipschitz}"
        )
      )
    else
      for
        _ <- validateInitial(initial, smooth.variableRows)
        _ <- validateStep(config.stepSafety / smooth.lipschitz)
        _ <- initialDual.fold[Either[FirstOrderError, Unit]](Right(()))(value =>
          if value.cols != initial.cols then Left(FirstOrderError.InvalidConfiguration("initial dual column mismatch"))
          else validateMatrix("initial dual", value, operator.linearOperator.rows)
        )
        mappedInitial <- forward(operator.linearOperator, initial, "smooth-composite forward", execution, false)
        initialObjective <- smoothCompositeValue(smooth, direct, functional, initial, mappedInitial)
        result <- iterateSmoothComposite(
          smooth,
          direct,
          functional,
          operator,
          initial,
          initialObjective,
          config,
          execution,
          initialDual
        )
      yield result

  private def iterateSingle(
      method: FirstOrderMethod,
      initial: DMat,
      initialObjective: Double,
      config: FirstOrderConfig,
      execution: OptimizationExecution,
      step: Double,
      update: (DMat, Double) => Either[FirstOrderError, (DMat, Double)]
  ): Either[FirstOrderError, FirstOrderSolution] =
    var current = initial
    var objective = initialObjective
    var objectiveChange = 0.0
    var iteration = 0
    val verificationStep = Math.min(1.0, step)
    var pending = update(current, step)
    def checkedMeasure(next: DMat): Either[FirstOrderError, (Double, Double)] =
      if verificationStep == step then residualMeasure(current, next, step)
      else
        update(current, verificationStep).flatMap((checked, _) => residualMeasure(current, checked, verificationStep))
    var measurement = pending.flatMap((next, _) => checkedMeasure(next))
    val initialMeasure = measurement
    val scaleValue = initialMeasure match
      case Left(error)          => return Left(error)
      case Right((residual, _)) => Math.max(1.0, residual)
    val threshold = config.tolerance.threshold(scaleValue)
    if !threshold.isFinite then return Left(FirstOrderError.NumericalFailure("residual tolerance overflow"))
    var running = true
    var residual = 0.0
    var resolution = 0.0
    var status = FirstOrderStoppingStatus.IterationLimit
    while running do
      pending match
        case Left(error)                  => return Left(error)
        case Right((next, nextObjective)) =>
          measurement match
            case Left(error)                    => return Left(error)
            case Right((measured, uncertainty)) =>
              residual = measured
              resolution = uncertainty
          execution.record(
            FirstOrderSolution(
              current,
              None,
              objective,
              FirstOrderStoppingStatus.IterationLimit,
              certificate(
                method,
                current,
                None,
                objective,
                residual,
                0.0,
                objectiveChange,
                iteration,
                config,
                scaleValue,
                1.0,
                resolution,
                0.0
              ),
              primalStep = step
            )
          )
          if !execution.continue then return Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled))
          if residual <= threshold && resolution <= threshold then
            status = FirstOrderStoppingStatus.Converged
            running = false
          else if iteration >= config.maxIterations then running = false
          else if sameValues(current, next) then
            status = FirstOrderStoppingStatus.NumericalStagnation
            running = false
          else
            objectiveChange = Math.abs(nextObjective - objective)
            if !objectiveChange.isFinite then return Left(FirstOrderError.NumericalFailure("objective change overflow"))
            current = next
            objective = nextObjective
            iteration += 1
            pending = update(current, step)
            measurement = pending.flatMap((next, _) => checkedMeasure(next))
    val currentCertificate = certificate(
      method,
      current,
      None,
      objective,
      residual,
      0.0,
      objectiveChange,
      iteration,
      config,
      scaleValue,
      1.0,
      resolution,
      0.0
    )
    Right(FirstOrderSolution(current, None, objective, status, currentCertificate, primalStep = step))

  private def iteratePrimalDual(
      primalObjective: ProximalObjective,
      functional: LinearCompositeFunctional,
      operator: BoundedLinearOperator,
      initial: DMat,
      initialObjective: Double,
      config: FirstOrderConfig,
      execution: OptimizationExecution,
      initialDual: Option[DMat]
  ): Either[FirstOrderError, FirstOrderSolution] =
    val denominator = Math.max(1.0, operator.normUpperBound)
    val primalStep = config.stepSafety / denominator
    val dualStep = config.stepSafety / denominator
    validateStep(primalStep) match
      case Left(error) => return Left(error)
      case _           => ()
    var primal = initial
    var extrapolated = initial
    var dual = initialDual.getOrElse(DMat.zeros(operator.linearOperator.rows, initial.cols))
    validateLike("initial dual", dual, DMat.zeros(operator.linearOperator.rows, initial.cols)) match
      case Left(error) => return Left(error)
      case _           => ()
    var objective = initialObjective
    var primalResidual = Double.MaxValue
    var dualResidual = Double.MaxValue
    var primalResolution = 0.0
    var dualResolution = 0.0
    var primalScale = 1.0
    var dualScale = 1.0
    var stagnated = false
    var objectiveChange = Double.MaxValue
    var iteration = 0
    var converged = false
    var error = Option.empty[FirstOrderError]
    while iteration < config.maxIterations && !converged && !stagnated && error.isEmpty do
      val updated =
        for
          mapped <- forward(operator.linearOperator, extrapolated, "linear-composite forward", execution, false)
          dualTrial <- affineInput(dual, mapped, dualStep)
          nextDual <- functional.proximalConjugate(dualTrial, dualStep)
          _ <- validateLike("dual proximal result", nextDual, dual)
          adjoint <- forward(operator.linearOperator.adjoint, nextDual, "linear-composite adjoint", execution, true)
          primalTrial <- affineInput(primal, adjoint, -primalStep)
          nextPrimal <- primalObjective.proximal(primalTrial, primalStep)
          _ <- validateLike("primal proximal result", nextPrimal, primal)
          mappedNext <- forward(operator.linearOperator, nextPrimal, "linear-composite forward", execution, false)
          nextObjective <- compositeValue(primalObjective, functional, nextPrimal, mappedNext)
          residuals <- fixedPointResiduals(
            primalObjective,
            functional,
            operator.linearOperator,
            nextPrimal,
            nextDual,
            Math.min(1.0, primalStep),
            Math.min(1.0, dualStep),
            execution,
            mappedNext
          )
        yield (nextPrimal, nextDual, nextObjective, residuals)
      updated match
        case Left(value)                                             => error = Some(value)
        case Right((nextPrimal, nextDual, nextObjective, residuals)) =>
          stagnated = sameValues(primal, nextPrimal) && sameValues(dual, nextDual)
          extrapolated = add(
            nextPrimal,
            scale(subtract(nextPrimal, primal), config.extrapolation)
          )
          primal = nextPrimal
          dual = nextDual
          primalResidual = residuals._1
          dualResidual = residuals._2
          primalResolution = residuals._3
          dualResolution = residuals._4
          objectiveChange = Math.abs(nextObjective - objective)
          objective = nextObjective
          iteration += 1
          if iteration == 1 then
            primalScale = Math.max(1.0, primalResidual)
            dualScale = Math.max(1.0, dualResidual)
          if !objectiveChange.isFinite || !config.tolerance.threshold(primalScale).isFinite ||
            !config.tolerance.threshold(dualScale).isFinite
          then error = Some(FirstOrderError.NumericalFailure("stopping diagnostic overflow"))
          converged = error.isEmpty && primalResidual <= config.tolerance.threshold(primalScale) &&
            dualResidual <= config.tolerance.threshold(dualScale) &&
            primalResolution <= config.tolerance.threshold(primalScale) &&
            dualResolution <= config.tolerance.threshold(dualScale) &&
            objectiveChange <= config.tolerance.threshold(Math.max(1.0, Math.abs(objective)))
          if error.isEmpty then
            execution.record(
              FirstOrderSolution(
                primal,
                Some(dual),
                objective,
                FirstOrderStoppingStatus.IterationLimit,
                certificate(
                  FirstOrderMethod.LinearCompositePrimalDual,
                  primal,
                  Some(dual),
                  objective,
                  primalResidual,
                  dualResidual,
                  objectiveChange,
                  iteration,
                  config,
                  primalScale,
                  dualScale,
                  primalResolution,
                  dualResolution
                ),
                primalStep = primalStep,
                dualStep = dualStep
              )
            )
            if !execution.continue then
              error = Some(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled))
    error match
      case Some(value) => Left(value)
      case None        =>
        val status =
          if converged then FirstOrderStoppingStatus.Converged
          else if stagnated then FirstOrderStoppingStatus.NumericalStagnation
          else FirstOrderStoppingStatus.IterationLimit
        val currentCertificate = certificate(
          FirstOrderMethod.LinearCompositePrimalDual,
          primal,
          Some(dual),
          objective,
          primalResidual,
          dualResidual,
          objectiveChange,
          iteration,
          config,
          primalScale,
          dualScale,
          primalResolution,
          dualResolution
        )
        Right(
          FirstOrderSolution(
            primal,
            Some(dual),
            objective,
            status,
            currentCertificate,
            primalStep = primalStep,
            dualStep = dualStep
          )
        )

  private def iterateSmoothComposite(
      smooth: SmoothObjective,
      direct: ProximalTerm,
      functional: LinearCompositeFunctional,
      operator: BoundedLinearOperator,
      initial: DMat,
      initialObjective: Double,
      config: FirstOrderConfig,
      execution: OptimizationExecution,
      initialDual: Option[DMat]
  ): Either[FirstOrderError, FirstOrderSolution] =
    val normBound = operator.normUpperBound
    val primalStep = config.stepSafety / (smooth.lipschitz + normBound)
    val dualStep = if normBound == 0.0 then 1.0 else config.stepSafety / normBound
    (validateStep(primalStep), validateStep(dualStep)) match
      case (Left(error), _) => return Left(error)
      case (_, Left(error)) => return Left(error)
      case _                => ()
    var primal = initial
    var extrapolated = initial
    var dual = initialDual.getOrElse(DMat.zeros(operator.linearOperator.rows, initial.cols))
    validateLike("initial dual", dual, DMat.zeros(operator.linearOperator.rows, initial.cols)) match
      case Left(error) => return Left(error)
      case _           => ()
    var objective = initialObjective
    var primalResidual = Double.MaxValue
    var dualResidual = Double.MaxValue
    var primalResolution = 0.0
    var dualResolution = 0.0
    var primalScale = 1.0
    var dualScale = 1.0
    var stagnated = false
    var objectiveChange = Double.MaxValue
    var iteration = 0
    var converged = false
    var error = Option.empty[FirstOrderError]
    while iteration < config.maxIterations && !converged && !stagnated && error.isEmpty do
      val updated =
        for
          mapped <- forward(operator.linearOperator, extrapolated, "smooth-composite forward", execution, false)
          dualTrial <- affineInput(dual, mapped, dualStep)
          nextDual <- functional.proximalConjugate(dualTrial, dualStep)
          _ <- validateLike("smooth-composite dual proximal result", nextDual, dual)
          gradient <- smooth.gradient(primal)
          _ <- validateLike("smooth-composite gradient", gradient, primal)
          adjoint <- forward(operator.linearOperator.adjoint, nextDual, "smooth-composite adjoint", execution, true)
          primalTrial <- affineInput(primal, add(gradient, adjoint), -primalStep)
          nextPrimal <- direct.proximal(primalTrial, primalStep)
          _ <- validateLike("smooth-composite primal proximal result", nextPrimal, primal)
          mappedNext <- forward(operator.linearOperator, nextPrimal, "smooth-composite forward", execution, false)
          nextObjective <- smoothCompositeValue(smooth, direct, functional, nextPrimal, mappedNext)
          residuals <- smoothCompositeResiduals(
            smooth,
            direct,
            functional,
            operator.linearOperator,
            nextPrimal,
            nextDual,
            Math.min(1.0, primalStep),
            Math.min(1.0, dualStep),
            execution,
            mappedNext
          )
        yield (nextPrimal, nextDual, nextObjective, residuals)
      updated match
        case Left(value)                                             => error = Some(value)
        case Right((nextPrimal, nextDual, nextObjective, residuals)) =>
          stagnated = sameValues(primal, nextPrimal) && sameValues(dual, nextDual)
          extrapolated = add(nextPrimal, scale(subtract(nextPrimal, primal), config.extrapolation))
          primal = nextPrimal
          dual = nextDual
          primalResidual = residuals._1
          dualResidual = residuals._2
          primalResolution = residuals._3
          dualResolution = residuals._4
          objectiveChange = Math.abs(nextObjective - objective)
          objective = nextObjective
          iteration += 1
          if iteration == 1 then
            primalScale = Math.max(1.0, primalResidual)
            dualScale = Math.max(1.0, dualResidual)
          if !objectiveChange.isFinite || !config.tolerance.threshold(primalScale).isFinite ||
            !config.tolerance.threshold(dualScale).isFinite
          then error = Some(FirstOrderError.NumericalFailure("stopping diagnostic overflow"))
          converged = error.isEmpty && primalResidual <= config.tolerance.threshold(primalScale) &&
            dualResidual <= config.tolerance.threshold(dualScale) &&
            primalResolution <= config.tolerance.threshold(primalScale) &&
            dualResolution <= config.tolerance.threshold(dualScale) &&
            objectiveChange <= config.tolerance.threshold(Math.max(1.0, Math.abs(objective)))
          if error.isEmpty then
            execution.record(
              FirstOrderSolution(
                primal,
                Some(dual),
                objective,
                FirstOrderStoppingStatus.IterationLimit,
                certificate(
                  FirstOrderMethod.SmoothCompositePrimalDual,
                  primal,
                  Some(dual),
                  objective,
                  primalResidual,
                  dualResidual,
                  objectiveChange,
                  iteration,
                  config,
                  primalScale,
                  dualScale,
                  primalResolution,
                  dualResolution
                ),
                primalStep = primalStep,
                dualStep = dualStep
              )
            )
            if !execution.continue then
              error = Some(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled))
    error match
      case Some(value) => Left(value)
      case None        =>
        val status =
          if converged then FirstOrderStoppingStatus.Converged
          else if stagnated then FirstOrderStoppingStatus.NumericalStagnation
          else FirstOrderStoppingStatus.IterationLimit
        val currentCertificate = certificate(
          FirstOrderMethod.SmoothCompositePrimalDual,
          primal,
          Some(dual),
          objective,
          primalResidual,
          dualResidual,
          objectiveChange,
          iteration,
          config,
          primalScale,
          dualScale,
          primalResolution,
          dualResolution
        )
        Right(
          FirstOrderSolution(
            primal,
            Some(dual),
            objective,
            status,
            currentCertificate,
            primalStep = primalStep,
            dualStep = dualStep
          )
        )

  private def smoothCompositeResiduals(
      smooth: SmoothObjective,
      direct: ProximalTerm,
      functional: LinearCompositeFunctional,
      operator: DoubleLinearOperator,
      primal: DMat,
      dual: DMat,
      primalStep: Double,
      dualStep: Double,
      execution: OptimizationExecution,
      mapped: DMat
  ): Either[FirstOrderError, (Double, Double, Double, Double)] =
    for
      gradient <- smooth.gradient(primal)
      _ <- validateLike("smooth-composite certificate gradient", gradient, primal)
      adjoint <- forward(operator.adjoint, dual, "smooth-composite certificate adjoint", execution, true)
      primalTrial <- affineInput(primal, add(gradient, adjoint), -primalStep)
      primalFixed <- direct.proximal(primalTrial, primalStep)
      _ <- validateLike("smooth-composite certificate primal proximal", primalFixed, primal)
      dualTrial <- affineInput(dual, mapped, dualStep)
      dualFixed <- functional.proximalConjugate(dualTrial, dualStep)
      _ <- validateLike("smooth-composite certificate dual proximal", dualFixed, dual)
      primalMeasure <- residualMeasure(primal, primalFixed, primalStep)
      dualMeasure <- residualMeasure(dual, dualFixed, dualStep)
    yield (primalMeasure._1, dualMeasure._1, primalMeasure._2, dualMeasure._2)

  private def fixedPointResiduals(
      primalObjective: ProximalObjective,
      functional: LinearCompositeFunctional,
      operator: DoubleLinearOperator,
      primal: DMat,
      dual: DMat,
      primalStep: Double,
      dualStep: Double,
      execution: OptimizationExecution,
      mapped: DMat
  ): Either[FirstOrderError, (Double, Double, Double, Double)] =
    for
      adjoint <- forward(operator.adjoint, dual, "certificate adjoint", execution, true)
      primalTrial <- affineInput(primal, adjoint, -primalStep)
      primalFixed <- primalObjective.proximal(primalTrial, primalStep)
      _ <- validateLike("certificate primal proximal", primalFixed, primal)
      dualTrial <- affineInput(dual, mapped, dualStep)
      dualFixed <- functional.proximalConjugate(dualTrial, dualStep)
      _ <- validateLike("certificate dual proximal", dualFixed, dual)
      primalMeasure <- residualMeasure(primal, primalFixed, primalStep)
      dualMeasure <- residualMeasure(dual, dualFixed, dualStep)
    yield (primalMeasure._1, dualMeasure._1, primalMeasure._2, dualMeasure._2)

  private def certificate(
      method: FirstOrderMethod,
      primal: DMat,
      dual: Option[DMat],
      objective: Double,
      primalResidual: Double,
      dualResidual: Double,
      objectiveChange: Double,
      iterations: Int,
      config: FirstOrderConfig,
      primalScale: Double = 1.0,
      dualScale: Double = 1.0,
      primalResolution: Double = 0.0,
      dualResolution: Double = 0.0
  ): FirstOrderCertificate =
    FirstOrderCertificate(
      ValueSummary.from(primal),
      dual.map(ValueSummary.from),
      objective,
      primalResidual,
      dualResidual,
      objectiveChange,
      iterations,
      FirstOrderSettings(
        method,
        config.maxIterations,
        config.tolerance,
        config.stepSafety,
        config.extrapolation,
        primalScale,
        dualScale,
        AlgorithmSettings.FixedStep,
        config.control.maxEvaluations
      ),
      primalResolution,
      dualResolution,
      Some(primal),
      dual
    )

  private def combinedValue(
      smooth: SmoothObjective,
      term: ProximalTerm,
      at: DMat
  ): Either[FirstOrderError, Double] =
    for
      left <- smooth.value(at).flatMap(validateScalar("smooth objective", _))
      right <- term.value(at).flatMap(validateScalar("proximal term", _))
      result <- validateScalar("combined objective", left + right)
    yield result

  private def compositeValue(
      primal: ProximalObjective,
      functional: LinearCompositeFunctional,
      at: DMat,
      mapped: DMat
  ): Either[FirstOrderError, Double] =
    for
      left <- primal.value(at).flatMap(validateScalar("primal objective", _))
      right <- functional.value(mapped).flatMap(validateScalar("composite functional", _))
      result <- validateScalar("composite objective", left + right)
    yield result

  private def smoothCompositeValue(
      smooth: SmoothObjective,
      direct: ProximalTerm,
      functional: LinearCompositeFunctional,
      at: DMat,
      mapped: DMat
  ): Either[FirstOrderError, Double] =
    for
      smoothValue <- smooth.value(at).flatMap(validateScalar("smooth objective", _))
      directValue <- direct.value(at).flatMap(validateScalar("direct proximal term", _))
      compositeValue <- functional.value(mapped).flatMap(validateScalar("linear-composite functional", _))
      result <- validateScalar("smooth-composite objective", smoothValue + directValue + compositeValue)
    yield result

  private def forward(
      operator: DoubleLinearOperator,
      value: DMat,
      context: String,
      execution: OptimizationExecution,
      adjoint: Boolean
  ): Either[FirstOrderError, DMat] =
    execution
      .invoke(if adjoint then EvaluationKind.Adjoint else EvaluationKind.Forward)(
        operator.applyTo(value).left.map(error => FirstOrderError.OperatorFailure(context, error))
      )
      .flatMap: result =>
        if result.cols != value.cols then
          Left(
            FirstOrderError.InvalidConfiguration(
              s"$context expected ${value.cols} columns, got ${result.cols}"
            )
          )
        else validateMatrix(context, result, operator.rows).map(_ => result)

/** Verify only that the image of `basis` is contained in the constraint null space, to the supplied tolerance. This
  * does not establish independence, rank, or completeness: even a zero basis satisfies containment. A caller using a
  * full parameterization must establish those properties separately.
  */
object ExactLinearReduction:
  def verify(
      basis: DoubleLinearOperator,
      constraint: DoubleLinearOperator,
      tolerance: FirstOrderTolerance
  ): Either[FirstOrderError, LinearReductionCertificate] =
    if basis.rows < 0 || basis.cols < 0 || constraint.rows < 0 || constraint.cols < 0 then
      Left(FirstOrderError.InvalidReduction("operator dimensions must be non-negative"))
    else if basis.rows != constraint.cols then
      Left(
        FirstOrderError.InvalidReduction(
          s"basis target dimension ${basis.rows} does not match constraint source dimension ${constraint.cols}"
        )
      )
    else
      for
        image <- basis
          .applyTo(DMat.eye(basis.cols))
          .left
          .map: error =>
            FirstOrderError.OperatorFailure("null-space basis", error)
        _ <-
          if image.cols != basis.cols then
            Left(FirstOrderError.InvalidReduction(s"basis image expected ${basis.cols} columns, got ${image.cols}"))
          else validateMatrix("null-space basis image", image, basis.rows)
        residualValue <- constraint
          .applyTo(image)
          .left
          .map: error =>
            FirstOrderError.OperatorFailure("null-space constraint", error)
        _ <-
          if residualValue.cols != basis.cols then
            Left(
              FirstOrderError.InvalidReduction(
                s"null-space residual expected ${basis.cols} columns, got ${residualValue.cols}"
              )
            )
          else validateMatrix("null-space residual", residualValue, constraint.rows)
        residual = maxAbs(residualValue)
        threshold = tolerance.threshold(Math.max(1.0, ValueSummary.from(image).maxAbs))
        _ <-
          if !threshold.isFinite then
            Left(FirstOrderError.InvalidReduction("null-space tolerance threshold is not finite"))
          else if residual <= threshold then Right(())
          else Left(FirstOrderError.InvalidReduction(s"null-space residual $residual exceeds threshold $threshold"))
      yield LinearReductionCertificate(
        basis.cols,
        basis.rows,
        constraint.rows,
        ValueSummary.from(image),
        ValueSummary.from(residualValue),
        residual,
        threshold,
        tolerance
      )

/** Evidence of null-space containment only; not rank or completeness. */
final case class LinearReductionCertificate(
    freeRows: Int,
    semanticRows: Int,
    constraintRows: Int,
    basisImage: ValueSummary,
    constraintImage: ValueSummary,
    residual: Double,
    threshold: Double,
    tolerance: FirstOrderTolerance
):
  require(residual.isFinite && residual >= 0.0, "reduction residual must be finite and non-negative")
  require(threshold.isFinite && threshold >= 0.0, "reduction threshold must be finite and non-negative")

private def validateMatrix(
    context: String,
    value: DMat,
    expectedRows: Int
): Either[FirstOrderError, Unit] =
  if value.rows != expectedRows then Left(FirstOrderError.ShapeMismatch(context, expectedRows, value.rows))
  else
    var row = 0
    var linearIndex = 0
    var error = Option.empty[FirstOrderError]
    while row < value.rows && error.isEmpty do
      var column = 0
      while column < value.cols && error.isEmpty do
        val current = value(row, column)
        if !current.isFinite then error = Some(FirstOrderError.NonFiniteValue(context, linearIndex, current))
        column += 1
        linearIndex += 1
      row += 1
    error.toLeft(())

private def validateLike(
    context: String,
    actual: DMat,
    expected: DMat
): Either[FirstOrderError, Unit] =
  if actual.cols != expected.cols then
    Left(FirstOrderError.InvalidConfiguration(s"$context expected ${expected.cols} columns, got ${actual.cols}"))
  else validateMatrix(context, actual, expected.rows)

private def validateScalar(context: String, value: Double): Either[FirstOrderError, Double] =
  if value.isFinite then Right(value)
  else Left(FirstOrderError.NonFiniteValue(context, 0, value))

private def add(left: DMat, right: DMat): DMat =
  DMat.tabulate(left.rows, left.cols): (row, column) =>
    left(row, column) + right(row, column)

private def subtract(left: DMat, right: DMat): DMat =
  DMat.tabulate(left.rows, left.cols): (row, column) =>
    left(row, column) - right(row, column)

private def scale(value: DMat, factor: Double): DMat =
  DMat.tabulate(value.rows, value.cols): (row, column) =>
    value(row, column) * factor

private def maxAbs(value: DMat): Double =
  var result = 0.0
  var row = 0
  while row < value.rows do
    var column = 0
    while column < value.cols do
      result = Math.max(result, Math.abs(value(row, column)))
      column += 1
    row += 1
  result

private def validateInitial(value: DMat, rows: Int): Either[FirstOrderError, Unit] =
  if rows <= 0 || value.cols <= 0 then
    Left(FirstOrderError.InvalidConfiguration("optimization variables must have positive rows and columns"))
  else validateMatrix("initial value", value, rows)

private def validateStep(step: Double): Either[FirstOrderError, Unit] =
  if step.isFinite && step > 0.0 then Right(())
  else Left(FirstOrderError.NumericalFailure(s"derived step must be finite and positive, got $step"))

private def affineInput(value: DMat, direction: DMat, factor: Double): Either[FirstOrderError, DMat] =
  val result = DMat.tabulate(value.rows, value.cols): (row, column) =>
    value(row, column) + factor * direction(row, column)
  validateLike("solver affine input", result, value).map(_ => result)

private def sameValues(left: DMat, right: DMat): Boolean =
  if left.rows != right.rows || left.cols != right.cols then false
  else if left eq right then true
  else
    var same = true
    var row = 0
    while row < left.rows && same do
      var column = 0
      while column < left.cols && same do
        same = left(row, column) == right(row, column)
        column += 1
      row += 1
    same

/** The second value is a resolution safeguard, not a bound on oracle error. */
private def residualMeasure(at: DMat, fixed: DMat, step: Double): Either[FirstOrderError, (Double, Double)] =
  var residual = 0.0
  var resolution = 0.0
  var row = 0
  while row < at.rows do
    var column = 0
    while column < at.cols do
      val a = at(row, column)
      val b = fixed(row, column)
      residual = Math.max(residual, Math.abs(a - b) / step)
      resolution = Math.max(resolution, (Math.ulp(Math.max(Math.abs(a), Math.abs(b))) / step) * 2.0)
      column += 1
    row += 1
  if !residual.isFinite || !resolution.isFinite then
    Left(FirstOrderError.NumericalFailure("fixed-point residual or its numerical resolution overflowed"))
  else Right((residual, resolution))
