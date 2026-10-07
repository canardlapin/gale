package gale.optim

enum BoxOptimizationError:
  case InvalidInput(detail: String)
  case OracleFailure(detail: String)
  case NonFiniteOracle

  def message: String = this match
    case InvalidInput(detail)  => detail
    case OracleFailure(detail) => detail
    case NonFiniteOracle       => "objective value and every gradient entry must be finite"

/** A finite, nonempty box. Coordinates may be fixed by equal lower/upper bounds. */
final case class BoxBounds private (lower: Vector[Double], upper: Vector[Double]):
  def dimension: Int = lower.length

object BoxBounds:
  def from(lower: Vector[Double], upper: Vector[Double]): Either[BoxOptimizationError, BoxBounds] =
    if lower.isEmpty || lower.length != upper.length ||
      !lower.indices.forall(i =>
        lower(i).isFinite && upper(i).isFinite && lower(i) <= upper(i) &&
          (upper(i) - lower(i)).isFinite
      )
    then
      Left(
        BoxOptimizationError.InvalidInput(
          "bounds must have equal positive dimensions and finite ordered endpoints/widths"
        )
      )
    else Right(BoxBounds(lower, upper))

/** Arrays are borrowed for this call only: read coordinates and fill every gradient entry. A refused call is charged
  * and stops that trajectory while preserving its last finite point.
  */
trait BoxDifferentiableObjective:
  def dimension: Int
  def evaluate(coordinates: Array[Double], gradient: Array[Double]): Either[BoxOptimizationError, Double]

final case class BoxQuasiNewtonConfig private (
    maxIterations: Int,
    maxEvaluations: Int,
    maxLineSearch: Int,
    projectedGradientTolerance: Double,
    armijo: Double,
    contraction: Double
)

object BoxQuasiNewtonConfig:
  def from(
      maxIterations: Int,
      maxEvaluations: Int,
      maxLineSearch: Int,
      projectedGradientTolerance: Double,
      armijo: Double = 1e-4,
      contraction: Double = 0.5
  ): Either[BoxOptimizationError, BoxQuasiNewtonConfig] =
    if maxIterations < 0 || maxEvaluations < 1 || maxLineSearch < 1 ||
      !projectedGradientTolerance.isFinite || projectedGradientTolerance <= 0.0 ||
      !armijo.isFinite || armijo <= 0.0 || armijo >= 1.0 ||
      !contraction.isFinite || contraction <= 0.0 || contraction >= 1.0
    then
      Left(
        BoxOptimizationError.InvalidInput("invalid iteration/evaluation limits or line-search/stationarity settings")
      )
    else
      Right(
        BoxQuasiNewtonConfig(
          maxIterations,
          maxEvaluations,
          maxLineSearch,
          projectedGradientTolerance,
          armijo,
          contraction
        )
      )

enum BoxOptimizationStatus:
  /** First-order stationarity only: no curvature or global-optimality claim. */
  case Stationary
  case IterationLimit
  case EvaluationLimit
  case LineSearchLimit
  case NoRepresentableStep
  case OracleRefused(error: BoxOptimizationError)

final case class BoxOptimizationPoint private[optim] (
    coordinates: Vector[Double],
    value: Double,
    gradient: Vector[Double],
    projectedGradientNorm: Double
)

final case class BoxOptimizationWork(evaluations: Int, iterations: Int, rejectedSteps: Int, metricResets: Int)

/** Point is absent only when the initial oracle call refuses. Counts include failed calls. */
final case class BoxOptimizationResult private[optim] (
    point: Option[BoxOptimizationPoint],
    status: BoxOptimizationStatus,
    work: BoxOptimizationWork
)

/** Full-memory projected BFGS with an Armijo search along the projected proposal.
  *
  * Active coordinates are fixed while updating the free inverse metric. Active-set changes reset that metric;
  * non-descent directions fall back to projected steepest descent. This is a small/moderate-dimensional local solver,
  * not L-BFGS-B or a global optimizer. No positive objective Hessian is assumed. Coordinates and objective scaling are
  * caller policy. Retained primitive storage is one n*n metric plus O(n) scratch, independent of iteration count. No
  * objective evaluation is made outside the box or beyond the declared evaluation cap.
  */
object BoxQuasiNewton:
  def minimize(
      objective: BoxDifferentiableObjective,
      bounds: BoxBounds,
      initial: Vector[Double],
      config: BoxQuasiNewtonConfig
  ): Either[BoxOptimizationError, BoxOptimizationResult] =
    val n = bounds.dimension
    if n > 1024 || objective.dimension != n || initial.length != n ||
      !initial.indices.forall(i =>
        initial(i).isFinite && initial(i) >= bounds.lower(i) && initial(i) <= bounds.upper(i)
      )
    then
      Left(
        BoxOptimizationError.InvalidInput(
          "initial point/objective must match a box of dimension 1..1024 and lie inside it"
        )
      )
    else Right(run(objective, bounds, initial, config))

  private def run(
      objective: BoxDifferentiableObjective,
      bounds: BoxBounds,
      initial: Vector[Double],
      config: BoxQuasiNewtonConfig
  ): BoxOptimizationResult =
    val n = bounds.dimension
    val x = initial.toArray
    val gradient = new Array[Double](n)
    val candidate = new Array[Double](n)
    val nextGradient = new Array[Double](n)
    val projected = new Array[Double](n)
    val direction = new Array[Double](n)
    val step = new Array[Double](n)
    val y = new Array[Double](n)
    val hy = new Array[Double](n)
    val active = new Array[Boolean](n)
    val nextActive = new Array[Boolean](n)
    val metric = new Array[Double](n * n)
    var evaluations = 0
    var iterations = 0
    var rejected = 0
    var resets = 0
    var value = Double.NaN
    var validPoint = false

    def resetMetric(): Unit =
      java.util.Arrays.fill(metric, 0.0)
      var i = 0
      while i < n do
        metric(i * n + i) = 1.0
        i += 1

    def boundActive(at: Array[Double], g: Array[Double], out: Array[Boolean]): Unit =
      var i = 0
      while i < n do
        out(i) = bounds.lower(i) == bounds.upper(i) ||
          (at(i) <= bounds.lower(i) && g(i) > 0.0) || (at(i) >= bounds.upper(i) && g(i) < 0.0)
        i += 1

    def projectedNorm(): Double =
      boundActive(x, gradient, active)
      var norm = 0.0
      var i = 0
      while i < n do
        projected(i) = if active(i) then 0.0 else gradient(i)
        norm = math.max(norm, math.abs(projected(i)))
        i += 1
      norm

    def finish(status: BoxOptimizationStatus): BoxOptimizationResult =
      val point =
        if validPoint then Some(BoxOptimizationPoint(x.toVector, value, gradient.toVector, projectedNorm())) else None
      BoxOptimizationResult(point, status, BoxOptimizationWork(evaluations, iterations, rejected, resets))

    def evaluate(at: Array[Double], out: Array[Double]): Either[BoxOptimizationError, Double] =
      evaluations += 1
      java.util.Arrays.fill(out, Double.NaN)
      objective
        .evaluate(at, out)
        .flatMap: energy =>
          if energy.isFinite && out.forall(_.isFinite) then Right(energy)
          else Left(BoxOptimizationError.NonFiniteOracle)

    evaluate(x, gradient) match
      case Left(error)         => return finish(BoxOptimizationStatus.OracleRefused(error))
      case Right(initialValue) => value = initialValue; validPoint = true
    resetMetric()
    while true do
      val norm = projectedNorm()
      if norm <= config.projectedGradientTolerance then return finish(BoxOptimizationStatus.Stationary)
      if iterations >= config.maxIterations then return finish(BoxOptimizationStatus.IterationLimit)
      if evaluations >= config.maxEvaluations then return finish(BoxOptimizationStatus.EvaluationLimit)
      var slope = 0.0
      var pgNorm2 = 0.0
      var dNorm2 = 0.0
      var i = 0
      while i < n do
        var sum = 0.0
        var j = 0
        while j < n do
          sum += metric(i * n + j) * projected(j)
          j += 1
        val proposed = -sum
        direction(i) =
          if active(i) || (x(i) <= bounds.lower(i) && proposed < 0.0) ||
            (x(i) >= bounds.upper(i) && proposed > 0.0)
          then 0.0
          else proposed
        slope += gradient(i) * direction(i)
        pgNorm2 += projected(i) * projected(i)
        dNorm2 += direction(i) * direction(i)
        i += 1
      if !slope.isFinite || slope >= -1e-12 * math.sqrt(pgNorm2) * math.sqrt(dNorm2) then
        resetMetric()
        resets += 1
        i = 0
        while i < n do
          direction(i) = -projected(i)
          i += 1
      var accepted = false
      var alpha = 1.0
      var attempt = 0
      var nextValue = Double.NaN
      var anyMoved = false
      while !accepted && attempt < config.maxLineSearch do
        slope = 0.0
        var moved = false
        i = 0
        while i < n do
          candidate(i) = math.max(bounds.lower(i), math.min(bounds.upper(i), x(i) + alpha * direction(i)))
          step(i) = candidate(i) - x(i)
          moved = moved || step(i) != 0.0
          slope += gradient(i) * step(i)
          i += 1
        anyMoved = anyMoved || moved
        if moved && slope.isFinite && slope < 0.0 then
          if evaluations >= config.maxEvaluations then return finish(BoxOptimizationStatus.EvaluationLimit)
          evaluate(candidate, nextGradient) match
            case Left(error)   => return finish(BoxOptimizationStatus.OracleRefused(error))
            case Right(energy) =>
              nextValue = energy
              accepted = energy <= value + config.armijo * slope
        if !accepted then rejected += 1
        alpha *= config.contraction
        attempt += 1
      if !accepted then
        return finish(
          if anyMoved then BoxOptimizationStatus.LineSearchLimit else BoxOptimizationStatus.NoRepresentableStep
        )
      boundActive(candidate, nextGradient, nextActive)
      var changed = false
      var ys = 0.0
      var sNorm2 = 0.0
      var yNorm2 = 0.0
      i = 0
      while i < n do
        changed = changed || active(i) != nextActive(i)
        y(i) = if active(i) then 0.0 else nextGradient(i) - gradient(i)
        ys += y(i) * step(i)
        sNorm2 += step(i) * step(i)
        yNorm2 += y(i) * y(i)
        i += 1
      if changed then
        resetMetric()
        resets += 1
      else if ys.isFinite && ys > 1e-12 * math.sqrt(sNorm2) * math.sqrt(yNorm2) then
        var yhy = 0.0
        i = 0
        while i < n do
          hy(i) = 0.0
          var j = 0
          while j < n do
            hy(i) += metric(i * n + j) * y(j)
            j += 1
          yhy += y(i) * hy(i)
          i += 1
        val coefficient = (1.0 + yhy / ys) / ys
        var metricFinite = coefficient.isFinite
        i = 0
        while i < n do
          var j = 0
          while j <= i do
            val updated = metric(i * n + j) + coefficient * step(i) * step(j) -
              (step(i) * hy(j) + hy(i) * step(j)) / ys
            metric(i * n + j) = updated
            metric(j * n + i) = updated
            metricFinite = metricFinite && updated.isFinite
            j += 1
          i += 1
        if !metricFinite then
          resetMetric()
          resets += 1
      System.arraycopy(candidate, 0, x, 0, n)
      System.arraycopy(nextGradient, 0, gradient, 0, n)
      value = nextValue
      iterations += 1
    finish(BoxOptimizationStatus.IterationLimit)
