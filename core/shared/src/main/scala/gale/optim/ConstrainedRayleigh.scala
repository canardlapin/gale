package gale.optim

import gale.linalg.{DMat, DVec}

/** Failure at the trust boundary of a constrained generalized-Rayleigh solve. */
enum ConstrainedRayleighError:
  case DimensionMismatch(numeratorRows: Int, numeratorCols: Int, denominatorRows: Int, denominatorCols: Int)
  case InvalidConfiguration(reason: String)
  case NonFiniteInput(label: String, row: Int, column: Int, value: Double)
  case NonPositiveDenominator(value: Double)
  case InfeasibleInitialPoint

  def message: String =
    this match
      case DimensionMismatch(nr, nc, dr, dc) =>
        s"generalized-Rayleigh matrices must be non-empty square matrices of equal shape, got ${nr}x${nc} and ${dr}x${dc}"
      case InvalidConfiguration(reason) => reason
      case NonFiniteInput(label, row, column, value) =>
        s"$label contains non-finite value $value at ($row, $column)"
      case NonPositiveDenominator(value) =>
        s"generalized-Rayleigh denominator must be positive, got $value"
      case InfeasibleInitialPoint =>
        "the projected initial point has zero denominator norm"

/** Closed convex cone used by the projected-Rayleigh numerical kernel.
  *
  * The catalog is intentionally numerical. Domain meanings and admissible
  * combinations belong to the consuming model compiler.
  */
enum RayleighCone:
  case NonnegativeOrthant

  private[optim] def project(value: DVec): DVec =
    this match
      case NonnegativeOrthant =>
        val out = new Array[Double](value.length)
        var index = 0
        while index < value.length do
          out(index) = Math.max(0.0, value(index))
          index += 1
        DVec.fromArray(out)

  private[optim] def violation(value: DVec): Double =
    this match
      case NonnegativeOrthant =>
        var largest = 0.0
        var index = 0
        while index < value.length do
          largest = Math.max(largest, Math.max(0.0, -value(index)))
          index += 1
        largest

  /** KKT residual after the equality-normalization multiplier has been
    * eliminated by the Rayleigh root.
    */
  private[optim] def stationarity(value: DVec, tangentGradient: DVec, activeTolerance: Double): Double =
    this match
      case NonnegativeOrthant =>
        var squared = 0.0
        var index = 0
        while index < value.length do
          val residual =
            if value(index) > activeTolerance then tangentGradient(index)
            else Math.max(0.0, tangentGradient(index))
          squared += residual * residual
          index += 1
        Math.sqrt(squared)

final case class ProjectedRayleighConfig(
    tolerance: Double = 1e-9,
    maxIterations: Int = 5000,
    initialStep: Double = 1.0,
    minimumStep: Double = 1e-12,
    backtrackingFactor: Double = 0.5
)

enum ProjectedRayleighTermination:
  case Converged
  case IterationLimit
  case StepUnderflow

final case class ProjectedRayleighCertificate(
    stationarityResidual: Double,
    constraintViolation: Double,
    normalizationError: Double,
    objectiveChange: Double
):
  require(stationarityResidual.isFinite && stationarityResidual >= 0.0)
  require(constraintViolation.isFinite && constraintViolation >= 0.0)
  require(normalizationError.isFinite && normalizationError >= 0.0)
  require(objectiveChange.isFinite && objectiveChange >= 0.0)

final case class ProjectedRayleighResult(
    direction: DVec,
    root: Double,
    iterations: Int,
    termination: ProjectedRayleighTermination,
    certificate: ProjectedRayleighCertificate
):
  require(root.isFinite)
  require(iterations >= 0)

  def converged: Boolean = termination == ProjectedRayleighTermination.Converged

/** Portable projected ascent for a cone-constrained generalized Rayleigh
  * quotient `x' A x / x' B x`.
  *
  * The kernel never forms `B^-1`. Every accepted iterate is projected onto the
  * declared cone and normalized directly in the `B` geometry. The certificate
  * reports KKT stationarity, feasibility, and normalization independently;
  * because the feasible normalized problem is non-convex, convergence certifies
  * a stationary point rather than a global optimum.
  */
object ProjectedRayleigh:
  def solve(
      numerator: DMat,
      denominator: DMat,
      cone: RayleighCone,
      initial: DVec,
      config: ProjectedRayleighConfig = ProjectedRayleighConfig()
  ): Either[ConstrainedRayleighError, ProjectedRayleighResult] =
    for
      _ <- validate(numerator, denominator, initial, config)
      normalized <- normalize(cone.project(initial), denominator)
    yield iterate(numerator, denominator, cone, normalized, config)

  /** Deterministic multi-start solve using the positive uniform vector and all
    * coordinate rays. The best converged stationary point is returned; if no
    * start converges, the best finite iterate is retained with its termination.
    */
  def solveNonnegative(
      numerator: DMat,
      denominator: DMat,
      config: ProjectedRayleighConfig = ProjectedRayleighConfig()
  ): Either[ConstrainedRayleighError, ProjectedRayleighResult] =
    val dimension = numerator.rows
    val starts = Vector.newBuilder[DVec]
    starts += DVec.fromArray(Array.fill(dimension)(1.0))
    var coordinate = 0
    while coordinate < dimension do
      val values = new Array[Double](dimension)
      values(coordinate) = 1.0
      starts += DVec.fromArray(values)
      coordinate += 1

    val candidates = starts.result()
    var best = Option.empty[ProjectedRayleighResult]
    var index = 0
    while index < candidates.length do
      solve(numerator, denominator, RayleighCone.NonnegativeOrthant, candidates(index), config) match
        case Left(error) => return Left(error)
        case Right(candidate) =>
          best match
            case None => best = Some(candidate)
            case Some(current) =>
              val betterTermination = candidate.converged && !current.converged
              val sameTermination = candidate.converged == current.converged
              if betterTermination || (sameTermination && candidate.root > current.root) then best = Some(candidate)
      index += 1
    best.toRight(ConstrainedRayleighError.InfeasibleInitialPoint)

  private def validate(
      numerator: DMat,
      denominator: DMat,
      initial: DVec,
      config: ProjectedRayleighConfig
  ): Either[ConstrainedRayleighError, Unit] =
    if numerator.rows <= 0 || numerator.rows != numerator.cols ||
        denominator.rows != denominator.cols || numerator.rows != denominator.rows ||
        initial.length != numerator.rows
    then
      Left(
        ConstrainedRayleighError.DimensionMismatch(
          numerator.rows,
          numerator.cols,
          denominator.rows,
          denominator.cols
        )
      )
    else if !config.tolerance.isFinite || config.tolerance <= 0.0 then
      Left(ConstrainedRayleighError.InvalidConfiguration("tolerance must be finite and positive"))
    else if config.maxIterations <= 0 then
      Left(ConstrainedRayleighError.InvalidConfiguration("maxIterations must be positive"))
    else if !config.initialStep.isFinite || config.initialStep <= 0.0 then
      Left(ConstrainedRayleighError.InvalidConfiguration("initialStep must be finite and positive"))
    else if !config.minimumStep.isFinite || config.minimumStep <= 0.0 || config.minimumStep > config.initialStep then
      Left(ConstrainedRayleighError.InvalidConfiguration("minimumStep must be finite, positive, and no larger than initialStep"))
    else if !config.backtrackingFactor.isFinite || config.backtrackingFactor <= 0.0 || config.backtrackingFactor >= 1.0 then
      Left(ConstrainedRayleighError.InvalidConfiguration("backtrackingFactor must lie strictly between zero and one"))
    else
      finite("numerator", numerator)
        .flatMap(_ => finite("denominator", denominator))
        .flatMap(_ => finiteInitial(initial))

  private def finite(label: String, matrix: DMat): Either[ConstrainedRayleighError, Unit] =
    var row = 0
    while row < matrix.rows do
      var column = 0
      while column < matrix.cols do
        val value = matrix(row, column)
        if !value.isFinite then return Left(ConstrainedRayleighError.NonFiniteInput(label, row, column, value))
        column += 1
      row += 1
    Right(())

  private def finiteInitial(initial: DVec): Either[ConstrainedRayleighError, Unit] =
    var failure = Option.empty[ConstrainedRayleighError]
    var index = 0
    while index < initial.length && failure.isEmpty do
      val value = initial(index)
      if !value.isFinite then failure = Some(ConstrainedRayleighError.NonFiniteInput("initial", index, 0, value))
      index += 1
    failure.toLeft(())

  private def iterate(
      numerator: DMat,
      denominator: DMat,
      cone: RayleighCone,
      initial: DVec,
      config: ProjectedRayleighConfig
  ): ProjectedRayleighResult =
    var current = initial
    var root = quotient(current, numerator, denominator)
    var previousRoot = root
    var iteration = 0
    var termination = ProjectedRayleighTermination.IterationLimit
    var running = true
    var objectiveChange = Double.PositiveInfinity

    while running && iteration < config.maxIterations do
      val tangent = tangentGradient(current, numerator, denominator, root)
      val scale = Math.max(1.0, matrixScale(numerator) + Math.abs(root) * matrixScale(denominator))
      val stationarity = cone.stationarity(current, tangent, config.tolerance)
      val normalizedError = Math.abs(quadratic(current, denominator) - 1.0)
      val violation = cone.violation(current)
      if stationarity <= config.tolerance * scale &&
          normalizedError <= config.tolerance && violation <= config.tolerance
      then
        termination = ProjectedRayleighTermination.Converged
        running = false
      else
        var step = config.initialStep / scale
        var accepted = Option.empty[(DVec, Double)]
        while accepted.isEmpty && step >= config.minimumStep do
          val trial = addScaled(current, tangent, step)
          normalize(cone.project(trial), denominator) match
            case Right(candidate) =>
              val candidateRoot = quotient(candidate, numerator, denominator)
              if candidateRoot + config.tolerance * scale >= root then accepted = Some(candidate -> candidateRoot)
              else step *= config.backtrackingFactor
            case Left(_) => step *= config.backtrackingFactor
        accepted match
          case None =>
            termination = ProjectedRayleighTermination.StepUnderflow
            running = false
          case Some((candidate, candidateRoot)) =>
            previousRoot = root
            current = candidate
            root = candidateRoot
            objectiveChange = Math.abs(root - previousRoot)
            iteration += 1

    val finalTangent = tangentGradient(current, numerator, denominator, root)
    val finalChange = if objectiveChange.isFinite then objectiveChange else 0.0
    ProjectedRayleighResult(
      current,
      root,
      iteration,
      termination,
      ProjectedRayleighCertificate(
        cone.stationarity(current, finalTangent, config.tolerance),
        cone.violation(current),
        Math.abs(quadratic(current, denominator) - 1.0),
        finalChange
      )
    )

  private def normalize(value: DVec, denominator: DMat): Either[ConstrainedRayleighError, DVec] =
    val normSquared = quadratic(value, denominator)
    if !normSquared.isFinite || normSquared <= 0.0 then
      Left(ConstrainedRayleighError.NonPositiveDenominator(normSquared))
    else
      val scale = 1.0 / Math.sqrt(normSquared)
      val out = new Array[Double](value.length)
      var index = 0
      while index < value.length do
        out(index) = scale * value(index)
        index += 1
      Right(DVec.fromArray(out))

  private def quotient(value: DVec, numerator: DMat, denominator: DMat): Double =
    quadratic(value, numerator) / quadratic(value, denominator)

  private def tangentGradient(value: DVec, numerator: DMat, denominator: DMat, root: Double): DVec =
    val av = multiply(numerator, value)
    val bv = multiply(denominator, value)
    val out = new Array[Double](value.length)
    var index = 0
    while index < value.length do
      out(index) = 2.0 * (av(index) - root * bv(index))
      index += 1
    DVec.fromArray(out)

  private def quadratic(value: DVec, matrix: DMat): Double =
    var result = 0.0
    var row = 0
    while row < matrix.rows do
      var projected = 0.0
      var column = 0
      while column < matrix.cols do
        projected += matrix(row, column) * value(column)
        column += 1
      result += value(row) * projected
      row += 1
    result

  private def multiply(matrix: DMat, value: DVec): DVec =
    val out = new Array[Double](matrix.rows)
    var row = 0
    while row < matrix.rows do
      var result = 0.0
      var column = 0
      while column < matrix.cols do
        result += matrix(row, column) * value(column)
        column += 1
      out(row) = result
      row += 1
    DVec.fromArray(out)

  private def addScaled(left: DVec, right: DVec, scale: Double): DVec =
    val out = new Array[Double](left.length)
    var index = 0
    while index < left.length do
      out(index) = left(index) + scale * right(index)
      index += 1
    DVec.fromArray(out)

  private def matrixScale(matrix: DMat): Double =
    var squared = 0.0
    var row = 0
    while row < matrix.rows do
      var column = 0
      while column < matrix.cols do
        val value = matrix(row, column)
        squared += value * value
        column += 1
      row += 1
    Math.sqrt(squared)
