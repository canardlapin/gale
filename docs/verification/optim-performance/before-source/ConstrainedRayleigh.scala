package gale.optim

import gale.linalg.{DMat, DVec, DenseDecompositions, LinAlgError}

/** Failure at the trust boundary of a constrained generalized-Rayleigh solve. */
enum ConstrainedRayleighError:
  case DimensionMismatch(numeratorRows: Int, numeratorCols: Int, denominatorRows: Int, denominatorCols: Int)
  case InvalidConfiguration(reason: String)
  case NonFiniteInput(label: String, row: Int, column: Int, value: Double)
  case AsymmetricInput(label: String, row: Int, column: Int, upper: Double, lower: Double)
  case DenominatorNotPositiveDefinite(index: Int)
  case NonFiniteComputation(label: String, value: Double)
  case NonPositiveDenominator(value: Double)
  case InfeasibleInitialPoint

  def message: String =
    this match
      case DimensionMismatch(nr, nc, dr, dc) =>
        s"generalized-Rayleigh matrices must be non-empty square matrices of equal shape, got ${nr}x${nc} and ${dr}x${dc}"
      case InvalidConfiguration(reason)              => reason
      case NonFiniteInput(label, row, column, value) =>
        s"$label contains non-finite value $value at ($row, $column)"
      case AsymmetricInput(label, row, column, upper, lower) =>
        s"$label must be symmetric: ($row, $column)=$upper differs from ($column, $row)=$lower"
      case DenominatorNotPositiveDefinite(index) =>
        s"generalized-Rayleigh denominator is not positive definite at leading minor $index"
      case NonFiniteComputation(label, value) =>
        s"$label produced non-finite value $value"
      case NonPositiveDenominator(value) =>
        s"generalized-Rayleigh denominator must be positive, got $value"
      case InfeasibleInitialPoint =>
        "the projected initial point has zero denominator norm"

/** Closed convex cone used by the projected-Rayleigh numerical kernel.
  *
  * The catalog is intentionally numerical. Domain meanings and admissible combinations belong to the consuming model
  * compiler.
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

  /** KKT residual after the equality-normalization multiplier has been eliminated by the Rayleigh root.
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

/** Portable projected ascent for a cone-constrained generalized Rayleigh quotient `x' A x / x' B x`.
  *
  * The kernel never forms `B^-1`. Every accepted iterate is projected onto the declared cone and normalized directly in
  * the `B` geometry. The certificate reports KKT stationarity, feasibility, and normalization independently; because
  * the feasible normalized problem is non-convex, convergence certifies a stationary point rather than a global
  * optimum.
  */
object ProjectedRayleigh:
  private final case class Geometry(numeratorScale: Double, denominatorScale: Double)

  def solve(
      numerator: DMat,
      denominator: DMat,
      cone: RayleighCone,
      initial: DVec,
      config: ProjectedRayleighConfig = ProjectedRayleighConfig()
  ): Either[ConstrainedRayleighError, ProjectedRayleighResult] =
    for
      _ <- validateInputs(numerator, denominator, initial, config)
      normalized <- normalize(cone.project(initial), denominator)
      geometry <- validateDenominator(denominator, numerator, config)
      result <- iterate(numerator, denominator, geometry, cone, normalized, config)
    yield result

  /** Deterministic multi-start solve using the positive uniform vector and all coordinate rays. The best converged
    * stationary point is returned; if no start converges, the best finite iterate is retained with its termination.
    */
  def solveNonnegative(
      numerator: DMat,
      denominator: DMat,
      config: ProjectedRayleighConfig = ProjectedRayleighConfig()
  ): Either[ConstrainedRayleighError, ProjectedRayleighResult] =
    validateDenominator(denominator, numerator, config)
      .flatMap(solveNonnegativeValidated(numerator, denominator, _, config))

  private def solveNonnegativeValidated(
      numerator: DMat,
      denominator: DMat,
      geometry: Geometry,
      config: ProjectedRayleighConfig
  ): Either[ConstrainedRayleighError, ProjectedRayleighResult] =
    val dimension = numerator.rows
    var best = Option.empty[ProjectedRayleighResult]
    var start = -1
    while start < dimension do
      val values = Array.fill(dimension)(if start == -1 then 1.0 else 0.0)
      if start >= 0 then values(start) = 1.0
      val candidate =
        normalize(DVec.fromArray(values), denominator)
          .flatMap(iterate(numerator, denominator, geometry, RayleighCone.NonnegativeOrthant, _, config)) match
          case Left(error)  => return Left(error)
          case Right(value) => value
      best match
        case None          => best = Some(candidate)
        case Some(current) =>
          val betterTermination = candidate.converged && !current.converged
          val sameTermination = candidate.converged == current.converged
          if betterTermination || (sameTermination && candidate.root > current.root) then best = Some(candidate)
      start += 1
    best.toRight(ConstrainedRayleighError.InfeasibleInitialPoint)

  private def validateInputs(
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
      Left(
        ConstrainedRayleighError.InvalidConfiguration(
          "minimumStep must be finite, positive, and no larger than initialStep"
        )
      )
    else if !config.backtrackingFactor.isFinite || config.backtrackingFactor <= 0.0 || config.backtrackingFactor >= 1.0
    then
      Left(ConstrainedRayleighError.InvalidConfiguration("backtrackingFactor must lie strictly between zero and one"))
    else
      finite("numerator", numerator)
        .flatMap(_ => finite("denominator", denominator))
        .flatMap(_ => symmetric("numerator", numerator))
        .flatMap(_ => symmetric("denominator", denominator))
        .flatMap(_ => finiteInitial(initial))

  private def validateDenominator(
      denominator: DMat,
      numerator: DMat,
      config: ProjectedRayleighConfig
  ): Either[ConstrainedRayleighError, Geometry] =
    validateMatrixConfiguration(numerator, denominator, config)
      .flatMap(_ => finite("numerator", numerator))
      .flatMap(_ => finite("denominator", denominator))
      .flatMap(_ => symmetric("numerator", numerator))
      .flatMap(_ => symmetric("denominator", denominator))
      .flatMap { _ =>
        DenseDecompositions.cholesky(denominator) match
          case Right(_) => Right(Geometry(matrixScale(numerator), matrixScale(denominator)))
          case Left(LinAlgError.NotPositiveDefinite(index)) =>
            Left(ConstrainedRayleighError.DenominatorNotPositiveDefinite(index))
          case Left(error) => Left(ConstrainedRayleighError.InvalidConfiguration(error.getMessage))
      }

  private def validateMatrixConfiguration(
      numerator: DMat,
      denominator: DMat,
      config: ProjectedRayleighConfig
  ): Either[ConstrainedRayleighError, Unit] =
    if numerator.rows <= 0 || numerator.rows != numerator.cols ||
      denominator.rows != denominator.cols || numerator.rows != denominator.rows
    then
      Left(
        ConstrainedRayleighError.DimensionMismatch(numerator.rows, numerator.cols, denominator.rows, denominator.cols)
      )
    else if !config.tolerance.isFinite || config.tolerance <= 0.0 then
      Left(ConstrainedRayleighError.InvalidConfiguration("tolerance must be finite and positive"))
    else if config.maxIterations <= 0 then
      Left(ConstrainedRayleighError.InvalidConfiguration("maxIterations must be positive"))
    else if !config.initialStep.isFinite || config.initialStep <= 0.0 then
      Left(ConstrainedRayleighError.InvalidConfiguration("initialStep must be finite and positive"))
    else if !config.minimumStep.isFinite || config.minimumStep <= 0.0 || config.minimumStep > config.initialStep then
      Left(
        ConstrainedRayleighError.InvalidConfiguration(
          "minimumStep must be finite, positive, and no larger than initialStep"
        )
      )
    else if !config.backtrackingFactor.isFinite || config.backtrackingFactor <= 0.0 || config.backtrackingFactor >= 1.0
    then
      Left(ConstrainedRayleighError.InvalidConfiguration("backtrackingFactor must lie strictly between zero and one"))
    else Right(())

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

  private def symmetric(label: String, matrix: DMat): Either[ConstrainedRayleighError, Unit] =
    var row = 1
    while row < matrix.rows do
      var column = 0
      while column < row do
        val lower = matrix(row, column)
        val upper = matrix(column, row)
        if lower != upper then return Left(ConstrainedRayleighError.AsymmetricInput(label, row, column, upper, lower))
        column += 1
      row += 1
    Right(())

  private def iterate(
      numerator: DMat,
      denominator: DMat,
      geometry: Geometry,
      cone: RayleighCone,
      initial: DVec,
      config: ProjectedRayleighConfig
  ): Either[ConstrainedRayleighError, ProjectedRayleighResult] =
    var current = initial
    var root = quotient(current, numerator, denominator) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    var previousRoot = root
    var iteration = 0
    var termination = ProjectedRayleighTermination.IterationLimit
    var running = true
    var objectiveChange = Double.PositiveInfinity

    while running && iteration < config.maxIterations do
      val tangent = tangentGradient(current, numerator, denominator, root) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      val scale = geometryScale(geometry, root) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      val residualThresholdScale = residualScale(scale, current) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      val stationarity = cone.stationarity(current, tangent, config.tolerance)
      if !stationarity.isFinite then
        return Left(ConstrainedRayleighError.NonFiniteComputation("stationarity residual", stationarity))
      val normalizedError = quadratic(current, denominator)
        .flatMap(finite("normalization quadratic", _))
        .map(value => Math.abs(value - 1.0)) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      val violation = cone.violation(current)
      if stationarity <= config.tolerance * residualThresholdScale &&
        normalizedError <= config.tolerance && violation <= config.tolerance
      then
        termination = ProjectedRayleighTermination.Converged
        running = false
      else
        var step = config.initialStep / scale
        var accepted = Option.empty[(DVec, Double)]
        while accepted.isEmpty && step * scale >= config.minimumStep do
          addScaled(current, tangent, step) match
            case Left(error)  => return Left(error)
            case Right(trial) =>
              normalize(cone.project(trial), denominator) match
                case Right(candidate) =>
                  quotient(candidate, numerator, denominator) match
                    case Left(error)          => return Left(error)
                    case Right(candidateRoot) =>
                      val improvement = candidateRoot + config.tolerance * Math.max(1.0, Math.abs(root))
                      if !improvement.isFinite then
                        return Left(
                          ConstrainedRayleighError.NonFiniteComputation("line-search improvement bound", improvement)
                        )
                      if sameVector(candidate, current) then step = 0.0
                      else if improvement >= root then accepted = Some(candidate -> candidateRoot)
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
            if !objectiveChange.isFinite then
              return Left(ConstrainedRayleighError.NonFiniteComputation("objective change", objectiveChange))
            iteration += 1

    val finalTangent = tangentGradient(current, numerator, denominator, root) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    val finalNormalization = quadratic(current, denominator).flatMap(finite("final normalization quadratic", _)) match
      case Left(error)  => return Left(error)
      case Right(value) => Math.abs(value - 1.0)
    val finalStationarity = cone.stationarity(current, finalTangent, config.tolerance)
    if !finalStationarity.isFinite then
      return Left(ConstrainedRayleighError.NonFiniteComputation("final stationarity residual", finalStationarity))
    val finalChange = if objectiveChange.isFinite then objectiveChange else 0.0
    Right(
      ProjectedRayleighResult(
        current,
        root,
        iteration,
        termination,
        ProjectedRayleighCertificate(
          finalStationarity,
          cone.violation(current),
          finalNormalization,
          finalChange
        )
      )
    )

  private def normalize(value: DVec, denominator: DMat): Either[ConstrainedRayleighError, DVec] =
    val normSquared = quadratic(value, denominator) match
      case Left(ConstrainedRayleighError.NonFiniteComputation(_, value)) =>
        return Left(ConstrainedRayleighError.NonPositiveDenominator(value))
      case Left(error)  => return Left(error)
      case Right(value) => value
    if normSquared <= 0.0 then Left(ConstrainedRayleighError.NonPositiveDenominator(normSquared))
    else
      val scale = 1.0 / Math.sqrt(normSquared)
      val out = new Array[Double](value.length)
      var index = 0
      while index < value.length do
        out(index) = scale * value(index)
        if !out(index).isFinite then
          return Left(ConstrainedRayleighError.NonFiniteComputation("normalization", out(index)))
        index += 1
      Right(DVec.fromArray(out))

  private def quotient(value: DVec, numerator: DMat, denominator: DMat): Either[ConstrainedRayleighError, Double] =
    for
      top <- quadratic(value, numerator)
      bottom <- quadratic(value, denominator)
      _ <- if bottom > 0.0 then Right(()) else Left(ConstrainedRayleighError.NonPositiveDenominator(bottom))
      result <- finite("Rayleigh quotient", top / bottom)
    yield result

  private def tangentGradient(
      value: DVec,
      numerator: DMat,
      denominator: DMat,
      root: Double
  ): Either[ConstrainedRayleighError, DVec] =
    val av = multiply(numerator, value, "numerator product") match
      case Left(error)   => return Left(error)
      case Right(result) => result
    val bv = multiply(denominator, value, "denominator product") match
      case Left(error)   => return Left(error)
      case Right(result) => result
    val out = new Array[Double](value.length)
    var index = 0
    while index < value.length do
      out(index) = 2.0 * (av(index) - root * bv(index))
      if !out(index).isFinite then
        return Left(ConstrainedRayleighError.NonFiniteComputation("tangent gradient", out(index)))
      index += 1
    Right(DVec.fromArray(out))

  private def quadratic(value: DVec, matrix: DMat): Either[ConstrainedRayleighError, Double] =
    var result = 0.0
    var row = 0
    while row < matrix.rows do
      var projected = 0.0
      var column = 0
      while column < matrix.cols do
        projected += matrix(row, column) * value(column)
        if !projected.isFinite then
          return Left(ConstrainedRayleighError.NonFiniteComputation("quadratic product", projected))
        column += 1
      result += value(row) * projected
      if !result.isFinite then return Left(ConstrainedRayleighError.NonFiniteComputation("quadratic form", result))
      row += 1
    Right(result)

  private def multiply(matrix: DMat, value: DVec, label: String): Either[ConstrainedRayleighError, DVec] =
    val out = new Array[Double](matrix.rows)
    var row = 0
    while row < matrix.rows do
      var result = 0.0
      var column = 0
      while column < matrix.cols do
        result += matrix(row, column) * value(column)
        if !result.isFinite then return Left(ConstrainedRayleighError.NonFiniteComputation(label, result))
        column += 1
      out(row) = result
      row += 1
    Right(DVec.fromArray(out))

  private def addScaled(left: DVec, right: DVec, scale: Double): Either[ConstrainedRayleighError, DVec] =
    val out = new Array[Double](left.length)
    var index = 0
    while index < left.length do
      out(index) = left(index) + scale * right(index)
      if !out(index).isFinite then
        return Left(ConstrainedRayleighError.NonFiniteComputation("line-search trial", out(index)))
      index += 1
    Right(DVec.fromArray(out))

  private def matrixScale(matrix: DMat): Double =
    var largest = 0.0
    var row = 0
    while row < matrix.rows do
      var column = 0
      while column < matrix.cols do
        val value = matrix(row, column)
        largest = Math.max(largest, Math.abs(value))
        column += 1
      row += 1
    largest

  private def geometryScale(geometry: Geometry, root: Double): Either[ConstrainedRayleighError, Double] =
    finite("gradient scale", Math.max(1.0, geometry.numeratorScale + Math.abs(root) * geometry.denominatorScale))

  private def residualScale(scale: Double, direction: DVec): Either[ConstrainedRayleighError, Double] =
    var directionScale = 0.0
    var index = 0
    while index < direction.length do
      directionScale = Math.max(directionScale, Math.abs(direction(index)))
      index += 1
    finite("stationarity scale", Math.max(1.0, scale * directionScale))

  private def finite(label: String, value: Double): Either[ConstrainedRayleighError, Double] =
    if value.isFinite then Right(value) else Left(ConstrainedRayleighError.NonFiniteComputation(label, value))

  private def sameVector(left: DVec, right: DVec): Boolean =
    var index = 0
    while index < left.length do
      if left(index) != right(index) then return false
      index += 1
    true
