package gale.numeric

import gale.linalg.OutwardInterval

enum PairedResidualInput:
  case Response, PreviousDesign, CandidateDesign, PreviousCoefficients, CandidateCoefficients

enum PairedResidualStage:
  case CoefficientDifference, ReturnedModels, NormalEquation, Gram, Gershgorin, OldGap, ProfileDifference

enum PairedResidualError:
  case InvalidDimensions(rows: Int, columns: Int)
  case MatrixSizeOverflow(rows: Int, columns: Int)
  case MissingInput(input: PairedResidualInput)
  case MissingScope
  case LengthMismatch(input: PairedResidualInput, expected: Int, actual: Int)
  case NonFiniteInput(input: PairedResidualInput, index: Int)
  case UnboundedArithmetic(stage: PairedResidualStage, row: Int, column: Int)

  def message: String = this match
    case InvalidDimensions(rows, columns) => s"paired residual comparison needs positive dimensions, got $rows x $columns"
    case MatrixSizeOverflow(rows, columns) => s"paired residual design storage overflows: $rows x $columns"
    case MissingInput(input) => s"paired residual input $input is absent"
    case MissingScope => "paired residual comparison scope is absent"
    case LengthMismatch(input, expected, actual) => s"paired residual input $input has length $actual, expected $expected"
    case NonFiniteInput(input, index) => s"paired residual input $input is nonfinite at $index"
    case UnboundedArithmetic(stage, row, column) => s"paired residual arithmetic is unbounded in $stage at ($row,$column)"

enum PairedResidualScope:
  case ReturnedModels, ProfileMinimum

enum ResidualDifferenceDirection:
  case Decrease, Increase, Unresolved

/** A bound on the previous coefficient error in squared-residual units.
  * For G0=D0' D0 and q0=D0' (z-D0 beta0), a positive certified Gram lower
  * eigenvalue implies RSS0(beta0)-min RSS0=q0' inverse(G0) q0, bounded by
  * ||q0||²/lambdaLower. Since min RSS1 <= RSS1(beta1), the stored-model
  * difference plus this old gap is an upper bound on min RSS1-min RSS0.
  * No candidate rank or candidate optimum needs to be assumed or computed.
  */
final class PreviousProfileMinimumBound private[numeric] (
    val gramEigenvalueLower: Double,
    val normalEquationSquaredUpper: Double,
    val suboptimalityUpper: Double,
    val differenceUpper: Double
)

enum ProfileMinimumUnresolved:
  case NonPositiveGramLowerBound(lower: Double)
  case DecreaseNotCertified(bound: PreviousProfileMinimumBound)

enum ProfileMinimumComparison:
  case NotRequested
  case CertifiedDecrease(bound: PreviousProfileMinimumBound)
  case Unresolved(reason: ProfileMinimumUnresolved)

/** Bounds concern the supplied, stored predictors and coefficients. No common
  * energy offset, prior, optimizer policy, refit, or altered coefficient enters
  * this comparison. Inputs must remain unchanged during the synchronous call.
  */
final class PairedResidualComparison private[numeric] (
    val modelDifference: OutwardInterval,
    val profileMinimum: ProfileMinimumComparison
):
  def direction: ResidualDifferenceDirection =
    if modelDifference.upper < 0.0 then ResidualDifferenceDirection.Decrease
    else if modelDifference.lower > 0.0 then ResidualDifferenceDirection.Increase
    else ResidualDifferenceDirection.Unresolved

  def certifiesProfileDecrease: Boolean = profileMinimum match
    case ProfileMinimumComparison.CertifiedDecrease(_) => true
    case _ => false

/** Thread-confined, reusable primitive scratch for paired row-major designs.
  * Array storage is exactly `2 * rows + 2 * columns` Double cells. Model work
  * is O(rows * columns); optional previous-Gram certification adds
  * O(rows * columns^2) work and no matrix-sized workspace. Outputs retain no
  * input or scratch aliases. The hot loops allocate no interval/tuple records.
  */
final class PairedResidualWorkspace private[numeric] (val rows: Int, val columns: Int):
  val arrayCells: Long = 2L * rows + 2L * columns
  private val betaDifferenceLower = new Array[Double](columns)
  private val betaDifferenceUpper = new Array[Double](columns)
  private val previousResidualLower = new Array[Double](rows)
  private val previousResidualUpper = new Array[Double](rows)
  private val arithmetic = new PairedResidualArithmetic

  def compare(
      response: Array[Double],
      previousDesign: Array[Double],
      candidateDesign: Array[Double],
      previousCoefficients: Array[Double],
      candidateCoefficients: Array[Double],
      scope: PairedResidualScope = PairedResidualScope.ProfileMinimum
  ): Either[PairedResidualError, PairedResidualComparison] =
    val cells = rows * columns // constructor admitted the Long product first
    for
      _ <- if scope == null then Left(PairedResidualError.MissingScope) else Right(())
      _ <- validate(response, rows, PairedResidualInput.Response)
      _ <- validate(previousDesign, cells, PairedResidualInput.PreviousDesign)
      _ <- validate(candidateDesign, cells, PairedResidualInput.CandidateDesign)
      _ <- validate(previousCoefficients, columns, PairedResidualInput.PreviousCoefficients)
      _ <- validate(candidateCoefficients, columns, PairedResidualInput.CandidateCoefficients)
      result <- compareFinite(response, previousDesign, candidateDesign, previousCoefficients, candidateCoefficients, scope)
    yield result

  private def validate(values: Array[Double], size: Int, input: PairedResidualInput): Either[PairedResidualError, Unit] =
    if values == null then Left(PairedResidualError.MissingInput(input))
    else if values.length != size then Left(PairedResidualError.LengthMismatch(input, size, values.length))
    else
      var i = 0
      while i < size do
        if !values(i).isFinite then return Left(PairedResidualError.NonFiniteInput(input, i))
        i += 1
      Right(())

  private def compareFinite(
      response: Array[Double], d0: Array[Double], d1: Array[Double],
      b0: Array[Double], b1: Array[Double], scope: PairedResidualScope
  ): Either[PairedResidualError, PairedResidualComparison] =
    arithmetic.reset()
    var column = 0
    while column < columns do
      arithmetic.subtract(b1(column), b1(column), b0(column), b0(column))
      if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.CoefficientDifference, -1, column))
      betaDifferenceLower(column) = arithmetic.lower
      betaDifferenceUpper(column) = arithmetic.upper
      column += 1
    var totalLower = 0.0
    var totalUpper = 0.0
    var row = 0
    while row < rows do
      var predictionLower = 0.0
      var predictionUpper = 0.0
      var deltaLower = 0.0
      var deltaUpper = 0.0
      var normalResidualLower = response(row)
      var normalResidualUpper = response(row)
      column = 0
      while column < columns && arithmetic.valid do
        val index = row * columns + column
        val old = d0(index)
        arithmetic.multiply(old, old, b0(column), b0(column))
        val termLower = arithmetic.lower
        val termUpper = arithmetic.upper
        arithmetic.add(predictionLower, predictionUpper, termLower, termUpper)
        predictionLower = arithmetic.lower
        predictionUpper = arithmetic.upper
        if scope == PairedResidualScope.ProfileMinimum then
          arithmetic.subtract(normalResidualLower, normalResidualUpper, termLower, termUpper)
          normalResidualLower = arithmetic.lower
          normalResidualUpper = arithmetic.upper
        arithmetic.subtract(d1(index), d1(index), old, old)
        val designLower = arithmetic.lower
        val designUpper = arithmetic.upper
        arithmetic.multiply(designLower, designUpper, b1(column), b1(column))
        arithmetic.add(deltaLower, deltaUpper, arithmetic.lower, arithmetic.upper)
        deltaLower = arithmetic.lower
        deltaUpper = arithmetic.upper
        arithmetic.multiply(old, old, betaDifferenceLower(column), betaDifferenceUpper(column))
        arithmetic.add(deltaLower, deltaUpper, arithmetic.lower, arithmetic.upper)
        deltaLower = arithmetic.lower
        deltaUpper = arithmetic.upper
        column += 1
      if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.ReturnedModels, row, column))
      previousResidualLower(row) = normalResidualLower
      previousResidualUpper(row) = normalResidualUpper
      arithmetic.subtract(response(row), response(row), predictionLower, predictionUpper)
      arithmetic.multiply(2.0, 2.0, arithmetic.lower, arithmetic.upper)
      arithmetic.subtract(deltaLower, deltaUpper, arithmetic.lower, arithmetic.upper)
      arithmetic.multiply(deltaLower, deltaUpper, arithmetic.lower, arithmetic.upper)
      arithmetic.add(totalLower, totalUpper, arithmetic.lower, arithmetic.upper)
      if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.ReturnedModels, row, -1))
      totalLower = arithmetic.lower
      totalUpper = arithmetic.upper
      row += 1
    val difference = OutwardInterval(totalLower, totalUpper).toOption.get // finite ordered arithmetic invariant
    scope match
      case PairedResidualScope.ReturnedModels =>
        Right(new PairedResidualComparison(difference, ProfileMinimumComparison.NotRequested))
      case PairedResidualScope.ProfileMinimum =>
        profileBound(d0, difference).map(profile => new PairedResidualComparison(difference, profile))

  /** Certify only an upper bound: min RSS1-min RSS0 <= deltaRSSUpper+oldGapUpper.
    * Gershgorin positivity is sufficient and incomplete. An unavailable lower
    * bound is unresolved, even when the actual Gram matrix is positive definite.
    */
  private def profileBound(d0: Array[Double], difference: OutwardInterval): Either[PairedResidualError, ProfileMinimumComparison] =
    var normalSquaredUpper = 0.0
    var left = 0
    while left < columns do
      var lower = 0.0
      var upper = 0.0
      var row = 0
      while row < rows && arithmetic.valid do
        val value = d0(row * columns + left)
        arithmetic.multiply(value, value, previousResidualLower(row), previousResidualUpper(row))
        arithmetic.add(lower, upper, arithmetic.lower, arithmetic.upper)
        lower = arithmetic.lower
        upper = arithmetic.upper
        row += 1
      val magnitude = math.max(math.abs(lower), math.abs(upper))
      arithmetic.multiply(magnitude, magnitude, magnitude, magnitude)
      arithmetic.add(0.0, normalSquaredUpper, arithmetic.lower, arithmetic.upper)
      if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.NormalEquation, -1, left))
      normalSquaredUpper = arithmetic.upper
      left += 1
    var eigenvalueLower = Double.PositiveInfinity
    left = 0
    while left < columns do
      var diagonalLower = 0.0
      var offDiagonalUpper = 0.0
      var right = 0
      while right < columns do
        var lower = 0.0
        var upper = 0.0
        var row = 0
        while row < rows && arithmetic.valid do
          val x = d0(row * columns + left)
          val y = d0(row * columns + right)
          arithmetic.multiply(x, x, y, y)
          arithmetic.add(lower, upper, arithmetic.lower, arithmetic.upper)
          lower = arithmetic.lower
          upper = arithmetic.upper
          row += 1
        if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.Gram, left, right))
        if left == right then diagonalLower = lower
        else
          val magnitude = math.max(math.abs(lower), math.abs(upper))
          arithmetic.add(0.0, offDiagonalUpper, magnitude, magnitude)
          if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.Gershgorin, left, right))
          offDiagonalUpper = arithmetic.upper
        right += 1
      arithmetic.subtract(diagonalLower, diagonalLower, offDiagonalUpper, offDiagonalUpper)
      if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.Gershgorin, left, -1))
      eigenvalueLower = math.min(eigenvalueLower, arithmetic.lower)
      left += 1
    if eigenvalueLower <= 0.0 then
      Right(ProfileMinimumComparison.Unresolved(ProfileMinimumUnresolved.NonPositiveGramLowerBound(eigenvalueLower)))
    else
      arithmetic.dividePositiveUpper(normalSquaredUpper, eigenvalueLower)
      if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.OldGap, -1, -1))
      val gapUpper = arithmetic.upper
      arithmetic.add(difference.upper, difference.upper, gapUpper, gapUpper)
      if !arithmetic.valid then return Left(PairedResidualError.UnboundedArithmetic(PairedResidualStage.ProfileDifference, -1, -1))
      val bound = new PreviousProfileMinimumBound(eigenvalueLower, normalSquaredUpper, gapUpper, arithmetic.upper)
      if bound.differenceUpper < 0.0 then Right(ProfileMinimumComparison.CertifiedDecrease(bound))
      else Right(ProfileMinimumComparison.Unresolved(ProfileMinimumUnresolved.DecreaseNotCertified(bound)))

object PairedResidualWorkspace:
  /** Allocation-free admission shape for callers enforcing their own budgets. */
  def requiredArrayCells(rows: Int, columns: Int): Either[PairedResidualError, Long] =
    if rows <= 0 || columns <= 0 then Left(PairedResidualError.InvalidDimensions(rows, columns))
    else if rows.toLong * columns.toLong > Int.MaxValue then Left(PairedResidualError.MatrixSizeOverflow(rows, columns))
    else Right(2L * rows + 2L * columns)

  def apply(rows: Int, columns: Int): Either[PairedResidualError, PairedResidualWorkspace] =
    requiredArrayCells(rows, columns).map(_ => new PairedResidualWorkspace(rows, columns))

/** Mutable primitive endpoints; invalid arithmetic never becomes a result. */
private final class PairedResidualArithmetic:
  var lower = 0.0
  var upper = 0.0
  var valid = true

  def reset(): Unit =
    lower = 0.0
    upper = 0.0
    valid = true

  private def rounded(lo: Double, hi: Double): Unit =
    val nextLower = java.lang.Math.nextDown(lo)
    val nextUpper = java.lang.Math.nextUp(hi)
    if !lo.isFinite || !hi.isFinite || !nextLower.isFinite || !nextUpper.isFinite then valid = false
    else
      lower = nextLower
      upper = nextUpper

  def add(aLower: Double, aUpper: Double, bLower: Double, bUpper: Double): Unit =
    rounded(aLower + bLower, aUpper + bUpper)

  def subtract(aLower: Double, aUpper: Double, bLower: Double, bUpper: Double): Unit =
    rounded(aLower - bUpper, aUpper - bLower)

  def multiply(aLower: Double, aUpper: Double, bLower: Double, bUpper: Double): Unit =
    val a = aLower * bLower
    val b = aLower * bUpper
    val c = aUpper * bLower
    val d = aUpper * bUpper
    if !a.isFinite || !b.isFinite || !c.isFinite || !d.isFinite then valid = false
    else rounded(math.min(math.min(a, b), math.min(c, d)), math.max(math.max(a, b), math.max(c, d)))

  def dividePositiveUpper(numerator: Double, denominator: Double): Unit =
    rounded(0.0, numerator / denominator)
