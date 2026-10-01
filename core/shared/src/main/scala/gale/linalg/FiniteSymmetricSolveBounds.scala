package gale.linalg

/** One entry in the canonical row-major upper triangle of a symmetric matrix.
  *
  * Every diagonal and off-diagonal entry must be present, including intervals containing exactly zero. An off-diagonal
  * entry represents both symmetric positions.
  */
final case class UpperIntervalEntry(row: Int, column: Int, value: OutwardInterval)

enum FiniteSolveBoundError:
  case NonPositiveSize(size: Int)
  case RhsLengthMismatch(expected: Int, actual: Int)
  case CandidateLengthMismatch(expected: Int, actual: Int)
  case NonFiniteCandidate(index: Int)
  case MissingUpperTriangleEntry(row: Int, column: Int)
  case NonCanonicalUpperTriangleEntry(expectedRow: Int, expectedColumn: Int, actualRow: Int, actualColumn: Int)
  case ExtraUpperTriangleEntry(row: Int, column: Int)
  case ArithmeticFailure(cause: OutwardIntervalError)
  case InsufficientGershgorinMargin(marginLower: Double)

/** An immutable certificate for precisely the matrix enclosure, right-hand side enclosure, and candidate consumed by
  * [[FiniteSymmetricSolveBounds]].
  */
final class FiniteSolveBound private[linalg] (
    val size: Int,
    val upperTriangleEntryCount: Long,
    val gershgorinMarginLower: Double,
    val residualL1Upper: Double,
    val solutionError2Upper: Double,
    val arithmeticPolicy: String
)

/** Streaming Gershgorin/residual certificate for a finite symmetric solve.
  *
  * The certificate accepts only strictly diagonally dominant lower bounds, which is sufficient for positive
  * definiteness but deliberately incomplete. It retains O(n) interval state and consumes the canonical upper triangle
  * exactly once.
  */
object FiniteSymmetricSolveBounds:
  val ArithmeticPolicy = "finite IEEE-754 outward interval primitives; one-pass canonical upper triangle"

  def certify(
      size: Int,
      upperTriangle: Iterator[UpperIntervalEntry],
      rhs: Vector[OutwardInterval],
      candidate: Vector[Double]
  ): Either[FiniteSolveBoundError, FiniteSolveBound] =
    if size <= 0 then Left(FiniteSolveBoundError.NonPositiveSize(size))
    else if rhs.length != size then Left(FiniteSolveBoundError.RhsLengthMismatch(size, rhs.length))
    else if candidate.length != size then Left(FiniteSolveBoundError.CandidateLengthMismatch(size, candidate.length))
    else
      var index = 0
      while index < size && candidate(index).isFinite do index += 1
      if index < size then Left(FiniteSolveBoundError.NonFiniteCandidate(index))
      else certifyFinite(size, upperTriangle, rhs, candidate)

  private def certifyFinite(
      size: Int,
      upperTriangle: Iterator[UpperIntervalEntry],
      rhs: Vector[OutwardInterval],
      candidate: Vector[Double]
  ): Either[FiniteSolveBoundError, FiniteSolveBound] =
    val zero = OutwardInterval.point(0.0).toOption.get
    val candidateIntervals = new Array[OutwardInterval](size)
    val rowOffDiagonalAbsoluteUpper = Array.fill(size)(zero)
    val residual = new Array[OutwardInterval](size)
    val diagonalLower = new Array[Double](size)

    var i = 0
    while i < size do
      candidateIntervals(i) = OutwardInterval.point(candidate(i)).toOption.get
      residual(i) = interval(zero - rhs(i)) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      i += 1

    var entries = 0L
    var row = 0
    while row < size do
      var column = row
      while column < size do
        if !upperTriangle.hasNext then return Left(FiniteSolveBoundError.MissingUpperTriangleEntry(row, column))
        val entry = upperTriangle.next()
        if entry.row != row || entry.column != column then
          return Left(FiniteSolveBoundError.NonCanonicalUpperTriangleEntry(row, column, entry.row, entry.column))
        entries += 1
        if row == column then diagonalLower(row) = entry.value.lower
        else
          val absoluteUpper = math.max(math.abs(entry.value.lower), math.abs(entry.value.upper))
          val absoluteInterval = OutwardInterval.point(absoluteUpper).toOption.get
          rowOffDiagonalAbsoluteUpper(row) = interval(rowOffDiagonalAbsoluteUpper(row) + absoluteInterval) match
            case Left(error)  => return Left(error)
            case Right(value) => value
          rowOffDiagonalAbsoluteUpper(column) = interval(rowOffDiagonalAbsoluteUpper(column) + absoluteInterval) match
            case Left(error)  => return Left(error)
            case Right(value) => value

        val leftContribution = interval(entry.value * candidateIntervals(column)) match
          case Left(error)  => return Left(error)
          case Right(value) => value
        residual(row) = interval(residual(row) + leftContribution) match
          case Left(error)  => return Left(error)
          case Right(value) => value
        if row != column then
          val rightContribution = interval(entry.value * candidateIntervals(row)) match
            case Left(error)  => return Left(error)
            case Right(value) => value
          residual(column) = interval(residual(column) + rightContribution) match
            case Left(error)  => return Left(error)
            case Right(value) => value
        column += 1
      row += 1

    if upperTriangle.hasNext then
      val extra = upperTriangle.next()
      Left(FiniteSolveBoundError.ExtraUpperTriangleEntry(extra.row, extra.column))
    else
      var marginLower = Double.PositiveInfinity
      i = 0
      while i < size do
        val margin = interval(
          OutwardInterval.point(diagonalLower(i)).toOption.get - rowOffDiagonalAbsoluteUpper(i)
        ) match
          case Left(error)  => return Left(error)
          case Right(value) => value.lower
        if margin < marginLower then marginLower = margin
        i += 1
      if !(marginLower > 0.0) then Left(FiniteSolveBoundError.InsufficientGershgorinMargin(marginLower))
      else
        var residualL1 = zero
        i = 0
        while i < size do
          val absoluteUpper = math.max(math.abs(residual(i).lower), math.abs(residual(i).upper))
          residualL1 = interval(residualL1 + OutwardInterval.point(absoluteUpper).toOption.get) match
            case Left(error)  => return Left(error)
            case Right(value) => value
          i += 1
        val errorUpper = interval(
          OutwardInterval.point(residualL1.upper).toOption.get / OutwardInterval.point(marginLower).toOption.get
        ) match
          case Left(error)  => return Left(error)
          case Right(value) => value.upper
        Right(new FiniteSolveBound(size, entries, marginLower, residualL1.upper, errorUpper, ArithmeticPolicy))

  private def interval[A](value: Either[OutwardIntervalError, A]): Either[FiniteSolveBoundError, A] =
    value.left.map(FiniteSolveBoundError.ArithmeticFailure.apply)
