package gale.linalg

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.sparse.Banded

/** Conditioning diagnostics for a banded SPD factor.
  *
  * The inverse norm bound comes from absolute-value triangular comparison
  * systems, not a pivot-ratio heuristic. It bounds the true inverse 1-norm in
  * exact arithmetic; floating-point evaluation is not a directed-rounding
  * certificate. Overflow gives infinity (and reciprocal bound zero).
  */
final case class BandedConditioning(
    matrixOneNorm: Double,
    inverseOneNormUpperBound: Double,
    conditionNumberUpperBound: Double,
    pivotRatio: Double
):
  def reciprocalConditionLowerBound: Double = 1.0 / conditionNumberUpperBound

private enum BandedSolve:
  case Symmetric, Lower, LowerTranspose

/** Immutable packed Cholesky factor of `A = L Lᵀ`.
  *
  * `lowerBands(i, d)` stores `L(i, i-d)` for `0 <= d <= bandwidth`;
  * positions with `d > i` are zero padding. Storage is `n * (bandwidth+1)`,
  * factorization costs O(n bandwidth²), and each RHS costs O(n bandwidth).
  * No operation expands the band to a dense square factor. Pure solves
  * allocate only their requested result.
  *
  * Pure solves preserve their RHS. In-place solves reuse a MutableDVec or an
  * open DMatBuilder without allocating numerical scratch. Dimension and
  * nonfinite-input errors precede writes; arithmetic overflow may leave a
  * partially solved destination. A closed builder follows its usual throwing
  * ownership contract. The factor can be shared between independent callers.
  */
final class BandedCholesky private (
    val lowerBands: DMat,
    val logDet: Double,
    val options: CholeskyOptions,
    private val matrixNorm: Double,
    private val relativePivot: Double
) extends ExactSolveFactor:
  val size: Int = lowerBands.rows
  val bandwidth: Int = lowerBands.cols - 1
  val diagnostics: FactorizationDiagnostics =
    FactorizationDiagnostics(rank = Some(size), rankTolerance = Some(options.pivotTolerance))

  private val data = lowerBands.data
  private val ld = lowerBands.cols

  def solve(b: DVec): Either[LinAlgError, DVec] = copySolve(b, BandedSolve.Symmetric)
  def solve(b: DMat): Either[LinAlgError, DMat] = copySolve(b, BandedSolve.Symmetric)
  /** A is symmetric, so its transpose solve is the same system. */
  def solveTranspose(b: DVec): Either[LinAlgError, DVec] = solve(b)
  def solveTranspose(b: DMat): Either[LinAlgError, DMat] = solve(b)
  def solveLower(b: DVec): Either[LinAlgError, DVec] = copySolve(b, BandedSolve.Lower)
  def solveLower(b: DMat): Either[LinAlgError, DMat] = copySolve(b, BandedSolve.Lower)
  def solveLowerTranspose(b: DVec): Either[LinAlgError, DVec] = copySolve(b, BandedSolve.LowerTranspose)
  def solveLowerTranspose(b: DMat): Either[LinAlgError, DMat] = copySolve(b, BandedSolve.LowerTranspose)

  def solveInPlace(b: MutableDVec): Either[LinAlgError, Unit] = inPlace(b, BandedSolve.Symmetric)
  def solveInPlace(b: DMatBuilder): Either[LinAlgError, Unit] = inPlace(b, BandedSolve.Symmetric)
  def solveTransposeInPlace(b: MutableDVec): Either[LinAlgError, Unit] = solveInPlace(b)
  def solveTransposeInPlace(b: DMatBuilder): Either[LinAlgError, Unit] = solveInPlace(b)
  def solveLowerInPlace(b: MutableDVec): Either[LinAlgError, Unit] = inPlace(b, BandedSolve.Lower)
  def solveLowerInPlace(b: DMatBuilder): Either[LinAlgError, Unit] = inPlace(b, BandedSolve.Lower)
  def solveLowerTransposeInPlace(b: MutableDVec): Either[LinAlgError, Unit] = inPlace(b, BandedSolve.LowerTranspose)
  def solveLowerTransposeInPlace(b: DMatBuilder): Either[LinAlgError, Unit] = inPlace(b, BandedSolve.LowerTranspose)

  private def copySolve(b: DVec, mode: BandedSolve): Either[LinAlgError, DVec] =
    if b.length != size then Left(LinAlgError.VectorLengthMismatch(size, b.length))
    else
      val out = b.toDoubleArrayOwnedCopy
      run(out, 0, 1, 1, 0, mode).map(_ => DVec.fromDoubleArrayOwned(out))

  private def copySolve(b: DMat, mode: BandedSolve): Either[LinAlgError, DMat] =
    if b.rows != size then Left(LinAlgError.DimensionMismatch(Shape(Rows(size), Cols(b.cols)), b.shape))
    else
      val out = b.toDoubleArrayCopyRowMajor
      run(out, 0, b.cols, b.cols, 1, mode).map(_ => DMat.fromDoubleArrayOwned(size, b.cols, out))

  private def inPlace(b: MutableDVec, mode: BandedSolve): Either[LinAlgError, Unit] =
    if b.length != size then Left(LinAlgError.VectorLengthMismatch(size, b.length))
    else run(b.data, b.offset.value, b.stride.value, 1, 0, mode)

  private def inPlace(b: DMatBuilder, mode: BandedSolve): Either[LinAlgError, Unit] =
    if b.rows != size then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(size), Cols(b.cols)), Shape(Rows(b.rows), Cols(b.cols))))
    else run(b.writableData, 0, b.cols, b.cols, 1, mode)

  private def run(
      rhs: DoubleArray,
      offset: Int,
      rowStride: Int,
      columns: Int,
      colStride: Int,
      mode: BandedSolve
  ): Either[LinAlgError, Unit] =
    var i = 0
    while i < size do
      var c = 0
      while c < columns do
        if !rhs(offset + i * rowStride + c * colStride).isFinite then
          return Left(LinAlgError.InvalidArgument(s"nonfinite banded solve RHS at ($i,$c)"))
        c += 1
      i += 1

    if mode != BandedSolve.LowerTranspose then
      i = 0
      while i < size do
        val dst = offset + i * rowStride
        var j = math.max(0, i - bandwidth)
        while j < i do
          val lij = data(i * ld + i - j)
          val src = offset + j * rowStride
          var c = 0
          while c < columns do
            rhs(dst + c * colStride) = rhs(dst + c * colStride) - lij * rhs(src + c * colStride)
            c += 1
          j += 1
        var c = 0
        while c < columns do
          val index = dst + c * colStride
          rhs(index) = rhs(index) / data(i * ld)
          if !rhs(index).isFinite then return Left(LinAlgError.InvalidArgument(s"nonfinite lower solve at ($i,$c)"))
          c += 1
        i += 1

    if mode != BandedSolve.Lower then
      i = size - 1
      while i >= 0 do
        val dst = offset + i * rowStride
        var j = i + 1
        val end = math.min(size - 1, i + bandwidth)
        while j <= end do
          val lji = data(j * ld + j - i)
          val src = offset + j * rowStride
          var c = 0
          while c < columns do
            rhs(dst + c * colStride) = rhs(dst + c * colStride) - lji * rhs(src + c * colStride)
            c += 1
          j += 1
        var c = 0
        while c < columns do
          val index = dst + c * colStride
          rhs(index) = rhs(index) / data(i * ld)
          if !rhs(index).isFinite then return Left(LinAlgError.InvalidArgument(s"nonfinite transpose solve at ($i,$c)"))
          c += 1
        i -= 1
    Right(())

  /** O(n bandwidth) comparison bound, computed once with O(n) temporary storage.
    * For diagonal A the bound is the exact condition number; cancellation in
    * triangular inverses can make it conservative for wider bands.
    */
  lazy val conditioning: BandedConditioning =
    if size == 0 then BandedConditioning(0.0, 0.0, 1.0, 1.0)
    else
      val bound = DoubleArray.alloc(size)
      var i = 0
      while i < size do
        var sum = 1.0
        var j = math.max(0, i - bandwidth)
        while j < i do
          val a = math.abs(data(i * ld + i - j))
          if a > 0.0 then sum += a * bound(j)
          j += 1
        bound(i) = sum / data(i * ld)
        i += 1
      var inverseNorm = 0.0
      i = size - 1
      while i >= 0 do
        var sum = bound(i)
        var j = i + 1
        val end = math.min(size - 1, i + bandwidth)
        while j <= end do
          val a = math.abs(data(j * ld + j - i))
          if a > 0.0 then sum += a * bound(j)
          j += 1
        bound(i) = sum / data(i * ld)
        inverseNorm = math.max(inverseNorm, bound(i))
        i -= 1
      BandedConditioning(matrixNorm, inverseNorm, math.max(1.0, matrixNorm * inverseNorm), relativePivot)

object BandedCholesky:
  /** Factor a lower band packed by row: `bands(i,d) = A(i,i-d)`, diagonal
    * in column zero. Padding at d > i is ignored, including nonfinite padding.
    * Active entries must be finite. The input may be strided and is preserved.
    * A zero-sized system uses a 0-by-1 band; otherwise 0 <= bandwidth < n.
    */
  def factorLower(bands: DMat, options: CholeskyOptions = CholeskyOptions.Default): Either[LinAlgError, BandedCholesky] =
    validateShape(bands.rows, bands.cols).flatMap: _ =>
      fromOwnedLower(bands.rows, bands.cols, bands.toDoubleArrayCopyRowMajor, options)

  /** Read the lower band of the existing general-band representation, with
    * the same lower-triangle symmetry convention as dense Cholesky.
    */
  def factorLower(matrix: Banded): Either[LinAlgError, BandedCholesky] = factorLower(matrix, CholeskyOptions.Default)

  def factorLower(matrix: Banded, options: CholeskyOptions): Either[LinAlgError, BandedCholesky] =
    if matrix.rows != matrix.cols then Left(LinAlgError.NonSquareMatrix(Shape(Rows(matrix.rows), Cols(matrix.cols))))
    else
      val n = matrix.rows
      val width = if n == 0 then 1 else math.min(matrix.kl, n - 1) + 1
      val bands = DMatBuilder.zeros(n, width)
      var i = 0
      while i < n do
        var d = 0
        while d < width && d <= i do
          bands(i, d) = matrix(i, i - d)
          d += 1
        i += 1
      bands.consumeBandedCholesky(options)

  private def validateShape(n: Int, width: Int): Either[LinAlgError, Unit] =
    if width < 1 || (n == 0 && width != 1) || (n > 0 && width > n) then
      Left(LinAlgError.InvalidArgument(s"lower band shape must be n-by-(b+1) with 0 <= b < n (0-by-1 if empty), got $n-by-$width"))
    else Right(())

  private[gale] def fromOwnedLower(
      n: Int,
      width: Int,
      bands: DoubleArray,
      options: CholeskyOptions
  ): Either[LinAlgError, BandedCholesky] =
    validateShape(n, width) match
      case Left(error) => return Left(error)
      case Right(_) => ()
    val bandwidth = width - 1
    val rowSums = DoubleArray.alloc(n)
    var i = 0
    while i < n do
      var d = 0
      while d < width do
        val index = i * width + d
        if d > i then bands(index) = 0.0
        else
          val value = bands(index)
          if !value.isFinite then return Left(LinAlgError.InvalidArgument(s"nonfinite lower band at ($i,$d)"))
          rowSums(i) = rowSums(i) + math.abs(value)
          if d > 0 then rowSums(i - d) = rowSums(i - d) + math.abs(value)
        d += 1
      i += 1
    var norm = 0.0
    i = 0
    while i < n do
      norm = math.max(norm, rowSums(i))
      i += 1
    var logDet = 0.0
    var minPivot = Double.PositiveInfinity
    var maxPivot = 0.0
    i = 0
    while i < n do
      val start = math.max(0, i - bandwidth)
      var j = start
      while j < i do
        var value = bands(i * width + i - j)
        var k = start
        while k < j do
          value -= bands(i * width + i - k) * bands(j * width + j - k)
          k += 1
        value /= bands(j * width)
        if !value.isFinite then return Left(LinAlgError.InvalidArgument(s"nonfinite banded factor at ($i,$j)"))
        bands(i * width + i - j) = value
        j += 1
      var pivot = bands(i * width)
      j = start
      while j < i do
        val value = bands(i * width + i - j)
        pivot -= value * value
        j += 1
      if !pivot.isFinite then return Left(LinAlgError.InvalidArgument(s"nonfinite banded pivot at $i"))
      if pivot <= options.pivotTolerance then return Left(LinAlgError.NotPositiveDefinite(i))
      minPivot = math.min(minPivot, pivot)
      maxPivot = math.max(maxPivot, pivot)
      bands(i * width) = math.sqrt(pivot)
      logDet += math.log(pivot)
      i += 1
    Right(new BandedCholesky(
      DMat.fromDoubleArrayOwned(n, width, bands),
      logDet,
      options,
      norm,
      if n == 0 then 1.0 else minPivot / maxPivot
    ))
