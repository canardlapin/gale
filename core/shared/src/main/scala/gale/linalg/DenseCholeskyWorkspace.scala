package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Allocation-free classification for repeated curvature checks. */
enum DensePositiveDefiniteness:
  case PositiveDefinite, NotPositiveDefinite, NonFiniteInput, InvalidSpan

/** Single-owner scratch for repeated dense SPD operations on row-major arrays.
  *
  * Successful calls allocate no objects or arrays after construction. Scratch is
  * owned here and never exposed; caller arrays are copied before computation.
  * Only the lower triangle participates. Failures leave caller arrays unchanged.
  * Instances are mutable and must not be shared between concurrent executions.
  * This is an execution resource, not an immutable factor or matrix value.
  */
final class DenseCholeskyWorkspace(val capacity: Int, val options: CholeskyOptions = CholeskyOptions.Default):
  require(capacity >= 0 && capacity.toLong * capacity <= Int.MaxValue, "invalid dense Cholesky capacity")
  private val lower = DoubleArray.alloc(capacity * capacity)
  private val rhs = DoubleArray.alloc(capacity)
  private val success: Either[LinAlgError, Unit] = Right(())

  /** Classify a lower-triangle SPD input without modifying it or allocating,
    * including the ordinary non-SPD path used by curvature checks.
    */
  def testPositiveDefinite(size: Int, values: Array[Double], offset: Int = 0): DensePositiveDefiniteness =
    if !matrixSpan(size, values, offset) then DensePositiveDefiniteness.InvalidSpan
    else
      var row = 0
      while row < size do
        var col = 0
        while col <= row do
          val value = values(offset + row * size + col)
          if !value.isFinite then return DensePositiveDefiniteness.NonFiniteInput
          lower(row * size + col) = value
          col += 1
        row += 1
      if DoubleKernels.dpotrfLower(size, lower, options.pivotTolerance) >= 0 then
        DensePositiveDefiniteness.NotPositiveDefinite
      else DensePositiveDefiniteness.PositiveDefinite

  /** Test positive definiteness without changing the input, using the configured
    * absolute pivot tolerance. Upper-triangle entries are ignored, including NaN.
    */
  def checkPositiveDefinite(size: Int, values: Array[Double], offset: Int = 0): Either[LinAlgError, Unit] =
    if !matrixSpan(size, values, offset) then invalidSpan()
    else factor(size, values, offset)

  /** Replace only the lower triangle with L, where A = L L^T. On failure no
    * element of `values` is changed; the upper triangle always remains unchanged.
    */
  def factorLowerInPlace(size: Int, values: Array[Double], offset: Int = 0): Either[LinAlgError, Unit] =
    if !matrixSpan(size, values, offset) then invalidSpan()
    else
      val result = factor(size, values, offset)
      if result.isRight then
        var row = 0
        while row < size do
          var col = 0
          while col <= row do
            values(offset + row * size + col) = lower(row * size + col)
            col += 1
          row += 1
      result

  /** Solve (L L^T) x = b using a stored row-major lower factor, replacing b
    * only on success. L must have finite entries and strictly positive diagonal.
    * Arbitrary valid array spans, including overlapping source/output, are safe:
    * all inputs are copied before any output is written.
    */
  def solveLowerInPlace(
      size: Int,
      factorValues: Array[Double],
      values: Array[Double],
      factorOffset: Int = 0,
      offset: Int = 0
  ): Either[LinAlgError, Unit] =
    if !matrixSpan(size, factorValues, factorOffset) || offset < 0 || offset.toLong + size > values.length then
      invalidSpan()
    else
      var row = 0
      while row < size do
        var col = 0
        while col <= row do
          val value = factorValues(factorOffset + row * size + col)
          if !value.isFinite then return Left(LinAlgError.InvalidArgument("non-finite lower factor"))
          if row == col && value <= 0.0 then return Left(LinAlgError.NotPositiveDefinite(row))
          lower(row * size + col) = value
          col += 1
        val value = values(offset + row)
        if !value.isFinite then return Left(LinAlgError.InvalidArgument("non-finite right-hand side"))
        rhs(row) = value
        row += 1
      // The transposed lower factor is upper triangular, with swapped strides.
      val forward = DoubleKernels.dtrsv(size, true, false, 0.0, lower, 0, size, 1, rhs, 0, 1)
      val backward = DoubleKernels.dtrsv(size, false, false, 0.0, lower, 0, 1, size, rhs, 0, 1)
      if forward >= 0 then Left(LinAlgError.SingularMatrix(forward))
      else if backward >= 0 then Left(LinAlgError.SingularMatrix(backward))
      else
        row = 0
        while row < size do
          if !rhs(row).isFinite then return Left(LinAlgError.InvalidArgument("non-finite Cholesky solution"))
          row += 1
        row = 0
        while row < size do
          values(offset + row) = rhs(row)
          row += 1
        success

  private def matrixSpan(size: Int, values: Array[Double], offset: Int): Boolean =
    size >= 0 && size <= capacity && offset >= 0 && offset.toLong + size.toLong * size <= values.length

  private def invalidSpan(): Either[LinAlgError, Unit] =
    Left(LinAlgError.InvalidArgument("dense Cholesky array span or workspace capacity mismatch"))

  private def factor(size: Int, values: Array[Double], offset: Int): Either[LinAlgError, Unit] =
    var row = 0
    while row < size do
      var col = 0
      while col <= row do
        val value = values(offset + row * size + col)
        if !value.isFinite then return Left(LinAlgError.InvalidArgument("non-finite SPD lower triangle"))
        lower(row * size + col) = value
        col += 1
      row += 1
    val pivot = DoubleKernels.dpotrfLower(size, lower, options.pivotTolerance)
    if pivot >= 0 then Left(LinAlgError.NotPositiveDefinite(pivot)) else success
