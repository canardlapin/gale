package gale.linalg

/** Zero-cost evidence that a value is symmetric.
  *
  * The wrapper is an opaque subtype of `A`: it allocates nothing and retains the
  * complete underlying API, while APIs that care about symmetry can require the
  * stronger type explicitly.
  */
opaque type Symmetric[A] <: A = A

object Symmetric:
  private[gale] inline def unsafe[A](value: A): Symmetric[A] = value

  extension [A](wrapped: Symmetric[A])
    inline def value: A = wrapped

/** Zero-cost evidence that a value is symmetric positive-definite.
  *
  * Positive-definite evidence refines [[Symmetric]], so it can be passed wherever
  * either the property wrapper or the underlying value is required.
  */
opaque type PositiveDefinite[A] <: Symmetric[A] = A

object PositiveDefinite:
  private[gale] inline def unsafe[A](value: A): PositiveDefinite[A] = value

  extension [A](wrapped: PositiveDefinite[A])
    inline def value: A = wrapped

/** Zero-cost evidence that a value is lower triangular. */
opaque type LowerTriangular[A] <: A = A

object LowerTriangular:
  private[gale] inline def unsafe[A](value: A): LowerTriangular[A] = value

  extension [A](wrapped: LowerTriangular[A])
    inline def value: A = wrapped

/** Zero-cost evidence that a value is upper triangular. */
opaque type UpperTriangular[A] <: A = A

object UpperTriangular:
  private[gale] inline def unsafe[A](value: A): UpperTriangular[A] = value

  extension [A](wrapped: UpperTriangular[A])
    inline def value: A = wrapped

extension (matrix: DMat)
  /** Assert symmetry without inspecting the matrix. */
  inline def assumeSymmetric: Symmetric[DMat] =
    Symmetric.unsafe(matrix)

  /** Assert positive-definiteness without inspecting the matrix. */
  inline def assumePositiveDefinite: PositiveDefinite[DMat] =
    PositiveDefinite.unsafe(matrix)

  /** Assert lower-triangular structure without inspecting the matrix. */
  inline def assumeLowerTriangular: LowerTriangular[DMat] =
    LowerTriangular.unsafe(matrix)

  /** Assert upper-triangular structure without inspecting the matrix. */
  inline def assumeUpperTriangular: UpperTriangular[DMat] =
    UpperTriangular.unsafe(matrix)

  /** Verify square shape, finite entries, and pairwise symmetry.
    *
    * A pair `(aᵢⱼ, aⱼᵢ)` matches when its absolute difference is at most
    * `tolerance * max(1, |aᵢⱼ|, |aⱼᵢ|)`. This combined absolute/relative rule
    * keeps the check meaningful for both tiny and large matrices.
    */
  def verifySymmetric(tolerance: Double = MatrixPropertyVerification.DefaultSymmetryTolerance)
      : Either[LinAlgError, Symmetric[DMat]] =
    MatrixPropertyVerification.verifySymmetric(matrix, tolerance)

  /** Verify symmetry with the default tolerance, then prove positive-definiteness
    * through the existing Cholesky factorization.
    */
  def verifyPositiveDefinite: Either[LinAlgError, PositiveDefinite[DMat]] =
    MatrixPropertyVerification.verifyPositiveDefinite(
      matrix,
      MatrixPropertyVerification.DefaultSymmetryTolerance
    )

  /** Verify symmetry with `tolerance`, then prove positive-definiteness through
    * the existing Cholesky factorization.
    */
  def verifyPositiveDefinite(tolerance: Double): Either[LinAlgError, PositiveDefinite[DMat]] =
    MatrixPropertyVerification.verifyPositiveDefinite(matrix, tolerance)

  /** Verify exact lower-triangular structure. */
  def verifyLowerTriangular: Either[LinAlgError, LowerTriangular[DMat]] =
    MatrixPropertyVerification.verifyLowerTriangular(matrix, 0.0)

  /** Verify lower-triangular structure, allowing forbidden entries with
    * magnitude at most `tolerance`.
    */
  def verifyLowerTriangular(tolerance: Double): Either[LinAlgError, LowerTriangular[DMat]] =
    MatrixPropertyVerification.verifyLowerTriangular(matrix, tolerance)

  /** Verify exact upper-triangular structure. */
  def verifyUpperTriangular: Either[LinAlgError, UpperTriangular[DMat]] =
    MatrixPropertyVerification.verifyUpperTriangular(matrix, 0.0)

  /** Verify upper-triangular structure, allowing forbidden entries with
    * magnitude at most `tolerance`.
    */
  def verifyUpperTriangular(tolerance: Double): Either[LinAlgError, UpperTriangular[DMat]] =
    MatrixPropertyVerification.verifyUpperTriangular(matrix, tolerance)

private object MatrixPropertyVerification:
  inline val DefaultSymmetryTolerance = 1.0e-12

  private def validateTolerance(tolerance: Double): Either[LinAlgError, Unit] =
    if tolerance < 0.0 || !tolerance.isFinite then
      Left(LinAlgError.InvalidArgument(s"property tolerance must be finite and non-negative, got $tolerance"))
    else Right(())

  private def validateSquare(matrix: DMat): Either[LinAlgError, Unit] =
    if matrix.rows != matrix.cols then Left(LinAlgError.NonSquareMatrix(matrix.shape))
    else Right(())

  def verifySymmetric(matrix: DMat, tolerance: Double): Either[LinAlgError, Symmetric[DMat]] =
    validateTolerance(tolerance)
      .flatMap(_ => validateSquare(matrix))
      .flatMap(_ => scanSymmetry(matrix, tolerance))

  private def scanSymmetry(matrix: DMat, tolerance: Double): Either[LinAlgError, Symmetric[DMat]] =
    var i = 0
    while i < matrix.rows do
      val diagonal = matrix(i, i)
      if !diagonal.isFinite then
        return Left(LinAlgError.InvalidArgument(s"matrix has a non-finite diagonal entry at ($i, $i)"))
      var j = 0
      while j < i do
        val lower = matrix(i, j)
        val upper = matrix(j, i)
        if !lower.isFinite || !upper.isFinite then
          return Left(LinAlgError.InvalidArgument(s"matrix has a non-finite symmetric pair at ($i, $j) / ($j, $i)"))
        val scale = math.max(1.0, math.max(math.abs(lower), math.abs(upper)))
        if math.abs(lower - upper) > tolerance * scale then
          return Left(
            LinAlgError.InvalidArgument(
              s"matrix is not symmetric at ($i, $j) / ($j, $i): $lower != $upper within tolerance $tolerance"
            )
          )
        j += 1
      i += 1
    Right(Symmetric.unsafe(matrix))

  def verifyPositiveDefinite(matrix: DMat, tolerance: Double): Either[LinAlgError, PositiveDefinite[DMat]] =
    verifySymmetric(matrix, tolerance).flatMap: _ =>
      matrix.cholesky.map(_ => PositiveDefinite.unsafe(matrix))

  def verifyLowerTriangular(matrix: DMat, tolerance: Double): Either[LinAlgError, LowerTriangular[DMat]] =
    verifyTriangular(matrix, tolerance, lower = true).map(_ => LowerTriangular.unsafe(matrix))

  def verifyUpperTriangular(matrix: DMat, tolerance: Double): Either[LinAlgError, UpperTriangular[DMat]] =
    verifyTriangular(matrix, tolerance, lower = false).map(_ => UpperTriangular.unsafe(matrix))

  private def verifyTriangular(matrix: DMat, tolerance: Double, lower: Boolean): Either[LinAlgError, Unit] =
    validateTolerance(tolerance)
      .flatMap(_ => validateSquare(matrix))
      .flatMap(_ => scanTriangle(matrix, tolerance, lower))

  private def scanTriangle(matrix: DMat, tolerance: Double, lower: Boolean): Either[LinAlgError, Unit] =
    var i = 0
    while i < matrix.rows do
      var j = 0
      while j < matrix.cols do
        val forbidden = if lower then j > i else i > j
        if forbidden then
          val value = matrix(i, j)
          if !value.isFinite || math.abs(value) > tolerance then
            val property = if lower then "lower" else "upper"
            return Left(
              LinAlgError.InvalidArgument(
                s"matrix is not $property triangular: entry ($i, $j) is $value with tolerance $tolerance"
              )
            )
        j += 1
      i += 1
    Right(())
