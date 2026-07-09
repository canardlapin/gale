package gale.linalg

import gale.kernel.DoubleKernels

/** Standalone triangular solves `L x = b` (lower) and `U x = b` (upper) for a
  * dense triangular matrix `A`.
  *
  * `A` is read as-is: only the triangle named by the entry point participates
  * (the opposite triangle is ignored, not required to be zero), and the diagonal
  * is used directly — a zero on the diagonal makes the system singular and yields
  * `Left(`[[LinAlgError.SingularMatrix]]`)` at that pivot. `A` may be any dense
  * view (strided, transposed, or a submatrix); the solve honours its strides.
  *
  * Both methods share the single strided substitution kernel
  * `gale.kernel.DoubleKernels.dtrsv`, which also backs the LU / Cholesky / QR
  * back-substitutions, so there is one triangular-solve implementation.
  */
object TriangularSolve:
  /** Solve `L x = b` for a lower-triangular `A` (forward substitution). */
  def lower(A: DMat, b: DVec): Either[LinAlgError, DVec] =
    solve(A, b, lower = true)

  /** Solve `U x = b` for an upper-triangular `A` (back substitution). */
  def upper(A: DMat, b: DVec): Either[LinAlgError, DVec] =
    solve(A, b, lower = false)

  private def solve(A: DMat, b: DVec, lower: Boolean): Either[LinAlgError, DVec] =
    val n = A.rows
    if A.cols != n then
      Left(LinAlgError.NonSquareMatrix(A.shape))
    else if b.length != n then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(n), Cols(1)), Shape(Rows(b.length), Cols(1))))
    else
      // Owned contiguous copy of b (offset 0, stride 1): mutated in place by the
      // kernel into the solution, then handed straight to the result vector.
      val x = b.toDoubleArrayOwnedCopy
      val info = DoubleKernels.dtrsv(
        n,
        lower,
        unit = false,
        0.0,
        A.data,
        A.offset.value,
        A.rowStride.value,
        A.colStride.value,
        x,
        0,
        1
      )
      if info >= 0 then Left(LinAlgError.SingularMatrix(info))
      else Right(DVec.fromDoubleArrayOwned(x))
