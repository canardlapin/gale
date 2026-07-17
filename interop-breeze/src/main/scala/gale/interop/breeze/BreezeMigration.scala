package gale.interop.breeze

import _root_.breeze.linalg.DenseMatrix
import _root_.breeze.linalg.DenseVector
import gale.linalg.*
import gale.spectral.Eigen
import gale.spectral.EigenSelection
import gale.spectral.EigenVectors

/** Migration '''sugar''' for porting Breeze call sites to gale: a few common Breeze
  * idioms, Breeze types in and out, computed by gale underneath (via
  * [[toBreezeCopy]] / [[fromBreezeCopy]]).
  *
  * This is deliberately '''minimal''' — a hop to ease a mechanical port, not a
  * Breeze compatibility layer. Each call copies across the boundary, so once a call
  * site is ported, prefer gale's native API (which returns `Either[LinAlgError, _]`
  * rather than throwing) and gale-native storage. Anything Breeze offers that gale
  * has no public equivalent for (e.g. a dense pseudo-inverse / `pinv`) is
  * intentionally absent; use [[leastSquares]] for overdetermined systems.
  *
  * Every method '''throws''' `gale.linalg.LinAlgError` on a structural failure
  * (singular matrix, non-SPD, …), mirroring Breeze's exception-throwing style.
  */
object BreezeMigration:

  private inline val BreezeSymmetryTolerance = 1.0e-7

  /** Breeze `A \ b` for a square system: solves `A x = b` via gale's LU. */
  def solve(a: DenseMatrix[Double], b: DenseVector[Double]): DenseVector[Double] =
    toBreezeCopy(fromBreezeCopy(a).solve(fromBreezeCopy(b)).orThrow)

  /** Breeze `A \ B` for several right-hand sides. Gale factors `A` once and
    * reuses the typed LU factor for every column of `B`.
    */
  def solve(a: DenseMatrix[Double], b: DenseMatrix[Double]): DenseMatrix[Double] =
    val ga = fromBreezeCopy(a)
    val gb = fromBreezeCopy(b)
    requireRhsRows(ga.rows, gb)
    val lu = ga.lu.orThrow
    val columns = Vector.tabulate(gb.cols)(j => lu.solve(gb.col(j)).orThrow)
    toBreezeCopy(Matrix.tabulate(ga.cols, gb.cols)((i, j) => columns(j)(i)))

  /** Breeze `det(A)`: determinant via gale's LU (sign included). */
  def det(a: DenseMatrix[Double]): Double =
    fromBreezeCopy(a).det.orThrow

  /** Breeze `cholesky(A)`: the lower Cholesky factor `L` with `A = L Lᵀ`. The shim
    * applies Breeze's non-empty, square, and `1e-7` scale-aware symmetry checks
    * before gale's lower-triangle kernel runs — so, exactly like Breeze (and unlike
    * gale's own `DMat.cholesky`), a matrix whose upper triangle disagrees with its
    * lower is '''rejected''', not silently read lower-triangle-only. Throws if `A`
    * fails those checks or is not positive-definite.
    */
  def cholesky(a: DenseMatrix[Double]): DenseMatrix[Double] =
    requireBreezeSymmetric(a)
    toBreezeCopy(fromBreezeCopy(a).cholesky.orThrow.lower)

  /** Breeze `eigSym(A)`: `(eigenvalues ascending, eigenvectors as columns)` for a
    * symmetric `A`. The shim applies Breeze's non-empty, square, and `1e-7`
    * scale-aware symmetry checks before Gale's lower-triangle eigensolver runs.
    */
  def eigSym(a: DenseMatrix[Double]): (DenseVector[Double], DenseMatrix[Double]) =
    requireBreezeSymmetric(a)
    val d = Eigen.eigSymmetric(fromBreezeCopy(a), EigenSelection.All, EigenVectors.Right).orThrow
    (toBreezeCopy(d.eigenvalues), toBreezeCopy(d.eigenvectors))

  /** Breeze `A \ b` for an overdetermined (tall) system: least-squares via gale's
    * QR. Throws on rank deficiency.
    */
  def leastSquares(a: DenseMatrix[Double], b: DenseVector[Double]): DenseVector[Double] =
    toBreezeCopy(fromBreezeCopy(a).leastSquares(fromBreezeCopy(b)).orThrow)

  /** Breeze `A \ B` for an overdetermined system with several right-hand sides.
    * Gale computes one QR factorization and reuses it for every column of `B`.
    */
  def leastSquares(a: DenseMatrix[Double], b: DenseMatrix[Double]): DenseMatrix[Double] =
    val ga = fromBreezeCopy(a)
    val gb = fromBreezeCopy(b)
    requireRhsRows(ga.rows, gb)
    val qr = ga.qr
    val columns = Vector.tabulate(gb.cols)(j => qr.solveLeastSquares(gb.col(j)).orThrow)
    toBreezeCopy(Matrix.tabulate(ga.cols, gb.cols)((i, j) => columns(j)(i)))

  private def requireRhsRows(expected: Int, b: DMat): Unit =
    if b.rows != expected then
      throw LinAlgError.DimensionMismatch(
        Shape(Rows(expected), Cols(b.cols)),
        Shape(Rows(b.rows), Cols(b.cols))
      )

  /** Match Breeze 2.1's public `cholesky` / `eigSym` input guard. Breeze checks
    * logical entries before copying the lower triangle into LAPACK; doing the same
    * here keeps this migration shim from silently accepting an asymmetric matrix
    * that the call site being migrated would have rejected.
    */
  private def requireBreezeSymmetric(a: DenseMatrix[Double]): Unit =
    if a.rows == 0 || a.cols == 0 then
      throw LinAlgError.InvalidArgument("matrix must be non-empty")
    if a.rows != a.cols then
      throw LinAlgError.NonSquareMatrix(Shape(Rows(a.rows), Cols(a.cols)))
    var i = 0
    while i < a.rows do
      var j = 0
      while j < i do
        val lower = a(i, j)
        val upper = a(j, i)
        val scale = math.max(1.0, math.max(math.abs(lower), math.abs(upper)))
        val close = lower == upper || math.abs(lower - upper) < scale * BreezeSymmetryTolerance
        if !close then
          throw LinAlgError.InvalidArgument(
            s"matrix is not symmetric at ($i, $j) / ($j, $i): $lower != $upper within Breeze tolerance $BreezeSymmetryTolerance"
          )
        j += 1
      i += 1
