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

  /** Breeze `A \ b` for a square system: solves `A x = b` via gale's LU. */
  def solve(a: DenseMatrix[Double], b: DenseVector[Double]): DenseVector[Double] =
    toBreezeCopy(fromBreezeCopy(a).solve(fromBreezeCopy(b)).orThrow)

  /** Breeze `det(A)`: determinant via gale's LU (sign included). */
  def det(a: DenseMatrix[Double]): Double =
    fromBreezeCopy(a).det.orThrow

  /** Breeze `cholesky(A)`: the lower Cholesky factor `L` with `A = L Lᵀ`. Throws if
    * `A` is not symmetric positive-definite.
    */
  def cholesky(a: DenseMatrix[Double]): DenseMatrix[Double] =
    toBreezeCopy(fromBreezeCopy(a).cholesky.orThrow.lower)

  /** Breeze `eigSym(A)`: `(eigenvalues ascending, eigenvectors as columns)` for a
    * symmetric `A` (only its lower triangle is read, as in Breeze/LAPACK).
    */
  def eigSym(a: DenseMatrix[Double]): (DenseVector[Double], DenseMatrix[Double]) =
    val d = Eigen.eigSymmetric(fromBreezeCopy(a), EigenSelection.All, EigenVectors.Right).orThrow
    (toBreezeCopy(d.eigenvalues), toBreezeCopy(d.eigenvectors))

  /** Breeze `A \ b` for an overdetermined (tall) system: least-squares via gale's
    * QR. Throws on rank deficiency.
    */
  def leastSquares(a: DenseMatrix[Double], b: DenseVector[Double]): DenseVector[Double] =
    toBreezeCopy(fromBreezeCopy(a).leastSquares(fromBreezeCopy(b)).orThrow)
