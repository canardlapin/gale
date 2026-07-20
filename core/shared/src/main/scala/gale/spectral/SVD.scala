package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError

/** The result of a (partial) singular value decomposition `A = U Σ Vᵀ`.
  *
  * Singular values are returned '''descending always''' (§ 8) — consistent with
  * dense SVD everywhere and with MATLAB, resolving the SciPy-`svds`-ascending
  * trap. `u` and `vt` hold the corresponding left and right singular vectors
  * (`vt` already transposed, i.e. `Vᵀ`); both are empty when only values were
  * requested. `rank` counts the '''returned''' singular values above the solve
  * tolerance relative to the largest '''returned''' one — for a
  * [[SingularOrder.Largest]] selection that is the numerical rank estimate at
  * that tolerance, but under [[SingularOrder.Smallest]] the reference value is
  * not `σ_max` of the matrix, so `rank` is a property of the returned set, not
  * of the whole matrix.
  */
final case class SVD private[spectral] (
    singularValues: DVec,
    u: DMat,
    vt: DMat,
    rank: Int,
    diagnostics: SpectralDiagnostics
):
  /** The number of returned singular values. */
  def size: Int =
    singularValues.length

  /** `Right(this)` when all requested triplets converged, else
    * `Left(`[[gale.linalg.LinAlgError.DidNotConverge]]`)`.
    */
  def requireConverged: Either[LinAlgError, SVD] =
    diagnostics.requireConverged(this)

  /** `Right(this)` only when the requested global singular-spectrum extreme is
    * independently certified. See
    * [[SpectralDiagnostics.requireExtremeCertified]].
    */
  def requireExtremeCertified: Either[LinAlgError, SVD] =
    diagnostics.requireExtremeCertified(this)
