package gale.spectral

import gale.linalg.DVec
import gale.linalg.LinAlgError

/** What a spectral solver has established about its requested result.
  *
  * `ResidualConverged` means every returned pair/triplet passed its residual
  * test inside the explored subspace. It does not prove that an iterative
  * solver reached the requested spectral extreme. `ExtremeCertified` adds an
  * independent membership certificate, currently a full-space reduction.
  */
enum SpectralConvergenceStatus:
  case NotConverged, ResidualConverged, ExtremeCertified

/** Convergence and degeneracy report attached to every spectral result.
  *
  * This is the contract that lets non-convergence be a `Right` rather than a
  * `Left` (§ Convergence & failure semantics of `docs/spectral-parity.md`),
  * mirroring [[gale.linalg.FactorizationDiagnostics]] and
  * [[gale.solvers.SolverResult]]: a structural/precondition violation is a
  * `Left(LinAlgError)`, but partial or zero convergence still returns a value
  * plus these diagnostics.
  *
  *   - `requested` — how many pairs/triplets the selection asked for.
  *   - `converged` — how many met the tolerance; the result contains exactly these.
  *   - `residuals` — the per-pair residual norms of the returned pairs.
  *   - `orthogonalityError` — `‖VᵀV − I‖` of the returned basis, a numerical
  *     quality signal.
  *   - `iterations` — iterations the solver took (`0` for a dense one-shot solve).
  *   - `rank` — numerical rank where meaningful (SVD/GSVD), else `None`.
  *   - `extremalityCertified` — whether the requested spectral membership was
  *     established independently of the returned residuals.
  */
final case class SpectralDiagnostics(
    requested: Int,
    converged: Int,
    residuals: DVec,
    orthogonalityError: Double,
    iterations: Int,
    rank: Option[Int] = None,
    extremalityCertified: Boolean = false
):
  /** True when every requested pair passed the solver's residual test
    * (`converged == requested`). For an iterative partial solver this is
    * convergence within the explored subspace, not proof that the requested
    * global spectral extreme was reached; inspect [[convergenceStatus]] when
    * that distinction matters.
    */
  def allConverged: Boolean =
    converged == requested

  /** Distinguishes incomplete residual convergence, residual convergence in an
    * explored subspace, and independently certified spectral membership.
    */
  def convergenceStatus: SpectralConvergenceStatus =
    if !allConverged then SpectralConvergenceStatus.NotConverged
    else if extremalityCertified then SpectralConvergenceStatus.ExtremeCertified
    else SpectralConvergenceStatus.ResidualConverged

  /** The largest returned residual (`0.0` when there are none), used as the
    * residual payload of a [[gale.linalg.LinAlgError.DidNotConverge]].
    */
  def worstResidual: Double =
    var worst = 0.0
    var i = 0
    while i < residuals.length do
      val r = residuals(i)
      if r > worst then worst = r
      i += 1
    worst

  /** Lift a residual-converged `result` to `Right`, or report residual
    * non-convergence as
    * `Left(`[[gale.linalg.LinAlgError.DidNotConverge]]`)`. This preserves the
    * historical residual-based contract; callers that need proof of requested
    * spectral membership must additionally require [[convergenceStatus]] to be
    * [[SpectralConvergenceStatus.ExtremeCertified]]. Each diagnostics-carrying
    * result exposes this as its own `requireConverged`; composing with the
    * existing `.orThrow` extension (`result.requireConverged.orThrow`) gives the
    * fail-fast residual form.
    */
  def requireConverged[A](result: A): Either[LinAlgError, A] =
    if allConverged then Right(result)
    else Left(LinAlgError.DidNotConverge(iterations, worstResidual))
