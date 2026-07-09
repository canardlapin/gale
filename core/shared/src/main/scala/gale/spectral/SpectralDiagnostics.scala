package gale.spectral

import gale.linalg.DVec
import gale.linalg.LinAlgError

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
  */
final case class SpectralDiagnostics(
    requested: Int,
    converged: Int,
    residuals: DVec,
    orthogonalityError: Double,
    iterations: Int,
    rank: Option[Int] = None
):
  /** True when every requested pair converged (`converged == requested`); the
    * spectral analogue of [[gale.linalg.FactorizationDiagnostics.isSuccess]].
    */
  def allConverged: Boolean =
    converged == requested

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

  /** Lift a converged `result` to `Right`, or report non-convergence as
    * `Left(`[[gale.linalg.LinAlgError.DidNotConverge]]`)`. Each diagnostics-
    * carrying result exposes this as its own `requireConverged`; composing with
    * the existing `.orThrow` extension (`result.requireConverged.orThrow`) gives
    * the fail-fast form.
    */
  def requireConverged[A](result: A): Either[LinAlgError, A] =
    if allConverged then Right(result)
    else Left(LinAlgError.DidNotConverge(iterations, worstResidual))
