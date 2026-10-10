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
  *   - `orthogonalityError` — the Gram error of the returned basis in the
  *     problem's natural inner product: `‖VᵀV − I‖` for ordinary symmetric
  *     problems, or `‖Vᵀ B V − I‖` for generalized symmetric-definite ones.
  *   - `iterations` — iterations the solver took (`0` for a dense one-shot solve).
  *   - `rank` — numerical rank where meaningful (SVD/GSVD), else `None`.
  *   - `extremalityCertified` — whether the requested spectral membership was
  *     established independently of the returned residuals.
  *   - `innerSolve` — aggregated work for algorithms that explicitly solve an
  *     inner linear system; `None` means no inner solves were performed.
  *
  * '''Deferred measurements.''' The dense one-shot facades (dense symmetric
  * eigen, ordinary and workspace routes, and dense SVD) decide nothing from
  * `residuals` or `orthogonalityError`, so they measure both '''on first
  * access''' and cache them: a caller that never reads them does not pay for
  * the `A·V` and `VᵀV` products. The measurement uses private snapshots of the
  * decomposed matrix and of the vectors, taken when the result was built, so
  * mutating the caller's input or the returned storage afterwards does not
  * change it. Every other solver supplies values it already computed.
  * Equality, `hashCode` and `toString` read both values (taking a deferred
  * measurement); [[copy]] keeps an unreplaced deferred value deferred.
  */
final class SpectralDiagnostics private (
    val requested: Int,
    val converged: Int,
    residualsSource: () => DVec,
    orthogonalitySource: () => Double,
    val iterations: Int,
    val rank: Option[Int],
    val extremalityCertified: Boolean,
    val innerSolve: Option[LinearSolveSummary]
):
  // Each source is dropped once measured, so a cached value does not keep the
  // snapshot matrices alive.
  private var pendingResiduals: () => DVec = residualsSource
  private var pendingOrthogonality: () => Double = orthogonalitySource

  /** Per-pair residual norms of the returned pairs. */
  lazy val residuals: DVec =
    val value = pendingResiduals()
    pendingResiduals = null
    value

  /** Gram error of the returned basis. */
  lazy val orthogonalityError: Double =
    val value = pendingOrthogonality()
    pendingOrthogonality = null
    value

  /** A copy with the given fields replaced. An unreplaced deferred measurement
    * stays deferred and is shared with this instance (measured once).
    */
  def copy(
      requested: Int = this.requested,
      converged: Int = this.converged,
      residuals: => DVec = this.residuals,
      orthogonalityError: => Double = this.orthogonalityError,
      iterations: Int = this.iterations,
      rank: Option[Int] = this.rank,
      extremalityCertified: Boolean = this.extremalityCertified,
      innerSolve: Option[LinearSolveSummary] = this.innerSolve
  ): SpectralDiagnostics =
    new SpectralDiagnostics(
      requested,
      converged,
      () => residuals,
      () => orthogonalityError,
      iterations,
      rank,
      extremalityCertified,
      innerSolve
    )

  override def equals(other: Any): Boolean =
    other match
      case that: SpectralDiagnostics =>
        requested == that.requested && converged == that.converged &&
        residuals == that.residuals && orthogonalityError == that.orthogonalityError &&
        iterations == that.iterations && rank == that.rank &&
        extremalityCertified == that.extremalityCertified && innerSolve == that.innerSolve
      case _ => false

  override def hashCode: Int =
    (requested, converged, residuals, orthogonalityError, iterations, rank, extremalityCertified, innerSolve).hashCode

  override def toString: String =
    s"SpectralDiagnostics($requested,$converged,$residuals,$orthogonalityError,$iterations,$rank," +
      s"$extremalityCertified,$innerSolve)"

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

  /** Require both residual convergence and an independent certificate that the
    * result belongs to the requested global spectral extreme.
    *
    * A residual failure returns [[gale.linalg.LinAlgError.DidNotConverge]]. A
    * residual-converged result without an extremality certificate returns the
    * distinct typed error
    * [[gale.linalg.LinAlgError.SpectralExtremeNotCertified]]. This method only
    * enforces an existing certificate; it never manufactures one.
    */
  def requireExtremeCertified[A](result: A): Either[LinAlgError, A] =
    convergenceStatus match
      case SpectralConvergenceStatus.NotConverged =>
        Left(LinAlgError.DidNotConverge(iterations, worstResidual))
      case SpectralConvergenceStatus.ResidualConverged =>
        Left(LinAlgError.SpectralExtremeNotCertified(iterations, worstResidual))
      case SpectralConvergenceStatus.ExtremeCertified =>
        Right(result)

object SpectralDiagnostics:
  /** Diagnostics with already-measured `residuals` and `orthogonalityError`. */
  def apply(
      requested: Int,
      converged: Int,
      residuals: DVec,
      orthogonalityError: Double,
      iterations: Int,
      rank: Option[Int] = None,
      extremalityCertified: Boolean = false,
      innerSolve: Option[LinearSolveSummary] = None
  ): SpectralDiagnostics =
    new SpectralDiagnostics(
      requested,
      converged,
      () => residuals,
      () => orthogonalityError,
      iterations,
      rank,
      extremalityCertified,
      innerSolve
    )

  /** Diagnostics whose `residuals` and `orthogonalityError` are measured on
    * first access and cached. Both sources must read only snapshots that no
    * caller can mutate.
    */
  private[spectral] def deferred(
      requested: Int,
      converged: Int,
      residuals: () => DVec,
      orthogonalityError: () => Double,
      iterations: Int,
      rank: Option[Int],
      extremalityCertified: Boolean
  ): SpectralDiagnostics =
    new SpectralDiagnostics(
      requested,
      converged,
      residuals,
      orthogonalityError,
      iterations,
      rank,
      extremalityCertified,
      None
    )
