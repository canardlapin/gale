package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.ExactSolveFactor
import gale.linalg.LinAlgError
import gale.linalg.PositiveDefinite
import gale.solvers.Preconditioner
import gale.solvers.SolverConfig
import gale.solvers.ToleranceMode

/** Which extreme of the spectrum a count-based eigen selection picks.
  *
  * Order legality depends on the problem type and is enforced at the solver
  * boundary (`Left(InvalidArgument)`), not here (§ Selection & ordering of
  * `docs/spectral-parity.md`): the algebraic orders and `BothEnds` are
  * symmetric-only, the real-part orders are nonsymmetric-only, and the magnitude
  * orders are legal for both. Imaginary-part orders are deferred out of v0.3.5.
  * The order chooses *membership*; the output layout is fixed separately
  * (symmetric ascending-algebraic, nonsymmetric by criterion).
  */
enum EigenOrder:
  case LargestMagnitude, SmallestMagnitude
  case LargestAlgebraic, SmallestAlgebraic
  case LargestRealPart, SmallestRealPart

  /** `eigsh`'s `'BE'`: ⌈k/2⌉ from the high end and ⌊k/2⌋ from the low end.
    * Symmetric-only.
    */
  case BothEnds

/** Which extreme of the singular spectrum a count-based SVD selection picks.
  * Singular values are real and nonnegative, so magnitude and algebraic order
  * coincide.
  */
enum SingularOrder:
  case Largest, Smallest

/** How many eigenpairs to compute and which ones.
  *
  * Supersedes the PRD's count-only `SpectralSelection`, which could not express
  * `eigh`'s subset-by-index / subset-by-value (§ 1, § 3 of the parity doc).
  * `IndexRange` and `ValueInterval` are symmetric-only. Every case carries
  * selection data only; all legality checks — `k` versus the dimension,
  * problem-type restrictions — belong to the solver and are returned as
  * `Left(InvalidArgument)`, never thrown here.
  */
enum EigenSelection:
  /** The full spectrum (dense decomposition). */
  case All

  /** The `k` eigenvalues at the extreme named by `order`. */
  case Count(k: Int, order: EigenOrder)

  /** The eigenvalues whose ascending-algebraic ranks lie in `[from, to]`
    * (0-based, inclusive). Symmetric-only.
    */
  case IndexRange(from: Int, to: Int)

  /** Every eigenvalue in the half-open interval `(lower, upper]`. The count of
    * results is data-dependent. Symmetric-only.
    */
  case ValueInterval(lower: Double, upper: Double)

/** How many singular triplets to compute and which ones. There is no index or
  * value subset — no ecosystem precedent worth cloning (§ 3, § 8).
  */
enum SingularSelection:
  case All
  case Count(k: Int, order: SingularOrder)

/** How the inner linear solve that shift-invert needs — `(A − σI) x = b`, or
  * `(A − σB) x = b` in the generalized case — is obtained.
  *
  * gale never auto-factorizes: a [[SpectralTarget.ShiftInvert]] either carries a
  * caller-prepared executable solve or explicitly requests the optional backend
  * capability. A plan is not a matrix operator and cannot be applied without
  * observing its solve diagnostics.
  */
enum LinearSolvePlan:
  /** Use this already-prepared, reusable solve capability. */
  case Use(solver: LinearSolveOperator)

  /** Ask the resolved backend for its advertised shift-invert solve. */
  case Backend

object LinearSolvePlan:
  /** Adapt a caller-created factor; this factory does not factor a matrix. */
  def direct(factor: ExactSolveFactor): LinearSolvePlan =
    LinearSolvePlan.Use(LinearSolveOperator.direct(factor))

  /** Build an explicit iterative solve for a caller-asserted SPD shifted system. */
  def iterative[A <: DoubleLinearOperator](
      operator: PositiveDefinite[A],
      config: SolverConfig = SolverConfig(),
      preconditioner: Preconditioner = Preconditioner.Identity,
      toleranceMode: ToleranceMode = ToleranceMode.RelativeToRhs
  ): Either[LinAlgError, LinearSolvePlan] =
    LinearSolveOperator
      .conjugateGradient(operator, config, preconditioner, toleranceMode)
      .map(LinearSolvePlan.Use.apply)

  /** Resolve a plan without implicit factorization. */
  def resolve(
      plan: LinearSolvePlan,
      a: DMat,
      b: Option[DMat],
      sigma: Double
  )(using backend: SpectralBackend): Either[LinAlgError, LinearSolveOperator] =
    def validateSize(solver: LinearSolveOperator): Either[LinAlgError, LinearSolveOperator] =
      if solver.size != a.rows then
        Left(LinAlgError.VectorLengthMismatch(a.rows, solver.size))
      else Right(solver)

    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else if b.exists(matrix => matrix.rows != a.rows || matrix.cols != a.cols) then
      Left(LinAlgError.InvalidArgument("shift-invert metric B must have the same square shape as A"))
    else if !sigma.isFinite then
      Left(LinAlgError.InvalidArgument(s"shift-invert sigma must be finite, got $sigma"))
    else
      plan match
        case LinearSolvePlan.Use(solver) =>
          validateSize(solver)
        case LinearSolvePlan.Backend =>
          if backend.capabilities.contains(SpectralCapability.ShiftInvertSolve) then
            backend.shiftInvertSolve(a, b, sigma).flatMap(validateSize)
          else
            Left(
              LinAlgError.UnsupportedOperation(
                s"shift-invert solve: ${backend.name} backend does not provide it"
              )
            )

/** A spectral transformation targeting eigenvalues near a real point `σ`.
  * Complex shifts (targeting off the real axis) are out of v0.3.5 — they would
  * need complex linear solves (§ Explicitly OUT).
  */
enum SpectralTarget:
  /** Prefer eigenvalues near `value` (a selection hint, no shift-invert). */
  case Around(value: Double)

  /** Shift-invert around a '''real''' `sigma`, using `plan` for the inner solve. */
  case ShiftInvert(sigma: Double, plan: LinearSolvePlan)

/** Which eigenvectors a nonsymmetric solve returns. The symmetric API restricts
  * itself to `ValuesOnly` versus vectors (left and right coincide for a
  * symmetric matrix), so it never exposes `Left`/`LeftAndRight`.
  */
enum EigenVectors:
  case ValuesOnly, Right, Left, LeftAndRight

/** Tuning knobs for the partial (iterative) spectral solvers, mirroring
  * `eigs`/`eigsh`/`svds` (§ 6–8).
  *
  *   - `tolerance` / `maxIterations` share [[gale.solvers.SolverConfig]]'s
  *     defaults.
  *   - `subspaceDimension` is the initial block-Krylov/Arnoldi basis size
  *     (`ncv`/`p`). `None` means the solver computes the default
  *     `min(n, max(2k+1, 20))` at the call site, where `n` and `k` are known.
  *     The symmetric solver grows this projection across thick restarts when
  *     needed; the option is not a hard memory cap.
  *   - `startVector` seeds the first block column (`v0`); `None` uses an internal
  *     default. The multiplicity-safe symmetric solver fills the other block
  *     columns with deterministic orthogonal probes.
  *   - `returnVectors` selects which vectors to compute.
  */
final case class SpectralOptions(
    tolerance: Double = 1e-10,
    maxIterations: Int = 1000,
    subspaceDimension: Option[Int] = None,
    startVector: Option[DVec] = None,
    returnVectors: EigenVectors = EigenVectors.Right
)

/** Tuning knobs for the partial generalized symmetric-definite operator solver
  * `A x = λ B x`, where `A` is symmetric and `B` is symmetric
  * positive-definite.
  *
  * This is deliberately distinct from [[SpectralOptions]]. LOBPCG evolves a
  * block of `k` vectors together, so reducing its initialization to one
  * `startVector` would hide the rank and multiplicity contract.
  *
  *   - `tolerance` is applied to the true generalized residual
  *     `‖A x - λ B x‖`.
  *   - `maxIterations` bounds outer block iterations; inner work performed by a
  *     preconditioner or later metric-solve provider reports separately.
  *   - `initialSubspace`, when present, must have shape `n × k` for a
  *     `Count(k, ...)` selection. The solver owns a snapshot and may
  *     B-orthonormalize, rank-reveal, and deterministically replenish it.
  *   - `returnVectors` is restricted to [[EigenVectors.ValuesOnly]] or
  *     [[EigenVectors.Right]] at the solver boundary.
  *
  * Constructors remain data-only. Shape, finiteness, rank, tolerance, and
  * vector-flag validation return typed `Left` values from the solver facade.
  */
final case class GeneralizedSpectralOptions(
    tolerance: Double = 1e-10,
    maxIterations: Int = 200,
    initialSubspace: Option[DMat] = None,
    returnVectors: EigenVectors = EigenVectors.Right
)
