package gale.spectral

import gale.linalg.DVec
import gale.solvers.SolverConfig

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
  * `(A − σB) x = b` in the generalized case — is carried out.
  *
  * gale never auto-factorizes (unlike MATLAB `eigs`, which factorizes `A − σB`
  * behind the scenes): a [[SpectralTarget.ShiftInvert]] must name its plan
  * explicitly. This is the minimal placeholder the types phase needs; the kernel
  * phase may extend it. It reuses the existing [[gale.solvers.SolverConfig]]
  * rather than introducing a parallel config type.
  */
enum LinearSolvePlan:
  /** Factor the shifted matrix once with a dense LU and reuse it across
    * iterations.
    */
  case Direct

  /** Solve the shifted system iteratively with the given configuration. */
  case Iterative(config: SolverConfig = SolverConfig())

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
