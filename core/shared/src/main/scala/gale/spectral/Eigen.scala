package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DenseDecompositions
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import gale.linalg.Rows
import gale.linalg.Shape
import gale.linalg.TriangularSolve

/** Public eigendecomposition entry points — symmetric (phase a) and nonsymmetric
  * (phase b) of `docs/spectral-parity.md`.
  *
  * '''Symmetric''' ([[EigenDecomposition]], § 1 / § 6). '''Eigenvalues are returned
  * ascending-algebraic always'''; a selection's `order` decides only *which*
  * eigenvalues are members, never their layout.
  *
  *   - [[eigSymmetric(a:gale\.linalg\.DMat*]] — dense full solve via the
  *     tridiagonal QL/QR kernel, the requested selection realized as a slice of the
  *     full ascending spectrum (the phase-a scope decision in § 1). Reads only the
  *     lower triangle of `a`.
  *   - [[eigSymmetric(op:gale\.linalg\.DoubleLinearOperator*]] — iterative partial
  *     solve (`eigsh`) via Lanczos with full reorthogonalization, for a matrix or
  *     a matrix-free [[gale.linalg.DoubleLinearOperator]].
  *
  * '''Nonsymmetric''' ([[NonsymmetricEigenDecomposition]], § 2 / § 7). A real input
  * can have complex eigenvalues in conjugate pairs; the output is ordered '''by the
  * selection criterion''' with conjugate pairs kept adjacent (positive-imaginary
  * first) and '''never split''' during selection.
  *
  *   - [[eigNonsymmetric(a:gale\.linalg\.DMat*]] — dense full solve via Hessenberg
  *     reduction + Francis double-shift QR, the selection realized as a permutation
  *     of the full spectrum. Right eigenvectors from the kernel; '''left'''
  *     eigenvectors (`wᴴ A = λ wᴴ`) via the biorthogonal `V⁻¹` route.
  *   - [[eigNonsymmetric(op:gale\.linalg\.DoubleLinearOperator*]] — iterative
  *     partial solve (`eigs`) via Arnoldi with full reorthogonalization; right
  *     eigenvectors only (left vectors need the dense path).
  *
  * '''Failure model''' (§ Convergence & failure semantics). Structural /
  * precondition violations are `Left(LinAlgError)`. For the '''dense''' paths,
  * kernel sweep exhaustion is also a `Left(DidNotConverge)` — there is no partial
  * dense result to hand back (in practice the QL/QR and Francis QR kernels always
  * converge within budget for a well-formed matrix). For the '''iterative''' paths,
  * non-convergence is never a `Left`: partial or zero convergence returns
  * `Right(result-with-only-the-converged-pairs + SpectralDiagnostics)`, and the
  * caller opts into fail-fast via `result.requireConverged`.
  */
object Eigen:

  // ===========================================================================
  // Dense symmetric eigendecomposition
  // ===========================================================================

  /** Dense symmetric eigendecomposition of `a` computing eigenvectors
    * ([[EigenVectors.Right]]). See the three-argument overload for the vector-flag
    * and failure details.
    */
  def eigSymmetric(a: DMat, selection: EigenSelection): Either[LinAlgError, EigenDecomposition] =
    eigSymmetric(a, selection, EigenVectors.Right)

  /** Dense symmetric eigendecomposition of `a`. Only the '''lower triangle''' of
    * `a` is read (the `Cholesky` precedent); the strict upper triangle is treated
    * as its mirror. `vectors` selects [[EigenVectors.ValuesOnly]] versus
    * [[EigenVectors.Right]] (left and right eigenvectors coincide for a symmetric
    * matrix, so [[EigenVectors.Left]]/[[EigenVectors.LeftAndRight]] are rejected
    * with `Left(InvalidArgument)`).
    *
    * The full ascending spectrum is computed once and the `selection` realized as
    * a slice of it: [[EigenSelection.Count]] takes the `k` extremes named by its
    * order, [[EigenSelection.IndexRange]] the ascending-rank window, and
    * [[EigenSelection.ValueInterval]] every eigenvalue in `(lower, upper]`. Output
    * is ascending regardless.
    *
    * `Left` on: non-square `a`; an [[EigenOrder]] illegal for a symmetric problem
    * ([[EigenOrder.LargestRealPart]]/[[EigenOrder.SmallestRealPart]]); `k` outside
    * `[1, n]`; an out-of-bounds `IndexRange`; an inverted `ValueInterval`; or (in
    * practice unreachable) kernel non-convergence.
    */
  def eigSymmetric(
      a: DMat,
      selection: EigenSelection,
      vectors: EigenVectors
  ): Either[LinAlgError, EigenDecomposition] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else
      val n = a.rows
      validateVectors(vectors) match
        case Left(error) => Left(error)
        case Right(wantVectors) =>
          validateDenseSelection(selection, n) match
            case Left(error) => Left(error)
            case Right(()) =>
              DenseSpectralKernels.symmetricEigen(a, wantVectors) match
                case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
                  // Dense path has no partial result to return, so exhaustion is a
                  // Left; residual is not measured on failure (reported as 0.0).
                  Left(LinAlgError.DidNotConverge(iters, 0.0))
                case Right(kernel) =>
                  val indices = denseSelectionIndices(selection, kernel.values, n)
                  Right(assembleDense(a, kernel.values, kernel.vectors, indices, wantVectors))

  // ===========================================================================
  // Iterative partial symmetric eigendecomposition (eigsh)
  // ===========================================================================

  /** Partial symmetric eigendecomposition of the operator `op` (an `n x n`
    * symmetric matrix or matrix-free operator), computing the `k` eigenpairs at
    * the extreme named by the selection's order via Lanczos with '''full
    * reorthogonalization''' over a '''growing''' Krylov subspace.
    *
    * Full reorthogonalization makes the method robust at `O(n·m²)` cost for a
    * subspace of size `m` (portable-correctness first; implicit restarting is a
    * later optimization). Starting from `m = subspaceDimension.getOrElse(min(n,
    * max(2k+1, 20)))` (clamped to `[k+1, n]`), each step builds an `m`-step
    * Lanczos factorization from a '''fixed''' start vector, solves the projected
    * tridiagonal with the shared QL/QR kernel, forms Ritz pairs `θ, y = V s`, and
    * measures the true residual `‖op·y − θ y‖`. A pair is converged when that
    * residual is `≤ tolerance · max(1, max|θ|)`. If the wanted pairs have not all
    * converged, the subspace '''grows''' (same start vector — a lossless
    * refinement whose extreme Ritz values improve by Cauchy interlacing and which
    * is exact at `m = n`; per-pair residuals are not guaranteed monotone under a
    * `maxIterations` truncation), up to `maxIterations` growth steps or a full
    * `m = n` subspace (whose Ritz values are the exact spectrum).
    * Well-separated extremes converge cheaply; tightly clustered or interior
    * targets (the classic slow case for plain Lanczos, e.g.
    * [[EigenOrder.SmallestMagnitude]] across zero) may only converge as `m`
    * approaches `n`, or return partial — shift-invert is the deferred remedy.
    *
    * '''Multiplicity caveat.''' Single-vector Lanczos resolves at most one
    * eigenvector per '''distinct''' eigenvalue: a repeated eigenvalue is returned
    * once, its multiplicity silently collapsed, and the result can still report
    * `allConverged = true` (each returned pair genuinely has a tiny residual).
    * With a degenerate extreme, the returned set is the top-`k` of the
    * '''distinct''' spectrum, not the top-`k` with multiplicity — plain `eigsh`
    * shares this limitation. Use the dense API when multiplicities matter.
    *
    * Non-convergence is a `Right`: the result holds only the converged pairs
    * (ascending), with `SpectralDiagnostics.converged < requested` and per-pair
    * residuals recorded. Eigenvalues that could not be resolved are simply absent;
    * `requireConverged` is the caller's opt-in to a `Left(DidNotConverge)`.
    *
    * `Left` on: a non-square operator (`NonSquareMatrix`) or a square one whose
    * shape disagrees with `n` (`DimensionMismatch`); a
    * non-positive `n`; `k ≤ 0` or `k ≥ n` (the message points at the dense API —
    * there is no silent dense fallback); an [[EigenOrder]] illegal for a symmetric
    * problem; a start vector of the wrong length or zero norm; or a `target`
    * (shift-invert / `Around`), which is '''deferred''' —
    * `Left(UnsupportedOperation)` until the `LinearSolvePlan` wiring lands.
    */
  def eigSymmetric(
      op: DoubleLinearOperator,
      n: Int,
      selection: EigenSelection,
      options: SpectralOptions = SpectralOptions(),
      target: Option[SpectralTarget] = None
  ): Either[LinAlgError, EigenDecomposition] =
    target match
      case Some(t) =>
        Left(LinAlgError.UnsupportedOperation(s"shift-invert / targeted selection ($t)"))
      case None =>
        selection match
          case count: EigenSelection.Count => eigSymmetricLanczos(op, n, count.k, count.order, options)
          case other =>
            Left(
              LinAlgError.InvalidArgument(
                s"the iterative solver requires EigenSelection.Count; ${other.getClass.getSimpleName} is dense-only"
              )
            )

  private def eigSymmetricLanczos(
      op: DoubleLinearOperator,
      n: Int,
      k: Int,
      order: EigenOrder,
      options: SpectralOptions
  ): Either[LinAlgError, EigenDecomposition] =
        if n <= 0 then Left(LinAlgError.InvalidArgument(s"dimension must be positive, got $n"))
        else if op.rows != op.cols then
          Left(LinAlgError.NonSquareMatrix(Shape(Rows(op.rows), Cols(op.cols))))
        else if op.rows != n then
          Left(LinAlgError.DimensionMismatch(Shape(Rows(n), Cols(n)), Shape(Rows(op.rows), Cols(op.cols))))
        else
          validateVectors(options.returnVectors) match
            case Left(error) => Left(error)
            case Right(wantVectors) =>
              if !symmetricOrderLegal(order) then
                Left(LinAlgError.InvalidArgument(s"$order is nonsymmetric-only; use an algebraic, magnitude, or both-ends order"))
              else if k <= 0 || k >= n then
                Left(
                  LinAlgError.InvalidArgument(
                    s"k=$k must be in [1, ${n - 1}] for the iterative solver; use the dense eigSymmetric for the full spectrum"
                  )
                )
              else
                startVectorFor(options.startVector, n) match
                  case Left(error) => Left(error)
                  case Right(v0)   => runLanczos(op, n, k, order, options, wantVectors, v0)

  // ===========================================================================
  // Dense nonsymmetric eigendecomposition (eig)
  // ===========================================================================

  /** Dense nonsymmetric eigendecomposition of `a` computing right eigenvectors
    * ([[EigenVectors.Right]]). See the three-argument overload for the vector-flag
    * and failure details.
    */
  def eigNonsymmetric(
      a: DMat,
      selection: EigenSelection
  ): Either[LinAlgError, NonsymmetricEigenDecomposition] =
    eigNonsymmetric(a, selection, EigenVectors.Right)

  /** Dense nonsymmetric eigendecomposition of `a` (phase b, § 2 of
    * `docs/spectral-parity.md`) via Hessenberg reduction + Francis double-shift QR
    * ([[gale.spectral.DenseSpectralKernels.nonsymmetricEigen]]), the requested
    * `selection` realized as a permutation of the full spectrum.
    *
    * A real input can have complex eigenvalues in conjugate pairs; the result
    * stores them structure-of-arrays with the real-Schur packing described on
    * [[NonsymmetricEigenDecomposition]] and read only through its typed accessors.
    *
    * '''Ordering (§ 2 guarantee).''' The output is ordered '''by the selection
    * criterion''' — the largest orders descending, the smallest orders ascending
    * (on magnitude or real part) — with ties broken by descending real part then
    * descending imaginary part, and '''conjugate pairs kept adjacent, positive-
    * imaginary member first''' (mandatory for the packing). [[EigenSelection.All]]
    * uses the [[EigenOrder.LargestMagnitude]] order as its canonical layout.
    *
    * '''Selection operates on pairs as units.''' A conjugate pair is selected and
    * moved as one two-column unit; it is '''never split'''. Consequently a
    * [[EigenSelection.Count]] whose `k`-boundary would fall between the two members
    * of a pair (they share a criterion value, so the boundary is exactly a tie)
    * returns the '''whole pair''' — i.e. `k+1` eigenvalues rather than `k`. This is
    * a dense one-shot solve, so `diagnostics.requested` and `converged` both report
    * the actual count returned and `allConverged` is true.
    *
    * `vectors` chooses [[EigenVectors.ValuesOnly]], [[EigenVectors.Right]],
    * [[EigenVectors.Left]], or [[EigenVectors.LeftAndRight]]. '''Left eigenvectors'''
    * follow the convention `wᴴ A = λ wᴴ`, '''unit 2-norm''' — a real `w` for a real
    * `λ`, and the same real-Schur SoA packing as the right vectors for a conjugate
    * pair — read through [[NonsymmetricEigenDecomposition.leftEigenvector]]. They are
    * recovered from the full right-eigenvector matrix `V` via the biorthogonality
    * identity (the left vectors are the conjugated rows of `V⁻¹`), the complex `V⁻¹`
    * formed with a real `2n×2n` embedded LU solve (no complex kernel tier), so the
    * pairing with the right vectors and eigenvalues is '''exact by construction'''.
    * A '''defective''' (non-diagonalizable) `a` has no full set of left eigenvectors;
    * gale detects that with a conditioning + residual guard on the recovery (see
    * `computeLeftVectors`) and returns `Left(SingularMatrix)` rather than a degenerate
    * basis. '''This diverges from LAPACK `dgeev`''', which returns (possibly
    * near-parallel, non-orthogonal) left vectors for a defective `A`; gale prefers the
    * explicit `Left`. Orthogonality error is reported as `0.0`: nonsymmetric
    * eigenvectors are not orthonormal, so `‖VᵀV − I‖` is not a meaningful signal.
    *
    * `Left` on: non-square `a`; an [[EigenOrder]] illegal for a nonsymmetric
    * problem ([[EigenOrder.LargestAlgebraic]]/[[EigenOrder.SmallestAlgebraic]]/
    * [[EigenOrder.BothEnds]] are symmetric-only); a `k` outside `[1, n]`;
    * [[EigenSelection.IndexRange]] or [[EigenSelection.ValueInterval]]
    * (symmetric-only); a defective `a` when left vectors are requested
    * (`SingularMatrix`); or (in practice unreachable) kernel non-convergence.
    */
  def eigNonsymmetric(
      a: DMat,
      selection: EigenSelection,
      vectors: EigenVectors
  ): Either[LinAlgError, NonsymmetricEigenDecomposition] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else
      val n = a.rows
      nonsymmetricVectorFlags(vectors) match
        case Left(error) => Left(error)
        case Right((wantRight, wantLeft)) =>
          validateNonsymmetricDenseSelection(selection, n) match
            case Left(error) => Left(error)
            case Right(()) =>
              // Left vectors are derived from the full right-eigenvector matrix, so
              // the kernel must produce right vectors whenever either side is wanted.
              DenseSpectralKernels.nonsymmetricEigen(a, wantRight || wantLeft) match
                case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
                  // Dense path has no partial result to return, so exhaustion is a
                  // Left; residual is not measured on failure (reported as 0.0).
                  Left(LinAlgError.DidNotConverge(iters, 0.0))
                case Right(kernel) =>
                  val leftFull: Either[LinAlgError, Option[DMat]] =
                    if wantLeft then computeLeftVectors(a, kernel).map(Some(_)) else Right(None)
                  leftFull match
                    case Left(error) => Left(error)
                    case Right(leftOpt) =>
                      val outIdx = nonsymDenseIndices(kernel.re, kernel.im, selection)
                      Right(assembleNonsymDense(a, kernel, leftOpt, outIdx, wantRight, wantLeft))

  // ===========================================================================
  // Iterative partial nonsymmetric eigendecomposition (eigs)
  // ===========================================================================

  /** Partial nonsymmetric eigendecomposition of the operator `op` (an `n x n`
    * matrix or matrix-free operator, § 7) — the `k` eigenpairs at the extreme named
    * by the selection's order — via '''Arnoldi with full reorthogonalization'''
    * (classical Gram–Schmidt twice) over a '''growing''' Krylov subspace from a
    * fixed start vector, mirroring the symmetric Lanczos path
    * ([[eigSymmetric(op:gale\.linalg\.DoubleLinearOperator*]]).
    *
    * Each build forms an `m`-step Arnoldi factorization, solves the projected
    * `m×m` upper-Hessenberg Rayleigh quotient `VᵀAV` with the shared Francis QR
    * kernel (which re-Hessenbergs it internally — a minor waste at the projected
    * size), lifts the projected Ritz pairs to the big space (`V·s`), and measures
    * the '''true''' residual `‖A(v_re + i·v_im) − λ(v_re + i·v_im)‖` in real
    * arithmetic (two matvecs per complex pair, one per real value). A pair is
    * converged when that residual is `≤ tolerance · max(1, max|λ|)`; '''pairs
    * converge or not as pairs'''. If the wanted set has not all converged the
    * subspace '''grows''' (same start vector; exact at `m = n`), up to
    * `maxIterations` growth steps or a full `m = n` subspace.
    *
    * '''Selection on pairs / boundary rule.''' As in the dense path a conjugate
    * pair is one unit and is never split, so a [[EigenSelection.Count]] whose
    * `k`-boundary falls inside a pair targets the whole pair (`k+1` values).
    * `diagnostics.requested` reports that never-split target (`k`, or `k+1` at a
    * boundary pair); `converged` counts how many of the targeted eigenvalues met
    * the tolerance, so `allConverged` stays a reliable success signal. (When a
    * happy breakdown yields an invariant subspace smaller than `k`, `requested`
    * stays `k` and `converged < k`, so `allConverged` is correctly false.)
    *
    * '''Left vectors deferred.''' Only right eigenvectors are available (kernel
    * limitation), so `options.returnVectors` must be [[EigenVectors.ValuesOnly]] or
    * [[EigenVectors.Right]]; [[EigenVectors.Left]]/[[EigenVectors.LeftAndRight]]
    * return `Left(UnsupportedOperation)`. Orthogonality error is `0.0` (the Ritz
    * vectors of a nonsymmetric matrix are not orthonormal).
    *
    * Non-convergence is a `Right`: the result holds only the converged pairs (in
    * the canonical order), with `converged < requested` and per-pair residuals
    * recorded; `requireConverged` is the caller's opt-in to a
    * `Left(DidNotConverge)`.
    *
    * `Left` on: a non-square operator (`NonSquareMatrix`) or a square one whose
    * shape disagrees with `n` (`DimensionMismatch`); a non-positive `n`; `k ≤ 0` or
    * `k ≥ n-1` (§ 7's `k < n-1`; the message points at the dense API); a
    * non-[[EigenSelection.Count]] selection; an [[EigenOrder]] illegal for a
    * nonsymmetric problem (algebraic / both-ends); a start vector of the wrong
    * length or zero norm; or a `target` (shift-invert / `Around`), which is
    * '''deferred''' — `Left(UnsupportedOperation)`.
    */
  def eigNonsymmetric(
      op: DoubleLinearOperator,
      n: Int,
      selection: EigenSelection,
      options: SpectralOptions = SpectralOptions(),
      target: Option[SpectralTarget] = None
  ): Either[LinAlgError, NonsymmetricEigenDecomposition] =
    target match
      case Some(t) =>
        Left(LinAlgError.UnsupportedOperation(s"shift-invert / targeted selection ($t)"))
      case None =>
        selection match
          case count: EigenSelection.Count =>
            eigNonsymmetricArnoldi(op, n, count.k, count.order, options)
          case other =>
            Left(
              LinAlgError.InvalidArgument(
                s"the iterative solver requires EigenSelection.Count; ${other.getClass.getSimpleName} is dense-only"
              )
            )

  private def eigNonsymmetricArnoldi(
      op: DoubleLinearOperator,
      n: Int,
      k: Int,
      order: EigenOrder,
      options: SpectralOptions
  ): Either[LinAlgError, NonsymmetricEigenDecomposition] =
    if n <= 0 then Left(LinAlgError.InvalidArgument(s"dimension must be positive, got $n"))
    else if op.rows != op.cols then
      Left(LinAlgError.NonSquareMatrix(Shape(Rows(op.rows), Cols(op.cols))))
    else if op.rows != n then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(n), Cols(n)), Shape(Rows(op.rows), Cols(op.cols))))
    else
      validateArnoldiVectors(options.returnVectors) match
        case Left(error) => Left(error)
        case Right(wantVectors) =>
          if !nonsymmetricOrderLegal(order) then
            Left(LinAlgError.InvalidArgument(s"$order is symmetric-only; use a magnitude or real-part order"))
          else if k <= 0 || k >= n - 1 then
            Left(
              LinAlgError.InvalidArgument(
                s"k=$k must be in [1, ${n - 2}] for the iterative nonsymmetric solver (k < n-1); use the dense eigNonsymmetric for the full spectrum"
              )
            )
          else
            startVectorFor(options.startVector, n) match
              case Left(error) => Left(error)
              case Right(v0)   => runArnoldi(op, n, k, order, options, wantVectors, v0)

  // ===========================================================================
  // Generalized symmetric-definite eigendecomposition (A x = λ B x, B SPD)
  // ===========================================================================

  /** Generalized symmetric-definite eigendecomposition `A x = λ B x` computing
    * `B`-orthonormal eigenvectors ([[EigenVectors.Right]]). See the four-argument
    * overload for the vector-flag and failure details.
    */
  def eigSymmetricGeneralized(
      a: DMat,
      b: DMat,
      selection: EigenSelection
  ): Either[LinAlgError, EigenDecomposition] =
    eigSymmetricGeneralized(a, b, selection, EigenVectors.Right)

  /** Generalized symmetric-definite eigendecomposition `A x = λ B x` with `A`
    * symmetric and `B` '''symmetric positive-definite''' (phase b, § 4 of
    * `docs/spectral-parity.md`; SciPy `eigh(A, B)` type 1 — types 2/3 are out).
    *
    * Reduces to a standard symmetric problem by Cholesky of `B` (`B = L Lᵀ`):
    * `C = L⁻¹ A L⁻ᵀ` is symmetric with the same real eigenvalues `λ`, solved by the
    * shared tridiagonal QL/QR kernel, and the eigenvectors back-transform as
    * `x = L⁻ᵀ y`. `C` is formed without inverses — solve `L Y = A`, then (since `C`
    * is symmetric) `L C = Yᵀ` — and symmetrized before the kernel. The result reuses
    * [[EigenDecomposition]]: real eigenvalues '''ascending-algebraic always''' (the
    * selection chooses membership, not layout), real eigenvectors.
    *
    * '''`B`-orthonormality.''' The returned eigenvectors are '''`B`-orthonormal'''
    * (`Xᵀ B X = I`), not Euclidean-orthonormal — the natural normalization for the
    * `B`-inner-product this problem lives in. Accordingly `diagnostics.residuals`
    * are the '''true generalized residuals''' `‖A x − λ B x‖` and
    * `diagnostics.orthogonalityError` is `‖Xᵀ B X − I‖_F` (both zero /
    * not-computed when values-only).
    *
    * '''Conditioning of `B`.''' The Cholesky reduction amplifies error roughly
    * with `κ(B)` (the standard `sygv`-family sensitivity): for an ill-conditioned
    * `B` the call still succeeds and `allConverged` is structurally `true` on
    * this dense one-shot path, so '''`diagnostics.residuals` is the honest
    * accuracy signal''' — check it when `B` is near-singular.
    *
    * Like [[eigSymmetric(a:gale\.linalg\.DMat*]], only the '''lower triangle''' of
    * `A` and of `B` is read (the `Cholesky` precedent); the strict upper triangles
    * are treated as their mirrors. `selection` is realized as a slice of the full
    * ascending spectrum, with the same legality as the symmetric dense path
    * ([[EigenSelection.Count]] with an algebraic/magnitude order,
    * [[EigenSelection.IndexRange]], [[EigenSelection.ValueInterval]] all legal;
    * real-part orders rejected). `vectors` selects [[EigenVectors.ValuesOnly]] vs
    * [[EigenVectors.Right]].
    *
    * '''Scope.''' Only the dense path ships here; the operator/iterative
    * generalized solver (§ 6's `B`-inner-product Lanczos) and the generalized
    * '''nonsymmetric''' pencil / QZ (§ 5, backend-scoped) are '''out''' of this
    * entry point.
    *
    * `Left` on: non-square `a` or `b` (`NonSquareMatrix`); disagreeing shapes
    * (`DimensionMismatch`); `B` not positive-definite (`NotPositiveDefinite`, the
    * same `Left` the dense `Cholesky` returns); an [[EigenOrder]] illegal for a
    * symmetric problem; `k` outside `[1, n]`; an out-of-bounds `IndexRange`; an
    * inverted `ValueInterval`; [[EigenVectors.Left]]/[[EigenVectors.LeftAndRight]];
    * or (in practice unreachable) kernel non-convergence.
    */
  def eigSymmetricGeneralized(
      a: DMat,
      b: DMat,
      selection: EigenSelection,
      vectors: EigenVectors
  ): Either[LinAlgError, EigenDecomposition] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else if b.rows != b.cols then Left(LinAlgError.NonSquareMatrix(b.shape))
    else if a.rows != b.rows then Left(LinAlgError.DimensionMismatch(a.shape, b.shape))
    else
      val n = a.rows
      validateVectors(vectors) match
        case Left(error) => Left(error)
        case Right(wantVectors) =>
          validateDenseSelection(selection, n) match
            case Left(error) => Left(error)
            case Right(()) =>
              DenseDecompositions.cholesky(b) match
                // B not SPD: reuse the Cholesky path's own Left (NotPositiveDefinite).
                case Left(error) => Left(error)
                case Right(chol) =>
                  val l = chol.lower
                  val aSym = mirrorLower(a)
                  val bSym = mirrorLower(b)
                  reduceToStandard(l, aSym, n) match
                    case Left(error) => Left(error)
                    case Right(c) =>
                      DenseSpectralKernels.symmetricEigen(c, wantVectors) match
                        case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
                          Left(LinAlgError.DidNotConverge(iters, 0.0))
                        case Right(kernel) =>
                          val indices = denseSelectionIndices(selection, kernel.values, n)
                          if wantVectors then
                            backTransformVectors(l, kernel.vectors.get, indices) match
                              case Left(error) => Left(error)
                              case Right(x) =>
                                Right(assembleGeneralized(aSym, bSym, kernel.values, indices, Some(x)))
                          else Right(assembleGeneralized(aSym, bSym, kernel.values, indices, None))

  // ===========================================================================
  // Generalized nonsymmetric eigendecomposition (QZ) — backend-scoped
  // ===========================================================================

  /** Generalized '''nonsymmetric''' eigendecomposition `A x = λ B x` of a general
    * pencil (parity § 5 / § 1.3 of `docs/spectral-backend-boundary.md`). This is a
    * '''backend-scoped''' seam: the pure core ships no QZ engine, so with the
    * default `given SpectralBackend` ([[SpectralBackend.none]]) — the only one in
    * scope unless an acceleration module is imported — it returns
    * `Left(UnsupportedOperation)`. A [[SpectralCapability.GeneralizedNonsymmetricEigen]]-capable
    * backend supplies the raw `(α, β)` spectrum and vectors; this facade validates
    * shape first, then '''canonicalizes''' the backend's output — the projective
    * ordering ([[generalizedIndices]]: infinities first, pairs adjacent), the
    * packing, and the re-derived homogeneous residuals `‖β A x − α B x‖` — and
    * assembles the sealed [[GeneralizedEigenDecomposition]]. The order is gale's,
    * never the engine's.
    *
    * Eigenvalues are '''projective''' `(α, β)`: `β = 0` marks an infinite eigenvalue
    * from a singular / rank-deficient `B`. Output order is the canonical
    * nonsymmetric order applied to `α/β` (descending magnitude, infinities first),
    * conjugate pairs adjacent (positive-imaginary first, equal `β`).
    *
    * `Left` on: non-square `a`/`b` (`NonSquareMatrix`); shape disagreement
    * (`DimensionMismatch`); or no capable backend (`UnsupportedOperation`).
    */
  def eigGeneralizedNonsymmetric(a: DMat, b: DMat, vectors: EigenVectors = EigenVectors.Right)(using
      backend: SpectralBackend
  ): Either[LinAlgError, GeneralizedEigenDecomposition] =
    validateSquarePencil(a, b) match
      case Left(error) => Left(error)
      case Right(_) =>
        if !backend.capabilities.contains(SpectralCapability.GeneralizedNonsymmetricEigen) then
          Left(
            LinAlgError.UnsupportedOperation(
              s"generalized nonsymmetric eigen (QZ): no spectral backend registered (backend: ${backend.name})"
            )
          )
        else
          backend.generalizedNonsymmetricEigen(a, b, vectors) match
            case Left(error) => Left(error)
            case Right(raw)  => assembleQz(a, b, raw)

  // ===========================================================================
  // Validation helpers
  // ===========================================================================

  /** Map the symmetric vector flag to "compute vectors?". `Left`/`LeftAndRight`
    * are nonsymmetric-only and rejected.
    */
  private def validateVectors(vectors: EigenVectors): Either[LinAlgError, Boolean] =
    vectors match
      case EigenVectors.ValuesOnly => Right(false)
      case EigenVectors.Right      => Right(true)
      case EigenVectors.Left | EigenVectors.LeftAndRight =>
        Left(
          LinAlgError.InvalidArgument(
            "symmetric eigenvectors: left and right coincide; use EigenVectors.Right or EigenVectors.ValuesOnly"
          )
        )

  private def symmetricOrderLegal(order: EigenOrder): Boolean =
    order match
      case EigenOrder.LargestRealPart | EigenOrder.SmallestRealPart => false
      case _                                                        => true

  private def validateDenseSelection(selection: EigenSelection, n: Int): Either[LinAlgError, Unit] =
    selection match
      case EigenSelection.All => Right(())
      case EigenSelection.Count(k, order) =>
        if !symmetricOrderLegal(order) then
          Left(LinAlgError.InvalidArgument(s"$order is nonsymmetric-only; use an algebraic, magnitude, or both-ends order"))
        else if k <= 0 || k > n then
          Left(LinAlgError.InvalidArgument(s"k=$k must be in [1, $n]"))
        else Right(())
      case EigenSelection.IndexRange(from, to) =>
        if from < 0 || to >= n || from > to then
          Left(LinAlgError.InvalidArgument(s"index range [$from, $to] must satisfy 0 <= from <= to <= ${n - 1}"))
        else Right(())
      case EigenSelection.ValueInterval(lower, upper) =>
        if lower > upper || lower.isNaN || upper.isNaN then
          Left(LinAlgError.InvalidArgument(s"value interval ($lower, $upper] must have lower <= upper"))
        else Right(())

  // ===========================================================================
  // Dense selection + assembly
  // ===========================================================================

  /** The ascending-order indices of `values` selected by `selection`. */
  private def denseSelectionIndices(selection: EigenSelection, values: DVec, n: Int): Array[Int] =
    selection match
      case EigenSelection.All             => (0 until n).toArray
      case EigenSelection.Count(k, order) => selectExtremeIndices(values, math.min(k, n), order)
      case EigenSelection.IndexRange(from, to) => (from to to).toArray
      case EigenSelection.ValueInterval(lower, upper) =>
        (0 until n).filter(i => values(i) > lower && values(i) <= upper).toArray

  /** The ascending indices of the `k` eigenvalues at the extreme named by `order`
    * (`values` is ascending). Ties break by index for determinism.
    */
  private def selectExtremeIndices(values: DVec, k: Int, order: EigenOrder): Array[Int] =
    val n = values.length
    val chosen: Array[Int] =
      order match
        case EigenOrder.SmallestAlgebraic => (0 until k).toArray
        case EigenOrder.LargestAlgebraic  => ((n - k) until n).toArray
        case EigenOrder.SmallestMagnitude =>
          (0 until n).sortBy(i => (math.abs(values(i)), i)).take(k).toArray
        case EigenOrder.LargestMagnitude =>
          (0 until n).sortBy(i => (-math.abs(values(i)), i)).take(k).toArray
        case EigenOrder.BothEnds =>
          val hi = (k + 1) / 2 // ceil(k/2) from the high end
          val lo = k / 2       // floor(k/2) from the low end
          ((0 until lo) ++ ((n - hi) until n)).toArray
        case EigenOrder.LargestRealPart | EigenOrder.SmallestRealPart =>
          Array.empty // unreachable: rejected in validation
    java.util.Arrays.sort(chosen)
    chosen

  /** Build the dense result: selected values (ascending), selected eigenvector
    * columns (or an empty matrix when values-only), and diagnostics with per-pair
    * residuals and orthogonality error computed from the actual matrix.
    */
  private def assembleDense(
      a: DMat,
      values: DVec,
      vectors: Option[DMat],
      indices: Array[Int],
      wantVectors: Boolean
  ): EigenDecomposition =
    val selValues = DVec.tabulate(indices.length)(i => values(indices(i)))
    val n = a.rows
    val (selVectors, residuals, orthoErr) =
      if wantVectors then
        val src = vectors.get
        val sel = DMat.tabulate(n, indices.length)((r, c) => src(r, indices(c)))
        // Residuals must be measured against the matrix actually decomposed: the
        // kernel reads only the lower triangle, so mirror it here too — otherwise
        // a non-mirror strict upper triangle yields spurious residuals against a
        // correct decomposition.
        val sym = DMat.tabulate(n, n)((i, j) => if i >= j then a(i, j) else a(j, i))
        (sel, perPairResiduals(v => sym * v, selValues, sel), orthogonalityError(sel))
      else (DMat.zeros(n, 0), DVec.zeros(indices.length), 0.0)
    val diagnostics =
      SpectralDiagnostics(
        requested = indices.length,
        converged = indices.length,
        residuals = residuals,
        orthogonalityError = orthoErr,
        iterations = 0,
        rank = None
      )
    EigenDecomposition(selValues, selVectors, diagnostics)

  // ===========================================================================
  // Lanczos
  // ===========================================================================

  private def runLanczos(
      op: DoubleLinearOperator,
      n: Int,
      k: Int,
      order: EigenOrder,
      options: SpectralOptions,
      wantVectors: Boolean,
      v0: DVec
  ): Either[LinAlgError, EigenDecomposition] =
    val ncv0 = math.min(n, math.max(options.subspaceDimension.getOrElse(math.max(2 * k + 1, 20)), k + 1))
    val maxSteps = math.max(options.maxIterations, 1)
    val tol = options.tolerance
    val growBy = math.max(k, 16)

    var m = ncv0
    var step = 0
    var done = false
    var failure: Option[LinAlgError] = None
    // The latest build's result. Growing from a fixed start vector improves the
    // extreme Ritz values (Cauchy interlacing) and is exact at m = n; per-pair
    // convergence is not strictly monotone, so under a maxIterations truncation
    // the final build is the reported one, not a guaranteed best-so-far.
    var result: EigenDecomposition = assembleRitz(n, k, Array.empty, Array.empty, Array.empty, wantVectors, 0)

    while step < maxSteps && !done && failure.isEmpty do
      val (basis, alpha, betaSub, mEff) = buildLanczos(op, n, v0, m)
      val diag = DVec.tabulate(mEff)(i => alpha(i))
      val off = DVec.tabulate(math.max(mEff - 1, 0))(i => betaSub(i))
      DenseSpectralKernels.symmetricTridiagonalEigen(diag, off, wantVectors = true) match
        case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
          // The projected tridiagonal is tiny; failing it is a genuine numerical
          // breakdown, not partial convergence — surface it rather than loop.
          failure = Some(LinAlgError.DidNotConverge(iters, 0.0))
        case Right(proj) =>
          val theta = proj.values
          val s = proj.vectors.get
          val anorm = math.max(1.0, maxAbs(theta))
          val wanted = selectExtremeIndices(theta, math.min(k, mEff), order)
          val ritz = wanted.map(col => ritzVector(basis, s, col, n, mEff))
          val residuals = wanted.indices.map(w => trueResidual(op, ritz(w), theta(wanted(w)))).toArray
          val convergedSlots = wanted.indices.filter(w => residuals(w) <= tol * anorm).toArray
          result = assembleRitz(
            n,
            k,
            convergedSlots.map(w => theta(wanted(w))),
            convergedSlots.map(w => ritz(w)),
            convergedSlots.map(w => residuals(w)),
            wantVectors,
            step + 1
          )
          val fullyConverged = convergedSlots.length == wanted.length && wanted.length == k
          // A happy breakdown (mEff < m) means the Krylov space is invariant — the
          // full information reachable from this start — so growing cannot help.
          if fullyConverged || mEff < m || m >= n then done = true
          else m = math.min(n, m + growBy)
      step += 1

    failure match
      case Some(error) => Left(error)
      case None        => Right(result)

  /** Build an `ncv`-step Lanczos factorization from the unit start vector, with
    * full reorthogonalization (classical Gram–Schmidt twice) against all prior
    * basis vectors. Returns `(basis, alpha, betaSub, m)` where `m <= ncv` is the
    * number of vectors actually built (`m < ncv` marks a happy breakdown — an
    * invariant subspace), `alpha(0..m-1)` the tridiagonal diagonal, and
    * `betaSub(0..m-2)` its subdiagonal.
    */
  private def buildLanczos(
      op: DoubleLinearOperator,
      n: Int,
      v0: DVec,
      ncv: Int
  ): (Array[DVec], Array[Double], Array[Double], Int) =
    val basis = new Array[DVec](ncv)
    val alpha = new Array[Double](ncv)
    val betaSub = new Array[Double](ncv)
    var current = v0
    var previous: DVec = null
    var betaPrev = 0.0
    var m = 0
    var j = 0
    var breakdown = false
    while j < ncv && !breakdown do
      basis(j) = current
      val w = MutableDVec.zeros(n)
      op.applyTo(current, w)
      val a = current.dot(w.asVec)
      alpha(j) = a
      w.axpyInPlace(-a, current)
      if previous != null then w.axpyInPlace(-betaPrev, previous)
      reorthogonalize(w, basis, j + 1)
      reorthogonalize(w, basis, j + 1)
      val beta = w.asVec.norm2
      m = j + 1
      if j < ncv - 1 then
        if beta <= 1e-12 * math.max(1.0, math.abs(a)) then breakdown = true
        else
          betaSub(j) = beta
          previous = current
          betaPrev = beta
          current = (w.asVec * (1.0 / beta)).copy
      j += 1
    (basis, alpha, betaSub, m)

  /** `w := (I − V Vᵀ) w` against the first `count` basis columns. */
  private def reorthogonalize(w: MutableDVec, basis: Array[DVec], count: Int): Unit =
    var i = 0
    while i < count do
      val c = basis(i).dot(w.asVec)
      w.axpyInPlace(-c, basis(i))
      i += 1

  /** The Ritz vector `V s(:,col)`, normalized to unit 2-norm. */
  private def ritzVector(basis: Array[DVec], s: DMat, col: Int, n: Int, m: Int): DVec =
    val y = MutableDVec.zeros(n)
    var j = 0
    while j < m do
      y.axpyInPlace(s(j, col), basis(j))
      j += 1
    val nrm = y.asVec.norm2
    if nrm > 0.0 then (y.asVec * (1.0 / nrm)).copy else y.asVec.copy

  /** True residual `‖op·y − θ y‖` of a Ritz pair. */
  private def trueResidual(op: DoubleLinearOperator, y: DVec, theta: Double): Double =
    val w = MutableDVec.zeros(y.length)
    op.applyTo(y, w)
    w.axpyInPlace(-theta, y)
    w.asVec.norm2

  /** Assemble the iterative result from the converged Ritz pairs, ascending. */
  private def assembleRitz(
      n: Int,
      requestedK: Int,
      values: Array[Double],
      vectors: Array[DVec],
      residuals: Array[Double],
      wantVectors: Boolean,
      iterations: Int
  ): EigenDecomposition =
    val order = values.indices.sortBy(i => (values(i), i)).toArray
    val sortedValues = DVec.tabulate(order.length)(i => values(order(i)))
    val sortedResiduals = DVec.tabulate(order.length)(i => residuals(order(i)))
    val (vectorsMat, orthoErr) =
      if wantVectors && order.nonEmpty then
        val mat = DMat.tabulate(n, order.length)((r, c) => vectors(order(c))(r))
        (mat, orthogonalityError(mat))
      else (DMat.zeros(n, 0), 0.0)
    val diagnostics =
      SpectralDiagnostics(
        requested = requestedK,
        converged = order.length,
        residuals = sortedResiduals,
        orthogonalityError = orthoErr,
        iterations = iterations,
        rank = None
      )
    EigenDecomposition(sortedValues, vectorsMat, diagnostics)

  // ===========================================================================
  // Shared numeric helpers
  // ===========================================================================

  /** Per-column residual `‖A v_i − λ_i v_i‖` for eigenvector columns. */
  private def perPairResiduals(applyOp: DVec => DVec, values: DVec, vectors: DMat): DVec =
    DVec.tabulate(vectors.cols): c =>
      val v = vectors.col(c)
      val diff = applyOp(v) - v * values(c)
      diff.norm2

  /** `‖VᵀV − I‖_F` of the eigenvector columns. */
  private def orthogonalityError(vectors: DMat): Double =
    val g = vectors.t * vectors
    var sum = 0.0
    var i = 0
    while i < g.rows do
      var j = 0
      while j < g.cols do
        val d = if i == j then g(i, j) - 1.0 else g(i, j)
        sum += d * d
        j += 1
      i += 1
    math.sqrt(sum)

  private def maxAbs(v: DVec): Double =
    var m = 0.0
    var i = 0
    while i < v.length do
      val a = math.abs(v(i))
      if a > m then m = a
      i += 1
    m

  /** A deterministic (LCG) start vector, or the caller's, normalized to unit
    * length. The seed sequence is bit-for-bit portable (32-bit `Int` arithmetic
    * wraps identically on JVM and Scala.js); the Krylov iteration built on it is
    * deterministic '''per platform''' but may differ between JVM and Scala.js in
    * the last bits, because the dense kernels use the platform's fused
    * multiply-add (`Math.fma` on the JVM, `a*b + c` on JS).
    */
  private def startVectorFor(provided: Option[DVec], n: Int): Either[LinAlgError, DVec] =
    provided match
      case Some(v) =>
        if v.length != n then Left(LinAlgError.VectorLengthMismatch(n, v.length))
        else
          val nrm = v.norm2
          if nrm == 0.0 then Left(LinAlgError.InvalidArgument("start vector must be nonzero"))
          else Right((v * (1.0 / nrm)).copy)
      case None =>
        val out = new Array[Double](n)
        var state = 123456789
        var i = 0
        while i < n do
          state = state * 1103515245 + 12345
          out(i) = ((state >>> 9) & 0x7fffff).toDouble / 0x800000.toDouble * 2.0 - 1.0
          i += 1
        val v = DVec.fromSeq(out.toIndexedSeq)
        val nrm = v.norm2
        Right(if nrm == 0.0 then v else (v * (1.0 / nrm)).copy)

  // ===========================================================================
  // Nonsymmetric validation
  // ===========================================================================

  /** Map the nonsymmetric vector flag to `(computeRight, computeLeft)` for the
    * '''dense''' path. All four cases are supported: left eigenvectors are derived
    * from the full right-eigenvector matrix, so the kernel still computes right
    * vectors internally whenever left vectors are requested.
    */
  private def nonsymmetricVectorFlags(vectors: EigenVectors): Either[LinAlgError, (Boolean, Boolean)] =
    vectors match
      case EigenVectors.ValuesOnly   => Right((false, false))
      case EigenVectors.Right        => Right((true, false))
      case EigenVectors.Left         => Right((false, true))
      case EigenVectors.LeftAndRight => Right((true, true))

  /** Map the vector flag to "compute right vectors?" for the '''iterative
    * (Arnoldi)''' path. Left / left-and-right are '''unsupported''' there — a Krylov
    * basis for `A` carries no left-eigenvector information, and two-sided Arnoldi is
    * out of scope — so they route to the dense API via `Left(UnsupportedOperation)`.
    */
  private def validateArnoldiVectors(vectors: EigenVectors): Either[LinAlgError, Boolean] =
    vectors match
      case EigenVectors.ValuesOnly => Right(false)
      case EigenVectors.Right      => Right(true)
      case EigenVectors.Left | EigenVectors.LeftAndRight =>
        Left(
          LinAlgError.UnsupportedOperation(
            "left eigenvectors from the iterative solver (a Krylov basis for A gives no left vectors; use the dense eigNonsymmetric)"
          )
        )

  /** Nonsymmetric order legality: magnitude and real-part orders only; the
    * algebraic orders and [[EigenOrder.BothEnds]] are symmetric-only (§ Selection).
    */
  private def nonsymmetricOrderLegal(order: EigenOrder): Boolean =
    order match
      case EigenOrder.LargestMagnitude | EigenOrder.SmallestMagnitude                     => true
      case EigenOrder.LargestRealPart | EigenOrder.SmallestRealPart                       => true
      case EigenOrder.LargestAlgebraic | EigenOrder.SmallestAlgebraic | EigenOrder.BothEnds => false

  private def validateNonsymmetricDenseSelection(selection: EigenSelection, n: Int): Either[LinAlgError, Unit] =
    selection match
      case EigenSelection.All => Right(())
      case EigenSelection.Count(k, order) =>
        if !nonsymmetricOrderLegal(order) then
          Left(LinAlgError.InvalidArgument(s"$order is symmetric-only; use a magnitude or real-part order"))
        else if k <= 0 || k > n then Left(LinAlgError.InvalidArgument(s"k=$k must be in [1, $n]"))
        else Right(())
      case _: EigenSelection.IndexRange =>
        Left(
          LinAlgError.InvalidArgument(
            "IndexRange selection is symmetric-only (ascending-algebraic rank has no nonsymmetric meaning); use Count or All"
          )
        )
      case _: EigenSelection.ValueInterval =>
        Left(LinAlgError.InvalidArgument("ValueInterval selection is symmetric-only; use Count or All"))

  // ===========================================================================
  // Nonsymmetric spectrum units, ordering, and selection
  // ===========================================================================

  /** The base indices of the spectrum's eigenvalue '''units''' in packed order: a
    * real eigenvalue is a size-1 unit (one column); a conjugate pair a size-2 unit
    * (two adjacent columns, positive-imaginary member first). Walks `im` once,
    * relying on the kernel's packing (a unit always starts on a nonnegative-
    * imaginary index).
    */
  private def spectrumUnitBases(im: DVec): Array[Int] =
    val bases = scala.collection.mutable.ArrayBuffer.empty[Int]
    var i = 0
    while i < im.length do
      bases += i
      if im(i) > 0.0 then i += 2 else i += 1
    bases.toArray

  /** Order unit base indices by the selection criterion (§ 2 ordering guarantee):
    * the primary key follows the order's direction (largest ⇒ descending, smallest
    * ⇒ ascending, on magnitude or real part), ties broken by descending real part
    * then descending imaginary part, then base index for a total, deterministic
    * order. The representative eigenvalue of a pair is its positive-imaginary
    * member (`re(base)`, `im(base) > 0`).
    */
  private def orderUnits(re: DVec, im: DVec, bases: Array[Int], order: EigenOrder): Array[Int] =
    bases.sortBy: base =>
      val r = re(base)
      val iimag = im(base)
      val mag = math.hypot(r, iimag)
      val primary =
        order match
          case EigenOrder.LargestMagnitude  => -mag
          case EigenOrder.SmallestMagnitude => mag
          case EigenOrder.LargestRealPart   => -r
          case EigenOrder.SmallestRealPart  => r
          case _                            => -mag // unreachable: legality validated
      (primary, -r, -iimag, base)

  /** Expand ordered unit base indices into the full column/eigenvalue index list,
    * appending both members of a pair (positive first) in order.
    */
  private def expandUnits(orderedBases: Array[Int], im: DVec): Array[Int] =
    val out = scala.collection.mutable.ArrayBuffer.empty[Int]
    var u = 0
    while u < orderedBases.length do
      val base = orderedBases(u)
      out += base
      if im(base) > 0.0 then out += (base + 1)
      u += 1
    out.toArray

  /** The kernel-order indices selected by `selection`, in canonical output order,
    * with conjugate pairs kept whole (never split). [[EigenSelection.Count]]
    * accumulates whole units until at least `k` eigenvalues are covered — so a
    * boundary pair yields `k+1`. [[EigenSelection.All]] returns every unit in the
    * [[EigenOrder.LargestMagnitude]] canonical order.
    */
  private def nonsymDenseIndices(re: DVec, im: DVec, selection: EigenSelection): Array[Int] =
    val bases = spectrumUnitBases(im)
    selection match
      case EigenSelection.All =>
        expandUnits(orderUnits(re, im, bases, EigenOrder.LargestMagnitude), im)
      case EigenSelection.Count(k, order) =>
        val sorted = orderUnits(re, im, bases, order)
        val out = scala.collection.mutable.ArrayBuffer.empty[Int]
        var count = 0
        var u = 0
        while u < sorted.length && count < k do
          val base = sorted(u)
          out += base
          count += 1
          if im(base) > 0.0 then
            out += (base + 1)
            count += 1
          u += 1
        out.toArray
      case _ => Array.empty[Int] // unreachable: IndexRange/ValueInterval validated away

  // ===========================================================================
  // Nonsymmetric assembly + residuals
  // ===========================================================================

  /** Per-eigenvalue residual `‖A v − λ v‖` in real arithmetic for a packed
    * nonsymmetric result. A real eigenvalue uses its single column; a conjugate
    * pair uses its two adjacent columns `(v_re, v_im)` and both members share the
    * one residual. Returns zeros when no vectors are present.
    */
  private def nonsymResiduals(applyOp: DVec => DVec, re: DVec, im: DVec, packed: DMat): DVec =
    val cnt = re.length
    if packed.cols == 0 then DVec.zeros(cnt)
    else
      val out = new Array[Double](cnt)
      var j = 0
      while j < cnt do
        val b = im(j)
        if b == 0.0 then
          val v = packed.col(j)
          out(j) = (applyOp(v) - v * re(j)).norm2
          j += 1
        else
          // Positive-imaginary member at j; its conjugate is at j+1 (shared residual).
          val vr = packed.col(j)
          val vi = packed.col(j + 1)
          val a = re(j)
          val realPart = applyOp(vr) - (vr * a - vi * b)
          val imagPart = applyOp(vi) - (vi * a + vr * b)
          val res = math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))
          out(j) = res
          out(j + 1) = res
          j += 2
      DVec.fromSeq(out.toIndexedSeq)

  /** Per-eigenvalue '''left''' residual `‖wᴴ A − λ wᴴ‖` in real arithmetic, using
    * `at = Aᵀ`. For `λ = a + b·i` and `w = w_re + i·w_im` the (transposed) residual
    * is `[Aᵀ w_re − a·w_re − b·w_im] + i·[a·w_im − Aᵀ w_im − b·w_re]`. A real
    * eigenvalue uses its single column; a conjugate pair shares one residual.
    */
  private def nonsymLeftResiduals(at: DMat, re: DVec, im: DVec, packed: DMat): DVec =
    val cnt = re.length
    if packed.cols == 0 then DVec.zeros(cnt)
    else
      val out = new Array[Double](cnt)
      var j = 0
      while j < cnt do
        val b = im(j)
        if b == 0.0 then
          val w = packed.col(j)
          out(j) = ((at * w) - w * re(j)).norm2
          j += 1
        else
          val wr = packed.col(j)
          val wi = packed.col(j + 1)
          val a = re(j)
          val realPart = (at * wr) - (wr * a) - (wi * b)
          val imagPart = (wi * a) - (at * wi) - (wr * b)
          val res = math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))
          out(j) = res
          out(j + 1) = res
          j += 2
      DVec.fromSeq(out.toIndexedSeq)

  /** Build the dense nonsymmetric result: permute the kernel's re/im, the right
    * packed columns, and (when requested) the left packed columns by `outIdx` — the
    * left columns move in the '''same lockstep''' as the right ones — compute
    * per-eigenvalue residuals, and report all-converged diagnostics (a dense
    * one-shot solve). For [[EigenVectors.LeftAndRight]] the reported residual is the
    * worse of the right (`‖A v − λ v‖`) and left (`‖wᴴ A − λ wᴴ‖`) sides.
    * Orthogonality error is `0.0` — nonsymmetric eigenvectors are not orthonormal.
    */
  private def assembleNonsymDense(
      a: DMat,
      kernel: DenseSpectralKernels.NonsymmetricEigen,
      leftFull: Option[DMat],
      outIdx: Array[Int],
      wantRight: Boolean,
      wantLeft: Boolean
  ): NonsymmetricEigenDecomposition =
    val cnt = outIdx.length
    val re = DVec.tabulate(cnt)(j => kernel.re(outIdx(j)))
    val im = DVec.tabulate(cnt)(j => kernel.im(outIdx(j)))
    val n = a.rows
    val rightSel =
      if wantRight then DMat.tabulate(n, cnt)((r, c) => kernel.vectors.get(r, outIdx(c)))
      else DMat.zeros(n, 0)
    val leftSel =
      if wantLeft then Some(DMat.tabulate(n, cnt)((r, c) => leftFull.get(r, outIdx(c))))
      else None
    val residuals =
      (wantRight, wantLeft) match
        case (true, true) =>
          val rr = nonsymResiduals(v => a * v, re, im, rightSel)
          val lr = nonsymLeftResiduals(a.t, re, im, leftSel.get)
          DVec.tabulate(cnt)(i => math.max(rr(i), lr(i)))
        case (true, false) => nonsymResiduals(v => a * v, re, im, rightSel)
        case (false, true) => nonsymLeftResiduals(a.t, re, im, leftSel.get)
        case (false, false) => DVec.zeros(cnt)
    val diagnostics =
      SpectralDiagnostics(
        requested = cnt,
        converged = cnt,
        residuals = residuals,
        orthogonalityError = 0.0,
        iterations = 0,
        rank = None
      )
    new NonsymmetricEigenDecomposition(re, im, rightSel, leftSel, diagnostics)

  /** Relative pivot floor for the left-vector conditioning guard: reject when the
    * embedded LU's smallest pivot is `≤ LeftVectorPivotFactor · 2n · largest pivot`.
    * At `128·ε` it sits far above the `~ε` ratio a degenerate (defective) eigenbasis
    * leaves and far below any diagonalizable basis' ratio.
    */
  private inline val LeftVectorPivotFactor = 128.0 * 2.220446049250313e-16

  /** Relative left-residual ceiling for the backstop guard: reject when the worst
    * `‖wᴴ A − λ wᴴ‖` exceeds `LeftVectorResidualTolerance · max(1, ‖A‖)`. Loose
    * enough to admit honestly ill-conditioned bases, tight enough to catch a garbage
    * (e.g. higher-order Jordan) basis whose residual is `O(1)`.
    */
  private inline val LeftVectorResidualTolerance = 1e-2

  /** Recover the '''full''' left eigenvectors of `a` from its full right-eigenvector
    * matrix `V` (the kernel's packed columns) via the biorthogonality identity: the
    * left eigenvectors (`wᴴ A = λ wᴴ`) are the '''conjugated rows of `V⁻¹`'''. The
    * complex `V⁻¹` is formed with the real `2n×2n` embedding
    * `[[Vr, −Vi], [Vi, Vr]]` and gale's LU (no complex kernel tier), so `V⁻¹`'s row
    * `i` — hence left vector `i` — pairs with right column `i` and eigenvalue `i`
    * '''exactly by construction''' (no eigenvalue matching). Result is `n × n` in the
    * real-Schur packing.
    *
    * '''Defective-input guards.''' A defective (or near-defective) `a` does '''not'''
    * make the embedding exactly singular — the Francis kernel perturbs the
    * repeated/parallel eigenvectors just enough that the LU sees only near-zero pivots
    * and would silently return a degenerate (duplicated or garbage) left basis. Two
    * guards reject those as `Left(SingularMatrix)` (a full set of left eigenvectors
    * does not exist for a defective `a`): (1) a '''conditioning''' check on the
    * embedded LU pivots ([[LeftVectorPivotFactor]]) — this catches
    * duplicated/nilpotent bases whose per-vector residual is a deceptive `0`; and
    * (2) a '''residual''' backstop ([[LeftVectorResidualTolerance]]) catching a
    * garbage basis that slips past the pivot check. Inputs that pass both may still be
    * ill-conditioned; their honest, larger residuals surface through `diagnostics`.
    */
  private def computeLeftVectors(
      a: DMat,
      kernel: DenseSpectralKernels.NonsymmetricEigen
  ): Either[LinAlgError, DMat] =
    val n = a.rows
    val p = kernel.vectors.get
    val wi = kernel.im
    // Reconstruct the complex right-eigenvector matrix V = Vr + i·Vi from the packed
    // real-Schur columns: a real column is real; a conjugate pair (k, k+1) is
    // (v_re + i·v_im) at k and its conjugate at k+1.
    val vr = Array.ofDim[Double](n, n)
    val vi = Array.ofDim[Double](n, n)
    var k = 0
    while k < n do
      if wi(k) == 0.0 then
        var r = 0
        while r < n do
          vr(r)(k) = p(r, k)
          vi(r)(k) = 0.0
          r += 1
        k += 1
      else
        var r = 0
        while r < n do
          val pk = p(r, k)
          val pk1 = p(r, k + 1)
          vr(r)(k) = pk
          vi(r)(k) = pk1
          vr(r)(k + 1) = pk
          vi(r)(k + 1) = -pk1
          r += 1
        k += 2
    val twoN = 2 * n
    val embedded = DMat.tabulate(twoN, twoN): (i, j) =>
      if i < n && j < n then vr(i)(j)
      else if i < n then -vi(i)(j - n)
      else if j < n then vi(i - n)(j)
      else vr(i - n)(j - n)
    DenseDecompositions.lu(embedded) match
      case Left(error) => Left(error)
      case Right(lu) =>
        // Guard 1 (conditioning): the LU pivots are the U diagonal. A degenerate
        // (defective) eigenbasis leaves a ~ε pivot ratio; reject well above it.
        var minPivot = Double.MaxValue
        var maxPivot = 0.0
        var pivotArg = 0
        var d = 0
        while d < twoN do
          val piv = math.abs(lu.packed(d, d))
          if piv < minPivot then
            minPivot = piv
            pivotArg = d
          if piv > maxPivot then maxPivot = piv
          d += 1
        if maxPivot == 0.0 || minPivot <= LeftVectorPivotFactor * twoN.toDouble * maxPivot then
          Left(LinAlgError.SingularMatrix(pivotArg))
        else
          // Solve embedded · [Xr; Xi] = [I; 0] one column at a time; V⁻¹ = Xr + i·Xi.
          val xr = Array.ofDim[Double](n, n)
          val xi = Array.ofDim[Double](n, n)
          var j = 0
          var failure: Option[LinAlgError] = None
          while j < n && failure.isEmpty do
            lu.solve(DVec.tabulate(twoN)(idx => if idx == j then 1.0 else 0.0)) match
              case Left(error) => failure = Some(error)
              case Right(sol) =>
                var i = 0
                while i < n do
                  xr(i)(j) = sol(i)
                  xi(i)(j) = sol(n + i)
                  i += 1
            j += 1
          failure match
            case Some(error) => Left(error)
            case None =>
              val packed = packLeftVectors(n, wi, xr, xi)
              // Guard 2 (residual backstop): a garbage basis that slipped past the
              // pivot check has a large left residual on the full spectrum.
              val residuals = nonsymLeftResiduals(a.t, kernel.re, kernel.im, packed)
              var worst = 0.0
              var i = 0
              while i < residuals.length do
                if residuals(i) > worst then worst = residuals(i)
                i += 1
              if worst > LeftVectorResidualTolerance * math.max(1.0, frobeniusNorm(a)) then
                Left(LinAlgError.SingularMatrix(0))
              else Right(packed)

  /** The Frobenius norm `‖A‖_F`. */
  private def frobeniusNorm(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        val v = a(i, j)
        sum += v * v
        j += 1
      i += 1
    math.sqrt(sum)

  /** Pack the left eigenvectors (conjugated rows of `V⁻¹ = Xr + i·Xi`) into the
    * real-Schur convention, unit 2-norm. Left vector `i` is `conj(row i of V⁻¹)`, so
    * its real part is row `i` of `Xr` and its imaginary part is `−(row i of Xi)`. A
    * real eigenvalue keeps only the (real) real part; a conjugate pair stores the
    * positive member's real part in column `i` and imaginary part in column `i+1`.
    */
  private def packLeftVectors(n: Int, wi: DVec, xr: Array[Array[Double]], xi: Array[Array[Double]]): DMat =
    val out = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      if wi(i) == 0.0 then
        var r = 0
        var nrm = 0.0
        while r < n do
          val v = xr(i)(r)
          out(r)(i) = v
          nrm += v * v
          r += 1
        nrm = math.sqrt(nrm)
        if nrm > 0.0 then
          r = 0
          while r < n do
            out(r)(i) = out(r)(i) / nrm
            r += 1
        i += 1
      else
        var r = 0
        var nrm = 0.0
        while r < n do
          val wRe = xr(i)(r)
          val wIm = -xi(i)(r)
          out(r)(i) = wRe
          out(r)(i + 1) = wIm
          nrm += wRe * wRe + wIm * wIm
          r += 1
        nrm = math.sqrt(nrm)
        if nrm > 0.0 then
          r = 0
          while r < n do
            out(r)(i) = out(r)(i) / nrm
            out(r)(i + 1) = out(r)(i + 1) / nrm
            r += 1
        i += 2
    DMat.tabulate(n, n)((r, c) => out(r)(c))

  // ===========================================================================
  // Arnoldi
  // ===========================================================================

  private def runArnoldi(
      op: DoubleLinearOperator,
      n: Int,
      k: Int,
      order: EigenOrder,
      options: SpectralOptions,
      wantVectors: Boolean,
      v0: DVec
  ): Either[LinAlgError, NonsymmetricEigenDecomposition] =
    // Floor the initial subspace at k+2 (room for a boundary conjugate pair beyond
    // the k-th value), clamped to n; k < n-1 guarantees k+2 <= n.
    val ncv0 = math.min(n, math.max(options.subspaceDimension.getOrElse(math.max(2 * k + 1, 20)), k + 2))
    val maxSteps = math.max(options.maxIterations, 1)
    val tol = options.tolerance
    val growBy = math.max(k, 16)

    var m = ncv0
    var step = 0
    var done = false
    var failure: Option[LinAlgError] = None
    // Latest build's result. Growing from a fixed start vector improves the extreme
    // Ritz values and is exact at m = n; per-pair convergence is not strictly
    // monotone, so under a maxIterations truncation the final build is reported.
    var result: NonsymmetricEigenDecomposition = emptyNonsym(n, k)

    while step < maxSteps && !done && failure.isEmpty do
      val (basis, h, mEff) = buildArnoldi(op, n, v0, m)
      val hMat = DMat.tabulate(mEff, mEff)((i, j) => h(i)(j))
      DenseSpectralKernels.nonsymmetricEigen(hMat, wantVectors = true) match
        case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
          // The projected Hessenberg is tiny; failing it is a genuine numerical
          // breakdown, not partial convergence — surface it rather than loop.
          failure = Some(LinAlgError.DidNotConverge(iters, 0.0))
        case Right(proj) =>
          val pRe = proj.re
          val pIm = proj.im
          val s = proj.vectors.get
          val anorm = math.max(1.0, maxMagnitude(pRe, pIm))
          val sorted = orderUnits(pRe, pIm, spectrumUnitBases(pIm), order)

          // Accumulate whole target units (never split) until >= k eigenvalues.
          val targetBases = scala.collection.mutable.ArrayBuffer.empty[Int]
          var targetSize = 0
          var u = 0
          while u < sorted.length && targetSize < k do
            val base = sorted(u)
            targetBases += base
            targetSize += (if pIm(base) > 0.0 then 2 else 1)
            u += 1

          // Lift each target unit to the big space, measure its true residual, and
          // keep it iff converged. Output stays in the sorted (canonical) order.
          val outRe = scala.collection.mutable.ArrayBuffer.empty[Double]
          val outIm = scala.collection.mutable.ArrayBuffer.empty[Double]
          val outCols = scala.collection.mutable.ArrayBuffer.empty[DVec]
          val outRes = scala.collection.mutable.ArrayBuffer.empty[Double]
          var allTargetsConverged = true
          var t = 0
          while t < targetBases.length do
            val base = targetBases(t)
            if pIm(base) > 0.0 then
              val a = pRe(base)
              val b = pIm(base)
              val (vr, vi) = ritzPair(basis, s, base, n, mEff)
              val realPart = (op * vr) - (vr * a - vi * b)
              val imagPart = (op * vi) - (vi * a + vr * b)
              val res = math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))
              if res <= tol * anorm then
                outRe += a
                outRe += a
                outIm += b
                outIm += (-b)
                outCols += vr
                outCols += vi
                outRes += res
                outRes += res
              else allTargetsConverged = false
            else
              val a = pRe(base)
              val vr = ritzReal(basis, s, base, n, mEff)
              val res = ((op * vr) - vr * a).norm2
              if res <= tol * anorm then
                outRe += a
                outIm += 0.0
                outCols += vr
                outRes += res
              else allTargetsConverged = false
            t += 1

          // requested = never-split target (k, or k+1 at a boundary pair); but when
          // an invariant subspace is too small to reach k, keep requested = k so
          // allConverged stays false.
          result = assembleArnoldi(
            n,
            math.max(k, targetSize),
            outRe.toArray,
            outIm.toArray,
            outCols.toArray,
            outRes.toArray,
            wantVectors,
            step + 1
          )
          val fullyConverged = targetSize >= k && allTargetsConverged
          // A happy breakdown (mEff < m) means the reachable Krylov space is
          // invariant, so growing cannot add information.
          if fullyConverged || mEff < m || m >= n then done = true
          else m = math.min(n, m + growBy)
      step += 1

    failure match
      case Some(error) => Left(error)
      case None        => Right(result)

  /** Build an `ncv`-step Arnoldi factorization from the unit start vector with full
    * reorthogonalization (classical Gram–Schmidt twice, accumulating the projection
    * coefficients into the Hessenberg column). Returns `(basis, h, mEff)` where
    * `mEff <= ncv` is the number of vectors built (`mEff < ncv` marks a happy
    * breakdown — an invariant subspace) and `h(i)(j)` is the `mEff x mEff` upper-
    * Hessenberg Rayleigh quotient `VᵀAV`; the coupling entry `h(mEff, mEff-1)`
    * beyond the square block is intentionally not stored (the true residual is
    * measured directly, so the cheap estimate is not needed).
    */
  private def buildArnoldi(
      op: DoubleLinearOperator,
      n: Int,
      v0: DVec,
      ncv: Int
  ): (Array[DVec], Array[Array[Double]], Int) =
    val basis = new Array[DVec](ncv)
    val h = Array.ofDim[Double](ncv, ncv)
    var current = v0
    var mEff = 0
    var j = 0
    var breakdown = false
    while j < ncv && !breakdown do
      basis(j) = current
      val w = MutableDVec.zeros(n)
      op.applyTo(current, w)
      val scale = w.asVec.norm2
      // First CGS pass: h(i,j) = v_i · (A v_j); deflate.
      var i = 0
      while i <= j do
        val c = basis(i).dot(w.asVec)
        h(i)(j) = c
        w.axpyInPlace(-c, basis(i))
        i += 1
      // Second pass (reorthogonalization): add the residual projections back into h.
      i = 0
      while i <= j do
        val c = basis(i).dot(w.asVec)
        h(i)(j) += c
        w.axpyInPlace(-c, basis(i))
        i += 1
      val beta = w.asVec.norm2
      mEff = j + 1
      if j < ncv - 1 then
        if beta <= 1e-12 * math.max(1.0, scale) then breakdown = true
        else
          h(j + 1)(j) = beta
          current = (w.asVec * (1.0 / beta)).copy
      j += 1
    (basis, h, mEff)

  /** The big-space Ritz vectors of a projected conjugate pair at packed column
    * `base`: `(V·s(:,base), V·s(:,base+1))` — the real and imaginary parts —
    * renormalized together to unit complex norm (`‖v_re‖² + ‖v_im‖² = 1`).
    */
  private def ritzPair(basis: Array[DVec], s: DMat, base: Int, n: Int, mEff: Int): (DVec, DVec) =
    val vr = combine(basis, s, base, n, mEff)
    val vi = combine(basis, s, base + 1, n, mEff)
    val nrm = math.sqrt(vr.dot(vr) + vi.dot(vi))
    if nrm > 0.0 then (vr * (1.0 / nrm), vi * (1.0 / nrm)) else (vr, vi)

  /** The big-space Ritz vector `V·s(:,base)` of a projected real eigenvalue,
    * normalized to unit 2-norm.
    */
  private def ritzReal(basis: Array[DVec], s: DMat, base: Int, n: Int, mEff: Int): DVec =
    val v = combine(basis, s, base, n, mEff)
    val nrm = v.norm2
    if nrm > 0.0 then v * (1.0 / nrm) else v

  /** `Σ_j s(j, col) · basis(j)` in `R^n`. */
  private def combine(basis: Array[DVec], s: DMat, col: Int, n: Int, mEff: Int): DVec =
    val y = MutableDVec.zeros(n)
    var j = 0
    while j < mEff do
      y.axpyInPlace(s(j, col), basis(j))
      j += 1
    y.asVec.copy

  /** The largest eigenvalue magnitude in a packed re/im spectrum. */
  private def maxMagnitude(re: DVec, im: DVec): Double =
    var m = 0.0
    var i = 0
    while i < re.length do
      val mag = math.hypot(re(i), im(i))
      if mag > m then m = mag
      i += 1
    m

  /** Assemble the iterative nonsymmetric result from the converged target units,
    * already in canonical order (pairs adjacent, positive-imaginary first).
    * `requested` is the never-split target count. Orthogonality error is `0.0`
    * (nonsymmetric Ritz vectors are not orthonormal).
    */
  private def assembleArnoldi(
      n: Int,
      requested: Int,
      outRe: Array[Double],
      outIm: Array[Double],
      outCols: Array[DVec],
      outRes: Array[Double],
      wantVectors: Boolean,
      iterations: Int
  ): NonsymmetricEigenDecomposition =
    val cnt = outRe.length
    val re = DVec.fromSeq(outRe.toIndexedSeq)
    val im = DVec.fromSeq(outIm.toIndexedSeq)
    val vectors =
      if wantVectors && cnt > 0 then DMat.tabulate(n, cnt)((r, c) => outCols(c)(r))
      else DMat.zeros(n, 0)
    val diagnostics =
      SpectralDiagnostics(
        requested = requested,
        converged = cnt,
        residuals = DVec.fromSeq(outRes.toIndexedSeq),
        orthogonalityError = 0.0,
        iterations = iterations,
        rank = None
      )
    new NonsymmetricEigenDecomposition(re, im, vectors, None, diagnostics)

  /** An empty iterative result (no pairs converged yet), used before the first
    * build; `requested` is the caller's `k`.
    */
  private def emptyNonsym(n: Int, requested: Int): NonsymmetricEigenDecomposition =
    new NonsymmetricEigenDecomposition(
      DVec.zeros(0),
      DVec.zeros(0),
      DMat.zeros(n, 0),
      None,
      SpectralDiagnostics(requested, 0, DVec.zeros(0), 0.0, 0, None)
    )

  // ===========================================================================
  // Generalized symmetric-definite reduction + assembly
  // ===========================================================================

  /** Mirror the lower triangle of `m` into a full symmetric matrix (the
    * lower-triangle-only read convention shared with `Cholesky` / the symmetric
    * kernel): `out(i,j) = m(i,j)` for `i >= j`, else `m(j,i)`.
    */
  private def mirrorLower(m: DMat): DMat =
    val n = m.rows
    DMat.tabulate(n, n)((i, j) => if i >= j then m(i, j) else m(j, i))

  /** Solve the triangular system `tri · X = M` column by column, `M` given by
    * `colOf`, returning `X` as an `n × cols` matrix. `lower` picks forward vs back
    * substitution. A `Left` (singular triangle) propagates — unreachable for a
    * factor `L` from a successful SPD Cholesky (strictly positive diagonal).
    */
  private def solveTriColumns(
      tri: DMat,
      lower: Boolean,
      cols: Int,
      colOf: Int => DVec
  ): Either[LinAlgError, DMat] =
    val n = tri.rows
    val xs = new Array[DVec](cols)
    var j = 0
    while j < cols do
      val solved = if lower then TriangularSolve.lower(tri, colOf(j)) else TriangularSolve.upper(tri, colOf(j))
      solved match
        case Left(error) => return Left(error)
        case Right(x)    => xs(j) = x
      j += 1
    Right(DMat.tabulate(n, cols)((i, jj) => xs(jj)(i)))

  /** Reduce `A x = λ B x` to the standard symmetric `C y = λ y` with
    * `C = L⁻¹ A L⁻ᵀ`, without forming inverses: solve `L Y = A` (`Y = L⁻¹A`), then
    * — since `C = L⁻¹ A L⁻ᵀ` is symmetric, `C = Cᵀ = L⁻¹ Yᵀ` — solve `L C = Yᵀ`.
    * The result is symmetrized `(C + Cᵀ)/2` so any rounding asymmetry does not
    * reach the kernel (which reads one triangle).
    */
  private def reduceToStandard(l: DMat, aSym: DMat, n: Int): Either[LinAlgError, DMat] =
    solveTriColumns(l, lower = true, n, j => aSym.col(j)) match
      case Left(error) => Left(error)
      case Right(y) =>
        // Column j of Yᵀ is row j of Y.
        solveTriColumns(l, lower = true, n, j => DVec.tabulate(n)(i => y(j, i))) match
          case Left(error) => Left(error)
          case Right(c) => Right(DMat.tabulate(n, n)((i, j) => 0.5 * (c(i, j) + c(j, i))))

  /** Back-transform the selected standard eigenvectors `y` to generalized
    * eigenvectors `x = L⁻ᵀ y` (solve `Lᵀ X = Y`), one selected column each.
    */
  private def backTransformVectors(
      l: DMat,
      kernelVectors: DMat,
      indices: Array[Int]
  ): Either[LinAlgError, DMat] =
    solveTriColumns(l.t, lower = false, indices.length, c => kernelVectors.col(indices(c)))

  /** Build the generalized result: selected eigenvalues (ascending), the
    * back-transformed `B`-orthonormal eigenvectors (or an empty matrix when
    * values-only), true generalized residuals `‖A x − λ B x‖`, and the
    * `B`-inner-product orthogonality error `‖Xᵀ B X − I‖_F`.
    */
  private def assembleGeneralized(
      aSym: DMat,
      bSym: DMat,
      fullValues: DVec,
      indices: Array[Int],
      vectors: Option[DMat]
  ): EigenDecomposition =
    val m = indices.length
    val selValues = DVec.tabulate(m)(i => fullValues(indices(i)))
    val n = aSym.rows
    val (x, residuals, orthoErr) =
      vectors match
        case Some(mat) =>
          val res = DVec.tabulate(m): c =>
            val xc = mat.col(c)
            ((aSym * xc) - (bSym * xc) * selValues(c)).norm2
          (mat, res, bOrthogonalityError(mat, bSym))
        case None =>
          (DMat.zeros(n, 0), DVec.zeros(m), 0.0)
    val diagnostics =
      SpectralDiagnostics(
        requested = m,
        converged = m,
        residuals = residuals,
        orthogonalityError = orthoErr,
        iterations = 0,
        rank = None
      )
    EigenDecomposition(selValues, x, diagnostics)

  /** `‖Xᵀ B X − I‖_F` — the `B`-inner-product orthogonality error of the
    * generalized eigenvectors (`B`-orthonormal ⇒ this is ~0).
    */
  private def bOrthogonalityError(x: DMat, bSym: DMat): Double =
    val g = x.t * (bSym * x)
    var sum = 0.0
    var i = 0
    while i < g.rows do
      var j = 0
      while j < g.cols do
        val d = if i == j then g(i, j) - 1.0 else g(i, j)
        sum += d * d
        j += 1
      i += 1
    math.sqrt(sum)

  // ===========================================================================
  // Generalized nonsymmetric (QZ) — canonicalization + assembly (facade-owned)
  // ===========================================================================

  private def validateSquarePencil(a: DMat, b: DMat): Either[LinAlgError, Int] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else if b.rows != b.cols then Left(LinAlgError.NonSquareMatrix(b.shape))
    else if a.rows != b.rows then Left(LinAlgError.DimensionMismatch(a.shape, b.shape))
    else Right(a.rows)

  /** The base indices of the projective spectrum's units in packed order (a real /
    * infinite singleton, or a finite conjugate pair, positive-imaginary first),
    * walking `alphaIm` once.
    */
  private def generalizedUnitBases(alphaIm: DVec): Array[Int] =
    val bases = scala.collection.mutable.ArrayBuffer.empty[Int]
    var i = 0
    while i < alphaIm.length do
      bases += i
      if alphaIm(i) > 0.0 then i += 2 else i += 1
    bases.toArray

  /** Canonical index order for a projective `(α, β)` spectrum (§ 1.3, D-b): sort
    * units by the criterion on `α/β` via '''cross-multiplication''' (never dividing
    * by `β`, so `β = 0` yields no `∞`/`NaN`); infinite eigenvalues (`β = 0`) are
    * partitioned to the criterion's extreme (largest magnitude ⇒ first) and
    * tie-broken among themselves by descending `|α|` then producer index; conjugate
    * pairs stay adjacent. Used by the facade with [[EigenOrder.LargestMagnitude]]
    * (the canonical layout); the other legal orders are supported for completeness.
    */
  private[spectral] def generalizedIndices(alphaRe: DVec, alphaIm: DVec, beta: DVec, order: EigenOrder): Array[Int] =
    val bases = generalizedUnitBases(alphaIm)
    val sorted = bases.sortWith((i, j) => generalizedUnitCompare(alphaRe, alphaIm, beta, i, j, order) < 0)
    val out = scala.collection.mutable.ArrayBuffer.empty[Int]
    var u = 0
    while u < sorted.length do
      val base = sorted(u)
      out += base
      if alphaIm(base) > 0.0 then out += (base + 1)
      u += 1
    out.toArray

  /** Total order on unit base indices `i`, `j` (negative ⇒ `i` before `j`):
    * partition by infinite class, then compare finite units by the criterion via
    * cross-multiplication, then deterministic tie-breaks.
    */
  private def generalizedUnitCompare(
      alphaRe: DVec,
      alphaIm: DVec,
      beta: DVec,
      i: Int,
      j: Int,
      order: EigenOrder
  ): Int =
    val classI = infiniteClass(alphaRe(i), beta(i), order)
    val classJ = infiniteClass(alphaRe(j), beta(j), order)
    if classI != classJ then Integer.compare(classI, classJ)
    else if classI != 0 then
      // Both infinite in the same bucket: descending |α|, then producer index.
      val c = java.lang.Double.compare(math.hypot(alphaRe(j), alphaIm(j)), math.hypot(alphaRe(i), alphaIm(i)))
      if c != 0 then c else Integer.compare(i, j)
    else
      val c = finiteCriterionCompare(alphaRe, alphaIm, beta, i, j, order)
      if c != 0 then c
      else
        // Tie-break: descending real part (αRe/β), then descending imag, then index.
        val re = java.lang.Double.compare(alphaRe(j) * beta(i), alphaRe(i) * beta(j))
        if re != 0 then re
        else
          val im = java.lang.Double.compare(alphaIm(j) * beta(i), alphaIm(i) * beta(j))
          if im != 0 then im else Integer.compare(i, j)

  /** Ordering bucket for a unit: `0` finite, `-1` infinite-at-front, `+1`
    * infinite-at-back — placing `β = 0` at the criterion's extreme without division.
    */
  private def infiniteClass(alphaReVal: Double, betaVal: Double, order: EigenOrder): Int =
    if betaVal != 0.0 then 0
    else
      order match
        case EigenOrder.LargestMagnitude  => -1
        case EigenOrder.SmallestMagnitude => 1
        case EigenOrder.LargestRealPart   => if alphaReVal >= 0.0 then -1 else 1
        case EigenOrder.SmallestRealPart  => if alphaReVal >= 0.0 then 1 else -1
        case _                            => -1

  /** Compare two '''finite''' units by the criterion on `α/β` via cross-
    * multiplication (`β > 0` here, so products keep the ratio order). Negative ⇒
    * `i` before `j`.
    */
  private def finiteCriterionCompare(
      alphaRe: DVec,
      alphaIm: DVec,
      beta: DVec,
      i: Int,
      j: Int,
      order: EigenOrder
  ): Int =
    order match
      case EigenOrder.LargestMagnitude =>
        java.lang.Double.compare(
          math.hypot(alphaRe(j), alphaIm(j)) * beta(i),
          math.hypot(alphaRe(i), alphaIm(i)) * beta(j)
        )
      case EigenOrder.SmallestMagnitude =>
        java.lang.Double.compare(
          math.hypot(alphaRe(i), alphaIm(i)) * beta(j),
          math.hypot(alphaRe(j), alphaIm(j)) * beta(i)
        )
      case EigenOrder.LargestRealPart =>
        java.lang.Double.compare(alphaRe(j) * beta(i), alphaRe(i) * beta(j))
      case EigenOrder.SmallestRealPart =>
        java.lang.Double.compare(alphaRe(i) * beta(j), alphaRe(j) * beta(i))
      case _ => 0

  /** Permute the columns of a packed vector matrix by `outIdx` (empty stays empty). */
  private def permuteQzCols(m: DMat, outIdx: Array[Int]): DMat =
    if m.cols == 0 then DMat.zeros(m.rows, 0)
    else DMat.tabulate(m.rows, outIdx.length)((r, c) => m(r, outIdx(c)))

  /** Assemble the sealed [[GeneralizedEigenDecomposition]] from a backend's raw QZ
    * output: guard the `β ≥ 0` backend contract, impose the canonical projective
    * order (moving each conjugate pair's two columns in lockstep), re-derive the
    * homogeneous residuals against `a`, `b` — including the '''left''' residual when
    * left vectors were returned (§ 1.5 honesty split) — and report a dense one-shot
    * diagnostics (rank = number of finite eigenvalues).
    */
  private def assembleQz(a: DMat, b: DMat, raw: RawGeneralizedEigen): Either[LinAlgError, GeneralizedEigenDecomposition] =
    // Guard the β ≥ 0 backend contract before sorting: a negative β would invert the
    // cross-multiplied criterion comparisons, making the sorter non-transitive (a
    // TimSort "violates its general contract" crash). Reject loudly instead.
    var bi = 0
    var negativeBeta = false
    while bi < raw.beta.length do
      if raw.beta(bi) < 0.0 then negativeBeta = true
      bi += 1
    if negativeBeta then
      Left(
        LinAlgError.InvalidArgument(
          "QZ backend contract violation: β must be ≥ 0 (LAPACK ggev convention); a negative β was returned"
        )
      )
    else
      val outIdx = generalizedIndices(raw.alphaRe, raw.alphaIm, raw.beta, EigenOrder.LargestMagnitude)
      val cnt = outIdx.length
      val aRe = DVec.tabulate(cnt)(j => raw.alphaRe(outIdx(j)))
      val aIm = DVec.tabulate(cnt)(j => raw.alphaIm(outIdx(j)))
      val bt = DVec.tabulate(cnt)(j => raw.beta(outIdx(j)))
      val rightSel = permuteQzCols(raw.rightPacked, outIdx)
      val leftSel = raw.leftPacked.map(m => permuteQzCols(m, outIdx))
      // Right residual (zeros when no right vectors); left residual (zeros when no
      // left vectors). Per eigenvalue the reported residual is the worse of the two,
      // so a left-only or left-and-right result reports an honest left residual
      // rather than a deceptive zero.
      val rightRes = generalizedQzResiduals(a, b, aRe, aIm, bt, rightSel)
      val leftRes = leftSel match
        case Some(l) => generalizedQzLeftResiduals(a.t, b.t, aRe, aIm, bt, l)
        case None    => DVec.zeros(cnt)
      val residuals = DVec.tabulate(cnt)(i => math.max(rightRes(i), leftRes(i)))
      var finite = 0
      var i = 0
      while i < cnt do
        if bt(i) != 0.0 then finite += 1
        i += 1
      val diagnostics =
        SpectralDiagnostics(
          requested = cnt,
          converged = cnt,
          residuals = residuals,
          orthogonalityError = 0.0,
          iterations = 0,
          rank = Some(finite)
        )
      Right(new GeneralizedEigenDecomposition(aRe, aIm, bt, rightSel, leftSel, diagnostics))

  /** Per-eigenvalue '''homogeneous''' residual `‖β A x − α B x‖` in real arithmetic
    * (finite ⇒ the pencil residual `A x = (α/β) B x`; infinite ⇒ `|α|·‖B x‖`,
    * testing `x ∈ null(B)`). Zeros when no vectors are present.
    */
  private def generalizedQzResiduals(a: DMat, b: DMat, aRe: DVec, aIm: DVec, bt: DVec, packed: DMat): DVec =
    val cnt = aRe.length
    if packed.cols == 0 then DVec.zeros(cnt)
    else
      val out = new Array[Double](cnt)
      var j = 0
      while j < cnt do
        val im = aIm(j)
        if im == 0.0 then
          val x = packed.col(j)
          out(j) = ((a * x) * bt(j) - (b * x) * aRe(j)).norm2
          j += 1
        else
          val xr = packed.col(j)
          val xi = packed.col(j + 1)
          val ar = aRe(j)
          val be = bt(j)
          val realPart = (a * xr) * be - (b * xr) * ar + (b * xi) * im
          val imagPart = (a * xi) * be - (b * xi) * ar - (b * xr) * im
          val res = math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))
          out(j) = res
          out(j + 1) = res
          j += 2
      DVec.fromSeq(out.toIndexedSeq)

  /** Per-eigenvalue homogeneous '''left''' residual `‖β wᴴ A − α wᴴ B‖` in real
    * arithmetic, using `at = Aᵀ`, `bt = Bᵀ` (finite ⇒ `wᴴ A = (α/β) wᴴ B`; infinite
    * ⇒ `|α|·‖wᴴ B‖`, testing `w` in the left null space of `B`). The transposed
    * residual is `[β Aᵀw_re − αRe Bᵀw_re − αIm Bᵀw_im] +
    * i[−β Aᵀw_im − αIm Bᵀw_re + αRe Bᵀw_im]`. Zeros when no left vectors.
    */
  private def generalizedQzLeftResiduals(at: DMat, bt: DMat, aRe: DVec, aIm: DVec, beta: DVec, packed: DMat): DVec =
    val cnt = aRe.length
    if packed.cols == 0 then DVec.zeros(cnt)
    else
      val out = new Array[Double](cnt)
      var j = 0
      while j < cnt do
        val im = aIm(j)
        if im == 0.0 then
          val w = packed.col(j)
          out(j) = ((at * w) * beta(j) - (bt * w) * aRe(j)).norm2
          j += 1
        else
          val wr = packed.col(j)
          val wi = packed.col(j + 1)
          val ar = aRe(j)
          val be = beta(j)
          val realPart = (at * wr) * be - (bt * wr) * ar - (bt * wi) * im
          val imagPart = (at * wi) * (-be) - (bt * wr) * im + (bt * wi) * ar
          val res = math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))
          out(j) = res
          out(j + 1) = res
          j += 2
      DVec.fromSeq(out.toIndexedSeq)
