package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import gale.linalg.Rows
import gale.linalg.Shape

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
  * first) and '''never split''' during selection. Only right eigenvectors are
  * computed (left vectors deferred).
  *
  *   - [[eigNonsymmetric(a:gale\.linalg\.DMat*]] — dense full solve via Hessenberg
  *     reduction + Francis double-shift QR, the selection realized as a permutation
  *     of the full spectrum.
  *   - [[eigNonsymmetric(op:gale\.linalg\.DoubleLinearOperator*]] — iterative
  *     partial solve (`eigs`) via Arnoldi with full reorthogonalization.
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
    * `vectors` chooses [[EigenVectors.ValuesOnly]] versus [[EigenVectors.Right]].
    * '''Left and left-and-right eigenvectors are deferred''' — the Francis QR
    * kernel computes right eigenvectors only — so [[EigenVectors.Left]] /
    * [[EigenVectors.LeftAndRight]] return `Left(UnsupportedOperation)` (the parity
    * doc lists left vectors in-b, but the kernel lacks them; deferring is the
    * honest scope). Orthogonality error is reported as `0.0`: nonsymmetric
    * eigenvectors are not orthonormal, so `‖VᵀV − I‖` is not a meaningful signal.
    *
    * `Left` on: non-square `a`; an [[EigenOrder]] illegal for a nonsymmetric
    * problem ([[EigenOrder.LargestAlgebraic]]/[[EigenOrder.SmallestAlgebraic]]/
    * [[EigenOrder.BothEnds]] are symmetric-only); a `k` outside `[1, n]`;
    * [[EigenSelection.IndexRange]] or [[EigenSelection.ValueInterval]]
    * (symmetric-only); [[EigenVectors.Left]]/[[EigenVectors.LeftAndRight]]; or (in
    * practice unreachable) kernel non-convergence.
    */
  def eigNonsymmetric(
      a: DMat,
      selection: EigenSelection,
      vectors: EigenVectors
  ): Either[LinAlgError, NonsymmetricEigenDecomposition] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else
      val n = a.rows
      validateNonsymmetricVectors(vectors) match
        case Left(error) => Left(error)
        case Right(wantVectors) =>
          validateNonsymmetricDenseSelection(selection, n) match
            case Left(error) => Left(error)
            case Right(()) =>
              DenseSpectralKernels.nonsymmetricEigen(a, wantVectors) match
                case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
                  // Dense path has no partial result to return, so exhaustion is a
                  // Left; residual is not measured on failure (reported as 0.0).
                  Left(LinAlgError.DidNotConverge(iters, 0.0))
                case Right(kernel) =>
                  val outIdx = nonsymDenseIndices(kernel.re, kernel.im, selection)
                  Right(assembleNonsymDense(a, kernel, outIdx, wantVectors))

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
      validateNonsymmetricVectors(options.returnVectors) match
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

  /** Map the nonsymmetric vector flag to "compute right vectors?". Left / and
    * left-and-right are '''deferred''' (the Francis QR kernel returns right vectors
    * only) — `Left(UnsupportedOperation)`, unlike the symmetric path where they are
    * meaningless (`InvalidArgument`).
    */
  private def validateNonsymmetricVectors(vectors: EigenVectors): Either[LinAlgError, Boolean] =
    vectors match
      case EigenVectors.ValuesOnly => Right(false)
      case EigenVectors.Right      => Right(true)
      case EigenVectors.Left | EigenVectors.LeftAndRight =>
        Left(
          LinAlgError.UnsupportedOperation(
            "left eigenvectors (the Francis QR kernel computes right eigenvectors only; left vectors are deferred)"
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

  /** Build the dense nonsymmetric result: permute the kernel's re/im and packed
    * vector columns by `outIdx`, compute per-eigenvalue residuals against `a`, and
    * report all-converged diagnostics (a dense one-shot solve). Orthogonality error
    * is `0.0` — nonsymmetric eigenvectors are not orthonormal, so that quality
    * metric does not apply.
    */
  private def assembleNonsymDense(
      a: DMat,
      kernel: DenseSpectralKernels.NonsymmetricEigen,
      outIdx: Array[Int],
      wantVectors: Boolean
  ): NonsymmetricEigenDecomposition =
    val cnt = outIdx.length
    val re = DVec.tabulate(cnt)(j => kernel.re(outIdx(j)))
    val im = DVec.tabulate(cnt)(j => kernel.im(outIdx(j)))
    val n = a.rows
    val (vectors, residuals) =
      if wantVectors then
        val packed = kernel.vectors.get
        val sel = DMat.tabulate(n, cnt)((r, c) => packed(r, outIdx(c)))
        (sel, nonsymResiduals(v => a * v, re, im, sel))
      else (DMat.zeros(n, 0), DVec.zeros(cnt))
    val diagnostics =
      SpectralDiagnostics(
        requested = cnt,
        converged = cnt,
        residuals = residuals,
        orthogonalityError = 0.0,
        iterations = 0,
        rank = None
      )
    new NonsymmetricEigenDecomposition(re, im, vectors, None, diagnostics)

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
