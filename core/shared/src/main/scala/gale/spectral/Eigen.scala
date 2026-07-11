package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import gale.linalg.Rows
import gale.linalg.Shape

/** Public symmetric eigendecomposition entry points (phase a of
  * `docs/spectral-parity.md`).
  *
  * Both entry points share the same result type ([[EigenDecomposition]]) and the
  * same ordering guarantee: '''eigenvalues are returned ascending-algebraic
  * always''' (§ 1, § 6). A selection's `order` decides only *which* eigenvalues
  * are members of the result, never their layout.
  *
  *   - [[eigSymmetric(a:gale\.linalg\.DMat*]] — dense full solve via the
  *     tridiagonal QL/QR kernel, with the requested selection realized as a slice
  *     of the full ascending spectrum (the phase-a scope decision in § 1). Reads
  *     only the lower triangle of `a`.
  *   - [[eigSymmetric(op:gale\.linalg\.DoubleLinearOperator*]] — iterative partial
  *     solve (`eigsh`) via Lanczos with full reorthogonalization, for a matrix or
  *     a matrix-free [[gale.linalg.DoubleLinearOperator]].
  *
  * '''Failure model''' (§ Convergence & failure semantics). Structural /
  * precondition violations are `Left(LinAlgError)`. For the '''dense''' path,
  * kernel sweep exhaustion is also a `Left(DidNotConverge)` — there is no partial
  * dense result to hand back (in practice the QL/QR kernel always converges within
  * its budget for a well-formed symmetric matrix). For the '''iterative''' path,
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
