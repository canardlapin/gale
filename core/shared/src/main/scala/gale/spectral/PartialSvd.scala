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

/** Public singular value decomposition — partial (`svds`, phase a of
  * `docs/spectral-parity.md` § 3 / § 8) and, on the dense entry points, the
  * '''full/economy''' SVD that phase a deferred.
  *
  * '''Full dense SVD''' ([[SingularSelection.All]], or `Count(k)` with
  * `k = min(m, n)`, on the `DMat` overloads). Served by the
  * [[DenseSvdKernel]] — Householder bidiagonalization plus Golub–Kahan–Reinsch
  * implicit-shift QR with the Wilkinson shift — in '''economy''' shapes: `U`
  * is `m×k`, `Vᵀ` is `k×n`, `k = min(m, n)`. This closes seam S7 of
  * `docs/spectral-backend-boundary.md`: the facade takes a
  * `using SpectralBackend` and routes to a
  * [[SpectralCapability.DenseSvd]]-capable backend when one is imported; the
  * pure kernel is the always-available default (so, unlike the other
  * backend-scoped seams, no-import is a `Right`, not a `Left`). The operator
  * (matrix-free) overload stays Count-only — a dense bidiagonal reduction
  * needs the materialized matrix.
  *
  * '''Partial SVD.''' The `k < min(m, n)` requested
  * triplets are computed by '''Golub–Kahan–Lanczos bidiagonalization''' of `A`
  * (two-sided, with full reorthogonalization of both the left and right vector
  * sequences), whose small bidiagonal `B` is diagonalized by running the existing
  * symmetric tridiagonal QL/QR kernel on the '''Jordan–Wielandt augmented
  * tridiagonal''' — the perfect-shuffle of `[[0, Bᵀ], [B, 0]]`, whose eigenvalues
  * are `±σ_i`. That augmented form (as opposed to the normal-equation `BᵀB`, whose
  * `σ²` squares the condition number) is what keeps the '''smallest''' singular
  * values usable. The generic symmetric QL/QR on the augmented form gives
  * '''absolute''' accuracy `~ε‖A‖`, not the high relative accuracy of a dedicated
  * bidiagonal kernel — acceptable for gale's portable-correctness scope, and
  * called out here so callers know the accuracy floor for tiny `σ`.
  *
  * '''Ordering.''' Singular values are returned '''descending always''' (§ 8),
  * regardless of whether [[SingularOrder.Largest]] or [[SingularOrder.Smallest]]
  * chose membership — resolving the SciPy-`svds`-ascending trap.
  *
  * '''Failure model''' (§ Convergence & failure semantics). Structural violations
  * are `Left`; iterative non-convergence is a `Right` carrying only the converged
  * triplets plus `SpectralDiagnostics`, never a `Left` — with one carve-out
  * mirroring [[Eigen]]: exhaustion of the '''inner''' augmented-tridiagonal QL/QR
  * kernel is a genuine numerical breakdown with no partial result and surfaces as
  * `Left(DidNotConverge)` (pathological in practice). `requireConverged` is the
  * caller's opt-in to fail-fast. Rank deficiency surfaces through the `rank`
  * field / diagnostic and through near-zero singular values in the output, not an
  * error.
  *
  * The two entry points mirror [[Eigen]]: a dense [[svd(a:gale\.linalg\.DMat*]]
  * convenience and a matrix-free
  * [[svd(op:gale\.linalg\.DoubleLinearOperator*]] operator form. The operator
  * '''must''' implement both `applyTo` (`A·v`) and `transposeApplyTo` (`Aᵀ·u`) —
  * a `DMat` does; a one-sided matrix-free operator does not and is unsupported:
  * its default `transposeApplyTo` '''throws''' `UnsupportedOperation` from the
  * primitive layer (the library's primitives-throw convention), it does not
  * return a `Left`.
  *
  * (Named `Svds` — the ecosystem's `svds` — rather than `Svd`, which would clash
  * with the [[SVD]] result type on a case-insensitive filesystem.)
  */
object Svds:

  // ===========================================================================
  // Dense SVD (full via the bidiagonal kernel; partial via Golub–Kahan–Lanczos)
  // ===========================================================================

  /** SVD of the dense `a`, computing singular vectors ([[EigenVectors.Right]]).
    * See the three-argument overload for the vector flag and failure details.
    */
  def svd(a: DMat, selection: SingularSelection)(using SpectralBackend): Either[LinAlgError, SVD] =
    svd(a, selection, EigenVectors.Right)

  /** SVD of the dense `a` — `A = U Σ Vᵀ` with the singular triplets named by
    * `selection`, `Σ` descending always. [[SingularSelection.All]] (and
    * equivalently `Count(k)` with `k = min(m, n)`, whichever order — membership
    * is the whole spectrum either way) computes the '''full/economy''' dense
    * SVD through the bidiagonal kernel: `U` `m×k`, `Vᵀ` `k×n`,
    * `k = min(m, n)`, routed to a [[SpectralCapability.DenseSvd]]-capable
    * `backend` when one is in scope (seam S7) and to the pure
    * [[DenseSvdKernel]] otherwise. A `Count` with `k < min(m, n)` runs the
    * partial Golub–Kahan–Lanczos path. `vectors` chooses
    * [[EigenVectors.ValuesOnly]] versus [[EigenVectors.Right]] (there is no
    * one-sided singular-vector mode, § 8: vectors means '''both''' `U` and
    * `V`, or neither; [[EigenVectors.Left]]/[[EigenVectors.LeftAndRight]] are
    * rejected). Uses the default [[SpectralOptions]]; for tolerance/subspace
    * control of the partial path drive the operator overload with `a` as the
    * operator.
    *
    * On the full path `rank` (and `diagnostics.rank`) counts the singular
    * values above the default solve tolerance (`1e-10`) relative to `σ_max`,
    * the same returned-set policy the [[SVD]] result type documents;
    * `diagnostics.iterations` is `0` (a dense one-shot solve).
    *
    * `Left` on: a non-positive dimension; an illegal vector flag; `k ≤ 0` or
    * `k > min(m, n)`; or (in practice unreachable) bidiagonal-QR sweep
    * exhaustion, `Left(DidNotConverge)` like the dense eigen paths.
    */
  def svd(a: DMat, selection: SingularSelection, vectors: EigenVectors)(using
      backend: SpectralBackend
  ): Either[LinAlgError, SVD] =
    val m = a.rows
    val n = a.cols
    if m <= 0 || n <= 0 then Left(LinAlgError.InvalidArgument(s"dimensions must be positive, got ${m}x$n"))
    else
      validateVectors(vectors) match
        case Left(error) => Left(error)
        case Right(wantVectors) =>
          val p = math.min(m, n)
          selection match
            case SingularSelection.All =>
              svdFullDense(a, wantVectors)
            case SingularSelection.Count(k, _) if k == p =>
              svdFullDense(a, wantVectors)
            case SingularSelection.Count(k, _) if k <= 0 || k > p =>
              Left(LinAlgError.InvalidArgument(s"k=$k must be in [1, $p]"))
            case _ =>
              svdPartial(a, m, n, selection, SpectralOptions(returnVectors = vectors))

  // ===========================================================================
  // Operator (matrix-free) partial SVD
  // ===========================================================================

  /** Partial SVD of the `rows × cols` operator `op`, which must apply both `A·v`
    * ([[gale.linalg.DoubleLinearOperator.applyTo]]) and `Aᵀ·u`
    * ([[gale.linalg.DoubleLinearOperator.transposeApplyTo]]). Rectangular is
    * legal (`rows ≠ cols`), so a shape disagreement with the operator is a
    * `DimensionMismatch`, never `NonSquareMatrix`.
    *
    * Selection is Count-only (§ 8): [[SingularSelection.Count]] with
    * [[SingularOrder.Largest]] or [[SingularOrder.Smallest]]; membership by the
    * order, layout always descending. `options.returnVectors` selects
    * values-only versus both `U` and `V`. A provided `options.startVector` seeds
    * the bidiagonalization in `R^min(m,n)` (the taller-oriented start space), so
    * its length must be `min(m, n)`. Convergence and failure semantics are as
    * documented on the object.
    */
  def svd(
      op: DoubleLinearOperator,
      rows: Int,
      cols: Int,
      selection: SingularSelection,
      options: SpectralOptions = SpectralOptions()
  ): Either[LinAlgError, SVD] =
    svdPartial(op, rows, cols, selection, options)

  // ===========================================================================
  // Full/economy dense SVD (seam S7) + pseudo-inverse
  // ===========================================================================

  /** Moore–Penrose pseudo-inverse of `a` (`n×m` for an `m×n` input) via the
    * full/economy dense SVD: `A⁺ = V Σ⁺ Uᵀ` with `Σ⁺` inverting exactly the
    * singular values above the cutoff and zeroing the rest.
    *
    * '''Cutoff convention''' (NumPy's `pinv` default semantics): a singular
    * value is treated as zero unless `σ_i > rcond · σ_max` with
    * `rcond = max(m, n) · ε` (`ε = 2⁻⁵² ≈ 2.22e-16`) — i.e. the cutoff is
    * `max(m, n) · ε · σ_max`, scale-aware per `docs/numerical-contract.md`. An
    * all-zero `a` therefore yields the all-zero `A⁺` (its correct
    * pseudo-inverse), never a division by zero.
    *
    * `Left` exactly when the underlying [[svd(a:gale\.linalg\.DMat*]] is: a
    * non-positive dimension, or (in practice unreachable) kernel
    * non-convergence.
    */
  def pinv(a: DMat)(using SpectralBackend): Either[LinAlgError, DMat] =
    svd(a, SingularSelection.All, EigenVectors.Right).map: s =>
      val m = a.rows
      val n = a.cols
      val p = s.size
      val sigmaMax = if p > 0 then s.singularValues(0) else 0.0
      val cutoff = math.max(m, n).toDouble * MachineEpsilon * sigmaMax
      // W = V·Σ⁺ (n×p): W(i, l) = Vᵀ(l, i) / σ_l above the cutoff, else 0.
      val w = DMat.tabulate(n, p): (i, l) =>
        val sigma = s.singularValues(l)
        if sigma > cutoff then s.vt(l, i) / sigma else 0.0
      w * s.u.t

  /** Route the full dense SVD (already validated: positive dims, legal vector
    * flag): a [[SpectralCapability.DenseSvd]]-capable backend computes the raw
    * factors, else the pure [[DenseSvdKernel]] does; either way the facade
    * canonicalizes and assembles ([[assembleFullSvd]]). Kernel sweep exhaustion
    * maps to `Left(DidNotConverge)` exactly like the dense eigen paths (no
    * partial result to hand back).
    */
  private def svdFullDense(a: DMat, wantVectors: Boolean)(using backend: SpectralBackend): Either[LinAlgError, SVD] =
    val raw: Either[LinAlgError, RawSvd] =
      if backend.capabilities.contains(SpectralCapability.DenseSvd) then backend.denseSvd(a, wantVectors)
      else
        DenseSvdKernel.svd(a, wantVectors) match
          case Left(DenseSvdKernel.SvdKernelFailure.DidNotConverge(iters)) =>
            // Dense path has no partial result to return, so exhaustion is a
            // Left; residual is not measured on failure (reported as 0.0).
            Left(LinAlgError.DidNotConverge(iters, 0.0))
          case Right(factors) => Right(factors)
    raw.map(assembleFullSvd(a, _, wantVectors))

  /** Canonicalize raw full-SVD factors into the sealed [[SVD]]: enforce
    * `σ ≥ 0` (a negative raw value flips with its `Vᵀ` row — legal because
    * `σ u vᵀ = (−σ) u (−v)ᵀ`), impose the '''descending''' order (§ 8, ties by
    * raw index), and re-derive what is checkable — per-triplet two-sided
    * residuals `max(‖A v − σ u‖, ‖Aᵀ u − σ v‖)`, the worse `U`/`V`
    * orthogonality error, and `rank` at the default solve tolerance relative
    * to `σ_max` (the returned-set policy [[SVD]] documents). Values-only
    * results carry empty factors, zero residuals, and zero orthogonality
    * error, mirroring [[Eigen]]'s dense assembly.
    */
  private def assembleFullSvd(a: DMat, raw: RawSvd, wantVectors: Boolean): SVD =
    val m = a.rows
    val n = a.cols
    val p = raw.sigma.length
    val sig = new Array[Double](p)
    val flip = new Array[Boolean](p)
    var i = 0
    while i < p do
      val s = raw.sigma(i)
      if s < 0.0 then
        sig(i) = -s
        flip(i) = true
      else sig(i) = s
      i += 1
    val order = sig.indices.sortBy(i => (-sig(i), i)).toArray
    val values = DVec.tabulate(p)(i => sig(order(i)))
    val (u, vt, residuals, orthoErr) =
      if wantVectors && raw.u.cols == p then
        val uMat = DMat.tabulate(m, p)((r, c) => raw.u(r, order(c)))
        val vtMat = DMat.tabulate(p, n): (r, c) =>
          val src = order(r)
          if flip(src) then -raw.vt(src, c) else raw.vt(src, c)
        val res = DVec.tabulate(p): c =>
          val uCol = uMat.col(c)
          val vRow = vtMat.row(c)
          val sigma = values(c)
          val rV = (a * vRow - uCol * sigma).norm2
          val rU = (a.t * uCol - vRow * sigma).norm2
          math.max(rV, rU)
        (uMat, vtMat, res, math.max(orthogonalityError(uMat), orthogonalityError(vtMat.t)))
      else (DMat.zeros(m, 0), DMat.zeros(0, n), DVec.zeros(p), 0.0)
    val tol = SpectralOptions().tolerance
    val sigmaMax = if p > 0 then values(0) else 0.0
    var rank = 0
    i = 0
    while i < p do
      if values(i) > tol * sigmaMax then rank += 1
      i += 1
    val diagnostics =
      SpectralDiagnostics(
        requested = p,
        converged = p,
        residuals = residuals,
        orthogonalityError = orthoErr,
        iterations = 0,
        rank = Some(rank)
      )
    SVD(values, u, vt, rank, diagnostics)

  /** IEEE machine epsilon for `Double` (2^-52), the `pinv` cutoff scale. */
  private val MachineEpsilon: Double = 2.220446049250313e-16

  // ===========================================================================
  // Core
  // ===========================================================================

  private def svdPartial(
      op: DoubleLinearOperator,
      m: Int,
      n: Int,
      selection: SingularSelection,
      options: SpectralOptions
  ): Either[LinAlgError, SVD] =
    if m <= 0 || n <= 0 then Left(LinAlgError.InvalidArgument(s"dimensions must be positive, got ${m}x$n"))
    else if op.rows != m || op.cols != n then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(m), Cols(n)), Shape(Rows(op.rows), Cols(op.cols))))
    else
      selection match
        case SingularSelection.All =>
          Left(
            LinAlgError.InvalidArgument(
              "the operator (matrix-free) SVD selects by Count(k, order); use the dense matrix overload (SingularSelection.All) for the full SVD"
            )
          )
        case SingularSelection.Count(k, order) =>
          validateVectors(options.returnVectors) match
            case Left(error) => Left(error)
            case Right(wantVectors) =>
              val p = math.min(m, n)
              if k <= 0 || k >= p then
                Left(
                  LinAlgError.InvalidArgument(
                    s"k=$k must be in [1, ${p - 1}] for partial SVD; use the dense matrix overload for the full SVD (k = min(m,n))"
                  )
                )
              else
                // Orient the bidiagonalization on the taller side: start in the
                // smaller space R^p so no Krylov dimension is wasted on the
                // |m−n|-dimensional null space (min(m,n) steps would otherwise miss
                // singular values). For m < n run on Aᵀ and swap left/right back.
                val transposed = m < n
                val effOp = if transposed then transposeOp(op, m, n) else op
                val effRows = if transposed then n else m
                startVectorFor(options.startVector, p) match
                  case Left(error) => Left(error)
                  case Right(v0) =>
                    runGolubKahan(op, effOp, m, n, effRows, p, transposed, k, order, options, wantVectors, v0)

  /** Map the SVD vector flag to "compute U and V?". Left/LeftAndRight have no
    * meaning for singular vectors (no left-vs-right eigenvector split) and there
    * is no one-sided mode in v0.3.5 (§ 8), so both are rejected.
    */
  private def validateVectors(vectors: EigenVectors): Either[LinAlgError, Boolean] =
    vectors match
      case EigenVectors.ValuesOnly => Right(false)
      case EigenVectors.Right      => Right(true)
      case EigenVectors.Left | EigenVectors.LeftAndRight =>
        Left(
          LinAlgError.InvalidArgument(
            "SVD returns both U and V or neither; use EigenVectors.Right or EigenVectors.ValuesOnly"
          )
        )

  private def runGolubKahan(
      op: DoubleLinearOperator,
      effOp: DoubleLinearOperator,
      m: Int,
      n: Int,
      effRows: Int,
      effCols: Int,
      transposed: Boolean,
      k: Int,
      order: SingularOrder,
      options: SpectralOptions,
      wantVectors: Boolean,
      v0: DVec
  ): Either[LinAlgError, SVD] =
    val p = effCols // = min(m, n)
    val ncv0 = math.min(p, math.max(options.subspaceDimension.getOrElse(math.max(2 * k + 1, 20)), k + 1))
    val maxSteps = math.max(options.maxIterations, 1)
    val tol = options.tolerance
    val growBy = math.max(k, 16)

    var mBuild = ncv0
    var step = 0
    var done = false
    var failure: Option[LinAlgError] = None
    // Latest build's result. The bidiagonalization grows from a fixed start
    // vector and is exact at mBuild = p; per-pair residuals are not strictly
    // monotone under a maxIterations truncation, so the final build is reported.
    var result: SVD = assembleSvd(m, n, k, Array.empty, Array.empty, Array.empty, Array.empty, wantVectors, tol, 0)

    while step < maxSteps && !done && failure.isEmpty do
      val (uB, vB, alpha, beta, mEff, pendingRight) = buildGolubKahan(effOp, effRows, effCols, v0, mBuild)
      if mEff == 0 then done = true // effOp·v0 ≈ 0: nothing to resolve.
      else
        // Jordan–Wielandt augmented tridiagonal: zero diagonal, off-diagonal the
        // perfect-shuffle interleave [α_0, β_0, α_1, β_1, …]. With a pending right
        // vector (α-breakdown after β_{mEff-1} cleared the threshold) the exact
        // projected bidiagonal is RECTANGULAR mEff×(mEff+1); its augmented form is
        // odd-dimensional (2·mEff+1) and keeps the trailing β_{mEff-1}. Otherwise
        // the square mEff×mEff form (2·mEff) applies and any computed trailing β
        // is the Rayleigh–Ritz residual coupling, correctly excluded.
        val vCount = if pendingRight then mEff + 1 else mEff
        val augDim = mEff + vCount
        val diag = DVec.zeros(augDim)
        val off = DVec.tabulate(augDim - 1)(idx => if idx % 2 == 0 then alpha(idx / 2) else beta(idx / 2))
        DenseSpectralKernels.symmetricTridiagonalEigen(diag, off, wantVectors = true) match
          case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
            failure = Some(LinAlgError.DidNotConverge(iters, 0.0))
          case Right(eig) =>
            // The spectrum is ±σ (plus one exact zero in the odd-dimensional
            // case, which belongs to the null space and is never selected): the
            // top mEff eigenvalues, ascending, are the σ's.
            val eigvecs = eig.vectors.get
            val sigmaMax = math.max(1.0, eig.values(augDim - 1))
            val kk = math.min(k, mEff)
            val wantedIdx = order match
              case SingularOrder.Largest  => (mEff - kk until mEff).toArray
              case SingularOrder.Smallest => (0 until kk).toArray
            val sigmas = new Array[Double](wantedIdx.length)
            val uVecs = new Array[DVec](wantedIdx.length)
            val vVecs = new Array[DVec](wantedIdx.length)
            val residuals = new Array[Double](wantedIdx.length)
            var w = 0
            while w < wantedIdx.length do
              val col = vCount + wantedIdx(w)
              val sigma = eig.values(col)
              val (uEff, vEff) = ritzTriplet(eigvecs, col, uB, vB, effRows, effCols, mEff, vCount)
              // Map effOp's singular vectors back to A's: when we ran on Aᵀ, its
              // left vector (R^n) is A's right vector and vice versa.
              val (aU, aV) = if transposed then (vEff, uEff) else (uEff, vEff)
              sigmas(w) = sigma
              uVecs(w) = aU
              vVecs(w) = aV
              residuals(w) = svdResidual(op, aU, aV, sigma, m, n)
              w += 1
            val convergedSlots = residuals.indices.filter(i => residuals(i) <= tol * sigmaMax).toArray
            result = assembleSvd(
              m,
              n,
              k,
              convergedSlots.map(sigmas),
              convergedSlots.map(uVecs),
              convergedSlots.map(vVecs),
              convergedSlots.map(residuals),
              wantVectors,
              tol,
              step + 1
            )
            val fullyConverged = convergedSlots.length == wantedIdx.length && wantedIdx.length == k
            if fullyConverged || mEff < mBuild || mBuild >= p then done = true
            else mBuild = math.min(p, mBuild + growBy)
      step += 1

    failure match
      case Some(error) => Left(error)
      case None        => Right(result)

  /** Golub–Kahan–Lanczos bidiagonalization from the unit right start vector `v0`.
    *
    * Returns `(uBasis, vBasis, alpha, beta, mEff)` where `mEff ≤ ncv` is the size
    * of the computed upper-bidiagonal `B` (`mEff < ncv` marks a breakdown — an
    * exhausted Krylov space), `alpha(0..mEff-1)` is `B`'s diagonal and
    * `beta(0..mEff-2)` its superdiagonal (`beta(mEff-1)` is the trailing residual
    * factor). Both the left (`R^m`) and right (`R^n`) sequences are fully
    * reorthogonalized (classical Gram–Schmidt twice), the discipline that keeps
    * the bidiagonalization numerically orthogonal.
    */
  private def buildGolubKahan(
      op: DoubleLinearOperator,
      m: Int,
      n: Int,
      v0: DVec,
      ncv: Int
  ): (Array[DVec], Array[DVec], Array[Double], Array[Double], Int, Boolean) =
    val uB = new Array[DVec](ncv + 1)
    val vB = new Array[DVec](ncv + 1)
    val alpha = new Array[Double](ncv)
    val beta = new Array[Double](ncv)
    vB(0) = v0
    var aMax = 0.0
    var j = 0
    var mEff = 0
    var stop = false
    var pendingRight = false
    while j < ncv && !stop do
      // Left vector: u_j from A v_j, deflating the previous u and reorthogonalizing.
      val wLeft = MutableDVec.zeros(m)
      op.applyTo(vB(j), wLeft)
      if j > 0 then wLeft.axpyInPlace(-beta(j - 1), uB(j - 1))
      reorthogonalize(wLeft, uB, j)
      reorthogonalize(wLeft, uB, j)
      val aj = wLeft.asVec.norm2
      if aj <= 1e-12 * math.max(1.0, aMax) then
        // α_j ≈ 0: left Krylov space exhausted. When a right vector v_j is
        // pending with a significant trailing β_{j-1} (always the case for j > 0,
        // since v_j was only created when β_{j-1} cleared the threshold), the
        // exact projected problem is the RECTANGULAR j×(j+1) bidiagonal that
        // keeps that coupling — truncating to the square j×j silently solves the
        // wrong matrix (this is the exact-rank-deficiency path: a start vector
        // with a null-space component always ends here).
        stop = true
        pendingRight = j > 0
      else
        alpha(j) = aj
        if aj > aMax then aMax = aj
        uB(j) = (wLeft.asVec * (1.0 / aj)).copy
        mEff = j + 1
        // Right vector: v_{j+1} from Aᵀ u_j, deflating v_j and reorthogonalizing.
        val wRight = MutableDVec.zeros(n)
        op.transposeApplyTo(uB(j), wRight)
        wRight.axpyInPlace(-aj, vB(j))
        reorthogonalize(wRight, vB, j + 1)
        reorthogonalize(wRight, vB, j + 1)
        val bj = wRight.asVec.norm2
        beta(j) = bj
        if j < ncv - 1 && bj > 1e-12 * math.max(1.0, aMax) then
          vB(j + 1) = (wRight.asVec * (1.0 / bj)).copy
        else
          stop = true // β_j ≈ 0 (right exhausted) or the subspace is full.
      j += 1
    (uB, vB, alpha, beta, mEff, pendingRight)

  /** Decode a Jordan–Wielandt eigenvector (column `col` of the augmented
    * eigenproblem) into an original-space singular triplet: the perfect-shuffle
    * even components are the projected right vector (combine with `vB`, `R^n`),
    * the odd components the projected left vector (combine with `uB`, `R^m`). Both
    * are normalized to unit length; `A v = σ u` holds because the augmented
    * eigenvector's two halves have equal norm.
    */
  private def ritzTriplet(
      eigvecs: DMat,
      col: Int,
      uB: Array[DVec],
      vB: Array[DVec],
      m: Int,
      n: Int,
      mEff: Int,
      vCount: Int
  ): (DVec, DVec) =
    val uRitz = MutableDVec.zeros(m)
    val vRitz = MutableDVec.zeros(n)
    var i = 0
    while i < vCount do
      vRitz.axpyInPlace(eigvecs(2 * i, col), vB(i))
      if i < mEff then uRitz.axpyInPlace(eigvecs(2 * i + 1, col), uB(i))
      i += 1
    (normalize(uRitz), normalize(vRitz))

  private def normalize(x: MutableDVec): DVec =
    val nrm = x.asVec.norm2
    if nrm > 0.0 then (x.asVec * (1.0 / nrm)).copy else x.asVec.copy

  /** A view of `op` transposed: `Aᵀ` as an `n × m` operator, reusing `op`'s
    * forward/transpose applies swapped. Lets the bidiagonalization always run on
    * the taller orientation.
    */
  private def transposeOp(op: DoubleLinearOperator, m: Int, n: Int): DoubleLinearOperator =
    new DoubleLinearOperator:
      def rows: Int = n
      def cols: Int = m
      def applyTo(x: DVec, into: MutableDVec): Unit = op.transposeApplyTo(x, into)
      override def transposeApplyTo(y: DVec, into: MutableDVec): Unit = op.applyTo(y, into)

  /** `max(‖A v − σ u‖, ‖Aᵀ u − σ v‖)` — the two-sided singular-triplet residual. */
  private def svdResidual(op: DoubleLinearOperator, u: DVec, v: DVec, sigma: Double, m: Int, n: Int): Double =
    val av = MutableDVec.zeros(m)
    op.applyTo(v, av)
    av.axpyInPlace(-sigma, u)
    val atu = MutableDVec.zeros(n)
    op.transposeApplyTo(u, atu)
    atu.axpyInPlace(-sigma, v)
    math.max(av.asVec.norm2, atu.asVec.norm2)

  /** Assemble the result: singular values '''descending''', `U` (`m×count`), `Vᵀ`
    * (`count×n`), numerical `rank` (count of returned `σ > tol·σ_max`), and
    * diagnostics whose residuals align with the descending values and whose
    * orthogonality error is the worse of `U` and `V`.
    */
  private def assembleSvd(
      m: Int,
      n: Int,
      requestedK: Int,
      sigmas: Array[Double],
      uVecs: Array[DVec],
      vVecs: Array[DVec],
      residuals: Array[Double],
      wantVectors: Boolean,
      tol: Double,
      iterations: Int
  ): SVD =
    val order = sigmas.indices.sortBy(i => (-sigmas(i), i)).toArray
    val sortedSigmas = DVec.tabulate(order.length)(i => sigmas(order(i)))
    val sortedResiduals = DVec.tabulate(order.length)(i => residuals(order(i)))
    val sigmaMax = if order.isEmpty then 0.0 else sigmas(order(0))
    var rank = 0
    var i = 0
    while i < order.length do
      if sigmas(order(i)) > tol * sigmaMax then rank += 1
      i += 1
    val (u, vt, orthoErr) =
      if wantVectors && order.nonEmpty then
        val uMat = DMat.tabulate(m, order.length)((r, c) => uVecs(order(c))(r))
        val vtMat = DMat.tabulate(order.length, n)((r, c) => vVecs(order(r))(c))
        (uMat, vtMat, math.max(orthogonalityError(uMat), orthogonalityError(vtMat.t)))
      else (DMat.zeros(m, 0), DMat.zeros(0, n), 0.0)
    val diagnostics =
      SpectralDiagnostics(
        requested = requestedK,
        converged = order.length,
        residuals = sortedResiduals,
        orthogonalityError = orthoErr,
        iterations = iterations,
        rank = Some(rank)
      )
    SVD(sortedSigmas, u, vt, rank, diagnostics)

  // ===========================================================================
  // Shared helpers (mirroring Eigen; kept local to keep visibility narrow)
  // ===========================================================================

  /** `w := (I − Q Qᵀ) w` against the first `count` columns of `basis`. */
  private def reorthogonalize(w: MutableDVec, basis: Array[DVec], count: Int): Unit =
    var i = 0
    while i < count do
      val c = basis(i).dot(w.asVec)
      w.axpyInPlace(-c, basis(i))
      i += 1

  /** `‖QᵀQ − I‖_F` of the columns of `mat`. */
  private def orthogonalityError(mat: DMat): Double =
    val g = mat.t * mat
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

  /** A deterministic (LCG) unit start vector in `R^n`, or the caller's normalized.
    * The seed sequence is bit-for-bit portable (32-bit `Int` wraps identically on
    * JVM and Scala.js); the bidiagonalization built on it is deterministic '''per
    * platform''' but may differ between JVM and Scala.js in the last bits — the
    * dense kernels use the platform's fused multiply-add.
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
  // Generalized SVD (gsvd) — pure, full-column-rank pencils only
  // ===========================================================================

  /** Below this the cosine/sine snaps to an exact `0`: `s ≤` ⇒
    * [[GeneralizedSingularValue.Infinite]], `c ≤` ⇒
    * [[GeneralizedSingularValue.Zero]]. The `c`/`s` here are '''normalized'''
    * (`c² + s² = 1`), so the tolerance is on the unit CS scale. It sits far above
    * the `~ε` a direct norm of a true-null direction produces, and far below any
    * genuinely finite ratio the tests exercise.
    */
  private val CsSnapTolerance: Double = 1e-9

  /** Generalized SVD of the pencil `(A, B)` (`A` `m×n`, `B` `p×n`) computing the
    * full factors ([[EigenVectors.Right]]). See the three-argument overload for the
    * vector flag, algorithm, and failure details.
    */
  def gsvd(a: DMat, b: DMat)(using backend: SpectralBackend): Either[LinAlgError, GeneralizedSVD] =
    gsvd(a, b, EigenVectors.Right)

  /** Generalized SVD of `(A, B)` — `A = U C Xᵀ`, `B = V S Xᵀ`, `CᵀC + SᵀS = I`
    * (§ 9 of `docs/spectral-parity.md`) — for a '''full-column-rank''' pencil
    * (`[A;B]` has rank `n`). Values are the generalized singular values `c_i/s_i`
    * ('''descending''', `Infinite` first / `Zero` last), typed as
    * [[GeneralizedSingularValue]].
    *
    * '''Algorithm''' (the pure route that avoids the not-yet-shipped full dense
    * SVD): QR of the stacked `[A;B]` (`(m+p)×n`) splits the thin `Q` into `Q1`
    * (top `m` rows) and `Q2` (bottom `p` rows) with `A = Q1 R`, `B = Q2 R`. The CS
    * decomposition comes from the symmetric-eigen kernel on the '''A-block Gram'''
    * `Q1ᵀQ1 = W C² Wᵀ` (the `n×n` block; a single documented choice — `Q2ᵀQ2`
    * shares `W` since `Q1ᵀQ1 + Q2ᵀQ2 = I`). Crucially the cosines and sines are
    * then taken as the '''direct column norms''' `c_i = ‖Q1 w_i‖`,
    * `s_i = ‖Q2 w_i‖` (not `s = √(1−c²)`), which keeps small `s` (and small `c`)
    * accurate rather than suffering the cancellation of `1−c²`. Finally
    * `U = Q1 W C⁻¹`, `V = Q2 W S⁻¹`, `Xᵀ = Wᵀ R`.
    *
    * '''Accuracy caveat.''' The Gram step delivers '''absolute''' accuracy `~ε` in
    * `c²`/`s²`; for `c` or `s` very near `0` the '''relative''' accuracy degrades
    * (those are exactly the near-`Infinite`/near-`Zero` values). The exact
    * `Infinite` (`s = 0`) / `Zero` (`c = 0`) classifications use the documented
    * snap tolerance `1e-9` on the normalized CS scale. A `Zero` leaves its `U`
    * column undetermined by this route (stored as zero); an `Infinite` its `V`
    * column — orthonormality is claimed only on the well-determined columns.
    *
    * `vectors` selects [[EigenVectors.ValuesOnly]] (only `c`/`s`/values) versus
    * [[EigenVectors.Right]] (the full `U`, `V`, `X`); there is no one-sided mode
    * ([[EigenVectors.Left]]/[[EigenVectors.LeftAndRight]] are rejected).
    *
    * '''Diagnostics.''' `diagnostics.residuals(i)` is the '''CS-identity defect'''
    * `|c_raw_i² + s_raw_i² − 1|` (the thin-`Q` column-orthonormality drift), '''not'''
    * an `A`/`B` reconstruction residual — despite `SpectralDiagnostics`'s generic
    * "per-pair residual" wording. `orthogonalityError` is the worst column-Gram
    * error `‖·ᵀ· − I‖_F` over the well-determined `U`/`V` columns.
    *
    * `Left` on: `A`/`B` with disagreeing column counts (`DimensionMismatch`); an
    * empty dimension (`InvalidArgument`); or a '''rank-deficient stacked pencil'''
    * (`rank([A;B]) < n`, including `m+p < n`) — `RankDeficient`, the honest scope
    * boundary (rank-deficient GSVD is deferred to a backend, § 9). Rank is judged
    * by the QR `R`-diagonal at gale's standard tolerance
    * (`2·max(m+p,n)·ε·max|R_ii|`).
    */
  def gsvd(a: DMat, b: DMat, vectors: EigenVectors)(using
      backend: SpectralBackend
  ): Either[LinAlgError, GeneralizedSVD] =
    if a.cols != b.cols then Left(LinAlgError.DimensionMismatch(a.shape, b.shape))
    else
      val m = a.rows
      val p = b.rows
      val n = a.cols
      if m <= 0 || p <= 0 || n <= 0 then
        Left(LinAlgError.InvalidArgument(s"GSVD requires nonempty A (${m}x$n) and B (${p}x$n)"))
      else
        validateVectors(vectors) match
          case Left(error) => Left(error)
          case Right(wantVectors) =>
            // m+p is an upper bound on rank([A;B]) (the QR has not run yet); when
            // it is below n the pencil cannot be full column rank. The reported
            // rank is therefore this bound, not a measured value.
            if m + p < n then rankDeficientRoute(a, b, wantVectors, m + p, n)
            else
              val stacked = DMat.tabulate(m + p, n)((i, j) => if i < m then a(i, j) else b(i - m, j))
              val qr = DenseDecompositions.qr(stacked)
              val rank = qr.diagnostics.rank.getOrElse(n)
              if rank < n then rankDeficientRoute(a, b, wantVectors, rank, n)
              else runGsvd(m, p, n, qr, wantVectors)

  /** Rank-deficient GSVD seam (S6 of `docs/spectral-backend-boundary.md`): a
    * [[SpectralCapability.RankDeficientGsvd]]-capable backend '''computes''' the
    * hard pencil (raw factors canonicalized here into a [[GeneralizedSVD]] with
    * `Infinite`/`Zero` typed values); with no such backend — the pure default — the
    * shipped `Left(RankDeficient)` stands.
    */
  private def rankDeficientRoute(a: DMat, b: DMat, wantVectors: Boolean, rank: Int, n: Int)(using
      backend: SpectralBackend
  ): Either[LinAlgError, GeneralizedSVD] =
    if backend.capabilities.contains(SpectralCapability.RankDeficientGsvd) then
      backend.rankDeficientGsvd(a, b, wantVectors).map(assembleRawGsvd)
    else Left(LinAlgError.RankDeficient(rank, n))

  private def runGsvd(
      m: Int,
      p: Int,
      n: Int,
      qr: gale.linalg.QR,
      wantVectors: Boolean
  ): Either[LinAlgError, GeneralizedSVD] =
    val q = qr.q
    val rFull = qr.r
    val q1 = DMat.tabulate(m, n)((i, j) => q(i, j))
    val q2 = DMat.tabulate(p, n)((i, j) => q(m + i, j))
    val rMat = DMat.tabulate(n, n)((i, j) => if i <= j then rFull(i, j) else 0.0)
    val gram = q1.t * q1
    DenseSpectralKernels.symmetricEigen(gram, wantVectors = true) match
      case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iters)) =>
        Left(LinAlgError.DidNotConverge(iters, 0.0))
      case Right(eig) =>
        Right(assembleGsvd(q1, q2, rMat, eig.vectors.get, m, p, n, wantVectors))

  private def classifyCs(c: Double, s: Double): GeneralizedSingularValue =
    if s <= CsSnapTolerance then GeneralizedSingularValue.Infinite
    else if c <= CsSnapTolerance then GeneralizedSingularValue.Zero
    else GeneralizedSingularValue.Finite(c / s)

  private def assembleGsvd(
      q1: DMat,
      q2: DMat,
      rMat: DMat,
      w: DMat,
      m: Int,
      p: Int,
      n: Int,
      wantVectors: Boolean
  ): GeneralizedSVD =
    val q1w = new Array[DVec](n)
    val q2w = new Array[DVec](n)
    val cRaw = new Array[Double](n)
    val sRaw = new Array[Double](n)
    val cN = new Array[Double](n)
    val sN = new Array[Double](n)
    val resid = new Array[Double](n)
    val typed = new Array[GeneralizedSingularValue](n)
    var i = 0
    while i < n do
      val wi = w.col(i)
      val a1 = q1 * wi
      val a2 = q2 * wi
      q1w(i) = a1
      q2w(i) = a2
      val cr = a1.norm2
      val sr = a2.norm2
      cRaw(i) = cr
      sRaw(i) = sr
      // Direct-norm CS values normalized so c²+s²=1 exactly (the raw pair already
      // sums to ~1 by thin-Q column orthonormality; normalizing pins the identity).
      val r = math.hypot(cr, sr)
      val cc = if r > 0.0 then cr / r else 0.0
      val ss = if r > 0.0 then sr / r else 0.0
      cN(i) = cc
      sN(i) = ss
      resid(i) = math.abs(cr * cr + sr * sr - 1.0)
      typed(i) = classifyCs(cc, ss)
      i += 1

    // Descending by ratio: Infinite (+∞) first, Zero (0) last; index breaks ties.
    val order = (0 until n).sortBy(idx => (-typed(idx).value, idx)).toArray
    val cOut = DVec.tabulate(n)(k => cN(order(k)))
    val sOut = DVec.tabulate(n)(k => sN(order(k)))
    val valuesOut = order.iterator.map(idx => typed(idx)).toIndexedSeq
    val residOut = DVec.tabulate(n)(k => resid(order(k)))

    val (u, v, x, orthoErr) =
      if wantVectors then
        // U = Q1 W C⁻¹ (unit columns from the raw norm), zeroed where c snaps to 0;
        // V = Q2 W S⁻¹, zeroed where s snaps to 0.
        val uMat = DMat.tabulate(m, n): (row, k) =>
          val idx = order(k)
          if cN(idx) > CsSnapTolerance then q1w(idx)(row) / cRaw(idx) else 0.0
        val vMat = DMat.tabulate(p, n): (row, k) =>
          val idx = order(k)
          if sN(idx) > CsSnapTolerance then q2w(idx)(row) / sRaw(idx) else 0.0
        // Xᵀ = Wᵀ R ⇒ X = Rᵀ W, columns permuted with the values.
        val xFull = rMat.t * w
        val xMat = DMat.tabulate(n, n)((row, k) => xFull(row, order(k)))
        val uErr = columnGramError((0 until n).filter(k => valuesOut(k) != GeneralizedSingularValue.Zero).map(uMat.col))
        val vErr = columnGramError((0 until n).filter(k => valuesOut(k) != GeneralizedSingularValue.Infinite).map(vMat.col))
        (uMat, vMat, xMat, math.max(uErr, vErr))
      else (DMat.zeros(m, 0), DMat.zeros(p, 0), DMat.zeros(n, 0), 0.0)

    val diagnostics =
      SpectralDiagnostics(
        requested = n,
        converged = n,
        residuals = residOut,
        orthogonalityError = orthoErr,
        iterations = 0,
        rank = Some(n)
      )
    GeneralizedSVD(u, v, x, cOut, sOut, valuesOut, diagnostics)

  /** `‖MᵀM − I‖_F` for the matrix whose columns are `cols` (`0` when empty). */
  private def columnGramError(cols: IndexedSeq[DVec]): Double =
    val k = cols.length
    if k == 0 then 0.0
    else
      var sum = 0.0
      var i = 0
      while i < k do
        var j = 0
        while j < k do
          val g = cols(i).dot(cols(j))
          val d = if i == j then g - 1.0 else g
          sum += d * d
          j += 1
        i += 1
      math.sqrt(sum)

  /** Canonicalize a backend's raw GSVD factors into the sealed [[GeneralizedSVD]]
    * (the S6 rank-deficient path): normalize each `(c, s)` to the unit CS scale,
    * classify the typed generalized singular value, impose the descending-ratio
    * order (`Infinite` first / `Zero` last), permute `U`/`V`/`X` columns in lockstep,
    * and report the CS-identity-defect residuals + well-determined orthogonality
    * error — the same layout the pure full-rank path yields.
    */
  private def assembleRawGsvd(raw: RawGsvd): GeneralizedSVD =
    val n = raw.c.length
    val cN = new Array[Double](n)
    val sN = new Array[Double](n)
    val resid = new Array[Double](n)
    val typed = new Array[GeneralizedSingularValue](n)
    var i = 0
    while i < n do
      val cr = raw.c(i)
      val sr = raw.s(i)
      val r = math.hypot(cr, sr)
      val cc = if r > 0.0 then cr / r else 0.0
      val ss = if r > 0.0 then sr / r else 0.0
      cN(i) = cc
      sN(i) = ss
      resid(i) = math.abs(cr * cr + sr * sr - 1.0)
      typed(i) = classifyCs(cc, ss)
      i += 1
    val order = (0 until n).sortBy(idx => (-typed(idx).value, idx)).toArray
    val cOut = DVec.tabulate(n)(k => cN(order(k)))
    val sOut = DVec.tabulate(n)(k => sN(order(k)))
    val valuesOut = order.iterator.map(idx => typed(idx)).toIndexedSeq
    val residOut = DVec.tabulate(n)(k => resid(order(k)))
    val wantVectors = raw.u.cols != 0
    val (u, v, x, orthoErr) =
      if wantVectors then
        val uMat = DMat.tabulate(raw.u.rows, n)((row, k) => raw.u(row, order(k)))
        val vMat = DMat.tabulate(raw.v.rows, n)((row, k) => raw.v(row, order(k)))
        val xMat = DMat.tabulate(raw.x.rows, n)((row, k) => raw.x(row, order(k)))
        val uErr = columnGramError((0 until n).filter(k => valuesOut(k) != GeneralizedSingularValue.Zero).map(uMat.col))
        val vErr = columnGramError((0 until n).filter(k => valuesOut(k) != GeneralizedSingularValue.Infinite).map(vMat.col))
        (uMat, vMat, xMat, math.max(uErr, vErr))
      else (DMat.zeros(raw.u.rows, 0), DMat.zeros(raw.v.rows, 0), DMat.zeros(raw.x.rows, 0), 0.0)
    val diagnostics =
      SpectralDiagnostics(
        requested = n,
        converged = n,
        residuals = residOut,
        orthogonalityError = orthoErr,
        iterations = 0,
        rank = Some(n)
      )
    GeneralizedSVD(u, v, x, cOut, sOut, valuesOut, diagnostics)
