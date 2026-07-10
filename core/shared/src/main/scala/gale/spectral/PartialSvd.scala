package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import gale.linalg.Rows
import gale.linalg.Shape

/** Public partial singular value decomposition (`svds`, phase a of
  * `docs/spectral-parity.md` § 3 / § 8).
  *
  * v0.3.5 ships '''partial''' SVD only — full/economy dense SVD is deferred (it
  * needs a bidiagonal-SVD kernel outside the phase-a plan, § 3). The `k` requested
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
  // Dense partial SVD
  // ===========================================================================

  /** Partial SVD of the dense `a`, computing singular vectors ([[EigenVectors.Right]]).
    * See the three-argument overload for the vector flag and failure details.
    */
  def svd(a: DMat, selection: SingularSelection): Either[LinAlgError, SVD] =
    svd(a, selection, EigenVectors.Right)

  /** Partial SVD of the dense `a` — `A = U Σ Vᵀ` with the `k` singular triplets
    * named by `selection`, `Σ` descending. `vectors` chooses
    * [[EigenVectors.ValuesOnly]] versus [[EigenVectors.Right]] (there is no
    * one-sided singular-vector mode in v0.3.5, § 8: vectors means '''both''' `U`
    * and `V`, or neither; [[EigenVectors.Left]]/[[EigenVectors.LeftAndRight]] are
    * rejected). Uses the default [[SpectralOptions]]; for tolerance/subspace
    * control drive the operator overload with `a` as the operator.
    *
    * `Left` on: a non-positive dimension; `k ≤ 0` or `k ≥ min(m, n)` (full SVD is
    * deferred — the message points there); [[SingularSelection.All]] (Count-only,
    * per § 8); or an illegal vector flag.
    */
  def svd(a: DMat, selection: SingularSelection, vectors: EigenVectors): Either[LinAlgError, SVD] =
    svdPartial(a, a.rows, a.cols, selection, SpectralOptions(returnVectors = vectors))

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
              "partial SVD selects by Count(k, order); full/All SVD is deferred (docs/spectral-parity.md §3)"
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
                    s"k=$k must be in [1, ${p - 1}] for partial SVD; full SVD (k = min(m,n)) is deferred"
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
    * 32-bit `Int` arithmetic wraps identically on JVM and Scala.js, so the seed
    * sequence — and the whole bidiagonalization — is bit-for-bit portable.
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
