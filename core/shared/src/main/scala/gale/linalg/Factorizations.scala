package gale.linalg

import gale.backend.Backend
import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.PlatformMath.fma

final case class FactorizationDiagnostics(
    info: Int = 0,
    rank: Option[Int] = None
):
  def isSuccess: Boolean =
    info == 0

/** The row permutation produced by LU pivoting: `apply(i)` is the original row
  * that ends up in position `i`.
  *
  * `toArray` / `toIndexSeq` are public `Int`-index accessors. P4 restricted
  * `Double` storage to the platform-specific interop doorway, but index arrays
  * are plain `Array[Int]` and safe to expose: both accessors hand back an
  * independent copy, so callers can read or mutate the result without touching
  * the factorization's internal state.
  */
final class PivotVector private (private val values: Array[Int]):
  def length: Int =
    values.length

  def apply(index: Int): Int =
    if index < 0 || index >= values.length then
      throw LinAlgError.IndexOutOfBounds(index, values.length)
    values(index)

  /** A fresh `Array[Int]` copy of the pivots; safe to mutate. */
  def toArray: Array[Int] =
    values.clone()

  /** An immutable snapshot of the pivots, for callers that want no copy to guard. */
  def toIndexSeq: IndexedSeq[Int] =
    values.toIndexedSeq

object PivotVector:
  def apply(values: Int*): PivotVector =
    new PivotVector(values.toArray)

  private[gale] def fromArray(values: Array[Int]): PivotVector =
    new PivotVector(values.clone())

final case class LU private[gale] (
    packed: DMat,
    pivots: PivotVector,
    parity: Int,
    diagnostics: FactorizationDiagnostics
):
  def solve(b: DVec): Either[LinAlgError, DVec] =
    DenseDecompositions.solve(this, b)

  def det: Either[LinAlgError, Double] =
    DenseDecompositions.det(this)

/** Cholesky factorization `A = L Lᵀ` of a symmetric positive-definite matrix,
  * holding the lower factor `L`.
  *
  * The factorization reads only the '''lower triangle''' of the input (including
  * the diagonal); the strict upper triangle is never inspected, so `A` is treated
  * as symmetric by assumption and any asymmetry in its upper triangle is ignored.
  * A non-positive pivot yields `Left(`[[LinAlgError.NotPositiveDefinite]]`)`.
  */
final case class Cholesky private[gale] (
    lower: DMat,
    diagnostics: FactorizationDiagnostics
):
  def solve(b: DVec): Either[LinAlgError, DVec] =
    DenseDecompositions.solve(this, b)

/** Householder QR of an `m x n` matrix, stored compactly.
  *
  * Rather than materialising the `m x m` orthogonal factor, the factorization
  * keeps the `min(m, n)` elementary reflectors (`reflectors`, column `k` holding
  * the length-`m` Householder vector `v_k`, zero above row `k`) and their scalars
  * (`tau`, where `tau(k) = 2 / (v_k · v_k)`), plus the upper-triangular `r`. The
  * dense `q` is rebuilt on demand; `solveLeastSquares` applies `Qᵀ` implicitly
  * through the reflectors, so it never forms `Q`.
  */
final case class QR private[gale] (
    reflectors: DMat,
    tau: DoubleArray,
    r: DMat,
    diagnostics: FactorizationDiagnostics
):
  /** The orthogonal factor `Q` (`m x m`), materialised from the stored
    * reflectors on first access and cached thereafter.
    */
  lazy val q: DMat =
    DenseDecompositions.materializeQ(this)

  def solveLeastSquares(b: DVec): Either[LinAlgError, DVec] =
    DenseDecompositions.solveLeastSquares(this, b)

object DenseDecompositions:
  def lu(A: DMat): Either[LinAlgError, LU] =
    if A.rows != A.cols then
      Left(LinAlgError.NonSquareMatrix(A.shape))
    else
      val n = A.rows
      val packed = A.toDoubleArrayCopyRowMajor
      val pivots = new Array[Int](n)
      var i = 0
      while i < n do
        pivots(i) = i
        i += 1

      var parity = 1
      var k = 0
      while k < n do
        var pivot = k
        var maxAbs = math.abs(packed(k * n + k))
        i = k + 1
        while i < n do
          val candidate = math.abs(packed(i * n + k))
          if candidate > maxAbs then
            maxAbs = candidate
            pivot = i
          i += 1

        if maxAbs == 0.0 || maxAbs.isNaN then
          return Left(LinAlgError.SingularMatrix(k))

        if pivot != k then
          swapRows(packed, n, k, pivot)
          val tmpPivot = pivots(k)
          pivots(k) = pivots(pivot)
          pivots(pivot) = tmpPivot
          parity = -parity

        val pivotValue = packed(k * n + k)
        i = k + 1
        while i < n do
          val ik = i * n + k
          packed(ik) = packed(ik) / pivotValue
          val multiplier = packed(ik)
          var j = k + 1
          while j < n do
            packed(i * n + j) = packed(i * n + j) - multiplier * packed(k * n + j)
            j += 1
          i += 1
        k += 1

      Right(
        LU(
          packed = DMat.fromDoubleArrayOwned(n, n, packed),
          pivots = PivotVector.fromArray(pivots),
          parity = parity,
          diagnostics = FactorizationDiagnostics(info = 0)
        )
      )

  def cholesky(A: DMat): Either[LinAlgError, Cholesky] =
    if A.rows != A.cols then
      Left(LinAlgError.NonSquareMatrix(A.shape))
    else
      val n = A.rows
      val lower = DoubleArray.alloc(n * n)
      var i = 0
      while i < n do
        var j = 0
        while j <= i do
          var sum = A(i, j)
          var k = 0
          while k < j do
            sum -= lower(i * n + k) * lower(j * n + k)
            k += 1
          if i == j then
            if sum <= 0.0 || sum.isNaN then
              return Left(LinAlgError.NotPositiveDefinite(i))
            lower(i * n + j) = math.sqrt(sum)
          else
            lower(i * n + j) = sum / lower(j * n + j)
          j += 1
        i += 1
      Right(
        Cholesky(
          lower = DMat.fromDoubleArrayOwned(n, n, lower),
          diagnostics = FactorizationDiagnostics(info = 0)
        )
      )

  def qr(A: DMat)(using Backend): QR =
    qr(A, DenseWorkspace.forQR(A.rows, A.cols))

  def qr(A: DMat, workspace: DenseWorkspace)(using backend: Backend): QR =
    val m = A.rows
    val n = A.cols
    val r = A.toDoubleArrayCopyRowMajor
    val limit = math.min(m, n)
    // Compact storage: reflectors column k holds v_k (length m, zero above row k),
    // tau(k) = 2 / (v_k · v_k). No m x m Q is formed here; it is rebuilt on demand.
    val reflectors = DoubleArray.alloc(m * limit)
    val tau = DoubleArray.alloc(limit)
    val scratch = workspace.work(DenseWorkspace.qrWorkSize(m, n))

    if DenseWorkspace.usesBlockedQR(m, n) then
      factorBlockedQR(r, m, n, reflectors, limit, tau, scratch)
    else
      factorUnblockedQR(r, m, n, reflectors, limit, tau, scratch)

    val rank = rankFromUpperTriangular(r, m, n)
    QR(
      reflectors = DMat.fromDoubleArrayOwned(m, limit, reflectors),
      tau = tau,
      r = DMat.fromDoubleArrayOwned(m, n, r),
      diagnostics = FactorizationDiagnostics(info = 0, rank = Some(rank))
    )

  /** Small-shape QR: scalar Householder generation with row-major rank-1
    * updates. Keeping this path avoids compact-WY setup overhead where Gale is
    * already ahead of Breeze.
    */
  private def factorUnblockedQR(
      r: DoubleArray,
      m: Int,
      n: Int,
      reflectors: DoubleArray,
      limit: Int,
      tau: DoubleArray,
      scratch: DoubleArray
  ): Unit =
    var k = 0
    while k < limit do
      factorHouseholder(r, m, n, reflectors, limit, tau, k)
      applyReflectorToColumns(r, m, n, reflectors, limit, tau(k), k, k + 1, n, scratch, 0)
      k += 1

  /** Large-shape QR using compact WY panels. Each panel is factored with the
    * unblocked kernel, but its reflectors are aggregated as
    * `H = I - V T Vᵀ`; the trailing matrix receives `Hᵀ` through two GEMMs.
    */
  private def factorBlockedQR(
      r: DoubleArray,
      m: Int,
      n: Int,
      reflectors: DoubleArray,
      limit: Int,
      tau: DoubleArray,
      scratch: DoubleArray
  )(using backend: Backend): Unit =
    val block = DenseWorkspace.QrBlockSize
    val vtOffset = 0
    val wOffset = block * m
    val tOffset = wOffset + block * n
    var panelStart = 0
    while panelStart < limit do
      val panelEnd = math.min(panelStart + block, limit)
      var k = panelStart
      while k < panelEnd do
        factorHouseholder(r, m, n, reflectors, limit, tau, k)
        applyReflectorToColumns(r, m, n, reflectors, limit, tau(k), k, k + 1, panelEnd, scratch, wOffset)
        k += 1

      if panelEnd < n then
        val panelWidth = panelEnd - panelStart
        formCompactWY(reflectors, m, limit, tau, panelStart, panelWidth, scratch, tOffset, wOffset, block)
        applyCompactWYTranspose(
          r, m, n, reflectors, limit, panelStart, panelEnd, panelWidth,
          scratch, vtOffset, wOffset, tOffset, block
        )
      panelStart = panelEnd

  /** Build one normalized Householder reflector with a scaled `dnrm2` norm.
    * The stored vector has `v(k)=1`; `tau=2/(vᵀv)` and the transformed diagonal
    * is written directly to `R`. This is the stable `dlarfg` convention.
    */
  private def factorHouseholder(
      r: DoubleArray,
      m: Int,
      n: Int,
      reflectors: DoubleArray,
      limit: Int,
      tau: DoubleArray,
      k: Int
  ): Unit =
    val diagIndex = k * n + k
    val x0 = r(diagIndex)
    val norm = DoubleKernels.dnrm2(m - k, r, diagIndex, n)
    reflectors(k * limit + k) = 1.0
    if norm > 0.0 then
      val beta = if x0 >= 0.0 then -norm else norm
      // Compute through ratios: x0-beta can overflow even when x0, beta, and
      // every input are finite (for example x0 ~= -beta ~= 1e308).
      val x0OverBeta = x0 / beta
      val denominatorOverBeta = x0OverBeta - 1.0
      val tauK = 1.0 - x0OverBeta
      tau(k) = tauK
      r(diagIndex) = beta
      var i = k + 1
      while i < m do
        val idx = i * n + k
        reflectors(i * limit + k) = (r(idx) / beta) / denominatorOverBeta
        r(idx) = 0.0
        i += 1
    else
      tau(k) = 0.0
      var i = k + 1
      while i < m do
        r(i * n + k) = 0.0
        i += 1

  /** Apply `I - tau*v*vᵀ` to a contiguous row-major column interval. The dot
    * vector is accumulated row-first so the matrix traversal stays unit-stride.
    */
  private def applyReflectorToColumns(
      r: DoubleArray,
      m: Int,
      n: Int,
      reflectors: DoubleArray,
      limit: Int,
      tauK: Double,
      k: Int,
      colFrom: Int,
      colUntil: Int,
      scratch: DoubleArray,
      scratchOffset: Int
  ): Unit =
    val width = colUntil - colFrom
    if tauK != 0.0 && width > 0 then
      var j = 0
      while j < width do
        scratch(scratchOffset + j) = 0.0
        j += 1
      var i = k
      while i < m do
        val vi = reflectors(i * limit + k)
        val rRow = i * n + colFrom
        j = 0
        while j < width do
          val wj = scratchOffset + j
          scratch(wj) = fma(vi, r(rRow + j), scratch(wj))
          j += 1
        i += 1
      j = 0
      while j < width do
        scratch(scratchOffset + j) = tauK * scratch(scratchOffset + j)
        j += 1
      i = k
      while i < m do
        val vi = reflectors(i * limit + k)
        val rRow = i * n + colFrom
        j = 0
        while j < width do
          r(rRow + j) = fma(-vi, scratch(scratchOffset + j), r(rRow + j))
          j += 1
        i += 1

  /** Form the upper-triangular compact-WY factor `T` for one reflector panel. */
  private def formCompactWY(
      reflectors: DoubleArray,
      m: Int,
      limit: Int,
      tau: DoubleArray,
      panelStart: Int,
      panelWidth: Int,
      scratch: DoubleArray,
      tOffset: Int,
      tempOffset: Int,
      tStride: Int
  )(using backend: Backend): Unit =
    var i = 0
    var j = 0
    while i < panelWidth do
      j = 0
      while j < panelWidth do
        scratch(tOffset + i * tStride + j) = 0.0
        j += 1
      i += 1

    i = 0
    while i < panelWidth do
      val globalI = panelStart + i
      val tauI = tau(globalI)
      if tauI != 0.0 then
        j = 0
        while j < i do
          var dot = 0.0
          var row = globalI
          while row < m do
            dot = fma(
              reflectors(row * limit + panelStart + j),
              reflectors(row * limit + globalI),
              dot
            )
            row += 1
          scratch(tempOffset + j) = -tauI * dot
          j += 1

        j = 0
        while j < i do
          var sum = 0.0
          var l = j
          while l < i do
            sum = fma(
              scratch(tOffset + j * tStride + l),
              scratch(tempOffset + l),
              sum
            )
            l += 1
          scratch(tOffset + j * tStride + i) = sum
          j += 1
      scratch(tOffset + i * tStride + i) = tauI
      i += 1

  /** Apply the panel product `Hᵀ = I - V Tᵀ Vᵀ` to the trailing matrix. */
  private def applyCompactWYTranspose(
      r: DoubleArray,
      m: Int,
      n: Int,
      reflectors: DoubleArray,
      limit: Int,
      panelStart: Int,
      panelEnd: Int,
      panelWidth: Int,
      scratch: DoubleArray,
      vtOffset: Int,
      wOffset: Int,
      tOffset: Int,
      tStride: Int
  )(using backend: Backend): Unit =
    val rowsRemaining = m - panelStart
    val trailingCols = n - panelEnd

    // Pack V^T row-major; V itself stays row-major in `reflectors` and feeds the
    // second GEMM without a copy.
    var j = 0
    while j < panelWidth do
      var row = 0
      while row < rowsRemaining do
        scratch(vtOffset + j * rowsRemaining + row) =
          reflectors((panelStart + row) * limit + panelStart + j)
        row += 1
      j += 1

    // W := V^T C.
    dispatchGemm(
      panelWidth, trailingCols, rowsRemaining,
      1.0,
      scratch, vtOffset, rowsRemaining, 1,
      r, panelStart * n + panelEnd, n, 1,
      0.0,
      scratch, wOffset, trailingCols, 1
    )

    // W := T^T W. Descending rows make the triangular multiply safe in place.
    var col = 0
    while col < trailingCols do
      var i = panelWidth - 1
      while i >= 0 do
        var sum = 0.0
        var l = 0
        while l <= i do
          sum = fma(
            scratch(tOffset + l * tStride + i),
            scratch(wOffset + l * trailingCols + col),
            sum
          )
          l += 1
        scratch(wOffset + i * trailingCols + col) = sum
        i -= 1
      col += 1

    // C := C - V W.
    dispatchGemm(
      rowsRemaining, trailingCols, panelWidth,
      -1.0,
      reflectors, panelStart * limit + panelStart, limit, 1,
      scratch, wOffset, trailingCols, 1,
      1.0,
      r, panelStart * n + panelEnd, n, 1
    )

  private def dispatchGemm(
      rows: Int,
      cols: Int,
      shared: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      aColStride: Int,
      b: DoubleArray,
      bOffset: Int,
      bRowStride: Int,
      bColStride: Int,
      beta: Double,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int,
      cColStride: Int
  )(using backend: Backend): Unit =
    if backend.routesGemm(rows, cols, shared) then
      backend.denseDouble.gemm(
        rows, cols, shared, alpha,
        a, aOffset, aRowStride, aColStride,
        b, bOffset, bRowStride, bColStride,
        beta, c, cOffset, cRowStride, cColStride
      )
    else
      DoubleKernels.dgemm(
        rows, cols, shared, alpha,
        a, aOffset, aRowStride, aColStride,
        b, bOffset, bRowStride, bColStride,
        beta, c, cOffset, cRowStride, cColStride
      )

  /** Rebuild the dense `m x m` orthogonal factor from the stored reflectors.
    *
    * `Q = H_0 H_1 ... H_{limit-1}` with `H_k = I - tau_k v_k v_kᵀ`. Starting from
    * the identity and right-multiplying by each `H_k` in order accumulates `Q`.
    */
  private[gale] def materializeQ(qr: QR): DMat =
    val reflectors = qr.reflectors
    val m = reflectors.rows
    val limit = reflectors.cols
    val tau = qr.tau
    // reflectors and q are both contiguous row-major with offset 0.
    val vData = reflectors.data
    val q = identityArray(m)
    var k = 0
    while k < limit do
      val beta = tau(k)
      if beta != 0.0 then
        var row = 0
        while row < m do
          var dot = 0.0
          var i = k
          while i < m do
            dot = fma(q(row * m + i), vData(i * limit + k), dot)
            i += 1
          dot *= beta
          i = k
          while i < m do
            q(row * m + i) = fma(-dot, vData(i * limit + k), q(row * m + i))
            i += 1
          row += 1
      k += 1
    DMat.fromDoubleArrayOwned(m, m, q)

  def solve(A: DMat, b: DVec): Either[LinAlgError, DVec] =
    lu(A).flatMap(_.solve(b))

  def rankEstimate(A: DMat)(using Backend): Int =
    // QR always populates rank, so a single factorization suffices. The ambient
    // backend is threaded through so the internal blocked-QR gemms route exactly
    // as they do under the `DMat.qr` facade's pure path — one dispatch policy,
    // one numerical answer.
    qr(A).diagnostics.rank.get

  /** 1-norm condition number estimate `||A||_1 * ||A^{-1}||_1`.
    *
    * `||A^{-1}||_1` is estimated by Hager's algorithm (Higham's refinement):
    * a handful of matvecs with `A^{-1}` and `A^{-T}` locate a large column of
    * the inverse without forming it. Each matvec is one LU back/forward solve,
    * so the total cost is two LU factorizations (of `A` and `Aᵀ`) plus O(5)
    * solves — not the `n` solves a full inverse would need. The estimate is a
    * lower bound on the true 1-norm, so the returned condition number never
    * exceeds the exact one.
    */
  def conditionEstimate(A: DMat): Either[LinAlgError, Double] =
    if A.rows != A.cols then
      Left(LinAlgError.NonSquareMatrix(A.shape))
    else
      lu(A) match
        case Left(_: LinAlgError.SingularMatrix) =>
          Right(Double.PositiveInfinity)
        case Left(error) =>
          Left(error)
        case Right(luA) =>
          lu(A.t) match
            case Left(_: LinAlgError.SingularMatrix) =>
              Right(Double.PositiveInfinity)
            case Left(error) =>
              Left(error)
            case Right(luAt) =>
              hagerInverseOneNorm(A.rows, luA, luAt) match
                case Left(error)        => Left(error)
                case Right(inverseNorm) => Right(norm1(A) * inverseNorm)

  /** Hager/Higham estimate of `||A^{-1}||_1` given LU factors of `A` and `Aᵀ`.
    *
    * `luA` solves `A y = x` (applies `A^{-1}`); `luAt` solves `Aᵀ z = ξ`
    * (applies `A^{-T}`). Iterates the standard column-search recurrence for at
    * most five steps, stopping when the search index repeats or the local
    * optimality test `||z||_inf <= zᵀx` holds.
    */
  private def hagerInverseOneNorm(n: Int, luA: LU, luAt: LU): Either[LinAlgError, Double] =
    if n == 0 then
      Right(0.0)
    else
      // Platform scratch reused across iterations: `xVec` / `xiVec` are stable
      // views over `x` / `xi`, re-solved after each in-place update, and the
      // solve results are read straight from their backing storage. No plain
      // Array or fromArray/toArray round-trips, so JS stays typed-array clean.
      val x = DoubleArray.alloc(n)
      val xi = DoubleArray.alloc(n)
      val xVec = new DVec(x, Offset.unsafe(0), Length.unsafe(n), Stride.unsafe(1))
      val xiVec = new DVec(xi, Offset.unsafe(0), Length.unsafe(n), Stride.unsafe(1))
      val start = 1.0 / n
      var i = 0
      while i < n do
        x(i) = start
        i += 1
      var est = 0.0
      var jOld = -1
      var iter = 0
      var converged = false
      while iter < 5 && !converged do
        solve(luA, xVec) match
          case Left(error) =>
            return Left(error)
          case Right(yVec) =>
            val yData = yVec.data
            val yOff = yVec.offset.value
            val yStep = yVec.stride.value
            est = 0.0
            i = 0
            while i < n do
              val yi = yData(yOff + i * yStep)
              est += math.abs(yi)
              xi(i) = if yi >= 0.0 then 1.0 else -1.0
              i += 1
            solve(luAt, xiVec) match
              case Left(error) =>
                return Left(error)
              case Right(zVec) =>
                val zData = zVec.data
                val zOff = zVec.offset.value
                val zStep = zVec.stride.value
                var zMax = 0.0
                var j = 0
                var zdotx = 0.0
                i = 0
                while i < n do
                  val zi = zData(zOff + i * zStep)
                  zdotx += zi * x(i)
                  val az = math.abs(zi)
                  if az > zMax then
                    zMax = az
                    j = i
                  i += 1
                if zMax <= zdotx || j == jOld then
                  converged = true
                else
                  jOld = j
                  i = 0
                  while i < n do
                    x(i) = 0.0
                    i += 1
                  x(j) = 1.0
        iter += 1
      Right(est)

  def solve(lu: LU, b: DVec): Either[LinAlgError, DVec] =
    val n = lu.packed.rows
    if lu.packed.cols != n then
      Left(LinAlgError.NonSquareMatrix(lu.packed.shape))
    else if b.length != n then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(n), Cols(1)), Shape(Rows(b.length), Cols(1))))
    else
      // Hoist the packed matrix's raw storage once and drive the shared
      // triangular-solve kernel: gather the permuted right-hand side, then a
      // unit-lower forward sweep and a non-unit-upper back sweep in place.
      val pData = lu.packed.data
      val pOff = lu.packed.offset.value
      val pRowStep = lu.packed.rowStride.value
      val pColStep = lu.packed.colStride.value
      val bData = b.data
      val bOff = b.offset.value
      val bStep = b.stride.value
      val pivots = lu.pivots
      // Platform scratch throughout, so JS keeps typed-array storage and the
      // result vector owns `x` directly (no Array round-trip).
      val x = DoubleArray.alloc(n)
      var i = 0
      while i < n do
        x(i) = bData(bOff + pivots(i) * bStep)
        i += 1
      // L y = Pb: unit lower diagonal, so the guard never trips here.
      DoubleKernels.dtrsv(n, lower = true, unit = true, 0.0, pData, pOff, pRowStep, pColStep, x, 0, 1)
      // U x = y: a successful LU has nonzero pivots, but guard exact zeros anyway.
      val info = DoubleKernels.dtrsv(n, lower = false, unit = false, 0.0, pData, pOff, pRowStep, pColStep, x, 0, 1)
      if info >= 0 then
        Left(LinAlgError.SingularMatrix(info))
      else
        Right(DVec.fromDoubleArrayOwned(x))

  def solve(cholesky: Cholesky, b: DVec): Either[LinAlgError, DVec] =
    val n = cholesky.lower.rows
    if cholesky.lower.cols != n then
      Left(LinAlgError.NonSquareMatrix(cholesky.lower.shape))
    else if b.length != n then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(n), Cols(1)), Shape(Rows(b.length), Cols(1))))
    else
      // Hoist the lower factor's raw storage once and drive the shared triangular
      // kernel twice: `L y = b` (non-unit lower), then `Lᵀ x = y`. The transpose
      // is the same storage with row/column strides swapped, so no copy is made.
      val lData = cholesky.lower.data
      val lOff = cholesky.lower.offset.value
      val lRowStep = cholesky.lower.rowStride.value
      val lColStep = cholesky.lower.colStride.value
      // Owned contiguous copy of b, mutated in place into the solution.
      val x = b.toDoubleArrayOwnedCopy
      val forward = DoubleKernels.dtrsv(n, lower = true, unit = false, 0.0, lData, lOff, lRowStep, lColStep, x, 0, 1)
      if forward >= 0 then
        Left(LinAlgError.NotPositiveDefinite(forward))
      else
        val back = DoubleKernels.dtrsv(n, lower = false, unit = false, 0.0, lData, lOff, lColStep, lRowStep, x, 0, 1)
        if back >= 0 then
          Left(LinAlgError.NotPositiveDefinite(back))
        else
          Right(DVec.fromDoubleArrayOwned(x))

  def solveLeastSquares(qr: QR, b: DVec): Either[LinAlgError, DVec] =
    val m = qr.reflectors.rows
    val n = qr.r.cols
    if qr.r.rows != m then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(m), Cols(n)), qr.r.shape))
    else if b.length != m then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(m), Cols(1)), Shape(Rows(b.length), Cols(1))))
    else if m < n then
      Left(LinAlgError.UnsupportedOperation("underdetermined least squares"))
    else
      val rank = qr.diagnostics.rank.getOrElse(rankFromMatrix(qr.r))
      if rank < n then
        Left(LinAlgError.RankDeficient(rank, n))
      else
        // y := Qᵀ b, applied implicitly as H_{limit-1} ... H_0 b through the
        // stored reflectors — no m x m Q is ever formed. Order matters: apply
        // H_0 first (the reflectors were applied to R in ascending k).
        val limit = qr.reflectors.cols
        val vData = qr.reflectors.data
        val tau = qr.tau
        // Owned platform copy of b: mutated in place as Qᵀ is applied.
        val y = b.toDoubleArrayOwnedCopy
        var k = 0
        while k < limit do
          val beta = tau(k)
          if beta != 0.0 then
            var dot = 0.0
            var i = k
            while i < m do
              dot = fma(vData(i * limit + k), y(i), dot)
              i += 1
            dot *= beta
            i = k
            while i < m do
              y(i) = fma(-dot, vData(i * limit + k), y(i))
              i += 1
          k += 1

        // Back-substitute the upper-triangular n x n leading block of R through
        // the shared kernel. `y` holds Qᵀb over its first n entries; copy those
        // into an owned length-n buffer and solve in place. A diagonal at or below
        // the rank tolerance marks rank deficiency at that column.
        val tolerance = rankToleranceFromMatrix(qr.r)
        val rData = qr.r.data
        val rOff = qr.r.offset.value
        val rRowStep = qr.r.rowStride.value
        val rColStep = qr.r.colStride.value
        val x = DoubleArray.alloc(n)
        DoubleKernels.dcopy(n, y, 0, 1, x, 0, 1)
        val info = DoubleKernels.dtrsv(n, lower = false, unit = false, tolerance, rData, rOff, rRowStep, rColStep, x, 0, 1)
        if info >= 0 then
          Left(LinAlgError.RankDeficient(info, n))
        else
          Right(DVec.fromDoubleArrayOwned(x))

  def det(lu: LU): Either[LinAlgError, Double] =
    val n = lu.packed.rows
    if lu.packed.cols != n then
      Left(LinAlgError.NonSquareMatrix(lu.packed.shape))
    else
      var out = lu.parity.toDouble
      var i = 0
      while i < n do
        out *= lu.packed(i, i)
        i += 1
      Right(out)

  private def swapRows(values: DoubleArray, cols: Int, r1: Int, r2: Int): Unit =
    var col = 0
    while col < cols do
      val i = r1 * cols + col
      val j = r2 * cols + col
      val tmp = values(i)
      values(i) = values(j)
      values(j) = tmp
      col += 1

  private def identityArray(size: Int): DoubleArray =
    val out = DoubleArray.alloc(size * size)
    var i = 0
    while i < size do
      out(i * size + i) = 1.0
      i += 1
    out

  private[gale] def rankFromMatrix(A: DMat): Int =
    rankFromUpperTriangular(A.toDoubleArrayCopyRowMajor, A.rows, A.cols)

  private def rankFromUpperTriangular(values: DoubleArray, rows: Int, cols: Int): Int =
    val tolerance = rankToleranceFromUpperTriangular(values, rows, cols)
    val limit = math.min(rows, cols)
    var rank = 0
    var i = 0
    while i < limit do
      if math.abs(values(i * cols + i)) > tolerance then
        rank += 1
      i += 1
    rank

  private[gale] def rankToleranceFromMatrix(A: DMat): Double =
    rankToleranceFromUpperTriangular(A.toDoubleArrayCopyRowMajor, A.rows, A.cols)

  private def rankToleranceFromUpperTriangular(values: DoubleArray, rows: Int, cols: Int): Double =
    val limit = math.min(rows, cols)
    var maxDiag = 0.0
    var i = 0
    while i < limit do
      maxDiag = math.max(maxDiag, math.abs(values(i * cols + i)))
      i += 1
    // Tolerance scales with the matrix magnitude; a hard floor at 1.0 would
    // wrongly rank uniformly tiny (but well-conditioned) matrices as deficient.
    // An all-zero matrix yields maxDiag == 0, hence tolerance 0, hence rank 0.
    // Two guards the final diagonal against the pair of rounded Householder
    // reductions that can leave an exactly dependent column a few ulps above the
    // textbook `max(m,n)*eps` cutoff.
    2.0 * math.max(rows, cols).toDouble * 2.220446049250313e-16 * maxDiag

  /** Matrix 1-norm: the maximum absolute column sum. */
  private def norm1(A: DMat): Double =
    val data = A.data
    val base = A.offset.value
    val rowStep = A.rowStride.value
    val colStep = A.colStride.value
    var out = 0.0
    var j = 0
    while j < A.cols do
      var sum = 0.0
      var i = 0
      var idx = base + j * colStep
      while i < A.rows do
        sum += math.abs(data(idx))
        idx += rowStep
        i += 1
      out = math.max(out, sum)
      j += 1
    out
