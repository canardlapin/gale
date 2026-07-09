package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Internal dense spectral kernels — the numerical foundation every v0.3.5
  * spectral feature builds on (`docs/spectral-parity.md`, fixed constraint 5).
  *
  * Two paths, all `private[gale]`:
  *
  *   - '''Symmetric.''' Householder tridiagonalization `A = Q T Qᵀ`
  *     ([[tridiagonalize]], Golub & Van Loan Alg. 8.3.1 / EISPACK `tred2`
  *     conventions) feeding an implicit-shift tridiagonal QL/QR eigensolver with
  *     Wilkinson shifts ([[symmetricTridiagonalEigen]], EISPACK `tql2`). The
  *     solver accumulates eigenvectors into a provided basis, so it composes both
  *     with the tridiagonalization `Q` (dense symmetric eigen, [[symmetricEigen]])
  *     and with an identity (the standalone tridiagonal problem the Lanczos
  *     projected problems need). Eigenvalues come out '''ascending-algebraic''',
  *     eigenvector columns permuted to match.
  *
  *   - '''Nonsymmetric.''' Householder reduction to upper Hessenberg form
  *     `A = Q H Qᵀ` ([[hessenberg]], EISPACK `orthes`/`ortran`) feeding Francis
  *     double-shift QR on the Hessenberg form ([[nonsymmetricEigen]], EISPACK
  *     `hqr2`). Eigenvalues are real/imaginary structure-of-arrays honouring the
  *     conjugate-pair convention (adjacent, positive-imaginary member first) with
  *     exact conjugate symmetry, so the [[NonsymmetricEigenDecomposition]]
  *     constructor accepts them directly; eigenvectors are returned in LAPACK
  *     real-Schur packed columns.
  *
  * '''Determinism.''' Every routine uses only `scala.math` (`sqrt`, `abs`,
  * comparisons) and a hand-rolled [[pythag]] built on IEEE-correctly-rounded
  * `sqrt`, so a given input yields bit-identical output on the JVM and Scala.js —
  * no `hypot`, no fused-multiply-add assumptions. Householder norms scale
  * explicitly (the `tred2`/`orthes` `scale` factor) to avoid intermediate
  * overflow/underflow, matching the dense QR kernel's discipline.
  *
  * '''Failure.''' The iterative solvers guard their sweep count. Reductions
  * ([[tridiagonalize]], [[hessenberg]]) are finite and cannot fail, so they
  * return their result directly; the QR iterations return
  * `Either[`[[SpectralKernelFailure]]`, _]`, signalling exhaustion as a typed
  * `Left` rather than looping forever. Callers map this onto the public
  * `SpectralDiagnostics` / `LinAlgError.DidNotConverge` idiom.
  */
private[gale] object DenseSpectralKernels:

  /** A spectral kernel exhausted its iteration budget without converging. Carries
    * the sweep count actually performed, which the caller surfaces as the
    * `iterations` of a `LinAlgError.DidNotConverge`.
    */
  enum SpectralKernelFailure:
    case DidNotConverge(iterations: Int)

  /** Householder tridiagonalization `A = Q T Qᵀ` of a symmetric matrix.
    *
    *   - `diagonal` — the length-`n` diagonal of `T`.
    *   - `offDiagonal` — the length-`max(n - 1, 0)` off-diagonal of `T`, with
    *     `offDiagonal(k) = T(k, k+1) = T(k+1, k)`.
    *   - `q` — the accumulated orthogonal factor (`n x n`) when requested, else
    *     `None` (the values-only reduction skips all accumulation work).
    */
  final case class Tridiagonalization(diagonal: DVec, offDiagonal: DVec, q: Option[DMat])

  /** Ascending-algebraic eigenvalues with, optionally, the matching orthonormal
    * eigenvectors as aligned columns (`vectors.get.col(i)` pairs with
    * `values(i)`).
    */
  final case class SymmetricEigen(values: DVec, vectors: Option[DMat])

  /** Householder reduction `A = Q H Qᵀ` to upper Hessenberg form. `h` has exact
    * zeros strictly below the subdiagonal; `q` is the accumulated orthogonal
    * factor when requested.
    */
  final case class Hessenberg(h: DMat, q: Option[DMat])

  /** Real/imaginary eigenvalue parts in real-Schur diagonal order — conjugate
    * pairs adjacent, positive-imaginary member first, with `im(j+1) == -im(j)`
    * and `re(j+1) == re(j)` holding '''exactly''' for a pair. `vectors`, when
    * present, holds the right eigenvectors in LAPACK real-Schur packed columns
    * (a real eigenvalue owns one column; a complex pair owns two, real part in
    * the lower-indexed column, imaginary part in the higher).
    */
  final case class NonsymmetricEigen(re: DVec, im: DVec, vectors: Option[DMat])

  // ---------------------------------------------------------------------------
  // Symmetric path
  // ---------------------------------------------------------------------------

  /** Tridiagonalize a symmetric `A` (`n x n`). Only the '''lower triangle''' of
    * `A` is read; the strict upper triangle is treated as its mirror image, so
    * any asymmetry there is ignored — matching the `Cholesky` precedent. When
    * `wantQ` is false the accumulation of `Q` is skipped entirely.
    */
  def tridiagonalize(a: DMat, wantQ: Boolean): Tridiagonalization =
    val n = a.rows
    require(a.cols == n, "tridiagonalize requires a square matrix")
    val work = symmetrizedLowerRowMajor(a, n)
    val d = DoubleArray.alloc(n)
    val e = DoubleArray.alloc(n) // EISPACK convention: e(0)=0, e(i)=T(i-1,i)
    tred2(n, work, d, e, wantQ)
    val offLen = math.max(n - 1, 0)
    val off = DoubleArray.alloc(offLen)
    var k = 0
    while k < offLen do
      off(k) = e(k + 1)
      k += 1
    val q =
      if wantQ then Some(DMat.fromDoubleArrayOwned(n, n, work)) else None
    Tridiagonalization(
      DVec.fromDoubleArrayOwned(d),
      DVec.fromDoubleArrayOwned(off),
      q
    )

  /** Eigenvalues (ascending) and optionally eigenvectors of a symmetric
    * tridiagonal matrix given by its `diagonal` (length `n`) and `offDiagonal`
    * (length `n - 1`, `offDiagonal(k) = T(k, k+1)`).
    *
    * This is the standalone tridiagonal solver: with `wantVectors` the
    * eigenvectors are accumulated onto an identity basis, so the returned columns
    * are the eigenvectors of `T` itself — exactly what a Lanczos projected
    * problem needs before transforming by its Krylov basis. `maxSweepsPerValue`
    * caps the QL/QR sweeps spent on any single eigenvalue (EISPACK's classic
    * bound is 30); exhaustion returns `Left(DidNotConverge)`.
    */
  def symmetricTridiagonalEigen(
      diagonal: DVec,
      offDiagonal: DVec,
      wantVectors: Boolean,
      maxSweepsPerValue: Int = 30
  ): Either[SpectralKernelFailure, SymmetricEigen] =
    val n = diagonal.length
    require(
      offDiagonal.length == math.max(n - 1, 0),
      s"offDiagonal length ${offDiagonal.length} must be ${math.max(n - 1, 0)} for diagonal length $n"
    )
    val d = diagonal.toDoubleArrayOwnedCopy
    // tql2's internal e(i) is the subdiagonal T(i-1,i) for i=1..n-1; e(0) unused.
    val e = DoubleArray.alloc(n)
    var i = 1
    while i < n do
      e(i) = offDiagonal(i - 1)
      i += 1
    val z = if wantVectors then Some(identityRowMajor(n)) else None
    solveTridiagonal(n, d, e, z, maxSweepsPerValue)

  /** Dense symmetric eigendecomposition `A V = V diag(λ)` with `λ` ascending and
    * `V` orthonormal (columns aligned with `values`). Composes
    * [[tridiagonalize]] with the tridiagonal QL/QR solver, accumulating
    * eigenvectors through the tridiagonalization `Q` when `wantVectors`. Reads
    * only the lower triangle of `A`.
    */
  def symmetricEigen(
      a: DMat,
      wantVectors: Boolean,
      maxSweepsPerValue: Int = 30
  ): Either[SpectralKernelFailure, SymmetricEigen] =
    val n = a.rows
    require(a.cols == n, "symmetricEigen requires a square matrix")
    val work = symmetrizedLowerRowMajor(a, n)
    val d = DoubleArray.alloc(n)
    val e = DoubleArray.alloc(n)
    tred2(n, work, d, e, wantVectors)
    // `work` now holds Q (when accumulating); tql2 rotates it into V.
    val z = if wantVectors then Some(work) else None
    solveTridiagonal(n, d, e, z, maxSweepsPerValue)

  /** Symmetric Householder tridiagonalization (EISPACK `tred2`), in place on the
    * `n x n` row-major `a`.
    *
    * On exit `d(0..n-1)` is the diagonal of `T` and `e(1..n-1)` its subdiagonal
    * (`e(0) = 0`). When `accumulate` is true `a` is overwritten with the
    * orthogonal `Q` such that `A = Q T Qᵀ`; when false the accumulation stores
    * and the final back-transform are skipped, and `a`'s contents are scratch.
    *
    * The reduction reads and updates only the lower triangle plus the row being
    * eliminated, so `A`'s upper triangle is never consulted. Householder norms
    * are formed on the explicitly `scale`-normalized row to avoid overflow.
    */
  private def tred2(n: Int, a: DoubleArray, d: DoubleArray, e: DoubleArray, accumulate: Boolean): Unit =
    // Reduce from the last row inward; row i is eliminated against columns 0..i-1.
    var i = n - 1
    while i >= 1 do
      val l = i - 1
      var h = 0.0
      var scale = 0.0
      if l > 0 then
        var k = 0
        while k <= l do
          scale += math.abs(a(i * n + k))
          k += 1
        if scale == 0.0 then
          // Row already reduced; its off-diagonal is the lone stored entry.
          e(i) = a(i * n + l)
        else
          k = 0
          while k <= l do
            val v = a(i * n + k) / scale
            a(i * n + k) = v
            h += v * v
            k += 1
          val f0 = a(i * n + l)
          val g0 = if f0 >= 0.0 then -math.sqrt(h) else math.sqrt(h)
          e(i) = scale * g0
          h -= f0 * g0
          a(i * n + l) = f0 - g0
          var f = 0.0
          var j = 0
          while j <= l do
            if accumulate then
              a(j * n + i) = a(i * n + j) / h
            // g = (A u)_j using only the lower triangle.
            var g = 0.0
            var kk = 0
            while kk <= j do
              g += a(j * n + kk) * a(i * n + kk)
              kk += 1
            kk = j + 1
            while kk <= l do
              g += a(kk * n + j) * a(i * n + kk)
              kk += 1
            e(j) = g / h
            f += e(j) * a(i * n + j)
            j += 1
          val hh = f / (h + h)
          j = 0
          while j <= l do
            val fj = a(i * n + j)
            val gj = e(j) - hh * fj
            e(j) = gj
            var kk = 0
            while kk <= j do
              val idx = j * n + kk
              a(idx) = a(idx) - (fj * e(kk) + gj * a(i * n + kk))
              kk += 1
            j += 1
      else
        e(i) = a(i * n + l)
      d(i) = h
      i -= 1

    if accumulate then d(0) = 0.0
    e(0) = 0.0
    i = 0
    while i < n do
      if accumulate then
        val l = i - 1
        if d(i) != 0.0 then
          var j = 0
          while j <= l do
            var g = 0.0
            var k = 0
            while k <= l do
              g += a(i * n + k) * a(k * n + j)
              k += 1
            k = 0
            while k <= l do
              val idx = k * n + j
              a(idx) = a(idx) - g * a(k * n + i)
              k += 1
            j += 1
        d(i) = a(i * n + i)
        a(i * n + i) = 1.0
        var j = 0
        while j <= l do
          a(j * n + i) = 0.0
          a(i * n + j) = 0.0
          j += 1
      else
        d(i) = a(i * n + i)
      i += 1

  /** Implicit-shift tridiagonal QL solver with Wilkinson shifts (EISPACK `tql2`).
    *
    * `d` holds the diagonal on entry, the eigenvalues on exit; `e(1..n-1)` holds
    * the subdiagonal on entry (`e(0)` arbitrary — the routine shifts it out).
    * When `z` is `Some(zData)` (an `n x n` row-major basis) its columns are
    * rotated in lockstep so that, starting from `Q` or the identity, they end as
    * the eigenvectors. Eigenvalues are sorted ascending afterwards with `z`'s
    * columns permuted to match. Returns `Left(DidNotConverge)` if any eigenvalue
    * needs more than `maxSweeps` QL sweeps.
    */
  private def solveTridiagonal(
      n: Int,
      d: DoubleArray,
      e: DoubleArray,
      z: Option[DoubleArray],
      maxSweeps: Int
  ): Either[SpectralKernelFailure, SymmetricEigen] =
    if n == 0 then
      return Right(SymmetricEigen(DVec.fromDoubleArrayOwned(d), z.map(_ => DMat.zeros(0, 0))))
    // Renumber the subdiagonal down by one (e(i-1) := e(i)), EISPACK convention.
    var i = 1
    while i < n do
      e(i - 1) = e(i)
      i += 1
    e(n - 1) = 0.0

    var l = 0
    while l < n do
      var iter = 0
      var continue = true
      while continue do
        // Find a small off-diagonal e(m) to split the problem at.
        var m = l
        var found = false
        while m < n - 1 && !found do
          val dd = math.abs(d(m)) + math.abs(d(m + 1))
          if math.abs(e(m)) <= Epsilon * dd then found = true
          else m += 1
        if m == l then
          continue = false
        else
          if iter == maxSweeps then
            return Left(SpectralKernelFailure.DidNotConverge(iter))
          iter += 1
          // Wilkinson-shifted implicit QL step on the block [l..m].
          var g = (d(l + 1) - d(l)) / (2.0 * e(l))
          var r = pythag(g, 1.0)
          g = d(m) - d(l) + e(l) / (g + sign(r, g))
          var s = 1.0
          var c = 1.0
          var p = 0.0
          var innerZero = false
          var iBt = m - 1
          while iBt >= l && !innerZero do
            var f = s * e(iBt)
            val b = c * e(iBt)
            r = pythag(f, g)
            e(iBt + 1) = r
            if r == 0.0 then
              // Recover from underflow: deflate and restart the sweep.
              d(iBt + 1) = d(iBt + 1) - p
              e(m) = 0.0
              innerZero = true
            else
              s = f / r
              c = g / r
              g = d(iBt + 1) - p
              r = (d(iBt) - g) * s + 2.0 * c * b
              p = s * r
              d(iBt + 1) = g + p
              g = c * r - b
              z.foreach: zData =>
                var k = 0
                while k < n do
                  val f2 = zData(k * n + iBt + 1)
                  zData(k * n + iBt + 1) = s * zData(k * n + iBt) + c * f2
                  zData(k * n + iBt) = c * zData(k * n + iBt) - s * f2
                  k += 1
              iBt -= 1
          if innerZero && iBt >= l then
            // The inner loop broke early (r == 0); retry without finalizing.
            ()
          else
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0
      l += 1

    sortAscending(n, d, z)
    val values = DVec.fromDoubleArrayOwned(d)
    val vectors = z.map(zData => DMat.fromDoubleArrayOwned(n, n, zData))
    Right(SymmetricEigen(values, vectors))

  /** Selection sort of `d` ascending, permuting the columns of the optional
    * `n x n` row-major `z` in lockstep. `n` is small (dense spectra), so the
    * `O(n²)` comparisons are irrelevant against the eigensolve.
    */
  private def sortAscending(n: Int, d: DoubleArray, z: Option[DoubleArray]): Unit =
    var i = 0
    while i < n - 1 do
      var k = i
      var p = d(i)
      var j = i + 1
      while j < n do
        if d(j) < p then
          k = j
          p = d(j)
        j += 1
      if k != i then
        d(k) = d(i)
        d(i) = p
        z.foreach: zData =>
          var row = 0
          while row < n do
            val tmp = zData(row * n + i)
            zData(row * n + i) = zData(row * n + k)
            zData(row * n + k) = tmp
            row += 1
      i += 1

  // ---------------------------------------------------------------------------
  // Nonsymmetric path
  // ---------------------------------------------------------------------------

  /** Householder reduction of a general `A` (`n x n`) to upper Hessenberg form
    * `A = Q H Qᵀ`. `H` has exact zeros strictly below the subdiagonal; `Q` is
    * accumulated only when `wantQ` is true. No balancing is applied (deferred per
    * `docs/spectral-parity.md` § 2).
    */
  def hessenberg(a: DMat, wantQ: Boolean): Hessenberg =
    val n = a.rows
    require(a.cols == n, "hessenberg requires a square matrix")
    val h = a.toDoubleArrayCopyRowMajor
    val ort = DoubleArray.alloc(n)
    orthes(n, h, ort)
    // Q must be accumulated from the reflectors stored below the subdiagonal
    // before those cells are cleared to expose a clean Hessenberg H.
    val q = if wantQ then Some(DMat.fromDoubleArrayOwned(n, n, ortran(n, h, ort))) else None
    zeroBelowSubdiagonal(n, h)
    Hessenberg(DMat.fromDoubleArrayOwned(n, n, h), q)

  /** Nonsymmetric eigendecomposition via Hessenberg reduction + Francis
    * double-shift QR. Eigenvalues come out in real-Schur diagonal order with the
    * conjugate-pair convention (adjacent, positive-imaginary first, exact
    * conjugate symmetry). With `wantVectors` the right eigenvectors are returned
    * in LAPACK real-Schur packed columns. `maxSweeps` (default `30·n`) bounds the
    * total QR steps; exhaustion returns `Left(DidNotConverge)`.
    */
  def nonsymmetricEigen(
      a: DMat,
      wantVectors: Boolean,
      maxSweeps: Int = -1
  ): Either[SpectralKernelFailure, NonsymmetricEigen] =
    val n = a.rows
    require(a.cols == n, "nonsymmetricEigen requires a square matrix")
    val budget = if maxSweeps < 0 then 30 * math.max(n, 1) else maxSweeps
    val h = a.toDoubleArrayCopyRowMajor
    val ort = DoubleArray.alloc(n)
    orthes(n, h, ort)
    val z = if wantVectors then ortran(n, h, ort) else DoubleArray.alloc(0)
    zeroBelowSubdiagonal(n, h)
    val wr = DoubleArray.alloc(n)
    val wi = DoubleArray.alloc(n)
    francisSchur(n, h, wr, wi, z, wantVectors, budget) match
      case Some(sweepsDone) =>
        Left(SpectralKernelFailure.DidNotConverge(sweepsDone))
      case None =>
        val vectors =
          if wantVectors then
            normalizeSchurVectors(n, wr, wi, z)
            Some(DMat.fromDoubleArrayOwned(n, n, z))
          else None
        Right(NonsymmetricEigen(DVec.fromDoubleArrayOwned(wr), DVec.fromDoubleArrayOwned(wi), vectors))

  /** Householder reduction to upper Hessenberg form (EISPACK `orthes`), in place
    * on the `n x n` row-major `h`, with `low = 0`, `high = n - 1`. The Householder
    * vectors are left in the strictly-below-subdiagonal cells (column `m-1`, rows
    * `m+1..high`) and in `ort(m)`, for [[ortran]] to accumulate `Q`. Norms are
    * formed on the explicitly `scale`-normalized column to avoid overflow.
    */
  private def orthes(n: Int, h: DoubleArray, ort: DoubleArray): Unit =
    val high = n - 1
    var m = 1
    while m <= high - 1 do
      var scale = 0.0
      ort(m) = 0.0
      var i = m
      while i <= high do
        scale += math.abs(h(i * n + (m - 1)))
        i += 1
      if scale != 0.0 then
        var hNorm = 0.0
        i = high
        while i >= m do
          val v = h(i * n + (m - 1)) / scale
          ort(i) = v
          hNorm += v * v
          i -= 1
        val g = if ort(m) >= 0.0 then -math.sqrt(hNorm) else math.sqrt(hNorm)
        hNorm -= ort(m) * g
        ort(m) = ort(m) - g
        // Left transform: H := (I - u uᵀ / hNorm) H over columns m..n-1.
        var j = m
        while j < n do
          var f = 0.0
          i = high
          while i >= m do
            f += ort(i) * h(i * n + j)
            i -= 1
          f = f / hNorm
          i = m
          while i <= high do
            val idx = i * n + j
            h(idx) = h(idx) - f * ort(i)
            i += 1
          j += 1
        // Right transform: H := H (I - u uᵀ / hNorm) over rows 0..high.
        i = 0
        while i <= high do
          var f = 0.0
          j = high
          while j >= m do
            f += ort(j) * h(i * n + j)
            j -= 1
          f = f / hNorm
          j = m
          while j <= high do
            val idx = i * n + j
            h(idx) = h(idx) - f * ort(j)
            j += 1
          i += 1
        ort(m) = scale * ort(m)
        h(m * n + (m - 1)) = scale * g
      m += 1

  /** Accumulate the orthogonal `Q` (EISPACK `ortran`) from the reflectors
    * [[orthes]] left in `h` and `ort`. Returns a fresh `n x n` row-major `Q`.
    * Must run before the reflector cells in `h` are cleared.
    */
  private def ortran(n: Int, h: DoubleArray, ort: DoubleArray): DoubleArray =
    val high = n - 1
    val z = identityRowMajor(n)
    var mp = high - 1
    while mp >= 1 do
      val sub = h(mp * n + (mp - 1))
      if sub != 0.0 then
        var i = mp + 1
        while i <= high do
          ort(i) = h(i * n + (mp - 1))
          i += 1
        var j = mp
        while j <= high do
          var g = 0.0
          i = mp
          while i <= high do
            g += ort(i) * z(i * n + j)
            i += 1
          // Double division (by ort(mp) then by the subdiagonal) avoids underflow.
          g = (g / ort(mp)) / sub
          i = mp
          while i <= high do
            val idx = i * n + j
            z(idx) = z(idx) + g * ort(i)
            i += 1
          j += 1
      mp -= 1
    z

  /** Zero the strictly-below-subdiagonal cells of an `n x n` row-major matrix,
    * turning the reflector-carrying reduction output into a clean Hessenberg.
    */
  private def zeroBelowSubdiagonal(n: Int, h: DoubleArray): Unit =
    var i = 2
    while i < n do
      var j = 0
      while j <= i - 2 do
        h(i * n + j) = 0.0
        j += 1
      i += 1

  /** Francis double-shift QR on the upper-Hessenberg `h` (EISPACK `hqr2`),
    * `low = 0`, `high = n - 1`.
    *
    * Reduces `h` to real-Schur (quasi-triangular) form and fills `wr`/`wi` with
    * the eigenvalues in diagonal order — conjugate pairs adjacent, positive-
    * imaginary member first, with exact conjugate symmetry (`wr(j+1) == wr(j)`,
    * `wi(j+1) == -wi(j)`). When `hasVectors` is true, `z` holds the accumulated
    * `Q` on entry and, on success, the real-Schur packed right eigenvectors on
    * exit (back-substitution through the quasi-triangular `T` followed by the `Q`
    * transform). When false, `z` is untouched and only the eigenvalues are formed.
    *
    * Returns `None` on success or `Some(sweeps)` if the total QR-step budget
    * `maxSweeps` is exhausted before every eigenvalue deflates.
    */
  private def francisSchur(
      n: Int,
      h: DoubleArray,
      wr: DoubleArray,
      wi: DoubleArray,
      z: DoubleArray,
      hasVectors: Boolean,
      maxSweeps: Int
  ): Option[Int] =
    if n == 0 then return None
    inline def hv(i: Int, j: Int): Double = h(i * n + j)
    inline def sh(i: Int, j: Int, v: Double): Unit = h(i * n + j) = v
    inline def zv(i: Int, j: Int): Double = z(i * n + j)
    inline def sz(i: Int, j: Int, v: Double): Unit = z(i * n + j) = v

    // Matrix norm over the (upper-Hessenberg) nonzeros, for the back-substitution
    // tolerances; computed on the input before the iteration destroys h.
    var norm = 0.0
    var kNorm = 0
    var iNorm = 0
    while iNorm < n do
      var j = kNorm
      while j < n do
        norm += math.abs(hv(iNorm, j))
        j += 1
      kNorm = iNorm
      iNorm += 1

    var en = n - 1
    var t = 0.0
    var itn = maxSweeps
    var ierr = -1
    while en >= 0 && ierr < 0 do
      var its = 0
      var deflated = false
      while !deflated && ierr < 0 do
        // Look for a single small sub-diagonal element to split at (-> l).
        var l = en
        var lFound = false
        while !lFound do
          if l == 0 then lFound = true
          else
            var s0 = math.abs(hv(l - 1, l - 1)) + math.abs(hv(l, l))
            if s0 == 0.0 then s0 = norm
            if s0 + math.abs(hv(l, l - 1)) == s0 then lFound = true
            else l -= 1

        val na = en - 1
        val enm2 = en - 2
        var x = hv(en, en)
        if l == en then
          // One real root found.
          sh(en, en, x + t)
          wr(en) = x + t
          wi(en) = 0.0
          en -= 1
          deflated = true
        else
          var y = hv(na, na)
          var w = hv(en, na) * hv(na, en)
          if l == na then
            // Two roots found (a real pair or a complex conjugate pair).
            var p = (y - x) / 2.0
            val q = p * p + w
            var zz = math.sqrt(math.abs(q))
            sh(en, en, x + t)
            val xEn = x + t
            sh(na, na, y + t)
            if q >= 0.0 then
              // Real pair: rotate to make the 2x2 block upper triangular.
              zz = p + sign(zz, p)
              wr(na) = xEn + zz
              wr(en) = wr(na)
              if zz != 0.0 then wr(en) = xEn - w / zz
              wi(na) = 0.0
              wi(en) = 0.0
              val xx = hv(en, na)
              val s2 = math.abs(xx) + math.abs(zz)
              var pRot = xx / s2
              var qRot = zz / s2
              val rRot = math.sqrt(pRot * pRot + qRot * qRot)
              pRot = pRot / rRot
              qRot = qRot / rRot
              var j = na
              while j < n do
                val zt = hv(na, j)
                sh(na, j, qRot * zt + pRot * hv(en, j))
                sh(en, j, qRot * hv(en, j) - pRot * zt)
                j += 1
              var i = 0
              while i <= en do
                val zt = hv(i, na)
                sh(i, na, qRot * zt + pRot * hv(i, en))
                sh(i, en, qRot * hv(i, en) - pRot * zt)
                i += 1
              if hasVectors then
                i = 0
                while i < n do
                  val zt = zv(i, na)
                  sz(i, na, qRot * zt + pRot * zv(i, en))
                  sz(i, en, qRot * zv(i, en) - pRot * zt)
                  i += 1
            else
              // Complex conjugate pair, positive-imaginary member first.
              wr(na) = xEn + p
              wr(en) = xEn + p
              wi(na) = zz
              wi(en) = -zz
            en -= 2
            deflated = true
          else if itn == 0 then
            ierr = en
          else
            if its == 10 || its == 20 then
              // Exceptional shift after stalling, to break periodic cycles.
              t += x
              var i = 0
              while i <= en do
                sh(i, i, hv(i, i) - x)
                i += 1
              val s2 = math.abs(hv(en, na)) + math.abs(hv(na, enm2))
              x = 0.75 * s2
              y = x
              w = -0.4375 * s2 * s2
            its += 1
            itn -= 1
            // Look for two consecutive small sub-diagonal elements (-> m).
            var m = enm2
            var pm = 0.0
            var qm = 0.0
            var rm = 0.0
            var mFound = false
            while !mFound do
              val zz0 = hv(m, m)
              val r0 = x - zz0
              val s0 = y - zz0
              pm = (r0 * s0 - w) / hv(m + 1, m) + hv(m, m + 1)
              qm = hv(m + 1, m + 1) - zz0 - r0 - s0
              rm = hv(m + 2, m + 1)
              val ss = math.abs(pm) + math.abs(qm) + math.abs(rm)
              pm = pm / ss
              qm = qm / ss
              rm = rm / ss
              if m == l then mFound = true
              else
                val tst1 = math.abs(pm) * (math.abs(hv(m - 1, m - 1)) + math.abs(zz0) + math.abs(hv(m + 1, m + 1)))
                val tst2 = tst1 + math.abs(hv(m, m - 1)) * (math.abs(qm) + math.abs(rm))
                if tst2 == tst1 then mFound = true
                else m -= 1

            // Clear the fill introduced two and three rows below the diagonal.
            var i2 = m + 2
            while i2 <= en do
              sh(i2, i2 - 2, 0.0)
              if i2 != m + 2 then sh(i2, i2 - 3, 0.0)
              i2 += 1

            // Double QR step over rows l..en and columns m..en.
            var k = m
            while k <= na do
              val notlas = k != na
              var pStep = 0.0
              var qStep = 0.0
              var rStep = 0.0
              var xAbs = 0.0
              var doStep = true
              if k == m then
                pStep = pm
                qStep = qm
                rStep = rm
              else
                pStep = hv(k, k - 1)
                qStep = hv(k + 1, k - 1)
                rStep = if notlas then hv(k + 2, k - 1) else 0.0
                xAbs = math.abs(pStep) + math.abs(qStep) + math.abs(rStep)
                if xAbs == 0.0 then doStep = false
                else
                  pStep = pStep / xAbs
                  qStep = qStep / xAbs
                  rStep = rStep / xAbs
              if doStep then
                val sStep = sign(math.sqrt(pStep * pStep + qStep * qStep + rStep * rStep), pStep)
                if k != m then sh(k, k - 1, -sStep * xAbs)
                else if l != m then sh(k, k - 1, -hv(k, k - 1))
                pStep = pStep + sStep
                val xr = pStep / sStep
                val yr = qStep / sStep
                val zzr = rStep / sStep
                qStep = qStep / pStep
                rStep = rStep / pStep
                if notlas then
                  var j = k
                  while j < n do
                    val pp = hv(k, j) + qStep * hv(k + 1, j) + rStep * hv(k + 2, j)
                    sh(k, j, hv(k, j) - pp * xr)
                    sh(k + 1, j, hv(k + 1, j) - pp * yr)
                    sh(k + 2, j, hv(k + 2, j) - pp * zzr)
                    j += 1
                  val jmax = math.min(en, k + 3)
                  var i = 0
                  while i <= jmax do
                    val pp = xr * hv(i, k) + yr * hv(i, k + 1) + zzr * hv(i, k + 2)
                    sh(i, k, hv(i, k) - pp)
                    sh(i, k + 1, hv(i, k + 1) - pp * qStep)
                    sh(i, k + 2, hv(i, k + 2) - pp * rStep)
                    i += 1
                  if hasVectors then
                    i = 0
                    while i < n do
                      val pp = xr * zv(i, k) + yr * zv(i, k + 1) + zzr * zv(i, k + 2)
                      sz(i, k, zv(i, k) - pp)
                      sz(i, k + 1, zv(i, k + 1) - pp * qStep)
                      sz(i, k + 2, zv(i, k + 2) - pp * rStep)
                      i += 1
                else
                  var j = k
                  while j < n do
                    val pp = hv(k, j) + qStep * hv(k + 1, j)
                    sh(k, j, hv(k, j) - pp * xr)
                    sh(k + 1, j, hv(k + 1, j) - pp * yr)
                    j += 1
                  val jmax = math.min(en, k + 3)
                  var i = 0
                  while i <= jmax do
                    val pp = xr * hv(i, k) + yr * hv(i, k + 1)
                    sh(i, k, hv(i, k) - pp)
                    sh(i, k + 1, hv(i, k + 1) - pp * qStep)
                    i += 1
                  if hasVectors then
                    i = 0
                    while i < n do
                      val pp = xr * zv(i, k) + yr * zv(i, k + 1)
                      sz(i, k, zv(i, k) - pp)
                      sz(i, k + 1, zv(i, k + 1) - pp * qStep)
                      i += 1
              k += 1
      // (inner while loops back to recompute l unless a root deflated)

    if ierr >= 0 then return Some(maxSweeps - itn)

    // Back-substitute eigenvectors of the quasi-triangular form, then transform
    // by the accumulated Q. Only when vectors were requested and the norm is
    // nonzero (a zero matrix leaves the identity Q as the eigenvectors).
    if hasVectors && norm != 0.0 then
      var en2 = n - 1
      while en2 >= 0 do
        val p = wr(en2)
        val q = wi(en2)
        val na = en2 - 1
        if q == 0.0 then
          // Real eigenvector.
          var m = en2
          sh(en2, en2, 1.0)
          if en2 >= 1 then
            var zzCarry = 0.0
            var sCarry = 0.0
            var i = en2 - 1
            while i >= 0 do
              val w = hv(i, i) - p
              var r = 0.0
              var j = m
              while j <= en2 do
                r += hv(i, j) * hv(j, en2)
                j += 1
              if wi(i) < 0.0 then
                zzCarry = w
                sCarry = r
              else
                m = i
                if wi(i) == 0.0 then
                  var tpiv = w
                  if tpiv == 0.0 then
                    tpiv = norm
                    var loop = true
                    while loop do
                      tpiv = 0.01 * tpiv
                      if !(norm + tpiv > norm) then loop = false
                  sh(i, en2, -r / tpiv)
                else
                  val x2 = hv(i, i + 1)
                  val y2 = hv(i + 1, i)
                  val q2 = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
                  val tval = (x2 * sCarry - zzCarry * r) / q2
                  sh(i, en2, tval)
                  if math.abs(x2) > math.abs(zzCarry) then
                    sh(i + 1, en2, (-r - w * tval) / x2)
                  else
                    sh(i + 1, en2, (-sCarry - y2 * tval) / zzCarry)
                val tabs = math.abs(hv(i, en2))
                if tabs != 0.0 && !(tabs + 1.0 / tabs > tabs) then
                  var j2 = i
                  while j2 <= en2 do
                    sh(j2, en2, hv(j2, en2) / tabs)
                    j2 += 1
              i -= 1
        else if q < 0.0 then
          // Complex eigenvector; the pair's positive-imaginary member is at na.
          var m = na
          if math.abs(hv(en2, na)) > math.abs(hv(na, en2)) then
            sh(na, na, q / hv(en2, na))
            sh(na, en2, -(hv(en2, en2) - p) / hv(en2, na))
          else
            val (cr, ci) = cdiv(0.0, -hv(na, en2), hv(na, na) - p, q)
            sh(na, na, cr)
            sh(na, en2, ci)
          sh(en2, na, 0.0)
          sh(en2, en2, 1.0)
          if en2 >= 2 then
            var zzCarry = 0.0
            var rCarry = 0.0
            var sCarry = 0.0
            var i = en2 - 2
            while i >= 0 do
              val w = hv(i, i) - p
              var ra = 0.0
              var sa = 0.0
              var j = m
              while j <= en2 do
                ra += hv(i, j) * hv(j, na)
                sa += hv(i, j) * hv(j, en2)
                j += 1
              if wi(i) < 0.0 then
                zzCarry = w
                rCarry = ra
                sCarry = sa
              else
                m = i
                if wi(i) == 0.0 then
                  val (cr, ci) = cdiv(-ra, -sa, w, q)
                  sh(i, na, cr)
                  sh(i, en2, ci)
                else
                  val x2 = hv(i, i + 1)
                  val y2 = hv(i + 1, i)
                  var vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
                  val vi = (wr(i) - p) * 2.0 * q
                  if vr == 0.0 && vi == 0.0 then
                    val tst1 =
                      norm * (math.abs(w) + math.abs(q) + math.abs(x2) + math.abs(y2) + math.abs(zzCarry))
                    vr = tst1
                    var loop = true
                    while loop do
                      vr = 0.01 * vr
                      if !(tst1 + vr > tst1) then loop = false
                  val (cr, ci) =
                    cdiv(x2 * rCarry - zzCarry * ra + q * sa, x2 * sCarry - zzCarry * sa - q * ra, vr, vi)
                  sh(i, na, cr)
                  sh(i, en2, ci)
                  if math.abs(x2) > math.abs(zzCarry) + math.abs(q) then
                    sh(i + 1, na, (-ra - w * hv(i, na) + q * hv(i, en2)) / x2)
                    sh(i + 1, en2, (-sa - w * hv(i, en2) - q * hv(i, na)) / x2)
                  else
                    val (cr2, ci2) =
                      cdiv(-rCarry - y2 * hv(i, na), -sCarry - y2 * hv(i, en2), zzCarry, q)
                    sh(i + 1, na, cr2)
                    sh(i + 1, en2, ci2)
                val tabs = math.max(math.abs(hv(i, na)), math.abs(hv(i, en2)))
                if tabs != 0.0 && !(tabs + 1.0 / tabs > tabs) then
                  var j2 = i
                  while j2 <= en2 do
                    sh(j2, na, hv(j2, na) / tabs)
                    sh(j2, en2, hv(j2, en2) / tabs)
                    j2 += 1
              i -= 1
        en2 -= 1

      // Multiply the accumulated Q by the triangular eigenvector matrix stored in
      // h's upper triangle, giving eigenvectors of the original matrix. Columns
      // are formed high-to-low so each read touches only untouched columns.
      var j = n - 1
      while j >= 0 do
        var i = 0
        while i < n do
          var zsum = 0.0
          var kf = 0
          while kf <= j do
            zsum += zv(i, kf) * hv(kf, j)
            kf += 1
          sz(i, j, zsum)
          i += 1
        j -= 1

    None

  /** Complex division `(ar + ai·i) / (br + bi·i)` via the EISPACK `cdiv`
    * scaled formula (scale by `|br| + |bi|`), avoiding intermediate overflow.
    */
  private def cdiv(ar: Double, ai: Double, br: Double, bi: Double): (Double, Double) =
    val s = math.abs(br) + math.abs(bi)
    val ars = ar / s
    val ais = ai / s
    val brs = br / s
    val bis = bi / s
    val d = brs * brs + bis * bis
    ((ars * brs + ais * bis) / d, (ais * brs - ars * bis) / d)

  /** Normalize each real-Schur packed eigenvector to unit 2-norm — a real column
    * to `‖v‖ = 1`, a conjugate pair's two columns so the complex vector
    * `v_re + i·v_im` has unit modulus. Scaling is real and shared across a pair,
    * so the packing invariant and the `A v = λ v` relation are preserved.
    */
  private def normalizeSchurVectors(n: Int, wr: DoubleArray, wi: DoubleArray, z: DoubleArray): Unit =
    var col = 0
    while col < n do
      if wi(col) > 0.0 then
        var nrm = 0.0
        var i = 0
        while i < n do
          val re = z(i * n + col)
          val im = z(i * n + col + 1)
          nrm += re * re + im * im
          i += 1
        nrm = math.sqrt(nrm)
        if nrm > 0.0 then
          i = 0
          while i < n do
            z(i * n + col) = z(i * n + col) / nrm
            z(i * n + col + 1) = z(i * n + col + 1) / nrm
            i += 1
        col += 2
      else
        // Real eigenvalue (wi == 0); wi < 0 cannot lead a column here because the
        // positive-imaginary member always precedes its conjugate.
        var nrm = 0.0
        var i = 0
        while i < n do
          val v = z(i * n + col)
          nrm += v * v
          i += 1
        nrm = math.sqrt(nrm)
        if nrm > 0.0 then
          i = 0
          while i < n do
            z(i * n + col) = z(i * n + col) / nrm
            i += 1
        col += 1

  // ---------------------------------------------------------------------------
  // Shared helpers
  // ---------------------------------------------------------------------------

  /** IEEE machine epsilon for `Double` (2^-52). */
  private inline val Epsilon = 2.220446049250313e-16

  /** `sqrt(a² + b²)` without forming the overflow-prone intermediate, using only
    * correctly-rounded `sqrt` so the result is identical on JVM and Scala.js
    * (unlike `math.hypot`, whose last bit can differ across platforms).
    */
  private def pythag(a: Double, b: Double): Double =
    val absa = math.abs(a)
    val absb = math.abs(b)
    if absa > absb then
      val ratio = absb / absa
      absa * math.sqrt(1.0 + ratio * ratio)
    else if absb == 0.0 then 0.0
    else
      val ratio = absa / absb
      absb * math.sqrt(1.0 + ratio * ratio)

  /** Magnitude of `a` carrying the sign of `b` (Fortran `SIGN`). */
  private inline def sign(a: Double, b: Double): Double =
    if b >= 0.0 then math.abs(a) else -math.abs(a)

  /** A fresh `n x n` row-major array whose lower triangle (and diagonal) is
    * copied from `a` and whose upper triangle mirrors it, so downstream reductions
    * see a genuinely symmetric matrix built from `a`'s lower triangle alone.
    */
  private def symmetrizedLowerRowMajor(a: DMat, n: Int): DoubleArray =
    val out = DoubleArray.alloc(n * n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        val v = a(i, j)
        out(i * n + j) = v
        out(j * n + i) = v
        j += 1
      i += 1
    out

  /** A fresh `n x n` row-major identity array. */
  private def identityRowMajor(n: Int): DoubleArray =
    val out = DoubleArray.alloc(n * n)
    var i = 0
    while i < n do
      out(i * n + i) = 1.0
      i += 1
    out
