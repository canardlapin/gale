package gale.spectral

import gale.linalg.DMat
import gale.platform.DoubleArray

/** The pure dense bidiagonal-SVD kernel — the piece whose absence deferred full
  * dense SVD out of v0.3.5 (parity § 3; seam S7 of
  * `docs/spectral-backend-boundary.md`).
  *
  * '''Algorithm''' (Golub–Kahan–Reinsch): Householder bidiagonalization of `A`
  * to an upper bidiagonal `B`, then implicit-shift QR on `B` with the Wilkinson
  * shift from the trailing 2×2, deflating on a negligible superdiagonal at the
  * scale-aware threshold `ε·‖B‖` (`‖B‖ = max_i(|d_i| + |e_i|)`, the standard
  * Golub–Reinsch `anorm` — consistent with `docs/numerical-contract.md`'s
  * "scale-aware tests" rule). A wide input (`m < n`) is handled by running on
  * `Aᵀ` and swapping `U`/`V` on output, so the working orientation is always
  * tall. Vectors are accumulated when requested; the values-only mode runs the
  * identical scalar recurrences and skips only the rotation applications.
  *
  * '''Output contract''' (raw, canonicalized by the `Svds` facade): `min(m, n)`
  * non-negative singular values in '''kernel order''' (not sorted — the facade
  * imposes the descending layout, exactly as it does for a backend's
  * [[RawSvd]]), with economy factors `U` (`m×k`) and `Vᵀ` (`k×n`),
  * `k = min(m, n)`; empty (`m×0` / `0×n`) factors when values-only.
  *
  * Deterministic per platform: fixed sweep order, no RNG, and only
  * correctly-rounded primitives ([[pythag]] instead of `math.hypot`, whose last
  * bit differs across platforms). Exhaustion of the per-value sweep budget is a
  * typed [[SvdKernelFailure.DidNotConverge]], mirroring
  * [[DenseSpectralKernels.SpectralKernelFailure]] (pathological in practice —
  * the classic bound is ~2 sweeps per singular value).
  */
private[gale] object DenseSvdKernel:

  /** Typed kernel failure, mirroring
    * [[DenseSpectralKernels.SpectralKernelFailure]]: the implicit-QR sweep
    * budget ran out. Carries the total sweeps performed.
    */
  enum SvdKernelFailure:
    case DidNotConverge(iterations: Int)

  /** IEEE machine epsilon for `Double` (2^-52). */
  private inline val Epsilon = 2.220446049250313e-16

  /** Implicit-QR sweeps allowed per singular value before the typed failure.
    * The classical Golub–Reinsch budget is 30; doubled for slack — convergence
    * is cubic once the shift locks on, so a well-formed matrix never gets near
    * either bound.
    */
  private inline val MaxSweepsPerValue = 60

  /** Full/economy dense SVD of `a` (any `m×n`, both dimensions positive —
    * validated by the facade, not here). Returns the raw factor carrier; the
    * facade sorts, measures residuals, and builds the sealed [[SVD]].
    */
  def svd(a: DMat, wantVectors: Boolean): Either[SvdKernelFailure, RawSvd] =
    if a.rows >= a.cols then svdTall(a, wantVectors)
    else
      // A = (Aᵀ)ᵀ: with Aᵀ = U' Σ V'ᵀ, A = V' Σ U'ᵀ — swap the factors back.
      svdTall(a.t, wantVectors).map: raw =>
        RawSvd(raw.sigma, raw.vt.t, raw.u.t)

  // ===========================================================================
  // Core (tall orientation, m >= n)
  // ===========================================================================

  private def svdTall(a: DMat, wantVectors: Boolean): Either[SvdKernelFailure, RawSvd] =
    val m = a.rows
    val n = a.cols
    val u = a.toDoubleArrayCopyRowMajor // m×n row-major; becomes economy U in place
    val w = DoubleArray.alloc(n)        // singular values
    val rv1 = DoubleArray.alloc(n)      // superdiagonal workspace
    val v = if wantVectors then DoubleArray.alloc(n * n) else DoubleArray.alloc(0)

    val anorm = bidiagonalize(m, n, u, w, rv1)
    if wantVectors then
      accumulateRight(n, u, v, rv1)
      accumulateLeft(m, n, u, w)
    diagonalize(m, n, u, w, rv1, v, anorm, wantVectors) match
      case Some(failure) => Left(failure)
      case None =>
        val sigma = gale.linalg.DVec.tabulate(n)(i => w(i))
        if wantVectors then
          Right(RawSvd(sigma, DMat.fromDoubleArrayOwned(m, n, u), DMat.fromDoubleArrayOwned(n, n, v).t))
        else Right(RawSvd(sigma, DMat.zeros(m, 0), DMat.zeros(0, n)))

  /** Householder bidiagonalization of the m×n row-major `u` in place: on return
    * `w` holds the diagonal, `rv1` the superdiagonal (`rv1(0) = 0`), and `u`
    * the scaled Householder reflectors the accumulation phases expand. Returns
    * `anorm = max_i(|w_i| + |rv1_i|)`, the deflation scale.
    */
  private def bidiagonalize(m: Int, n: Int, u: DoubleArray, w: DoubleArray, rv1: DoubleArray): Double =
    var g = 0.0
    var scale = 0.0
    var anorm = 0.0
    var i = 0
    while i < n do
      val l = i + 1
      rv1(i) = scale * g
      g = 0.0
      scale = 0.0
      if i < m then
        var k = i
        while k < m do
          scale += math.abs(u(k * n + i))
          k += 1
        if scale != 0.0 then
          var s = 0.0
          k = i
          while k < m do
            u(k * n + i) = u(k * n + i) / scale
            s += u(k * n + i) * u(k * n + i)
            k += 1
          val f = u(i * n + i)
          g = -sign(math.sqrt(s), f)
          val h = f * g - s
          u(i * n + i) = f - g
          var j = l
          while j < n do
            var dot = 0.0
            k = i
            while k < m do
              dot += u(k * n + i) * u(k * n + j)
              k += 1
            val fj = dot / h
            k = i
            while k < m do
              u(k * n + j) = u(k * n + j) + fj * u(k * n + i)
              k += 1
            j += 1
          k = i
          while k < m do
            u(k * n + i) = u(k * n + i) * scale
            k += 1
      w(i) = scale * g
      g = 0.0
      scale = 0.0
      if i < m && i != n - 1 then
        var k = l
        while k < n do
          scale += math.abs(u(i * n + k))
          k += 1
        if scale != 0.0 then
          var s = 0.0
          k = l
          while k < n do
            u(i * n + k) = u(i * n + k) / scale
            s += u(i * n + k) * u(i * n + k)
            k += 1
          val f = u(i * n + l)
          g = -sign(math.sqrt(s), f)
          val h = f * g - s
          u(i * n + l) = f - g
          k = l
          while k < n do
            rv1(k) = u(i * n + k) / h
            k += 1
          var j = l
          while j < m do
            var dot = 0.0
            k = l
            while k < n do
              dot += u(j * n + k) * u(i * n + k)
              k += 1
            k = l
            while k < n do
              u(j * n + k) = u(j * n + k) + dot * rv1(k)
              k += 1
            j += 1
          k = l
          while k < n do
            u(i * n + k) = u(i * n + k) * scale
            k += 1
      anorm = math.max(anorm, math.abs(w(i)) + math.abs(rv1(i)))
      i += 1
    anorm

  /** Expand the right-hand Householder reflectors stored in `u`'s rows into the
    * n×n row-major `v` (columns are right singular vectors of the bidiagonal
    * reduction), walking i from n−1 down so each reflector is applied to the
    * already-accumulated trailing block.
    */
  private def accumulateRight(n: Int, u: DoubleArray, v: DoubleArray, rv1: DoubleArray): Unit =
    var g = 0.0
    var l = 0
    var i = n - 1
    while i >= 0 do
      if i < n - 1 then
        if g != 0.0 then
          // Double division (by u(i,l) then g) avoids underflow, per Golub–Reinsch.
          var j = l
          while j < n do
            v(j * n + i) = (u(i * n + j) / u(i * n + l)) / g
            j += 1
          j = l
          while j < n do
            var s = 0.0
            var k = l
            while k < n do
              s += u(i * n + k) * v(k * n + j)
              k += 1
            k = l
            while k < n do
              v(k * n + j) = v(k * n + j) + s * v(k * n + i)
              k += 1
            j += 1
        var j = l
        while j < n do
          v(i * n + j) = 0.0
          v(j * n + i) = 0.0
          j += 1
      v(i * n + i) = 1.0
      g = rv1(i)
      l = i
      i -= 1

  /** Expand the left-hand Householder reflectors into the economy `U` (m×n, in
    * place over the reflector storage), walking i from n−1 down.
    */
  private def accumulateLeft(m: Int, n: Int, u: DoubleArray, w: DoubleArray): Unit =
    var i = n - 1
    while i >= 0 do
      val l = i + 1
      var g = w(i)
      var j = l
      while j < n do
        u(i * n + j) = 0.0
        j += 1
      if g != 0.0 then
        g = 1.0 / g
        j = l
        while j < n do
          var s = 0.0
          var k = l
          while k < m do
            s += u(k * n + i) * u(k * n + j)
            k += 1
          val f = (s / u(i * n + i)) * g
          k = i
          while k < m do
            u(k * n + j) = u(k * n + j) + f * u(k * n + i)
            k += 1
          j += 1
        var k = i
        while k < m do
          u(k * n + i) = u(k * n + i) * g
          k += 1
      else
        var k = i
        while k < m do
          u(k * n + i) = 0.0
          k += 1
      u(i * n + i) = u(i * n + i) + 1.0
      i -= 1

  /** Implicit-shift QR on the bidiagonal `(w, rv1)`, rotations accumulated into
    * `u`/`v` when `wantVectors`. Deflation and cancellation both test against
    * the scale-aware `ε·anorm`. On success every `w(i) ≥ 0` (a converged
    * negative value flips sign along with its `v` column). Returns the typed
    * failure when a value exhausts its sweep budget.
    */
  private def diagonalize(
      m: Int,
      n: Int,
      u: DoubleArray,
      w: DoubleArray,
      rv1: DoubleArray,
      v: DoubleArray,
      anorm: Double,
      wantVectors: Boolean
  ): Option[SvdKernelFailure] =
    val tol = Epsilon * anorm
    var totalSweeps = 0
    var k = n - 1
    while k >= 0 do
      var its = 0
      var converged = false
      while !converged do
        // Split search: l is the top of the active block; flag marks whether
        // rv1(l) still couples it to w(l-1) (then w(l-1) is negligible and the
        // coupling must be rotated away). rv1(0) is exactly 0, so l never
        // underruns and nm = l-1 is only read when l >= 1.
        var flag = true
        var l = k
        var nm = l - 1
        var searching = true
        while searching do
          nm = l - 1
          if math.abs(rv1(l)) <= tol then
            flag = false
            searching = false
          else if math.abs(w(nm)) <= tol then searching = false
          else l -= 1
        if flag then
          // Cancellation: chase rv1(l) to zero with Givens rotations against
          // the negligible w(nm), applied to U's columns (nm, l..k).
          var c = 0.0
          var s = 1.0
          var i = l
          var stop = false
          while i <= k && !stop do
            val f = s * rv1(i)
            rv1(i) = c * rv1(i)
            if math.abs(f) <= tol then stop = true
            else
              val g = w(i)
              var h = pythag(f, g)
              w(i) = h
              h = 1.0 / h
              c = g * h
              s = -f * h
              if wantVectors then
                var j = 0
                while j < m do
                  val y = u(j * n + nm)
                  val z = u(j * n + i)
                  u(j * n + nm) = y * c + z * s
                  u(j * n + i) = z * c - y * s
                  j += 1
            i += 1
        val z0 = w(k)
        if l == k then
          // Converged: enforce σ ≥ 0 by flipping the paired right vector.
          if z0 < 0.0 then
            w(k) = -z0
            if wantVectors then
              var j = 0
              while j < n do
                v(j * n + k) = -v(j * n + k)
                j += 1
          converged = true
        else if its >= MaxSweepsPerValue then return Some(SvdKernelFailure.DidNotConverge(totalSweeps))
        else
          its += 1
          totalSweeps += 1
          // Wilkinson shift from the trailing 2×2 of Bᵀ B, in shifted form.
          var x = w(l)
          val nm2 = k - 1
          var y = w(nm2)
          var g = rv1(nm2)
          var h = rv1(k)
          var f = ((y - z0) * (y + z0) + (g - h) * (g + h)) / (2.0 * h * y)
          g = pythag(f, 1.0)
          f = ((x - z0) * (x + z0) + h * ((y / (f + sign(g, f))) - h)) / x
          // One implicit QR sweep l..k: chase the bulge with paired rotations.
          var c = 1.0
          var s = 1.0
          var j = l
          while j <= nm2 do
            val i = j + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            var z = pythag(f, h)
            rv1(j) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = g * c - x * s
            h = y * s
            y = y * c
            if wantVectors then
              var jj = 0
              while jj < n do
                val xv = v(jj * n + j)
                val zv = v(jj * n + i)
                v(jj * n + j) = xv * c + zv * s
                v(jj * n + i) = zv * c - xv * s
                jj += 1
            z = pythag(f, h)
            w(j) = z
            if z != 0.0 then
              val zi = 1.0 / z
              c = f * zi
              s = h * zi
            f = c * g + s * y
            x = c * y - s * g
            if wantVectors then
              var jj = 0
              while jj < m do
                val yv = u(jj * n + j)
                val zv = u(jj * n + i)
                u(jj * n + j) = yv * c + zv * s
                u(jj * n + i) = zv * c - yv * s
                jj += 1
            j += 1
          rv1(l) = 0.0
          rv1(k) = f
          w(k) = x
      k -= 1
    None

  // ===========================================================================
  // Shared numeric helpers (duplicated from DenseSpectralKernels, which keeps
  // them private; both stay tiny and platform-identical)
  // ===========================================================================

  /** `sqrt(a² + b²)` without the overflow-prone intermediate, using only
    * correctly-rounded operations so the result is identical on JVM and
    * Scala.js (unlike `math.hypot`).
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
