package gale.spectral

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*

/** Cuppen divide and conquer for the eigenvectors of a symmetric tridiagonal
  * matrix (the LAPACK `dstedc`/`dlaed0`–`dlaed4` scheme), used by the dense
  * symmetric eigensolver when eigenvectors are wanted and the order is at
  * least [[MinOrder]].
  *
  * `T` is torn at the middle by a rank-one modification, the halves are solved
  * recursively (implicit QL at blocks of at most [[Leaf]] rows), and each merge
  * solves `diag(D) + ρ z zᵀ`:
  *
  *   - '''Deflation''' (`dlaed2`): a negligible `ρ|zᵢ|` or a pair of poles close
  *     enough that a Givens rotation zeroes one weight keeps that pair as it is.
  *   - '''Secular equation''' (`dlaed4`): each remaining root is found in its
  *     pole interval, relative to the nearer pole, by the middle-way rational
  *     iteration with a bisection safeguard, so every `dᵢ − λⱼ` is accurate.
  *   - '''Eigenvectors''' (Gu–Eisenstat): the weights are recomputed from the
  *     roots (the Löwner formula) before forming `ẑᵢ / (dᵢ − λⱼ)`, which keeps
  *     the merged vectors numerically orthogonal.
  *   - '''Update''': the new vectors are products of the merge coefficients
  *     with the old vectors, formed as two panel products that skip the zero
  *     blocks of `diag(Q₁, Q₂)`.
  *
  * Vectors are stored transposed: row `j` of the `n × n` region `zt` is the
  * eigenvector of eigenvalue `d(j)`. All storage comes from the caller: the
  * `n²` vector region, a [[doubleScratch]]-sized double region and an
  * [[indexScratch]]-sized index region. The routine is deterministic: it uses
  * only correctly rounded `+ − × ÷ √` (no fused multiply-add), so the JVM and
  * Scala.js produce the same bits.
  */
private[spectral] object TridiagonalDivideConquer:
  /** Smallest order solved by divide and conquer (QL below). Measured
    * crossover of the dense kernel: at 48 divide and conquer is 2.6% faster
    * on Zen 5 and level on an M3; at 32-40 QL is level or faster.
    */
  inline val MinOrder = 48

  /** Largest block solved directly by QL inside the recursion. */
  private inline val Leaf = 25

  /** Rows of new eigenvectors formed per panel product in a merge. */
  private inline val Panel = 32

  /** Length-`n` double vectors in the scratch region besides the two panels. */
  private inline val VectorRegions = 9

  private inline val Epsilon = 2.220446049250313e-16

  /** Secular iterations allowed per root before the solve is declared failed. */
  private inline val MaxSecularIterations = 400

  /** Double scratch, in elements, for an order-`n` solve: an `n²` packing
    * region plus `(9 + 3·32)·n` for the per-merge vectors and panels.
    */
  def doubleScratch(n: Int): Long =
    n.toLong * n.toLong + (VectorRegions + 3 * Panel).toLong * n.toLong

  /** Index scratch, in elements, for an order-`n` solve. */
  def indexScratch(n: Int): Long = 7L * n.toLong

  /** Eigen-decompose `T` with diagonal `d(0..n-1)` and off-diagonal
    * `off(offOffset + i) = T(i, i+1)`, `i < n-1` (overwritten). On success `d`
    * holds the eigenvalues ascending and `zt(ztOffset + j·n + c)` component `c`
    * of eigenvector `j`; returns `0`. Returns the failing iteration count
    * (positive) when a leaf QL or a secular root fails to converge. Requires
    * every entry of `T` to be finite.
    */
  def solve(
      n: Int,
      d: DoubleArray,
      off: DoubleArray,
      offOffset: Int,
      zt: DoubleArray,
      ztOffset: Int,
      work: DoubleArray,
      workOffset: Int,
      iwork: IndexArray,
      iworkOffset: Int,
      maxSweeps: Int
  ): Int =
    var i = 0
    while i < n * n do
      zt(ztOffset + i) = 0.0
      i += 1
    // Scale to unit max-norm (`dstedc` scales by the norm too); the
    // eigenvectors are scale-invariant and the eigenvalues are scaled back.
    var norm = 0.0
    i = 0
    while i < n do
      norm = math.max(norm, math.abs(d(i)))
      if i < n - 1 then norm = math.max(norm, math.abs(off(offOffset + i)))
      i += 1
    if norm == 0.0 then
      i = 0
      while i < n do
        zt(ztOffset + i * n + i) = 1.0
        i += 1
      return 0
    i = 0
    while i < n do
      d(i) = d(i) / norm
      if i < n - 1 then off(offOffset + i) = off(offOffset + i) / norm
      i += 1
    val solver = new Solver(n, d, off, offOffset, zt, ztOffset, work, workOffset, iwork, iworkOffset, maxSweeps)
    val failure = solver.split(0, n)
    if failure == 0 then
      i = 0
      while i < n do
        d(i) = d(i) * norm
        i += 1
    failure

  private final class Solver(
      n: Int,
      d: DoubleArray,
      off: DoubleArray,
      offOffset: Int,
      zt: DoubleArray,
      ztOffset: Int,
      work: DoubleArray,
      workOffset: Int,
      iwork: IndexArray,
      iworkOffset: Int,
      maxSweeps: Int
  ):
    // Double regions (each length n unless noted).
    private val pack = workOffset // n² packed old vectors
    private val vdd = pack + n * n // merged poles, deflated values
    private val vz = vdd + n // rank-one weights
    private val vdl = vz + n // non-deflated poles
    private val vw = vdl + n // non-deflated weights
    private val vwhat = vw + n // Gu–Eisenstat weights
    private val vtau = vwhat + n // root offset from its origin pole
    private val vlam = vtau + n // roots, then deflated values
    private val vtmp = vlam + n // secular deltas, leaf off-diagonal
    private val vtop = vtmp + n // Panel x n, top coefficient panel
    private val vbot = vtop + Panel * n // Panel x n, bottom coefficient panel
    private val vout = vbot + Panel * n // Panel x n, new vectors before placement
    // Index regions (each length n).
    private val iindx = iworkOffset // pole order
    private val ityp = iindx + n // column type 1 (top), 3 (bottom), 2 (both), 4 (deflated)
    private val ikeep = ityp + n // non-deflated, ascending pole order
    private val idefl = ikeep + n // deflated, ascending value order
    private val iorg = idefl + n // origin pole of each root
    private val isel = iorg + n // top/bottom panel membership
    private val ipos = isel + n // final position of each merged pair

    /** Solve the block `[lo, lo + size)`; 0 or a failing iteration count. */
    def split(lo: Int, size: Int): Int =
      if size <= Leaf then leaf(lo, size)
      else
        val n1 = size / 2
        val beta = off(offOffset + lo + n1 - 1)
        val rho = math.abs(beta)
        d(lo + n1 - 1) = d(lo + n1 - 1) - rho
        d(lo + n1) = d(lo + n1) - rho
        val left = split(lo, n1)
        if left != 0 then left
        else
          val right = split(lo + n1, size - n1)
          if right != 0 then right
          else merge(lo, n1, size - n1, rho, if beta < 0.0 then -1.0 else 1.0)

    /** Implicit QL (EISPACK `tql2`, the dense solver's sweep) on a leaf block,
      * rotating rows of the identity-initialized block of `zt`, then a
      * selection sort into ascending order.
      */
    private def leaf(lo: Int, size: Int): Int =
      val base = ztOffset + lo * n + lo
      var i = 0
      while i < size do
        zt(base + i * n + i) = 1.0
        work(vtmp + i) = if i < size - 1 then off(offOffset + lo + i) else 0.0
        i += 1
      var l = 0
      while l < size do
        var iter = 0
        var continue = true
        while continue do
          var m = l
          var found = false
          while m < size - 1 && !found do
            val dd = math.abs(d(lo + m)) + math.abs(d(lo + m + 1))
            // `T` is scaled to unit max-norm, so 1 is the norm-scaled test.
            if math.abs(work(vtmp + m)) <= Epsilon * math.max(dd, 1.0) then found = true
            else m += 1
          if m == l then continue = false
          else
            if iter == maxSweeps then return math.max(iter, 1)
            iter += 1
            var g = (d(lo + l + 1) - d(lo + l)) / (2.0 * work(vtmp + l))
            var r = pythag(g, 1.0)
            g = d(lo + m) - d(lo + l) + work(vtmp + l) / (g + (if g >= 0.0 then math.abs(r) else -math.abs(r)))
            var s = 1.0
            var c = 1.0
            var p = 0.0
            var innerZero = false
            var iBt = m - 1
            while iBt >= l && !innerZero do
              val f = s * work(vtmp + iBt)
              val b = c * work(vtmp + iBt)
              r = pythag(f, g)
              work(vtmp + iBt + 1) = r
              if r == 0.0 then
                d(lo + iBt + 1) = d(lo + iBt + 1) - p
                work(vtmp + m) = 0.0
                innerZero = true
              else
                s = f / r
                c = g / r
                g = d(lo + iBt + 1) - p
                r = (d(lo + iBt) - g) * s + 2.0 * c * b
                p = s * r
                d(lo + iBt + 1) = g + p
                g = c * r - b
                val row0 = base + iBt * n
                val row1 = row0 + n
                var k = 0
                while k < size do
                  val f2 = zt(row1 + k)
                  zt(row1 + k) = s * zt(row0 + k) + c * f2
                  zt(row0 + k) = c * zt(row0 + k) - s * f2
                  k += 1
                iBt -= 1
            if !(innerZero && iBt >= l) then
              d(lo + l) = d(lo + l) - p
              work(vtmp + l) = g
              work(vtmp + m) = 0.0
        l += 1
      i = 0
      while i < size - 1 do
        var k = i
        var p = d(lo + i)
        var j = i + 1
        while j < size do
          if d(lo + j) < p then
            k = j
            p = d(lo + j)
          j += 1
        if k != i then
          d(lo + k) = d(lo + i)
          d(lo + i) = p
          swapRows(base + i * n, base + k * n, size)
        i += 1
      0

    private def swapRows(row0: Int, row1: Int, length: Int): Unit =
      var k = 0
      while k < length do
        val t = zt(row0 + k)
        zt(row0 + k) = zt(row1 + k)
        zt(row1 + k) = t
        k += 1

    /** Merge the solved halves `[lo, lo+n1)` and `[lo+n1, lo+n1+n2)`: the
      * eigenpairs of `diag(D₁, D₂) + ρ z zᵀ` with `z` the last components of
      * the first half's vectors and (sign-adjusted) the first components of
      * the second half's.
      */
    private def merge(lo: Int, n1: Int, n2: Int, rhoTear: Double, sgn: Double): Int =
      val kk = n1 + n2
      val base = ztOffset + lo * n + lo
      // Weights scaled by 1/√2 (both halves are unit vectors), ρ by 2.
      val invSqrt2 = 1.0 / math.sqrt(2.0)
      var c = 0
      while c < kk do
        work(vdd + c) = d(lo + c)
        work(vz + c) =
          if c < n1 then zt(base + c * n + n1 - 1) * invSqrt2
          else sgn * zt(base + c * n + n1) * invSqrt2
        iwork(ityp + c) = if c < n1 then 1 else 3
        c += 1
      val rho = 2.0 * rhoTear
      // Merge the two ascending halves into pole order.
      var a = 0
      var b = n1
      var t = 0
      while t < kk do
        if b >= kk || (a < n1 && work(vdd + a) <= work(vdd + b)) then
          iwork(iindx + t) = a
          a += 1
        else
          iwork(iindx + t) = b
          b += 1
        t += 1
      var dmax = 0.0
      var zmax = 0.0
      c = 0
      while c < kk do
        dmax = math.max(dmax, math.abs(work(vdd + c)))
        zmax = math.max(zmax, math.abs(work(vz + c)))
        c += 1
      val tol = 8.0 * Epsilon * math.max(dmax, zmax)

      // Deflation (dlaed2): `keep` gets the surviving pairs in pole order,
      // `defl` the deflated ones.
      var k = 0
      var nd = 0
      if rho * zmax <= tol then
        t = 0
        while t < kk do
          iwork(ityp + iwork(iindx + t)) = 4
          iwork(idefl + nd) = iwork(iindx + t)
          nd += 1
          t += 1
      else
        var pj = -1
        t = 0
        while t < kk do
          val nj = iwork(iindx + t)
          if rho * math.abs(work(vz + nj)) <= tol then
            iwork(ityp + nj) = 4
            iwork(idefl + nd) = nj
            nd += 1
          else if pj < 0 then pj = nj
          else
            var s = work(vz + pj)
            var cc = work(vz + nj)
            val tau = pythag(cc, s)
            val gap = work(vdd + nj) - work(vdd + pj)
            cc = cc / tau
            s = -s / tau
            if math.abs(gap * cc * s) <= tol then
              // Rotate the pair so the weight of `pj` vanishes; `pj` deflates.
              work(vz + nj) = tau
              work(vz + pj) = 0.0
              if iwork(ityp + nj) != iwork(ityp + pj) then iwork(ityp + nj) = 2
              iwork(ityp + pj) = 4
              val rowP = base + pj * n
              val rowN = base + nj * n
              var col = 0
              while col < kk do
                val x = zt(rowP + col)
                val y = zt(rowN + col)
                zt(rowP + col) = cc * x + s * y
                zt(rowN + col) = cc * y - s * x
                col += 1
              val dp = work(vdd + pj)
              val dn = work(vdd + nj)
              work(vdd + pj) = dp * cc * cc + dn * s * s
              work(vdd + nj) = dp * s * s + dn * cc * cc
              iwork(idefl + nd) = pj
              nd += 1
            else
              iwork(ikeep + k) = pj
              k += 1
            pj = nj
          t += 1
        if pj >= 0 then
          iwork(ikeep + k) = pj
          k += 1
      // Deflated values ascending (insertion sort; rotations may perturb order).
      t = 1
      while t < nd do
        val v = iwork(idefl + t)
        val key = work(vdd + v)
        var s = t - 1
        while s >= 0 && work(vdd + iwork(idefl + s)) > key do
          iwork(idefl + s + 1) = iwork(idefl + s)
          s -= 1
        iwork(idefl + s + 1) = v
        t += 1

      // Pack the old vectors' nonzero parts so the block can be overwritten:
      // top panel rows (type 1/2, columns [0,n1)), bottom panel rows (type 2/3,
      // columns [n1,kk)), then the deflated rows in full.
      var m1 = 0
      var m2 = 0
      t = 0
      while t < k do
        val ty = iwork(ityp + iwork(ikeep + t))
        if ty != 3 then m1 += 1
        if ty != 1 then m2 += 1
        t += 1
      val packTop = pack
      val packBottom = packTop + m1 * n1
      val packDeflated = packBottom + m2 * n2
      var r1 = 0
      var r2 = 0
      t = 0
      while t < k do
        val src = iwork(ikeep + t)
        val ty = iwork(ityp + src)
        // isel: bit 1 = in top panel, bit 2 = in bottom panel.
        var membership = 0
        if ty != 3 then
          copy(zt, base + src * n, work, packTop + r1 * n1, n1)
          r1 += 1
          membership |= 1
        if ty != 1 then
          copy(zt, base + src * n + n1, work, packBottom + r2 * n2, n2)
          r2 += 1
          membership |= 2
        iwork(isel + t) = membership
        t += 1
      t = 0
      while t < nd do
        copy(zt, base + iwork(idefl + t) * n, work, packDeflated + t * kk, kk)
        t += 1

      // Secular roots, accumulating the Gu–Eisenstat products as they come.
      t = 0
      while t < k do
        work(vdl + t) = work(vdd + iwork(ikeep + t))
        work(vw + t) = work(vz + iwork(ikeep + t))
        work(vwhat + t) = 1.0
        t += 1
      var j = 0
      while j < k do
        val failed = secular(k, j, rho)
        if failed != 0 then return failed
        val org = iwork(iorg + j)
        val tauJ = work(vtau + j)
        work(vlam + j) = work(vdl + org) + tauJ
        var i = 0
        while i < k do
          val delta = (work(vdl + i) - work(vdl + org)) - tauJ
          work(vwhat + i) =
            if i == j then work(vwhat + i) * delta
            else work(vwhat + i) * (delta / (work(vdl + i) - work(vdl + j)))
          i += 1
        j += 1
      var i = 0
      while i < k do
        work(vwhat + i) = math.copySign(math.sqrt(-work(vwhat + i)), work(vw + i))
        i += 1

      // Final ascending positions: merge roots (ascending) with deflated values.
      a = 0
      b = 0
      t = 0
      while t < kk do
        if b >= nd || (a < k && work(vlam + a) <= work(vdd + iwork(idefl + b))) then
          iwork(ipos + a) = t
          a += 1
        else
          iwork(ipos + k + b) = t
          b += 1
        t += 1
      // Deflated rows and values go straight to their positions.
      t = 0
      while t < nd do
        val p = iwork(ipos + k + t)
        copy(work, packDeflated + t * kk, zt, base + p * n, kk)
        d(lo + p) = work(vdd + iwork(idefl + t))
        t += 1
      // New vectors, Panel rows at a time: coefficients ẑᵢ/(dᵢ − λⱼ),
      // normalized, gathered into top/bottom panels, times the packed rows.
      j = 0
      while j < k do
        val rows = math.min(Panel, k - j)
        var q = 0
        while q < rows do
          val root = j + q
          val org = iwork(iorg + root)
          val tauJ = work(vtau + root)
          var norm2 = 0.0
          i = 0
          while i < k do
            val v = work(vwhat + i) / ((work(vdl + i) - work(vdl + org)) - tauJ)
            work(vtmp + i) = v
            norm2 += v * v
            i += 1
          val nrm = math.sqrt(norm2)
          var top = 0
          var bottom = 0
          i = 0
          while i < k do
            val v = work(vtmp + i) / nrm
            val membership = iwork(isel + i)
            if (membership & 1) != 0 then
              work(vtop + q * m1 + top) = v
              top += 1
            if (membership & 2) != 0 then
              work(vbot + q * m2 + bottom) = v
              bottom += 1
            i += 1
          d(lo + iwork(ipos + root)) = work(vlam + root)
          q += 1
        // Both halves of the new rows at once, then each row to pos(j + q).
        if m1 > 0 then
          DenseSpectralKernels.multiplyRowMajor(rows, n1, m1, work, vtop, m1, work, packTop, n1, work, vout, kk)
        if m2 > 0 then
          DenseSpectralKernels.multiplyRowMajor(rows, n2, m2, work, vbot, m2, work, packBottom, n2, work, vout + n1, kk)
        q = 0
        while q < rows do
          val dest = base + iwork(ipos + j + q) * n
          if m1 > 0 then copy(work, vout + q * kk, zt, dest, n1) else zero(dest, n1)
          if m2 > 0 then copy(work, vout + q * kk + n1, zt, dest + n1, n2) else zero(dest + n1, n2)
          q += 1
        j += rows
      0

    private def zero(at: Int, length: Int): Unit =
      var c = 0
      while c < length do
        zt(at + c) = 0.0
        c += 1

    private def copy(src: DoubleArray, from: Int, dst: DoubleArray, to: Int, length: Int): Unit =
      var c = 0
      while c < length do
        dst(to + c) = src(from + c)
        c += 1

    /** Root `j` of `1/ρ + Σᵢ wᵢ² / (dlᵢ − λ) = 0` (`dl` strictly ascending,
      * `ρ > 0`), as `λ = dl(org) + τ` relative to the nearer pole `org`
      * (`dlaed4`'s middle-way iteration on the two bracketing poles, safeguarded
      * by bisection). Stores `org` and `τ`; returns 0, or the iteration count
      * on failure.
      */
    private def secular(k: Int, j: Int, rho: Double): Int =
      val rhoinv = 1.0 / rho
      if k == 1 then
        iwork(iorg + j) = 0
        work(vtau + j) = rho * work(vw) * work(vw)
        return 0
      val last = j == k - 1
      var org = j
      var lo = 0.0
      var hi = 0.0
      if last then
        var zz = 0.0
        var i = 0
        while i < k do
          zz += work(vw + i) * work(vw + i)
          i += 1
        org = k - 1
        hi = rho * zz
      else
        val mid = 0.5 * (work(vdl + j + 1) - work(vdl + j))
        var f = rhoinv
        var i = 0
        while i < k do
          f += work(vw + i) * work(vw + i) / ((work(vdl + i) - work(vdl + j)) - mid)
          i += 1
        if f >= 0.0 then hi = mid
        else
          org = j + 1
          lo = -mid
      val dorg = work(vdl + org)
      val pa = if last then k - 2 else j
      val pb = pa + 1
      var tau = 0.5 * (lo + hi)
      var iter = 0
      while true do
        var psi = 0.0
        var dpsi = 0.0
        var phi = 0.0
        var dphi = 0.0
        var erretm = 0.0
        var i = 0
        while i < k do
          val di = (work(vdl + i) - dorg) - tau
          work(vtmp + i) = di
          val ratio = work(vw + i) / di
          val term = work(vw + i) * ratio
          if i <= pa then
            psi += term
            dpsi += ratio * ratio
          else
            phi += term
            dphi += ratio * ratio
          erretm += math.abs(term)
          i += 1
        val fw = rhoinv + psi + phi
        val dw = dpsi + dphi
        erretm = 8.0 * erretm + 2.0 * rhoinv + 3.0 * math.abs(fw) + math.abs(tau) * dw
        if math.abs(fw) <= Epsilon * erretm || hi - lo <= 2.0 * Epsilon * math.max(math.abs(lo), math.abs(hi))
        then
          iwork(iorg + j) = org
          work(vtau + j) = tau
          return 0
        if iter == MaxSecularIterations then return iter
        if fw <= 0.0 then lo = math.max(lo, tau) else hi = math.min(hi, tau)
        val da = work(vtmp + pa)
        val db = work(vtmp + pb)
        val a = (da + db) * fw - da * db * dw
        val b = da * db * fw
        var eta = 0.0
        if last then
          val cc = math.abs(fw - da * dpsi - db * dphi)
          eta =
            if cc == 0.0 then hi - tau
            else if a >= 0.0 then (a + math.sqrt(math.abs(a * a - 4.0 * b * cc))) / (2.0 * cc)
            else 2.0 * b / (a - math.sqrt(math.abs(a * a - 4.0 * b * cc)))
          if fw * eta > 0.0 then eta = -fw / dw
        else
          val wa = work(vw + pa) / da
          val wb = work(vw + pb) / db
          val cc =
            if org == j then fw - db * dw - (work(vdl + pa) - work(vdl + pb)) * wa * wa
            else fw - da * dw - (work(vdl + pb) - work(vdl + pa)) * wb * wb
          eta =
            if cc == 0.0 then (if a == 0.0 then 0.0 else b / a)
            else if a <= 0.0 then (a - math.sqrt(math.abs(a * a - 4.0 * b * cc))) / (2.0 * cc)
            else 2.0 * b / (a + math.sqrt(math.abs(a * a - 4.0 * b * cc)))
          if fw * eta >= 0.0 then eta = -fw / dw
        val next = tau + eta
        tau = if next >= hi || next <= lo || next.isNaN then 0.5 * (lo + hi) else next
        iter += 1
      0

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
