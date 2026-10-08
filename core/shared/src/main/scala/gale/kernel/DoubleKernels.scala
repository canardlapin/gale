package gale.kernel

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.PlatformMath.fma

private[gale] object DoubleKernels:
  def ddot(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Double =
    if xStride == 1 && yStride == 1 then
      // Contiguous fast path: four independent accumulators break the reduction's
      // dependency chain so the JIT can pipeline/vectorize the multiply-adds
      // (the F2J trick). Reassociates the sum versus the scalar loop — fine within
      // the library's tolerances, and identical on JVM and Scala.js (shared code).
      var acc0 = 0.0
      var acc1 = 0.0
      var acc2 = 0.0
      var acc3 = 0.0
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < limit do
        acc0 = fma(x(xi), y(yi), acc0)
        acc1 = fma(x(xi + 1), y(yi + 1), acc1)
        acc2 = fma(x(xi + 2), y(yi + 2), acc2)
        acc3 = fma(x(xi + 3), y(yi + 3), acc3)
        xi += 4
        yi += 4
        i += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while i < n do
        acc = fma(x(xi), y(yi), acc)
        xi += 1
        yi += 1
        i += 1
      acc
    else
      var i = 0
      var xi = xOffset
      var yi = yOffset
      var acc = 0.0
      while i < n do
        acc = fma(x(xi), y(yi), acc)
        xi += xStride
        yi += yStride
        i += 1
      acc

  /** Euclidean norm via scaled accumulation (the LAPACK `dnrm2` recurrence).
    *
    * Tracks the running maximum magnitude `scale` and the scaled sum of squares
    * `ssq`, so `sqrt(sum x_i^2)` never forms the intermediate `sum x_i^2` that
    * would overflow for large elements (e.g. 1e155) or underflow to zero for
    * tiny ones (e.g. 1e-170). Ordinary inputs agree with `sqrt(dot(x, x))` to
    * full relative precision.
    */
  def dnrm2(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int
  ): Double =
    if n < 1 then 0.0
    else if n == 1 then math.abs(x(xOffset))
    else
      var scale = 0.0
      var ssq = 1.0
      var i = 0
      var xi = xOffset
      while i < n do
        val value = x(xi)
        if value != 0.0 then
          val abs = math.abs(value)
          if scale < abs then
            val ratio = scale / abs
            ssq = 1.0 + ssq * ratio * ratio
            scale = abs
          else
            val ratio = abs / scale
            ssq += ratio * ratio
        xi += xStride
        i += 1
      val norm = scale * math.sqrt(ssq)
      // The recurrence forms `Inf / Inf` once two infinite entries meet; a NaN
      // here therefore means "some NaN" or "several infinities". Rescan only on
      // that rare path so the result is `NaN` iff an entry is NaN, else `+Inf`.
      if norm.isNaN && !containsNaN(n, x, xOffset, xStride) then Double.PositiveInfinity
      else norm

  /** True when any of the `n` strided elements is NaN. */
  def containsNaN(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Boolean =
    var i = 0
    var xi = xOffset
    while i < n do
      if x(xi).isNaN then return true
      xi += xStride
      i += 1
    false

  def dcopy(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    var i = 0
    var xi = xOffset
    var yi = yOffset
    while i < n do
      y(yi) = x(xi)
      xi += xStride
      yi += yStride
      i += 1

  def daxpy(
      n: Int,
      alpha: Double,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    if xStride == 1 && yStride == 1 then
      // Contiguous fast path, unrolled 4x so independent lanes vectorize.
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < limit do
        y(yi) = fma(alpha, x(xi), y(yi))
        y(yi + 1) = fma(alpha, x(xi + 1), y(yi + 1))
        y(yi + 2) = fma(alpha, x(xi + 2), y(yi + 2))
        y(yi + 3) = fma(alpha, x(xi + 3), y(yi + 3))
        xi += 4
        yi += 4
        i += 4
      while i < n do
        y(yi) = fma(alpha, x(xi), y(yi))
        xi += 1
        yi += 1
        i += 1
    else
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < n do
        y(yi) = fma(alpha, x(xi), y(yi))
        xi += xStride
        yi += yStride
        i += 1

  def dscal(
      n: Int,
      alpha: Double,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int
  ): Unit =
    var i = 0
    var xi = xOffset
    while i < n do
      x(xi) = alpha * x(xi)
      xi += xStride
      i += 1

  def dadd(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int,
      out: DoubleArray,
      outOffset: Int,
      outStride: Int
  ): Unit =
    var i = 0
    var xi = xOffset
    var yi = yOffset
    var oi = outOffset
    while i < n do
      out(oi) = x(xi) + y(yi)
      xi += xStride
      yi += yStride
      oi += outStride
      i += 1

  def dsub(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int,
      out: DoubleArray,
      outOffset: Int,
      outStride: Int
  ): Unit =
    var i = 0
    var xi = xOffset
    var yi = yOffset
    var oi = outOffset
    while i < n do
      out(oi) = x(xi) - y(yi)
      xi += xStride
      yi += yStride
      oi += outStride
      i += 1

  def dgemv(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      rowStride: Int,
      colStride: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      beta: Double,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    val betaIsZero = beta == 0.0
    var row = 0
    var ai = aOffset
    var yi = yOffset
    while row < rows do
      var col = 0
      var aij = ai
      var xj = xOffset
      var acc = 0.0
      while col < cols do
        acc = fma(a(aij), x(xj), acc)
        aij += colStride
        xj += xStride
        col += 1
      y(yi) = if betaIsZero then alpha * acc else fma(alpha, acc, beta * y(yi))
      ai += rowStride
      yi += yStride
      row += 1

  def dgemvRowMajor(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      rowStride: Int,
      x: DoubleArray,
      xOffset: Int,
      beta: Double,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    val betaIsZero = beta == 0.0
    val limit = cols - (cols & 3)
    var row = 0
    var aRow = aOffset
    var yi = yOffset
    while row < rows do
      // Unroll the contiguous inner dot 4x with independent accumulators: both the
      // matrix row and x are unit-stride here, so the lanes vectorize.
      var acc0 = 0.0
      var acc1 = 0.0
      var acc2 = 0.0
      var acc3 = 0.0
      var col = 0
      var ai = aRow
      var xi = xOffset
      while col < limit do
        acc0 = fma(a(ai), x(xi), acc0)
        acc1 = fma(a(ai + 1), x(xi + 1), acc1)
        acc2 = fma(a(ai + 2), x(xi + 2), acc2)
        acc3 = fma(a(ai + 3), x(xi + 3), acc3)
        ai += 4
        xi += 4
        col += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while col < cols do
        acc = fma(a(ai), x(xi), acc)
        ai += 1
        xi += 1
        col += 1
      y(yi) = if betaIsZero then alpha * acc else fma(alpha, acc, beta * y(yi))
      aRow += rowStride
      yi += yStride
      row += 1

  def dgemvColMajor(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      colStride: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      beta: Double,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    var row = 0
    var yi = yOffset
    if beta == 0.0 then
      while row < rows do
        y(yi) = 0.0
        yi += yStride
        row += 1
    else
      while row < rows do
        y(yi) = beta * y(yi)
        yi += yStride
        row += 1

    var col = 0
    var aCol = aOffset
    var xj = xOffset
    while col < cols do
      val scale = alpha * x(xj)
      row = 0
      var ai = aCol
      yi = yOffset
      while row < rows do
        y(yi) = fma(scale, a(ai), y(yi))
        ai += 1
        yi += yStride
        row += 1
      aCol += colStride
      xj += xStride
      col += 1

  /** In-place row-major lower Cholesky. Returns -1 on success or the first
    * failed pivot. Only the lower triangle is read or written.
    */
  def dpotrfLower(n: Int, a: DoubleArray, tolerance: Double): Int =
    var row = 0
    while row < n do
      var col = 0
      while col <= row do
        var sum = a(row * n + col)
        var k = 0
        while k < col do
          sum -= a(row * n + k) * a(col * n + k)
          k += 1
        if !sum.isFinite then return col
        if row == col then
          if sum <= tolerance then return row
          a(row * n + col) = math.sqrt(sum)
        else a(row * n + col) = sum / a(col * n + col)
        col += 1
      row += 1
    -1

  /** In-place triangular solve `T x = b` (`x` holds `b` on entry, the solution on
    * exit). `lower` selects forward vs back substitution; `unit` skips the diagonal
    * division (implicit unit diagonal, so the stored diagonal is never read).
    *
    * Returns the index of the first diagonal whose magnitude is `<= tol` — a
    * singular / rank-deficient pivot — or `-1` on success. With `tol == 0` only an
    * exact zero diagonal fails; `unit == true` never inspects the diagonal and
    * always succeeds. Storage is fully strided, so transposed and submatrix views
    * (e.g. `Lᵀ` via swapped row/column strides) drive the same loop.
    */
  def dtrsv(
      n: Int,
      lower: Boolean,
      unit: Boolean,
      tol: Double,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      aColStride: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int
  ): Int =
    if lower then
      var i = 0
      var xi = xOffset
      var aRow = aOffset
      while i < n do
        var sum = x(xi)
        var aij = aRow
        var xj = xOffset
        var j = 0
        while j < i do
          sum = fma(-a(aij), x(xj), sum)
          aij += aColStride
          xj += xStride
          j += 1
        if unit then x(xi) = sum
        else
          val diag = a(aRow + i * aColStride)
          if math.abs(diag) <= tol then return i
          x(xi) = sum / diag
        aRow += aRowStride
        xi += xStride
        i += 1
      -1
    else
      var i = n - 1
      var xi = xOffset + (n - 1) * xStride
      var aRow = aOffset + (n - 1) * aRowStride
      while i >= 0 do
        var sum = x(xi)
        var aij = aRow + (i + 1) * aColStride
        var xj = xOffset + (i + 1) * xStride
        var j = i + 1
        while j < n do
          sum = fma(-a(aij), x(xj), sum)
          aij += aColStride
          xj += xStride
          j += 1
        if unit then x(xi) = sum
        else
          val diag = a(aRow + i * aColStride)
          if math.abs(diag) <= tol then return i
          x(xi) = sum / diag
        aRow -= aRowStride
        xi -= xStride
        i -= 1
      -1

  /** Above this element count (`64^3`) the row-major path blocks for cache reuse. */
  private inline val GemmBlockThreshold = 262144L
  // 128 measured best with the 4x4 register panel (n=256 gemm: 267 -> 356 ops/s
  // over 64; 256 gained ~3% more but streams B from L2 on larger inputs).
  private inline val GemmBlock = 128

  def dgemm(
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
  ): Unit =
    // Row-major operands (unit column stride) admit a cache-friendly i-k-j loop
    // that streams whole rows of B and C; a blocked variant kicks in for large
    // products. Any strided or transposed layout falls back to the general
    // i-j-k dot-product loop below, which honours arbitrary strides.
    if aColStride == 1 && bColStride == 1 && cColStride == 1 then
      if rows.toLong * cols.toLong * shared.toLong >= GemmBlockThreshold then
        dgemmBlockedRowMajor(
          rows, cols, shared, alpha, a, aOffset, aRowStride, b, bOffset, bRowStride, beta, c, cOffset, cRowStride
        )
      else
        dgemmRowMajor(
          rows, cols, shared, alpha, a, aOffset, aRowStride, b, bOffset, bRowStride, beta, c, cOffset, cRowStride
        )
    else
      dgemmStrided(
        rows, cols, shared, alpha, a, aOffset, aRowStride, aColStride, b, bOffset, bRowStride, bColStride, beta, c,
        cOffset, cRowStride, cColStride
      )

  private def dgemmStrided(
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
  ): Unit =
    val betaIsZero = beta == 0.0
    var row = 0
    var cRow = cOffset
    var aRow = aOffset
    while row < rows do
      var col = 0
      var cij = cRow
      var bCol = bOffset
      while col < cols do
        var k = 0
        var aik = aRow
        var bkj = bCol
        var acc = 0.0
        while k < shared do
          acc = fma(a(aik), b(bkj), acc)
          aik += aColStride
          bkj += bRowStride
          k += 1
        c(cij) = if betaIsZero then alpha * acc else fma(alpha, acc, beta * c(cij))
        cij += cColStride
        bCol += bColStride
        col += 1
      cRow += cRowStride
      aRow += aRowStride
      row += 1

  /** Scale C in place by `beta` (or zero it when `beta == 0`) before an i-k-j
    * accumulation. Row-major with unit column stride, so each row is contiguous.
    */
  private def scaleRowMajor(
      rows: Int,
      cols: Int,
      beta: Double,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    if beta == 1.0 then ()
    else
      val zero = beta == 0.0
      var row = 0
      var cRow = cOffset
      while row < rows do
        var col = 0
        var cij = cRow
        while col < cols do
          c(cij) = if zero then 0.0 else beta * c(cij)
          cij += 1
          col += 1
        cRow += cRowStride
        row += 1

  private def dgemmRowMajor(
      rows: Int,
      cols: Int,
      shared: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      b: DoubleArray,
      bOffset: Int,
      bRowStride: Int,
      beta: Double,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    val assign = beta == 0.0 && shared > 0
    if !assign then scaleRowMajor(rows, cols, beta, c, cOffset, cRowStride)
    gemmPanel(
      0, rows, 0, cols, 0, shared, alpha, a, aOffset, aRowStride, b, bOffset,
      bRowStride, c, cOffset, cRowStride, assign
    )

  /** Register-blocked accumulation micro-kernel over the tile `C[iStart:iEnd,
    * jStart:jEnd] += alpha·A[·, kStart:kEnd]·B[kStart:kEnd, ·]` (row-major, unit
    * column stride).
    *
    * The `4×4` interior holds a `C` tile in '''16 accumulators across the whole
    * k-loop''': each k-step loads 4 values of `A` and 4 of `B` and issues 16 fused
    * multiply-adds, so `C` is read/written once per tile rather than once per
    * k-step — cutting the `C` memory traffic that limited the plain unroll-and-jam
    * by a factor of the k-extent. The `0–3` leftover columns (still 4 rows at a
    * time) and `0–3` leftover rows fall to unrolled/scalar tails.
    */
  private def gemmPanel(
      iStart: Int,
      iEnd: Int,
      jStart: Int,
      jEnd: Int,
      kStart: Int,
      kEnd: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      b: DoubleArray,
      bOffset: Int,
      bRowStride: Int,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int,
      assign: Boolean
  ): Unit =
    val iMain = iStart + ((iEnd - iStart) & ~3)
    val jMain = jStart + ((jEnd - jStart) & ~3)
    var i = iStart
    while i < iMain do
      val aRow0 = aOffset + i * aRowStride
      val aRow1 = aRow0 + aRowStride
      val aRow2 = aRow1 + aRowStride
      val aRow3 = aRow2 + aRowStride
      val cRow0 = cOffset + i * cRowStride
      val cRow1 = cRow0 + cRowStride
      val cRow2 = cRow1 + cRowStride
      val cRow3 = cRow2 + cRowStride
      var j = jStart
      while j < jMain do
        var c00 = 0.0; var c01 = 0.0; var c02 = 0.0; var c03 = 0.0
        var c10 = 0.0; var c11 = 0.0; var c12 = 0.0; var c13 = 0.0
        var c20 = 0.0; var c21 = 0.0; var c22 = 0.0; var c23 = 0.0
        var c30 = 0.0; var c31 = 0.0; var c32 = 0.0; var c33 = 0.0
        var k = kStart
        var bRow = bOffset + kStart * bRowStride
        while k < kEnd do
          val a0 = a(aRow0 + k)
          val a1 = a(aRow1 + k)
          val a2 = a(aRow2 + k)
          val a3 = a(aRow3 + k)
          val b0 = b(bRow + j)
          val b1 = b(bRow + j + 1)
          val b2 = b(bRow + j + 2)
          val b3 = b(bRow + j + 3)
          c00 = fma(a0, b0, c00); c01 = fma(a0, b1, c01); c02 = fma(a0, b2, c02); c03 = fma(a0, b3, c03)
          c10 = fma(a1, b0, c10); c11 = fma(a1, b1, c11); c12 = fma(a1, b2, c12); c13 = fma(a1, b3, c13)
          c20 = fma(a2, b0, c20); c21 = fma(a2, b1, c21); c22 = fma(a2, b2, c22); c23 = fma(a2, b3, c23)
          c30 = fma(a3, b0, c30); c31 = fma(a3, b1, c31); c32 = fma(a3, b2, c32); c33 = fma(a3, b3, c33)
          bRow += bRowStride
          k += 1
        storeTile4(c, cRow0, j, alpha, c00, c01, c02, c03, assign)
        storeTile4(c, cRow1, j, alpha, c10, c11, c12, c13, assign)
        storeTile4(c, cRow2, j, alpha, c20, c21, c22, c23, assign)
        storeTile4(c, cRow3, j, alpha, c30, c31, c32, c33, assign)
        j += 4
      // Leftover 0–3 columns, still four rows at a time (one B load, four FMAs).
      while j < jEnd do
        var s0 = 0.0
        var s1 = 0.0
        var s2 = 0.0
        var s3 = 0.0
        var k = kStart
        var bRow = bOffset + kStart * bRowStride
        while k < kEnd do
          val bv = b(bRow + j)
          s0 = fma(a(aRow0 + k), bv, s0)
          s1 = fma(a(aRow1 + k), bv, s1)
          s2 = fma(a(aRow2 + k), bv, s2)
          s3 = fma(a(aRow3 + k), bv, s3)
          bRow += bRowStride
          k += 1
        val x0 = cRow0 + j; c(x0) = if assign then alpha * s0 else fma(alpha, s0, c(x0))
        val x1 = cRow1 + j; c(x1) = if assign then alpha * s1 else fma(alpha, s1, c(x1))
        val x2 = cRow2 + j; c(x2) = if assign then alpha * s2 else fma(alpha, s2, c(x2))
        val x3 = cRow3 + j; c(x3) = if assign then alpha * s3 else fma(alpha, s3, c(x3))
        j += 1
      i += 4
    // Leftover 0–3 rows: scalar dot over the k-extent.
    while i < iEnd do
      val aRow = aOffset + i * aRowStride
      val cRow = cOffset + i * cRowStride
      var j = jStart
      while j < jEnd do
        var s = 0.0
        var k = kStart
        var bRow = bOffset + kStart * bRowStride
        while k < kEnd do
          s = fma(a(aRow + k), b(bRow + j), s)
          bRow += bRowStride
          k += 1
        val idx = cRow + j
        c(idx) = if assign then alpha * s else fma(alpha, s, c(idx))
        j += 1
      i += 1

  /** Store one row of a register tile: `C[row, j..j+3] += alpha·(t0..t3)`. */
  private inline def storeTile4(
      c: DoubleArray,
      cRow: Int,
      j: Int,
      alpha: Double,
      t0: Double,
      t1: Double,
      t2: Double,
      t3: Double,
      assign: Boolean
  ): Unit =
    val i0 = cRow + j
    c(i0) = if assign then alpha * t0 else fma(alpha, t0, c(i0))
    val i1 = cRow + j + 1
    c(i1) = if assign then alpha * t1 else fma(alpha, t1, c(i1))
    val i2 = cRow + j + 2
    c(i2) = if assign then alpha * t2 else fma(alpha, t2, c(i2))
    val i3 = cRow + j + 3
    c(i3) = if assign then alpha * t3 else fma(alpha, t3, c(i3))

  private def dgemmBlockedRowMajor(
      rows: Int,
      cols: Int,
      shared: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      b: DoubleArray,
      bOffset: Int,
      bRowStride: Int,
      beta: Double,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    if beta != 0.0 || shared == 0 then
      scaleRowMajor(rows, cols, beta, c, cOffset, cRowStride)
    var ii = 0
    while ii < rows do
      val iMax = math.min(ii + GemmBlock, rows)
      var kk = 0
      while kk < shared do
        val kMax = math.min(kk + GemmBlock, shared)
        var jj = 0
        while jj < cols do
          val jMax = math.min(jj + GemmBlock, cols)
          gemmPanel(
            ii, iMax, jj, jMax, kk, kMax, alpha, a, aOffset, aRowStride, b,
            bOffset, bRowStride, c, cOffset, cRowStride, assign = beta == 0.0 && kk == 0
          )
          jj += GemmBlock
        kk += GemmBlock
      ii += GemmBlock

  /** Symmetric rank-k product `C := AᵀA` for a '''row-major''' `A` (`m × k`, unit
    * column stride), writing the full symmetric `k × k` result into `c` (row-major,
    * unit column stride). '''Assign-only:''' there is no `alpha`/`beta`; every
    * output cell is overwritten, so the input contents of `C` are ignored.
    *
    * A general gemm computing `Aᵀ·A` sees `Aᵀ` column-strided. This kernel instead
    * visits the upper triangle in `4×4` output tiles, holds each tile in registers
    * across the full `m` reduction, and mirrors it once. Compared with a row-wise
    * rank-1 update, each `C` cell is written once rather than `m` times. Full tiles
    * use 16 independent FMA accumulators; boundary tiles use the same scalar
    * reduction with exact triangular bounds.
    */
  def dsyrkRowMajor(
      m: Int,
      k: Int,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    var ib = 0
    while ib < k do
      val iEnd = math.min(ib + 4, k)
      var jb = ib
      while jb < k do
        val jEnd = math.min(jb + 4, k)
        if iEnd - ib == 4 && jEnd - jb == 4 then
          var c00 = 0.0; var c01 = 0.0; var c02 = 0.0; var c03 = 0.0
          var c10 = 0.0; var c11 = 0.0; var c12 = 0.0; var c13 = 0.0
          var c20 = 0.0; var c21 = 0.0; var c22 = 0.0; var c23 = 0.0
          var c30 = 0.0; var c31 = 0.0; var c32 = 0.0; var c33 = 0.0
          var l = 0
          var aRow = aOffset
          while l < m do
            val a0 = a(aRow + ib)
            val a1 = a(aRow + ib + 1)
            val a2 = a(aRow + ib + 2)
            val a3 = a(aRow + ib + 3)
            val b0 = a(aRow + jb)
            val b1 = a(aRow + jb + 1)
            val b2 = a(aRow + jb + 2)
            val b3 = a(aRow + jb + 3)
            c00 = fma(a0, b0, c00); c01 = fma(a0, b1, c01); c02 = fma(a0, b2, c02); c03 = fma(a0, b3, c03)
            c10 = fma(a1, b0, c10); c11 = fma(a1, b1, c11); c12 = fma(a1, b2, c12); c13 = fma(a1, b3, c13)
            c20 = fma(a2, b0, c20); c21 = fma(a2, b1, c21); c22 = fma(a2, b2, c22); c23 = fma(a2, b3, c23)
            c30 = fma(a3, b0, c30); c31 = fma(a3, b1, c31); c32 = fma(a3, b2, c32); c33 = fma(a3, b3, c33)
            aRow += aRowStride
            l += 1

          val r0 = cOffset + ib * cRowStride + jb
          val r1 = r0 + cRowStride
          val r2 = r1 + cRowStride
          val r3 = r2 + cRowStride
          c(r0) = c00; c(r0 + 1) = c01; c(r0 + 2) = c02; c(r0 + 3) = c03
          if jb == ib then
            c(r1 + 1) = c11; c(r1 + 2) = c12; c(r1 + 3) = c13
            c(r2 + 2) = c22; c(r2 + 3) = c23
            c(r3 + 3) = c33
          else
            c(r1) = c10; c(r1 + 1) = c11; c(r1 + 2) = c12; c(r1 + 3) = c13
            c(r2) = c20; c(r2 + 1) = c21; c(r2 + 2) = c22; c(r2 + 3) = c23
            c(r3) = c30; c(r3 + 1) = c31; c(r3 + 2) = c32; c(r3 + 3) = c33
        else
          var i = ib
          while i < iEnd do
            var j = math.max(jb, i)
            while j < jEnd do
              var acc = 0.0
              var l = 0
              var aRow = aOffset
              while l < m do
                acc = fma(a(aRow + i), a(aRow + j), acc)
                aRow += aRowStride
                l += 1
              c(cOffset + i * cRowStride + j) = acc
              j += 1
            i += 1
        jb += 4
      ib += 4
    // Mirror the computed upper triangle into the lower.
    var i = 1
    while i < k do
      val cRowI = cOffset + i * cRowStride
      var j = 0
      while j < i do
        c(cRowI + j) = c(cOffset + j * cRowStride + i)
        j += 1
      i += 1

  // ---------------------------------------------------------------------------
  // Reductions and elementwise numerics (public facade: DVec/DMat members and
  // gale.linalg.Numerics). Every kernel has a unit-stride fast path and a
  // strided fallback; NaN and infinities follow IEEE arithmetic unless a kernel
  // documents a different propagation rule.
  // ---------------------------------------------------------------------------

  /** Sum of `n` elements. The contiguous path uses four accumulators, so it
    * reassociates relative to a left-to-right loop (shared JVM/JS code).
    */
  def dsum(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    if xStride == 1 then
      var acc0 = 0.0
      var acc1 = 0.0
      var acc2 = 0.0
      var acc3 = 0.0
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      while i < limit do
        acc0 += x(xi)
        acc1 += x(xi + 1)
        acc2 += x(xi + 2)
        acc3 += x(xi + 3)
        xi += 4
        i += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while i < n do
        acc += x(xi)
        xi += 1
        i += 1
      acc
    else
      var acc = 0.0
      var i = 0
      var xi = xOffset
      while i < n do
        acc += x(xi)
        xi += xStride
        i += 1
      acc

  /** `sum x_i / divisor`, each term divided before it is added: the
    * overflow-safe fallback for a mean whose plain sum is not finite. With
    * `divisor = n` every term is at most `max|x_i| / n`, so infinities and NaN
    * keep their IEEE results while finite entries cannot overflow except by
    * rounding at `±Double.MaxValue` (see [[clampMeanOverflow]]).
    */
  def dsumDivided(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, divisor: Double): Double =
    var acc = 0.0
    var i = 0
    var xi = xOffset
    while i < n do
      acc += x(xi) / divisor
      xi += xStride
      i += 1
    acc

  /** Whether any of the `n` entries is `±Infinity`. */
  def dcontainsInfinity(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Boolean =
    var i = 0
    var xi = xOffset
    while i < n do
      if x(xi).isInfinite then return true
      xi += xStride
      i += 1
    false

  /** A mean lies between the smallest and largest entry, so a `±Inf` from
    * [[dsumDivided]] over finite entries is rounding at the edge of the range
    * (for example `n` copies of `Double.MaxValue`): return `±Double.MaxValue`.
    */
  inline def clampMeanOverflow(mean: Double, hasInfiniteEntry: => Boolean): Double =
    if mean.isInfinite && !hasInfiniteEntry then math.copySign(Double.MaxValue, mean) else mean

  /** Overflow-safe mean of one strided line whose plain sum was not finite. */
  def dmeanOfNonFiniteSum(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    clampMeanOverflow(dsumDivided(n, x, xOffset, xStride, n.toDouble), dcontainsInfinity(n, x, xOffset, xStride))

  /** Sum of absolute values (the vector 1-norm); `0.0` when `n == 0`. */
  def dasum(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    if xStride == 1 then
      var acc0 = 0.0
      var acc1 = 0.0
      var acc2 = 0.0
      var acc3 = 0.0
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      while i < limit do
        acc0 += math.abs(x(xi))
        acc1 += math.abs(x(xi + 1))
        acc2 += math.abs(x(xi + 2))
        acc3 += math.abs(x(xi + 3))
        xi += 4
        i += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while i < n do
        acc += math.abs(x(xi))
        xi += 1
        i += 1
      acc
    else
      var acc = 0.0
      var i = 0
      var xi = xOffset
      while i < n do
        acc += math.abs(x(xi))
        xi += xStride
        i += 1
      acc

  /** Largest absolute value (the vector ∞-norm); `0.0` when `n == 0`. A NaN
    * element makes the result NaN (`math.max` propagates NaN on JVM and JS).
    */
  def damax(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    if xStride == 1 then
      var m0 = 0.0
      var m1 = 0.0
      var m2 = 0.0
      var m3 = 0.0
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      while i < limit do
        m0 = math.max(m0, math.abs(x(xi)))
        m1 = math.max(m1, math.abs(x(xi + 1)))
        m2 = math.max(m2, math.abs(x(xi + 2)))
        m3 = math.max(m3, math.abs(x(xi + 3)))
        xi += 4
        i += 4
      var m = math.max(math.max(m0, m1), math.max(m2, m3))
      while i < n do
        m = math.max(m, math.abs(x(xi)))
        xi += 1
        i += 1
      m
    else
      var m = 0.0
      var i = 0
      var xi = xOffset
      while i < n do
        m = math.max(m, math.abs(x(xi)))
        xi += xStride
        i += 1
      m

  /** Index of the first maximum, `-1` when `n == 0`. The first NaN wins: its
    * index is returned as soon as it is seen, so a NaN anywhere propagates to
    * the value read back at the returned index.
    */
  def dmaxIndex(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Int =
    extremeIndex(n, x, xOffset, xStride, Double.NegativeInfinity)(_ > _)

  /** Index of the first minimum, `-1` when `n == 0`; the first NaN wins. */
  def dminIndex(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Int =
    extremeIndex(n, x, xOffset, xStride, Double.PositiveInfinity)(_ < _)

  /** Shared first-occurrence search. `better(v, best)` must be a strict IEEE
    * comparison (false for NaN). The contiguous path keeps four lanes, each the
    * first occurrence within its residue class; the combine step breaks value
    * ties toward the smaller index, so the result equals a sequential scan.
    */
  private inline def extremeIndex(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, worst: Double)(
      inline better: (Double, Double) => Boolean
  ): Int =
    if n <= 0 then -1
    else if xStride == 1 then
      var b0 = worst
      var b1 = worst
      var b2 = worst
      var b3 = worst
      var i0 = Int.MaxValue
      var i1 = Int.MaxValue
      var i2 = Int.MaxValue
      var i3 = Int.MaxValue
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      var nan = -1
      while i < limit && nan < 0 do
        val v0 = x(xi)
        val v1 = x(xi + 1)
        val v2 = x(xi + 2)
        val v3 = x(xi + 3)
        if better(v0, b0) then
          b0 = v0; i0 = i
        else if v0 != v0 then nan = i
        if nan < 0 then
          if better(v1, b1) then
            b1 = v1; i1 = i + 1
          else if v1 != v1 then nan = i + 1
        if nan < 0 then
          if better(v2, b2) then
            b2 = v2; i2 = i + 2
          else if v2 != v2 then nan = i + 2
        if nan < 0 then
          if better(v3, b3) then
            b3 = v3; i3 = i + 3
          else if v3 != v3 then nan = i + 3
        xi += 4
        i += 4
      if nan >= 0 then nan
      else
        var best = b0
        var bestIndex = i0
        if better(b1, best) || (b1 == best && i1 < bestIndex) then
          best = b1; bestIndex = i1
        if better(b2, best) || (b2 == best && i2 < bestIndex) then
          best = b2; bestIndex = i2
        if better(b3, best) || (b3 == best && i3 < bestIndex) then
          best = b3; bestIndex = i3
        while i < n && nan < 0 do
          val v = x(xi)
          if better(v, best) then
            best = v; bestIndex = i
          else if v != v then nan = i
          xi += 1
          i += 1
        if nan >= 0 then nan
        // Every element equalled `worst` (e.g. all -Inf for a max): the first wins.
        else if bestIndex == Int.MaxValue then 0
        else bestIndex
    else
      var best = x(xOffset)
      if best != best then 0
      else
        var bestIndex = 0
        var nan = -1
        var i = 1
        var xi = xOffset + xStride
        while i < n && nan < 0 do
          val v = x(xi)
          if better(v, best) then
            best = v; bestIndex = i
          else if v != v then nan = i
          xi += xStride
          i += 1
        if nan >= 0 then nan else bestIndex

  /** Sum of squares, `sum x_i^2` (four fma accumulators when contiguous). */
  def dsumsq(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    ddot(n, x, xOffset, xStride, x, xOffset, xStride)

  /** Scaled sum of squares `sum (x_i / scale)^2` for a positive finite `scale`. */
  def dsumsqScaled(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, scale: Double): Double =
    var acc = 0.0
    var i = 0
    var xi = xOffset
    while i < n do
      val v = x(xi) / scale
      acc = fma(v, v, acc)
      xi += xStride
      i += 1
    acc

  // The optimistic unscaled sum of squares is trusted when it is finite and at
  // least this large. Each square that underflows loses at most 2^-1075, so even
  // Int.MaxValue such losses stay ~1e-35 relative to the total: negligible.
  // Smaller totals (tiny entries, or exactly zero) and overflow rescan scaled.
  private val FrobeniusTrustedMin = 1e-280

  /** Frobenius norm of a `rows×cols` strided block, overflow- and underflow-safe.
    *
    * One optimistic pass forms the plain fma sum of squares. Only when that
    * total is non-finite, zero, or below a safe floor does a second pair of
    * passes find the largest magnitude `m` and sum `(a/m)^2`, giving
    * `m * sqrt(...)`. NaN anywhere gives NaN; otherwise an infinite entry gives
    * `+Inf`; an empty block gives `0.0`. A block whose storage is one
    * contiguous run (row- or column-major), or a single row or column, is
    * reduced line by line in one kernel call.
    */
  def dnrmFrobenius(
      rows: Int,
      cols: Int,
      x: DoubleArray,
      xOffset: Int,
      rowStride: Int,
      colStride: Int
  ): Double =
    if rows == 0 || cols == 0 then 0.0
    else
      val contiguous = (colStride == 1 && rowStride == cols) || (rowStride == 1 && colStride == rows)
      val byRows = if rows == 1 then true else if cols == 1 then false else colStride <= rowStride
      val lines = if contiguous then 1 else if byRows then rows else cols
      val lineLength = if contiguous then rows * cols else if byRows then cols else rows
      val lineStep = if byRows then rowStride else colStride
      val elementStep = if contiguous then 1 else if byRows then colStride else rowStride
      var ssq = 0.0
      var line = 0
      while line < lines do
        ssq += dsumsq(lineLength, x, xOffset + line * lineStep, elementStep)
        line += 1
      if ssq.isFinite && ssq >= FrobeniusTrustedMin then math.sqrt(ssq)
      else
        var m = 0.0
        line = 0
        while line < lines do
          m = math.max(m, damax(lineLength, x, xOffset + line * lineStep, elementStep))
          line += 1
        if m == 0.0 || m.isNaN || m.isInfinite then m
        else
          var scaled = 0.0
          line = 0
          while line < lines do
            scaled += dsumsqScaled(lineLength, x, xOffset + line * lineStep, elementStep, m)
            line += 1
          m * math.sqrt(scaled)

  /** `y_i := f(x_i)`; `y` may be `x` itself (same offset and stride). */
  inline def dmapInto(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  )(inline f: Double => Double): Unit =
    if xStride == 1 && yStride == 1 then
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < limit do
        y(yi) = f(x(xi))
        y(yi + 1) = f(x(xi + 1))
        y(yi + 2) = f(x(xi + 2))
        y(yi + 3) = f(x(xi + 3))
        xi += 4
        yi += 4
        i += 4
      while i < n do
        y(yi) = f(x(xi))
        xi += 1
        yi += 1
        i += 1
    else
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < n do
        y(yi) = f(x(xi))
        xi += xStride
        yi += yStride
        i += 1

  def dexpInto(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    dmapInto(n, x, xOffset, xStride, y, yOffset, yStride)(math.exp)

  def dlogInto(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    dmapInto(n, x, xOffset, xStride, y, yOffset, yStride)(math.log)

  def dlog1pInto(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    dmapInto(n, x, xOffset, xStride, y, yOffset, yStride)(math.log1p)

  def dexpm1Into(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    dmapInto(n, x, xOffset, xStride, y, yOffset, yStride)(math.expm1)

  /** Logistic sigmoid `1 / (1 + exp(-x))`, evaluated without overflow: with
    * `t = exp(-|x|) <= 1` it is `1/(1+t)` for `x >= 0` and `t/(1+t)` otherwise.
    * `sigmoid(+Inf) = 1`, `sigmoid(-Inf) = 0`, `sigmoid(NaN) = NaN`.
    */
  inline def sigmoid(v: Double): Double =
    val t = math.exp(-math.abs(v))
    if v >= 0.0 then 1.0 / (1.0 + t) else t / (1.0 + t)

  def dsigmoidInto(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    dmapInto(n, x, xOffset, xStride, y, yOffset, yStride)(v => sigmoid(v))

  /** `sum exp(x_i - shift)` with four accumulators when contiguous. */
  def dsumExpShifted(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, shift: Double): Double =
    if xStride == 1 then
      var acc0 = 0.0
      var acc1 = 0.0
      var acc2 = 0.0
      var acc3 = 0.0
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      while i < limit do
        acc0 += math.exp(x(xi) - shift)
        acc1 += math.exp(x(xi + 1) - shift)
        acc2 += math.exp(x(xi + 2) - shift)
        acc3 += math.exp(x(xi + 3) - shift)
        xi += 4
        i += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while i < n do
        acc += math.exp(x(xi) - shift)
        xi += 1
        i += 1
      acc
    else
      var acc = 0.0
      var i = 0
      var xi = xOffset
      while i < n do
        acc += math.exp(x(xi) - shift)
        xi += xStride
        i += 1
      acc

  /** The value at [[dmaxIndex]]: the first NaN if any, else the maximum. */
  private def maxValue(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    x(xOffset + dmaxIndex(n, x, xOffset, xStride) * xStride)

  /** `log(sum exp(x_i))` by the two-pass max shift `m + log(sum exp(x_i - m))`.
    *
    * Empty gives `-Inf` (the log of an empty sum). Otherwise NaN anywhere gives
    * NaN; else any `+Inf` gives `+Inf`; else all `-Inf` gives `-Inf`. Finite
    * inputs never overflow: the shifted maximum contributes exactly `exp(0) = 1`.
    */
  def dlogSumExp(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    if n == 0 then Double.NegativeInfinity
    else
      val m = maxValue(n, x, xOffset, xStride)
      if m.isNaN || m.isInfinite then m
      else m + math.log(dsumExpShifted(n, x, xOffset, xStride, m))

  /** `y := softmax(x) = exp(x - m) / sum exp(x - m)`. If the maximum `m` is not
    * finite (NaN anywhere, any `+Inf`, or all `-Inf`) every output is NaN.
    */
  def dsoftmaxInto(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    if n > 0 then
      val m = maxValue(n, x, xOffset, xStride)
      if m.isNaN || m.isInfinite then fillStrided(n, Double.NaN, y, yOffset, yStride)
      else
        val total = dexpShiftedSumInto(n, x, xOffset, xStride, m, y, yOffset, yStride)
        // Divide rather than multiply by `1 / total`: one rounding per entry.
        dmapInto(n, y, yOffset, yStride, y, yOffset, yStride)(v => v / total)

  /** `y_i := exp(x_i - shift)`, returning `sum y_i` from the same pass (four
    * accumulators when both operands are contiguous).
    */
  def dexpShiftedSumInto(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      shift: Double,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Double =
    if xStride == 1 && yStride == 1 then
      var acc0 = 0.0
      var acc1 = 0.0
      var acc2 = 0.0
      var acc3 = 0.0
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < limit do
        val e0 = math.exp(x(xi) - shift)
        val e1 = math.exp(x(xi + 1) - shift)
        val e2 = math.exp(x(xi + 2) - shift)
        val e3 = math.exp(x(xi + 3) - shift)
        y(yi) = e0
        y(yi + 1) = e1
        y(yi + 2) = e2
        y(yi + 3) = e3
        acc0 += e0
        acc1 += e1
        acc2 += e2
        acc3 += e3
        xi += 4
        yi += 4
        i += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while i < n do
        val e = math.exp(x(xi) - shift)
        y(yi) = e
        acc += e
        xi += 1
        yi += 1
        i += 1
      acc
    else
      var acc = 0.0
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < n do
        val e = math.exp(x(xi) - shift)
        y(yi) = e
        acc += e
        xi += xStride
        yi += yStride
        i += 1
      acc

  /** `y := x - logSumExp(x)`, evaluated as `(x - m) - log(sum exp(x - m))`.
    * Non-finite maxima give an all-NaN result, as in [[dsoftmaxInto]].
    */
  def dlogSoftmaxInto(
      n: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    if n > 0 then
      val m = maxValue(n, x, xOffset, xStride)
      if m.isNaN || m.isInfinite then fillStrided(n, Double.NaN, y, yOffset, yStride)
      else
        val logTotal = math.log(dsumExpShifted(n, x, xOffset, xStride, m))
        dmapInto(n, x, xOffset, xStride, y, yOffset, yStride)(v => (v - m) - logTotal)

  private def fillStrided(n: Int, value: Double, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    var i = 0
    var yi = yOffset
    while i < n do
      y(yi) = value
      yi += yStride
      i += 1
