package gale.kernel

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

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
        acc0 += x(xi) * y(yi)
        acc1 += x(xi + 1) * y(yi + 1)
        acc2 += x(xi + 2) * y(yi + 2)
        acc3 += x(xi + 3) * y(yi + 3)
        xi += 4
        yi += 4
        i += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while i < n do
        acc += x(xi) * y(yi)
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
        acc += x(xi) * y(yi)
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
      scale * math.sqrt(ssq)

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
        y(yi) = alpha * x(xi) + y(yi)
        y(yi + 1) = alpha * x(xi + 1) + y(yi + 1)
        y(yi + 2) = alpha * x(xi + 2) + y(yi + 2)
        y(yi + 3) = alpha * x(xi + 3) + y(yi + 3)
        xi += 4
        yi += 4
        i += 4
      while i < n do
        y(yi) = alpha * x(xi) + y(yi)
        xi += 1
        yi += 1
        i += 1
    else
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < n do
        y(yi) = alpha * x(xi) + y(yi)
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
        acc += a(aij) * x(xj)
        aij += colStride
        xj += xStride
        col += 1
      y(yi) = if betaIsZero then alpha * acc else alpha * acc + beta * y(yi)
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
        acc0 += a(ai) * x(xi)
        acc1 += a(ai + 1) * x(xi + 1)
        acc2 += a(ai + 2) * x(xi + 2)
        acc3 += a(ai + 3) * x(xi + 3)
        ai += 4
        xi += 4
        col += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while col < cols do
        acc += a(ai) * x(xi)
        ai += 1
        xi += 1
        col += 1
      y(yi) = if betaIsZero then alpha * acc else alpha * acc + beta * y(yi)
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
        y(yi) = y(yi) + scale * a(ai)
        ai += 1
        yi += yStride
        row += 1
      aCol += colStride
      xj += xStride
      col += 1

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
          sum -= a(aij) * x(xj)
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
          sum -= a(aij) * x(xj)
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
          acc += a(aik) * b(bkj)
          aik += aColStride
          bkj += bRowStride
          k += 1
        c(cij) = if betaIsZero then alpha * acc else alpha * acc + beta * c(cij)
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
    scaleRowMajor(rows, cols, beta, c, cOffset, cRowStride)
    gemmPanel(0, rows, 0, cols, 0, shared, alpha, a, aOffset, aRowStride, b, bOffset, bRowStride, c, cOffset, cRowStride)

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
      cRowStride: Int
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
          c00 += a0 * b0; c01 += a0 * b1; c02 += a0 * b2; c03 += a0 * b3
          c10 += a1 * b0; c11 += a1 * b1; c12 += a1 * b2; c13 += a1 * b3
          c20 += a2 * b0; c21 += a2 * b1; c22 += a2 * b2; c23 += a2 * b3
          c30 += a3 * b0; c31 += a3 * b1; c32 += a3 * b2; c33 += a3 * b3
          bRow += bRowStride
          k += 1
        storeTile4(c, cRow0, j, alpha, c00, c01, c02, c03)
        storeTile4(c, cRow1, j, alpha, c10, c11, c12, c13)
        storeTile4(c, cRow2, j, alpha, c20, c21, c22, c23)
        storeTile4(c, cRow3, j, alpha, c30, c31, c32, c33)
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
          s0 += a(aRow0 + k) * bv
          s1 += a(aRow1 + k) * bv
          s2 += a(aRow2 + k) * bv
          s3 += a(aRow3 + k) * bv
          bRow += bRowStride
          k += 1
        val x0 = cRow0 + j; c(x0) = c(x0) + alpha * s0
        val x1 = cRow1 + j; c(x1) = c(x1) + alpha * s1
        val x2 = cRow2 + j; c(x2) = c(x2) + alpha * s2
        val x3 = cRow3 + j; c(x3) = c(x3) + alpha * s3
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
          s += a(aRow + k) * b(bRow + j)
          bRow += bRowStride
          k += 1
        val idx = cRow + j
        c(idx) = c(idx) + alpha * s
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
      t3: Double
  ): Unit =
    val i0 = cRow + j
    c(i0) = c(i0) + alpha * t0
    val i1 = cRow + j + 1
    c(i1) = c(i1) + alpha * t1
    val i2 = cRow + j + 2
    c(i2) = c(i2) + alpha * t2
    val i3 = cRow + j + 3
    c(i3) = c(i3) + alpha * t3

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
          gemmPanel(ii, iMax, jj, jMax, kk, kMax, alpha, a, aOffset, aRowStride, b, bOffset, bRowStride, c, cOffset, cRowStride)
          jj += GemmBlock
        kk += GemmBlock
      ii += GemmBlock

  /** Symmetric rank-k product `C := AᵀA` for a '''row-major''' `A` (`m × k`, unit
    * column stride), writing the full symmetric `k × k` result into `c` (row-major,
    * unit column stride, assumed zero on entry — `DMat.zeros` guarantees it).
    * '''Assign-only:''' there is no `alpha`/`beta` and `C` must be zero on entry —
    * do NOT wire this into an accumulating (`C += …`) or scaled path without
    * adding those semantics first.
    *
    * A general gemm computing `Aᵀ·A` sees `Aᵀ` column-strided (stride `k`), so both
    * operands stream cache-hostilely; this kernel instead accumulates each row of
    * `A` as a '''triangular outer product''' into the upper triangle, then mirrors.
    * That halves the flops and keeps every inner access unit-stride: for row `l`,
    * `C[i, i..k-1] += A[l,i] · A[l, i..k-1]` streams the contiguous tail of `A`'s
    * row against the contiguous tail of `C`'s row. The inner loop is unrolled 4x.
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
    var l = 0
    var aRow = aOffset
    while l < m do
      var i = 0
      while i < k do
        val ali = a(aRow + i)
        val cRow = cOffset + i * cRowStride
        val count = k - i
        val jLimit = i + (count - (count & 3))
        var j = i
        while j < jLimit do
          val i0 = cRow + j
          c(i0) = c(i0) + ali * a(aRow + j)
          val i1 = cRow + j + 1
          c(i1) = c(i1) + ali * a(aRow + j + 1)
          val i2 = cRow + j + 2
          c(i2) = c(i2) + ali * a(aRow + j + 2)
          val i3 = cRow + j + 3
          c(i3) = c(i3) + ali * a(aRow + j + 3)
          j += 4
        while j < k do
          val idx = cRow + j
          c(idx) = c(idx) + ali * a(aRow + j)
          j += 1
        i += 1
      aRow += aRowStride
      l += 1
    // Mirror the computed upper triangle into the lower.
    var i = 1
    while i < k do
      val cRowI = cOffset + i * cRowStride
      var j = 0
      while j < i do
        c(cRowI + j) = c(cOffset + j * cRowStride + i)
        j += 1
      i += 1
