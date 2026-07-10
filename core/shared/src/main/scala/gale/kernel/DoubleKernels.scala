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
