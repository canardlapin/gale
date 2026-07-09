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
    var row = 0
    var aRow = aOffset
    var yi = yOffset
    while row < rows do
      var col = 0
      var ai = aRow
      var xi = xOffset
      var acc = 0.0
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
  private inline val GemmBlock = 64

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
    var row = 0
    var cRow = cOffset
    var aRow = aOffset
    while row < rows do
      var k = 0
      var aik = aRow
      var bRow = bOffset
      while k < shared do
        val scaled = alpha * a(aik)
        var col = 0
        var cij = cRow
        var bkj = bRow
        while col < cols do
          c(cij) = c(cij) + scaled * b(bkj)
          cij += 1
          bkj += 1
          col += 1
        aik += 1
        bRow += bRowStride
        k += 1
      cRow += cRowStride
      aRow += aRowStride
      row += 1

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
          var i = ii
          while i < iMax do
            val cRow = cOffset + i * cRowStride
            val aRow = aOffset + i * aRowStride
            var k = kk
            while k < kMax do
              val scaled = alpha * a(aRow + k)
              val bRow = bOffset + k * bRowStride
              var col = jj
              var cij = cRow + jj
              var bkj = bRow + jj
              while col < jMax do
                c(cij) = c(cij) + scaled * b(bkj)
                cij += 1
                bkj += 1
                col += 1
              k += 1
            i += 1
          jj += GemmBlock
        kk += GemmBlock
      ii += GemmBlock
