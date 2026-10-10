package gale.bench

import gale.kernel.DoubleKernels
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** The dense LU factorization and multi-right-hand-side LU and Cholesky solves
  * as they stood before the paired-column LU and the blocked `dtrsmLeft` kernel
  * (one-column elimination; one `dtrsv` per column for LU; an element-accessor
  * triple loop for Cholesky), kept verbatim so JVM and Scala.js harnesses can
  * pair them against the live code in one process.
  */
object SolveRef:
  /** One-column-at-a-time right-looking LU, wrapped exactly as the live `lu`. */
  def lu(a: DMat): LU =
    val n = a.rows
    val packed = a.toDoubleArrayCopyRowMajor
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
      if maxAbs == 0.0 || maxAbs.isNaN then throw LinAlgError.SingularMatrix(k)
      if pivot != k then
        var col = 0
        while col < n do
          val t = packed(k * n + col)
          packed(k * n + col) = packed(pivot * n + col)
          packed(pivot * n + col) = t
          col += 1
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
    LU(
      packed = DMat.fromDoubleArrayOwned(n, n, packed),
      pivots = PivotVector.fromArray(pivots),
      parity = parity,
      diagnostics = FactorizationDiagnostics(info = 0)
    )

  def luSolve(lu: LU, b: DMat): DMat =
    val n = lu.packed.rows
    val rhsCols = b.cols
    val packed = lu.packed
    val packedData = packed.data
    val packedOffset = packed.offset.value
    val packedRowStride = packed.rowStride.value
    val packedColStride = packed.colStride.value
    val values = DoubleArray.alloc(n * rhsCols)
    var row = 0
    while row < n do
      val sourceRow = lu.pivots(row)
      var rhs = 0
      while rhs < rhsCols do
        values(row * rhsCols + rhs) = b(sourceRow, rhs)
        rhs += 1
      row += 1
    var rhs = 0
    while rhs < rhsCols do
      DoubleKernels.dtrsv(
        n, lower = true, unit = true, 0.0, packedData, packedOffset, packedRowStride, packedColStride, values, rhs,
        rhsCols
      )
      val info = DoubleKernels.dtrsv(
        n, lower = false, unit = false, 0.0, packedData, packedOffset, packedRowStride, packedColStride, values, rhs,
        rhsCols
      )
      if info >= 0 then throw LinAlgError.SingularMatrix(info)
      rhs += 1
    DMat.fromDoubleArrayOwned(n, rhsCols, values)

  def choleskySolve(cholesky: Cholesky, b: DMat): DMat =
    val n = cholesky.lower.rows
    val rhsCols = b.cols
    val x = b.toDoubleArrayCopyRowMajor
    if !finite(x) then throw LinAlgError.InvalidArgument("non-finite Cholesky right-hand side")
    val lower = cholesky.lower
    var row = 0
    while row < n do
      val diagonal = lower(row, row)
      var rhs = 0
      while rhs < rhsCols do
        var value = x(row * rhsCols + rhs)
        var k = 0
        while k < row do
          value -= lower(row, k) * x(k * rhsCols + rhs)
          k += 1
        x(row * rhsCols + rhs) = value / diagonal
        rhs += 1
      row += 1
    row = n - 1
    while row >= 0 do
      val diagonal = lower(row, row)
      var rhs = 0
      while rhs < rhsCols do
        var value = x(row * rhsCols + rhs)
        var k = row + 1
        while k < n do
          value -= lower(k, row) * x(k * rhsCols + rhs)
          k += 1
        x(row * rhsCols + rhs) = value / diagonal
        rhs += 1
      row -= 1
    if !finite(x) then throw LinAlgError.InvalidArgument("non-finite Cholesky solution")
    DMat.fromDoubleArrayOwned(n, rhsCols, x)

  private def finite(values: DoubleArray): Boolean =
    var index = 0
    while index < values.length do
      if !values(index).isFinite then return false
      index += 1
    true
