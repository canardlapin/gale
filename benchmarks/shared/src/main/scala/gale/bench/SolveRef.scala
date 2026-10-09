package gale.bench

import gale.kernel.DoubleKernels
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** The multi-right-hand-side LU and Cholesky solves as they stood before the
  * blocked `dtrsmLeft` kernel (one `dtrsv` per column for LU; an element-accessor
  * triple loop for Cholesky), kept verbatim so JVM and Scala.js harnesses can
  * pair them against the live `solve` in one process.
  */
object SolveRef:
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
