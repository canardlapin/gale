package gale.linalg

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Single-owner, row-major construction buffer for an immutable [[DMat]].
  *
  * `DMatBuilder` is the allocation-conscious shared construction seam for code
  * that fills matrices in primitive loops. It exposes logical element writes,
  * never the platform backing array. [[result]] transfers ownership of that
  * backing storage to the returned matrix without copying and permanently
  * closes the builder, preventing mutable aliasing through the public API.
  */
final class DMatBuilder private (val rows: Int, val cols: Int, private[gale] val data: DoubleArray):
  private var open = true

  def size: Int =
    rows * cols

  def apply(row: Int, col: Int): Double =
    requireOpen()
    checkRow(row)
    checkCol(col)
    data(row * cols + col)

  def update(row: Int, col: Int, value: Double): Unit =
    requireOpen()
    checkRow(row)
    checkCol(col)
    data(row * cols + col) = value

  /** Write by contiguous row-major logical index. This is useful for hot
    * sequential fill loops while keeping platform storage encapsulated.
    */
  def updateRowMajor(index: Int, value: Double): Unit =
    requireOpen()
    if index < 0 || index >= size then
      throw LinAlgError.IndexOutOfBounds(index, size)
    data(index) = value

  def fill(value: Double): Unit =
    requireOpen()
    var i = 0
    while i < size do
      data(i) = value
      i += 1

  /** Transfer this builders storage to an immutable matrix without copying.
    * Every subsequent read, write, fill, or second `result()` call fails.
    */
  def result(): DMat =
    requireOpen()
    open = false
    DMat.fromDoubleArrayOwned(rows, cols, data)

  private[gale] def mutableColumn(index: Int): MutableDVec =
    requireOpen()
    checkCol(index)
    new MutableDVec(
      data,
      Offset.unsafe(index),
      Length.unsafe(rows),
      Stride.unsafe(if cols == 0 then 1 else cols)
    )

  private[gale] def mutableRow(index: Int): MutableDVec =
    requireOpen()
    checkRow(index)
    new MutableDVec(
      data,
      Offset.unsafe(index * cols),
      Length.unsafe(cols),
      Stride.unsafe(1)
    )

  private def requireOpen(): Unit =
    if !open then
      throw LinAlgError.UnsupportedOperation("DMatBuilder is closed after result()")

  private def checkRow(row: Int): Unit =
    if row < 0 || row >= rows then
      throw LinAlgError.IndexOutOfBounds(row, rows)

  private def checkCol(col: Int): Unit =
    if col < 0 || col >= cols then
      throw LinAlgError.IndexOutOfBounds(col, cols)

object DMatBuilder:
  def zeros(rows: Int, cols: Int): DMatBuilder =
    DMat.requireStorable(rows, cols)
    new DMatBuilder(rows, cols, DoubleArray.alloc(rows * cols))

  /** Initialize a builder from exactly one owned logical row-major copy.
    * Strided and transposed inputs are normalized without exposing their backing
    * storage; [[DMatBuilder.result]] then transfers this copy without recopying.
    */
  def from(matrix: DMat): DMatBuilder =
    new DMatBuilder(matrix.rows, matrix.cols, matrix.toDoubleArrayCopyRowMajor)
