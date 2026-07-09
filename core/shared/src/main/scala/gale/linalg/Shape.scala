package gale.linalg

import scala.annotation.targetName

opaque type Rows = Int
object Rows:
  def apply(value: Int): Rows =
    require(value >= 0, "rows must be non-negative")
    value

  private[gale] def unsafe(value: Int): Rows =
    value

extension (rows: Rows)
  @targetName("rowsValue")
  inline def value: Int = rows

opaque type Cols = Int
object Cols:
  def apply(value: Int): Cols =
    require(value >= 0, "cols must be non-negative")
    value

  private[gale] def unsafe(value: Int): Cols =
    value

extension (cols: Cols)
  @targetName("colsValue")
  inline def value: Int = cols

opaque type Length = Int
object Length:
  def apply(value: Int): Length =
    require(value >= 0, "length must be non-negative")
    value

  private[gale] def unsafe(value: Int): Length =
    value

extension (length: Length)
  @targetName("lengthValue")
  inline def value: Int = length

opaque type Offset = Int
object Offset:
  def apply(value: Int): Offset =
    require(value >= 0, "offset must be non-negative")
    value

  private[gale] def unsafe(value: Int): Offset =
    value

extension (offset: Offset)
  @targetName("offsetValue")
  inline def value: Int = offset

opaque type Stride = Int
object Stride:
  def apply(value: Int): Stride =
    require(value > 0, "stride must be positive")
    value

  private[gale] def unsafe(value: Int): Stride =
    value

extension (stride: Stride)
  @targetName("strideValue")
  inline def value: Int = stride

opaque type NonZeroCount = Int
object NonZeroCount:
  def apply(value: Int): NonZeroCount =
    require(value >= 0, "non-zero count must be non-negative")
    value

extension (count: NonZeroCount)
  @targetName("nonZeroCountValue")
  inline def value: Int = count

final case class Shape(rows: Rows, cols: Cols):
  def size: Long = rows.value.toLong * cols.value.toLong
