package gale.backend.jvm.`native`

import gale.linalg.{DMat, Matrix}

import java.lang.foreign.{Arena, MemorySegment, ValueLayout}

enum Layout:
  case RowMajor, ColMajor

/** Dense off-heap matrix with an explicit FFM lifetime and leading dimension.
  *
  * Instances created by [[NativeDMat.allocate]] own their arena and release it
  * from `close()`. Instances created by [[NativeDMat.allocateIn]] or `toNative`
  * borrow the caller's arena; their `close()` is intentionally a no-op and the
  * caller closes the arena. In either case FFM checks use-after-close.
  */
final class NativeDMat private (
    private[gale] val memory: MemorySegment,
    val rows: Int,
    val cols: Int,
    val leadingDimension: Int,
    val layout: Layout,
    private val ownedArena: Option[Arena]
) extends AutoCloseable:
  private val elementCount =
    if rows == 0 || cols == 0 then 0L
    else layout match
      case Layout.RowMajor => (rows - 1L) * leadingDimension + cols
      case Layout.ColMajor => (cols - 1L) * leadingDimension + rows

  require(rows >= 0 && cols >= 0, s"negative native matrix shape: $rows x $cols")
  require(
    leadingDimension >= (layout match
      case Layout.RowMajor => math.max(1, cols)
      case Layout.ColMajor => math.max(1, rows)),
    s"leading dimension $leadingDimension is too small for $rows x $cols $layout"
  )
  require(memory.byteSize() >= elementCount * java.lang.Double.BYTES,
    s"memory segment is too small for $rows x $cols $layout")

  def isAlive: Boolean = memory.scope().isAlive()

  def apply(row: Int, col: Int): Double =
    memory.get(ValueLayout.JAVA_DOUBLE, byteOffset(row, col))

  def update(row: Int, col: Int, value: Double): Unit =
    memory.set(ValueLayout.JAVA_DOUBLE, byteOffset(row, col), value)

  def toHeap: DMat =
    val total = Math.multiplyExact(rows, cols)
    val values = new Array[Double](total)
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        values(i * cols + j) = apply(i, j)
        j += 1
      i += 1
    DMat.fromArrayRowMajor(rows, cols, values)

  override def close(): Unit =
    ownedArena.foreach { arena =>
      if arena.scope().isAlive() then arena.close()
    }

  private def byteOffset(row: Int, col: Int): Long =
    if row < 0 || row >= rows then throw IndexOutOfBoundsException(s"row $row outside [0, $rows)")
    if col < 0 || col >= cols then throw IndexOutOfBoundsException(s"column $col outside [0, $cols)")
    val element = layout match
      case Layout.RowMajor => row.toLong * leadingDimension + col
      case Layout.ColMajor => col.toLong * leadingDimension + row
    element * java.lang.Double.BYTES

object NativeDMat:
  def allocate(rows: Int, cols: Int, layout: Layout = Layout.ColMajor): NativeDMat =
    val arena = Arena.ofConfined()
    try allocateImpl(rows, cols, layout, arena, Some(arena))
    catch
      case error: Throwable =>
        arena.close()
        throw error

  def allocateIn(rows: Int, cols: Int, layout: Layout = Layout.ColMajor)(using arena: Arena): NativeDMat =
    allocateImpl(rows, cols, layout, arena, None)

  extension (matrix: DMat)
    def toNative(layout: Layout = Layout.ColMajor)(using arena: Arena): NativeDMat =
      val native = allocateIn(matrix.rows, matrix.cols, layout)
      var i = 0
      while i < matrix.rows do
        var j = 0
        while j < matrix.cols do
          native(i, j) = matrix(i, j)
          j += 1
        i += 1
      native

  private def allocateImpl(
      rows: Int,
      cols: Int,
      layout: Layout,
      arena: Arena,
      owner: Option[Arena]
  ): NativeDMat =
    require(rows >= 0 && cols >= 0, s"negative native matrix shape: $rows x $cols")
    val leadingDimension = layout match
      case Layout.RowMajor => math.max(1, cols)
      case Layout.ColMajor => math.max(1, rows)
    val elements = Math.multiplyExact(rows.toLong, cols.toLong)
    val bytes = Math.multiplyExact(elements, java.lang.Double.BYTES.toLong)
    // Arena.allocate(0, alignment) is not portable across all supported JDKs.
    val segment = arena.allocate(math.max(1L, bytes), java.lang.Double.BYTES)
    new NativeDMat(segment, rows, cols, leadingDimension, layout, owner)
