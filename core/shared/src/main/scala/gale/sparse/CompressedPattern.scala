package gale.sparse

import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*

/** Allocation-free callback for one coordinate in a compressed pattern. */
trait SparsePositionConsumer:
  def apply(row: Int, col: Int): Unit

/** Allocation-free callback for an entry within one compressed lane.
  * `storedIndex` addresses the corresponding numeric value in storage order.
  */
trait SparseLaneEntryConsumer:
  def apply(minorIndex: Int, storedIndex: Int): Unit

/** Immutable, validated CSR structure independent of numeric values.
  *
  * Public factories defensively copy offsets and indices. The backing storage is
  * never exposed; traversal and lane queries are the interop surface. A pattern
  * may be shared by any number of immutable matrices.
  */
final class CSRPattern private[sparse] (
    private[sparse] val storage: CompressedPatternStorage,
    val rows: Int,
    val cols: Int
):
  def nnz: Int = storage.nnz

  /** True when columns are strictly increasing within every row. */
  def hasCanonicalFormat: Boolean = storage.hasCanonicalFormat

  def rowNnz(row: Int): Int =
    checkRow(row)
    storage.laneEnd(row) - storage.laneStart(row)

  def columnIndexAt(storedIndex: Int): Int =
    storage.minorAt(storedIndex)

  def foreachRow(row: Int)(consumer: SparseLaneEntryConsumer): Unit =
    checkRow(row)
    storage.foreachLane(row, consumer)

  def foreachStoredPosition(consumer: SparsePositionConsumer): Unit =
    var row = 0
    while row < rows do
      var p = storage.laneStart(row)
      val end = storage.laneEnd(row)
      while p < end do
        consumer(row, storage.minorAtUnchecked(p))
        p += 1
      row += 1

  /** Zero-copy structural transpose. Immutable offset/index storage is shared. */
  def t: CSCPattern =
    new CSCPattern(storage, cols, rows)

  /** Bind values through exactly one owned logical copy. */
  def bind(values: DVec): Either[LinAlgError, CSR] =
    if values.length != nnz then Left(LinAlgError.VectorLengthMismatch(nnz, values.length))
    else Right(bindOwned(values.toDoubleArrayOwnedCopy))

  /** Bind a caller array through exactly one defensive copy. */
  def bind(values: Array[Double]): Either[LinAlgError, CSR] =
    if values.length != nnz then Left(LinAlgError.VectorLengthMismatch(nnz, values.length))
    else Right(bindOwned(DoubleArray.fromArray(values)))

  /** Allocate a single-owner numeric builder for this pattern. */
  def valuesBuilder(initialValue: Double = 0.0): CSRValuesBuilder =
    val builder = new CSRValuesBuilder(this, DoubleArray.alloc(nnz))
    if initialValue != 0.0 then builder.fill(initialValue)
    builder

  /** Allocate reusable numeric destination storage for symbolic sparse plans.
    * Unlike [[valuesBuilder]], this destination remains open across evaluations;
    * [[CSRValuesDestination.snapshot]] makes the explicit owned immutable copy.
    */
  def valuesDestination(initialValue: Double = 0.0): CSRValuesDestination =
    val destination = new CSRValuesDestination(this, DoubleArray.alloc(nnz))
    if initialValue != 0.0 then destination.fill(initialValue)
    destination

  private[sparse] def bindOwned(values: DoubleArray): CSR =
    new CSR(rows, cols, storage.majorOffsets, storage.minorIndices, values)

  private[sparse] def sharesStorageWith(that: CSRPattern): Boolean =
    rows == that.rows && cols == that.cols && storage.sharesStorageWith(that.storage)

  private def checkRow(row: Int): Unit =
    if row < 0 || row >= rows then throw LinAlgError.IndexOutOfBounds(row, rows)

  override def equals(other: Any): Boolean =
    other match
      case that: CSRPattern =>
        rows == that.rows && cols == that.cols && storage.sameStructure(that.storage)
      case _ => false

  override def hashCode(): Int =
    var hash = 31 * rows + cols
    hash = 31 * hash + storage.structuralHash
    hash

  override def toString: String =
    s"CSRPattern($rows,$cols,nnz=$nnz)"

object CSRPattern:
  /** Validate and defensively copy CSR row offsets and column indices. */
  def checked(
      rows: Int,
      cols: Int,
      rowOffsets: Array[Int],
      columnIndices: Array[Int]
  ): Either[LinAlgError, CSRPattern] =
    CompressedPatternStorage
      .checked(rows, cols, rowOffsets, columnIndices, "CSR row offsets", "CSR column")
      .map(new CSRPattern(_, rows, cols))

  /** Collection-oriented checked factory; inputs are copied into owned storage. */
  def checked(
      rows: Int,
      cols: Int,
      rowOffsets: IndexedSeq[Int],
      columnIndices: IndexedSeq[Int]
  ): Either[LinAlgError, CSRPattern] =
    checked(rows, cols, rowOffsets.toArray, columnIndices.toArray)

/** Immutable, validated CSC structure independent of numeric values. */
final class CSCPattern private[sparse] (
    private[sparse] val storage: CompressedPatternStorage,
    val rows: Int,
    val cols: Int
):
  def nnz: Int = storage.nnz

  /** True when rows are strictly increasing within every column. */
  def hasCanonicalFormat: Boolean = storage.hasCanonicalFormat

  def columnNnz(col: Int): Int =
    checkCol(col)
    storage.laneEnd(col) - storage.laneStart(col)

  def rowIndexAt(storedIndex: Int): Int =
    storage.minorAt(storedIndex)

  def foreachColumn(col: Int)(consumer: SparseLaneEntryConsumer): Unit =
    checkCol(col)
    storage.foreachLane(col, consumer)

  def foreachStoredPosition(consumer: SparsePositionConsumer): Unit =
    var col = 0
    while col < cols do
      var p = storage.laneStart(col)
      val end = storage.laneEnd(col)
      while p < end do
        consumer(storage.minorAtUnchecked(p), col)
        p += 1
      col += 1

  /** Zero-copy structural transpose. Immutable offset/index storage is shared. */
  def t: CSRPattern =
    new CSRPattern(storage, cols, rows)

  /** Bind values through exactly one owned logical copy. */
  def bind(values: DVec): Either[LinAlgError, CSC] =
    if values.length != nnz then Left(LinAlgError.VectorLengthMismatch(nnz, values.length))
    else Right(bindOwned(values.toDoubleArrayOwnedCopy))

  /** Bind a caller array through exactly one defensive copy. */
  def bind(values: Array[Double]): Either[LinAlgError, CSC] =
    if values.length != nnz then Left(LinAlgError.VectorLengthMismatch(nnz, values.length))
    else Right(bindOwned(DoubleArray.fromArray(values)))

  /** Allocate a single-owner numeric builder for this pattern. */
  def valuesBuilder(initialValue: Double = 0.0): CSCValuesBuilder =
    val builder = new CSCValuesBuilder(this, DoubleArray.alloc(nnz))
    if initialValue != 0.0 then builder.fill(initialValue)
    builder

  private[sparse] def bindOwned(values: DoubleArray): CSC =
    new CSC(rows, cols, storage.majorOffsets, storage.minorIndices, values)

  private def checkCol(col: Int): Unit =
    if col < 0 || col >= cols then throw LinAlgError.IndexOutOfBounds(col, cols)

  override def equals(other: Any): Boolean =
    other match
      case that: CSCPattern =>
        rows == that.rows && cols == that.cols && storage.sameStructure(that.storage)
      case _ => false

  override def hashCode(): Int =
    var hash = 31 * rows + cols
    hash = 31 * hash + storage.structuralHash
    hash

  override def toString: String =
    s"CSCPattern($rows,$cols,nnz=$nnz)"

object CSCPattern:
  /** Validate and defensively copy CSC column offsets and row indices. */
  def checked(
      rows: Int,
      cols: Int,
      columnOffsets: Array[Int],
      rowIndices: Array[Int]
  ): Either[LinAlgError, CSCPattern] =
    CompressedPatternStorage
      .checked(cols, rows, columnOffsets, rowIndices, "CSC column offsets", "CSC row")
      .map(new CSCPattern(_, rows, cols))

  /** Collection-oriented checked factory; inputs are copied into owned storage. */
  def checked(
      rows: Int,
      cols: Int,
      columnOffsets: IndexedSeq[Int],
      rowIndices: IndexedSeq[Int]
  ): Either[LinAlgError, CSCPattern] =
    checked(rows, cols, columnOffsets.toArray, rowIndices.toArray)

/** Single-owner numeric construction buffer for a [[CSRPattern]]. */
final class CSRValuesBuilder private[sparse] (
    val pattern: CSRPattern,
    private val values: DoubleArray
):
  private var open = true

  def length: Int = pattern.nnz

  def apply(storedIndex: Int): Double =
    requireOpen()
    checkIndex(storedIndex)
    values(storedIndex)

  def update(storedIndex: Int, value: Double): Unit =
    requireOpen()
    checkIndex(storedIndex)
    values(storedIndex) = value

  def fill(value: Double): Unit =
    requireOpen()
    var i = 0
    while i < length do
      values(i) = value
      i += 1

  /** Transfer numeric storage without copying and permanently close the builder. */
  def result(): CSR =
    requireOpen()
    open = false
    pattern.bindOwned(values)

  private def requireOpen(): Unit =
    if !open then
      throw LinAlgError.UnsupportedOperation("CSRValuesBuilder is closed after result()")

  private def checkIndex(index: Int): Unit =
    if index < 0 || index >= length then throw LinAlgError.IndexOutOfBounds(index, length)

/** Single-owner numeric construction buffer for a [[CSCPattern]]. */
final class CSCValuesBuilder private[sparse] (
    val pattern: CSCPattern,
    private val values: DoubleArray
):
  private var open = true

  def length: Int = pattern.nnz

  def apply(storedIndex: Int): Double =
    requireOpen()
    checkIndex(storedIndex)
    values(storedIndex)

  def update(storedIndex: Int, value: Double): Unit =
    requireOpen()
    checkIndex(storedIndex)
    values(storedIndex) = value

  def fill(value: Double): Unit =
    requireOpen()
    var i = 0
    while i < length do
      values(i) = value
      i += 1

  /** Transfer numeric storage without copying and permanently close the builder. */
  def result(): CSC =
    requireOpen()
    open = false
    pattern.bindOwned(values)

  private def requireOpen(): Unit =
    if !open then
      throw LinAlgError.UnsupportedOperation("CSCValuesBuilder is closed after result()")

  private def checkIndex(index: Int): Unit =
    if index < 0 || index >= length then throw LinAlgError.IndexOutOfBounds(index, length)

/** Shared major-offset/minor-index invariant model for CSR rows and CSC columns. */
private[sparse] final class CompressedPatternStorage private[sparse] (
    val majorOffsets: IndexArray,
    val minorIndices: IndexArray
):
  def nnz: Int = minorIndices.length

  def laneStart(lane: Int): Int = majorOffsets(lane)

  def laneEnd(lane: Int): Int = majorOffsets(lane + 1)

  def minorAt(storedIndex: Int): Int =
    if storedIndex < 0 || storedIndex >= nnz then
      throw LinAlgError.IndexOutOfBounds(storedIndex, nnz)
    minorIndices(storedIndex)

  def minorAtUnchecked(storedIndex: Int): Int = minorIndices(storedIndex)

  def foreachLane(lane: Int, consumer: SparseLaneEntryConsumer): Unit =
    var p = laneStart(lane)
    val end = laneEnd(lane)
    while p < end do
      consumer(minorIndices(p), p)
      p += 1

  lazy val hasCanonicalFormat: Boolean =
    var lane = 0
    var canonical = true
    while lane + 1 < majorOffsets.length && canonical do
      var p = majorOffsets(lane)
      val end = majorOffsets(lane + 1)
      var previous = -1
      while p < end && canonical do
        val current = minorIndices(p)
        if current <= previous then canonical = false
        previous = current
        p += 1
      lane += 1
    canonical

  def sameStructure(that: CompressedPatternStorage): Boolean =
    sameArray(majorOffsets, that.majorOffsets) && sameArray(minorIndices, that.minorIndices)

  def sharesStorageWith(that: CompressedPatternStorage): Boolean =
    sharesArrays(that.majorOffsets, that.minorIndices)

  def sharesArrays(offsets: IndexArray, indices: IndexArray): Boolean =
    (majorOffsets.asInstanceOf[AnyRef] eq offsets.asInstanceOf[AnyRef]) &&
      (minorIndices.asInstanceOf[AnyRef] eq indices.asInstanceOf[AnyRef])

  lazy val structuralHash: Int =
    var hash = 1
    var i = 0
    while i < majorOffsets.length do
      hash = 31 * hash + majorOffsets(i)
      i += 1
    i = 0
    while i < minorIndices.length do
      hash = 31 * hash + minorIndices(i)
      i += 1
    hash

  private def sameArray(left: IndexArray, right: IndexArray): Boolean =
    if left.length != right.length then false
    else
      var i = 0
      while i < left.length do
        if left(i) != right(i) then return false
        i += 1
      true

private[sparse] object CompressedPatternStorage:
  def checked(
      majorSize: Int,
      minorSize: Int,
      offsets: Array[Int],
      indices: Array[Int],
      offsetsName: String,
      minorName: String
  ): Either[LinAlgError, CompressedPatternStorage] =
    validate(majorSize, minorSize, offsets, indices, offsetsName, minorName).map: _ =>
      new CompressedPatternStorage(IndexArray.fromArray(offsets), IndexArray.fromArray(indices))

  private def validate(
      majorSize: Int,
      minorSize: Int,
      offsets: Array[Int],
      indices: Array[Int],
      offsetsName: String,
      minorName: String
  ): Either[LinAlgError, Unit] =
    if majorSize < 0 || minorSize < 0 then
      Left(LinAlgError.InvalidArgument(s"compressed pattern shape must be non-negative, got ${majorSize}x${minorSize}"))
    else if majorSize == Int.MaxValue || offsets.length != majorSize + 1 then
      Left(
        LinAlgError.InvalidArgument(
          s"$offsetsName length must be major dimension + 1 (${majorSize.toLong + 1L}), got ${offsets.length}"
        )
      )
    else if offsets(0) != 0 then
      Left(LinAlgError.InvalidArgument(s"$offsetsName must start at 0, got ${offsets(0)}"))
    else
      var lane = 0
      while lane < majorSize do
        val start = offsets(lane)
        val end = offsets(lane + 1)
        if start < 0 || end < start || end > indices.length then
          return Left(
            LinAlgError.InvalidArgument(
              s"invalid $offsetsName lane $lane: [$start, $end) for nnz ${indices.length}"
            )
          )
        lane += 1
      if offsets(majorSize) != indices.length then
        Left(
          LinAlgError.InvalidArgument(
            s"final $offsetsName value ${offsets(majorSize)} must equal nnz ${indices.length}"
          )
        )
      else
        var p = 0
        while p < indices.length do
          val index = indices(p)
          if index < 0 || index >= minorSize then
            return Left(
              LinAlgError.InvalidArgument(
                s"$minorName index $index at stored position $p is outside [0, $minorSize)"
              )
            )
          p += 1
        Right(())
