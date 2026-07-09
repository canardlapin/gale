package gale.sparse

import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*
import scala.collection.mutable.ArrayBuffer

enum DuplicatePolicy:
  case Sum, Last, Error

trait SparseMatrix[A] extends Matrix[A]:
  def nnz: Int
  def density: Double =
    if rows == 0 || cols == 0 then 0.0 else nnz.toDouble / (rows.toDouble * cols.toDouble)

final case class COOEntry(row: Int, col: Int, value: Double)

final class COO private[gale] (
    val rows: Int,
    val cols: Int,
    private val rowIndices: IndexArray,
    private val colIndices: IndexArray,
    private val entryValues: DoubleArray,
    val hasCanonicalFormat: Boolean
) extends SparseMatrix[Double]
    with DoubleLinearOperator:
  def nnz: Int =
    entryValues.length

  def entries: Seq[COOEntry] =
    val builder = Vector.newBuilder[COOEntry]
    builder.sizeHint(nnz)
    var i = 0
    while i < nnz do
      builder += COOEntry(rowIndices(i), colIndices(i), entryValues(i))
      i += 1
    builder.result()

  /** Value at `(row, col)`, summing every stored entry at that coordinate.
    *
    * On a canonical COO (no duplicate coordinates) this is the single stored value
    * (or `0.0`), which agrees exactly with [[CSR.apply]] on the same matrix. On a
    * non-canonical COO it sums duplicates — the construction-oriented, sum-on-
    * assembly semantics COO carries until it is canonicalized or converted.
    */
  def apply(row: Int, col: Int): Double =
    checkRow(row)
    checkCol(col)
    var i = 0
    var out = 0.0
    while i < nnz do
      if rowIndices(i) == row && colIndices(i) == col then
        out += entryValues(i)
      i += 1
    out

  def row(index: Int): DVec =
    checkRow(index)
    val out = DVec.zeros(cols)
    var p = 0
    while p < nnz do
      if rowIndices(p) == index then
        out.data(colIndices(p)) = out.data(colIndices(p)) + entryValues(p)
      p += 1
    out

  def col(index: Int): DVec =
    checkCol(index)
    val out = DVec.zeros(rows)
    var p = 0
    while p < nnz do
      if colIndices(p) == index then
        out.data(rowIndices(p)) = out.data(rowIndices(p)) + entryValues(p)
      p += 1
    out

  def t: COO =
    // Swapping (row, col) turns row-major order into column-major order, which is
    // not canonical in general; do not carry the flag through the transpose.
    COO(cols, rows, IndexArray.copy(colIndices), IndexArray.copy(rowIndices), DoubleArray.copy(entryValues), canonical = false)

  // COO is construction-oriented: matvec scatters directly over the entry
  // triples (summing duplicates) without building a CSR first. For repeated
  // products, convert once via `toCSR`/`toCSC` and reuse that.
  override def applyTo(x: DVec, into: MutableDVec): Unit =
    if x.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if into.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(into.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, into.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val yData = into.data
    val yOff = into.offset.value
    val yStep = into.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val rIdx = rowIndices
    val cIdx = colIndices
    val vData = entryValues
    var row = 0
    while row < rows do
      yData(yOff + row * yStep) = 0.0
      row += 1
    var p = 0
    val count = vData.length
    while p < count do
      val r = rIdx(p)
      yData(yOff + r * yStep) = yData(yOff + r * yStep) + vData(p) * xData(xOff + cIdx(p) * xStep)
      p += 1

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    if x.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if into.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(into.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, into.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val yData = into.data
    val yOff = into.offset.value
    val yStep = into.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val rIdx = rowIndices
    val cIdx = colIndices
    val vData = entryValues
    var col = 0
    while col < cols do
      yData(yOff + col * yStep) = 0.0
      col += 1
    var p = 0
    val count = vData.length
    while p < count do
      val c = cIdx(p)
      yData(yOff + c * yStep) = yData(yOff + c * yStep) + vData(p) * xData(xOff + rIdx(p) * xStep)
      p += 1

  def toCSR: CSR =
    CSR.fromCOO(this)

  def toCSC: CSC =
    CSC.fromCOO(this)

  /** Canonical form of this COO: entries sorted row-major, duplicate coordinates
    * summed, and entries that are (or sum to) zero dropped. Reuses the builder's
    * sort/sum, then prunes zeros so cancellations disappear too.
    */
  def canonicalize: COO =
    val builder = new COOBuilder(rows, cols)
    var i = 0
    while i < nnz do
      builder.add(rowIndices(i), colIndices(i), entryValues(i))
      i += 1
    builder.canonicalize(DuplicatePolicy.Sum).pruneZeros.toCOO()

  private[gale] def rawRows: IndexArray =
    IndexArray.copy(rowIndices)

  private[gale] def rawCols: IndexArray =
    IndexArray.copy(colIndices)

  private[gale] def rawValues: DoubleArray =
    DoubleArray.copy(entryValues)

  private def checkRow(row: Int): Unit =
    if row < 0 || row >= rows then
      throw LinAlgError.IndexOutOfBounds(row, rows)

  private def checkCol(col: Int): Unit =
    if col < 0 || col >= cols then
      throw LinAlgError.IndexOutOfBounds(col, cols)

object COO:
  private[gale] def apply(
      rows: Int,
      cols: Int,
      rowIndices: IndexArray,
      colIndices: IndexArray,
      values: DoubleArray,
      canonical: Boolean
  ): COO =
    new COO(rows, cols, rowIndices, colIndices, values, canonical)

final class COOBuilder private[gale] (val rows: Int, val cols: Int):
  require(rows >= 0 && cols >= 0, "sparse matrix shape must be non-negative")

  private val rowIndices = ArrayBuffer.empty[Int]
  private val colIndices = ArrayBuffer.empty[Int]
  private val entryValues = ArrayBuffer.empty[Double]

  def nnz: Int =
    entryValues.length

  def add(row: Int, col: Int, value: Double): COOBuilder =
    if row < 0 || row >= rows then
      throw LinAlgError.IndexOutOfBounds(row, rows)
    if col < 0 || col >= cols then
      throw LinAlgError.IndexOutOfBounds(col, cols)
    rowIndices += row
    colIndices += col
    entryValues += value
    this

  def canonicalize(duplicates: DuplicatePolicy = DuplicatePolicy.Sum): COOBuilder =
    val coo = toCOO(duplicates)
    rowIndices.clear()
    colIndices.clear()
    entryValues.clear()
    coo.entries.foreach { entry =>
      rowIndices += entry.row
      colIndices += entry.col
      entryValues += entry.value
    }
    this

  def pruneZeros: COOBuilder =
    prune(absBelow = 0.0)

  def prune(absBelow: Double): COOBuilder =
    require(absBelow >= 0.0, "prune threshold must be non-negative")
    var write = 0
    var read = 0
    while read < entryValues.length do
      if math.abs(entryValues(read)) > absBelow then
        rowIndices(write) = rowIndices(read)
        colIndices(write) = colIndices(read)
        entryValues(write) = entryValues(read)
        write += 1
      read += 1
    rowIndices.dropRightInPlace(rowIndices.length - write)
    colIndices.dropRightInPlace(colIndices.length - write)
    entryValues.dropRightInPlace(entryValues.length - write)
    this

  def sortedIndices: Boolean =
    var i = 1
    while i < rowIndices.length do
      val prevRow = rowIndices(i - 1)
      val prevCol = colIndices(i - 1)
      val row = rowIndices(i)
      val col = colIndices(i)
      if row < prevRow || (row == prevRow && col < prevCol) then
        return false
      i += 1
    true

  def hasCanonicalFormat: Boolean =
    sortedIndices && !hasDuplicates

  def toCOO(duplicates: DuplicatePolicy = DuplicatePolicy.Sum): COO =
    val entries =
      rowIndices.indices
        .map(i => COOEntry(rowIndices(i), colIndices(i), entryValues(i)))
        .toArray
        .sortWith { (a, b) =>
          a.row < b.row || (a.row == b.row && a.col < b.col)
        }
    val rowsOut = ArrayBuffer.empty[Int]
    val colsOut = ArrayBuffer.empty[Int]
    val valuesOut = ArrayBuffer.empty[Double]
    var i = 0
    while i < entries.length do
      val row = entries(i).row
      val col = entries(i).col
      var value = entries(i).value
      var j = i + 1
      while j < entries.length && entries(j).row == row && entries(j).col == col do
        duplicates match
          case DuplicatePolicy.Sum =>
            value += entries(j).value
          case DuplicatePolicy.Last =>
            value = entries(j).value
          case DuplicatePolicy.Error =>
            throw LinAlgError.InvalidArgument(s"duplicate sparse entry at ($row, $col)")
        j += 1
      rowsOut += row
      colsOut += col
      valuesOut += value
      i = j
    COO(
      rows,
      cols,
      IndexArray.fromArray(rowsOut.toArray),
      IndexArray.fromArray(colsOut.toArray),
      DoubleArray.fromArray(valuesOut.toArray),
      canonical = true
    )

  def toCSR(duplicates: DuplicatePolicy = DuplicatePolicy.Sum): CSR =
    toCOO(duplicates).toCSR

  def toCSC(duplicates: DuplicatePolicy = DuplicatePolicy.Sum): CSC =
    toCOO(duplicates).toCSC

  private def hasDuplicates: Boolean =
    var i = 1
    while i < rowIndices.length do
      if rowIndices(i) == rowIndices(i - 1) && colIndices(i) == colIndices(i - 1) then
        return true
      i += 1
    false

final class CSR private[gale] (
    val rows: Int,
    val cols: Int,
    private[gale] val rowPtr: IndexArray,
    private[gale] val colIdx: IndexArray,
    private[gale] val values: DoubleArray
) extends SparseMatrix[Double]
    with DoubleLinearOperator:
  def nnz: Int =
    values.length

  /** Value at `(row, col)`: the stored value whose column matches, or `0.0`.
    *
    * On a canonical CSR (sorted, no duplicate columns) there is at most one match,
    * so this is an exact lookup agreeing with [[COO.apply]]. A non-canonical CSR
    * with duplicate columns is never produced by the public API; if one is built
    * directly, `apply` returns the first stored match while the matvec sums them,
    * so canonicalize first for a well-defined value.
    */
  def apply(row: Int, col: Int): Double =
    checkRow(row)
    checkCol(col)
    var p = rowPtr(row)
    val end = rowPtr(row + 1)
    while p < end do
      if colIdx(p) == col then return values(p)
      p += 1
    0.0

  def row(index: Int): DVec =
    checkRow(index)
    val out = DVec.zeros(cols)
    var p = rowPtr(index)
    val end = rowPtr(index + 1)
    while p < end do
      out.data(colIdx(p)) = values(p)
      p += 1
    out

  def col(index: Int): DVec =
    checkCol(index)
    // Single pass over the row ranges: for each row scan its (sorted) column
    // slice for `index`, rather than re-entering apply() row by row.
    val out = DVec.zeros(rows)
    val outData = out.data
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    var row = 0
    while row < rows do
      var p = rPtr(row)
      val end = rPtr(row + 1)
      var found = false
      while p < end && !found do
        if cIdx(p) == index then
          outData(row) = vData(p)
          found = true
        p += 1
      row += 1
    out

  def t: CSC =
    // Zero-copy transpose: an m x n CSR reinterpreted with rowPtr as column
    // pointers and colIdx as row indices is exactly the n x m matrix in CSC
    // form. Arrays are never mutated after construction, so they are shared.
    new CSC(cols, rows, rowPtr, colIdx, values)

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    mulInto(x, into)

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    tMulInto(x, into)

  override def *(x: DVec): DVec =
    val out = MutableDVec.zeros(rows)
    mulInto(x, out)
    out.asVec

  def mulInto(x: DVec, y: MutableDVec): Unit =
    if x.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if y.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(y.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, y.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val yData = y.data
    val yOff = y.offset.value
    val yStep = y.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    var row = 0
    while row < rows do
      var acc = 0.0
      var p = rPtr(row)
      val end = rPtr(row + 1)
      while p < end do
        acc += vData(p) * xData(xOff + cIdx(p) * xStep)
        p += 1
      yData(yOff + row * yStep) = acc
      row += 1

  def tMulInto(x: DVec, y: MutableDVec): Unit =
    if x.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if y.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(y.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, y.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val yData = y.data
    val yOff = y.offset.value
    val yStep = y.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    var col = 0
    while col < cols do
      yData(yOff + col * yStep) = 0.0
      col += 1
    var row = 0
    while row < rows do
      var p = rPtr(row)
      val end = rPtr(row + 1)
      val scale = xData(xOff + row * xStep)
      while p < end do
        val c = cIdx(p)
        yData(yOff + c * yStep) = yData(yOff + c * yStep) + vData(p) * scale
        p += 1
      row += 1

  def *(B: DMat): DMat =
    if cols != B.rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(B.cols)), Shape(Rows(B.rows), Cols(B.cols)))
    val out = DMat.zeros(rows, B.cols)
    // out is freshly allocated, hence contiguous row-major with offset 0; the
    // flat `row * outCols + col` indexing below relies on exactly that.
    assert(out.isContiguousRowMajor, "CSR*DMat output must be contiguous row-major")
    val outData = out.data
    val outCols = out.cols
    val bData = B.data
    val bOff = B.offset.value
    val bRowStep = B.rowStride.value
    val bColStep = B.colStride.value
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    val bCols = B.cols
    var row = 0
    while row < rows do
      var p = rPtr(row)
      val end = rPtr(row + 1)
      val outRow = row * outCols
      while p < end do
        val k = cIdx(p)
        val scale = vData(p)
        val bRow = bOff + k * bRowStep
        var col = 0
        var bIdx = bRow
        var outIdx = outRow
        while col < bCols do
          outData(outIdx) = outData(outIdx) + scale * bData(bIdx)
          bIdx += bColStep
          outIdx += 1
          col += 1
        p += 1
      row += 1
    out

  def +(that: CSR): CSR =
    zipValues(that)(_ + _)

  def -(that: CSR): CSR =
    zipValues(that)(_ - _)

  def *(alpha: Double): CSR =
    mapValues(_ * alpha)

  def mapValues(f: Double => Double): CSR =
    // Structure is preserved: map the stored values in one pass and compact out
    // any that map to zero, keeping the existing (sorted) column order. No COO
    // rebuild and no re-sort.
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    val n = vData.length
    val outColIdx = new Array[Int](n)
    val outValues = new Array[Double](n)
    val outRowPtr = new Array[Int](rows + 1)
    var write = 0
    var row = 0
    while row < rows do
      outRowPtr(row) = write
      var p = rPtr(row)
      val end = rPtr(row + 1)
      while p < end do
        val value = f(vData(p))
        if value != 0.0 then
          outColIdx(write) = cIdx(p)
          outValues(write) = value
          write += 1
        p += 1
      row += 1
    outRowPtr(rows) = write
    new CSR(
      rows,
      cols,
      toIndexArray(outRowPtr, rows + 1),
      toIndexArray(outColIdx, write),
      toDoubleArray(outValues, write)
    )

  def zipValues(that: CSR)(f: (Double, Double) => Double): CSR =
    requireSameShape(that)
    // Two-pointer merge of the two canonical (sorted, unique) column ranges per
    // row. Missing entries contribute 0.0 to `f`; results that vanish are
    // pruned. Output stays canonical without any HashMap or re-sort.
    val aPtr = rowPtr
    val aIdx = colIdx
    val aVal = values
    val bPtr = that.rowPtr
    val bIdx = that.colIdx
    val bVal = that.values
    val maxNnz = aVal.length + bVal.length
    val outColIdx = new Array[Int](maxNnz)
    val outValues = new Array[Double](maxNnz)
    val outRowPtr = new Array[Int](rows + 1)
    var write = 0
    var row = 0
    while row < rows do
      outRowPtr(row) = write
      var pa = aPtr(row)
      val aEnd = aPtr(row + 1)
      var pb = bPtr(row)
      val bEnd = bPtr(row + 1)
      while pa < aEnd && pb < bEnd do
        val ca = aIdx(pa)
        val cb = bIdx(pb)
        if ca < cb then
          val v = f(aVal(pa), 0.0)
          if v != 0.0 then
            outColIdx(write) = ca
            outValues(write) = v
            write += 1
          pa += 1
        else if cb < ca then
          val v = f(0.0, bVal(pb))
          if v != 0.0 then
            outColIdx(write) = cb
            outValues(write) = v
            write += 1
          pb += 1
        else
          val v = f(aVal(pa), bVal(pb))
          if v != 0.0 then
            outColIdx(write) = ca
            outValues(write) = v
            write += 1
          pa += 1
          pb += 1
      while pa < aEnd do
        val v = f(aVal(pa), 0.0)
        if v != 0.0 then
          outColIdx(write) = aIdx(pa)
          outValues(write) = v
          write += 1
        pa += 1
      while pb < bEnd do
        val v = f(0.0, bVal(pb))
        if v != 0.0 then
          outColIdx(write) = bIdx(pb)
          outValues(write) = v
          write += 1
        pb += 1
      row += 1
    outRowPtr(rows) = write
    new CSR(
      rows,
      cols,
      toIndexArray(outRowPtr, rows + 1),
      toIndexArray(outColIdx, write),
      toDoubleArray(outValues, write)
    )

  def toDense(maxEntries: Int = Int.MaxValue): DMat =
    val entries = rows.toLong * cols.toLong
    if entries > maxEntries.toLong then
      throw LinAlgError.UnsupportedOperation(s"dense conversion would allocate $entries entries")
    val out = DMat.zeros(rows, cols)
    var row = 0
    while row < rows do
      var p = rowPtr(row)
      val end = rowPtr(row + 1)
      while p < end do
        out.data(row * cols + colIdx(p)) = values(p)
        p += 1
      row += 1
    out

  def toCSR: CSR =
    this

  def toCSC: CSC =
    val counts = new Array[Int](cols)
    var p = 0
    while p < nnz do
      counts(colIdx(p)) += 1
      p += 1
    val colPtr = prefixCounts(counts)
    val next = colPtr.clone()
    // Fill the nnz-sized outputs directly in platform storage: the fill is the
    // only write, so the CSC adopts them without a second fromArray copy.
    val rowIdxOut = IndexArray.alloc(nnz)
    val valuesOut = DoubleArray.alloc(nnz)
    var row = 0
    while row < rows do
      p = rowPtr(row)
      val end = rowPtr(row + 1)
      while p < end do
        val col = colIdx(p)
        val dest = next(col)
        rowIdxOut(dest) = row
        valuesOut(dest) = values(p)
        next(col) += 1
        p += 1
      row += 1
    new CSC(rows, cols, toIndexArray(colPtr, cols + 1), rowIdxOut, valuesOut)

  def toCOO: COO =
    val builder = Sparse.coo(rows, cols)
    var row = 0
    while row < rows do
      var p = rowPtr(row)
      val end = rowPtr(row + 1)
      while p < end do
        builder.add(row, colIdx(p), values(p))
        p += 1
      row += 1
    builder.toCOO()

  def diagonal: DVec =
    val n = math.min(rows, cols)
    val out = DVec.zeros(n)
    var i = 0
    while i < n do
      out.data(i) = apply(i, i)
      i += 1
    out

  def trace: Double =
    val n = math.min(rows, cols)
    var i = 0
    var out = 0.0
    while i < n do
      out += apply(i, i)
      i += 1
    out

  /** True iff column indices are strictly increasing within every row (sorted,
    * no duplicates) and no stored value is an explicit zero — the canonical CSR
    * contract. Computed once on first access and cached.
    */
  lazy val hasCanonicalFormat: Boolean =
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    var ok = true
    var row = 0
    while row < rows && ok do
      var p = rPtr(row)
      val end = rPtr(row + 1)
      var prev = -1
      while p < end && ok do
        val c = cIdx(p)
        if c <= prev || vData(p) == 0.0 then ok = false
        prev = c
        p += 1
      row += 1
    ok

  /** Order the column indices ascending within each row, carrying values along.
    * Duplicate columns and explicit zeros are preserved (use [[canonicalize]] to
    * also sum duplicates and drop zeros). The row structure is unchanged, so the
    * existing row pointers are shared.
    */
  def sortedIndices: CSR =
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    val n = vData.length
    val outColIdx = new Array[Int](n)
    val outValues = new Array[Double](n)
    var row = 0
    while row < rows do
      val start = rPtr(row)
      val end = rPtr(row + 1)
      var k = start
      while k < end do
        outColIdx(k) = cIdx(k)
        outValues(k) = vData(k)
        k += 1
      insertionSortRange(outColIdx, outValues, start, end)
      row += 1
    new CSR(rows, cols, rowPtr, toIndexArray(outColIdx, n), toDoubleArray(outValues, n))

  /** Drop entries whose value is an exact zero, preserving order and structure. */
  def pruneZeros: CSR =
    prune(absBelow = 0.0)

  /** Drop entries whose magnitude is `<= absBelow`, preserving the order of the
    * survivors. Structure-only: no sort or duplicate summation is performed.
    */
  def prune(absBelow: Double): CSR =
    require(absBelow >= 0.0, "prune threshold must be non-negative")
    // Reuse mapValues' compaction: mapping a dropped entry to 0.0 removes it.
    mapValues(v => if math.abs(v) > absBelow then v else 0.0)

  /** Fully canonicalize: sort columns within each row, sum duplicate columns, and
    * drop entries that are (or sum to) zero. The result satisfies
    * [[hasCanonicalFormat]] and preserves the sum-semantics matrix.
    */
  def canonicalize: CSR =
    val rPtr = rowPtr
    val cIdx = colIdx
    val vData = values
    val n = vData.length
    val outColIdx = new Array[Int](n)
    val outValues = new Array[Double](n)
    val outRowPtr = new Array[Int](rows + 1)
    // Per-row scratch, sized to the widest possible row (the whole nnz).
    val colBuf = new Array[Int](n)
    val valBuf = new Array[Double](n)
    var write = 0
    var row = 0
    while row < rows do
      outRowPtr(row) = write
      val start = rPtr(row)
      val end = rPtr(row + 1)
      val width = end - start
      var k = 0
      while k < width do
        colBuf(k) = cIdx(start + k)
        valBuf(k) = vData(start + k)
        k += 1
      insertionSortRange(colBuf, valBuf, 0, width)
      var p = 0
      while p < width do
        val c = colBuf(p)
        var sum = valBuf(p)
        var q = p + 1
        while q < width && colBuf(q) == c do
          sum += valBuf(q)
          q += 1
        if sum != 0.0 then
          outColIdx(write) = c
          outValues(write) = sum
          write += 1
        p = q
      row += 1
    outRowPtr(rows) = write
    new CSR(rows, cols, toIndexArray(outRowPtr, rows + 1), toIndexArray(outColIdx, write), toDoubleArray(outValues, write))

  private def checkRow(row: Int): Unit =
    if row < 0 || row >= rows then throw LinAlgError.IndexOutOfBounds(row, rows)

  private def checkCol(col: Int): Unit =
    if col < 0 || col >= cols then throw LinAlgError.IndexOutOfBounds(col, cols)

  private def requireSameShape(that: CSR): Unit =
    if rows != that.rows || cols != that.cols then
      throw LinAlgError.DimensionMismatch(shape, that.shape)

final class CSC private[gale] (
    val rows: Int,
    val cols: Int,
    private[gale] val colPtr: IndexArray,
    private[gale] val rowIdx: IndexArray,
    private[gale] val values: DoubleArray
) extends SparseMatrix[Double]
    with DoubleLinearOperator:
  def nnz: Int =
    values.length

  def apply(row: Int, col: Int): Double =
    toCSR(row, col)

  def row(index: Int): DVec =
    toCSR.row(index)

  def col(index: Int): DVec =
    if index < 0 || index >= cols then throw LinAlgError.IndexOutOfBounds(index, cols)
    val out = DVec.zeros(rows)
    var p = colPtr(index)
    val end = colPtr(index + 1)
    while p < end do
      out.data(rowIdx(p)) = values(p)
      p += 1
    out

  def t: CSR =
    // Zero-copy transpose: an m x n CSC reinterpreted with colPtr as row
    // pointers and rowIdx as column indices is exactly the n x m matrix in CSR
    // form. Arrays are shared (never mutated post-construction).
    new CSR(cols, rows, colPtr, rowIdx, values)

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    mulInto(x, into)

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    tMulInto(x, into)

  /** `y := A x` by column scatter: each column contributes `x(col)` scaled into
    * the rows it touches. No CSR conversion.
    */
  def mulInto(x: DVec, y: MutableDVec): Unit =
    if x.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if y.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(y.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, y.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val yData = y.data
    val yOff = y.offset.value
    val yStep = y.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val cPtr = colPtr
    val rIdx = rowIdx
    val vData = values
    var row = 0
    while row < rows do
      yData(yOff + row * yStep) = 0.0
      row += 1
    var col = 0
    while col < cols do
      var p = cPtr(col)
      val end = cPtr(col + 1)
      val scale = xData(xOff + col * xStep)
      while p < end do
        val r = rIdx(p)
        yData(yOff + r * yStep) = yData(yOff + r * yStep) + vData(p) * scale
        p += 1
      col += 1

  /** `y := Aᵀ x` by per-column gather (each CSC column is a row of `Aᵀ`). */
  def tMulInto(x: DVec, y: MutableDVec): Unit =
    if x.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if y.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(y.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, y.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val yData = y.data
    val yOff = y.offset.value
    val yStep = y.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val cPtr = colPtr
    val rIdx = rowIdx
    val vData = values
    var col = 0
    while col < cols do
      var p = cPtr(col)
      val end = cPtr(col + 1)
      var acc = 0.0
      while p < end do
        acc += vData(p) * xData(xOff + rIdx(p) * xStep)
        p += 1
      yData(yOff + col * yStep) = acc
      col += 1

  override def *(x: DVec): DVec =
    val out = MutableDVec.zeros(rows)
    mulInto(x, out)
    out.asVec

  def toCSR: CSR =
    val counts = new Array[Int](rows)
    var p = 0
    while p < nnz do
      counts(rowIdx(p)) += 1
      p += 1
    val rowPtr = prefixCounts(counts)
    val next = rowPtr.clone()
    // Platform outputs filled in place; the CSR adopts them (no fromArray copy).
    val colIdxOut = IndexArray.alloc(nnz)
    val valuesOut = DoubleArray.alloc(nnz)
    var col = 0
    while col < cols do
      p = colPtr(col)
      val end = colPtr(col + 1)
      while p < end do
        val row = rowIdx(p)
        val dest = next(row)
        colIdxOut(dest) = col
        valuesOut(dest) = values(p)
        next(row) += 1
        p += 1
      col += 1
    new CSR(rows, cols, toIndexArray(rowPtr, rows + 1), colIdxOut, valuesOut)

  def toCSC: CSC =
    this

  // Canonicalization mirrors CSR through the zero-copy transpose: `t` is this CSC
  // reinterpreted as a CSR whose rows are this matrix's columns, so sorting,
  // summing, and pruning within CSR rows is exactly the CSC column contract. Each
  // result is transposed back to CSC; the CSR work builds fresh arrays and never
  // mutates shared storage, so the reinterpretation stays safe.

  /** True iff row indices are strictly increasing within every column (sorted, no
    * duplicates) and no stored value is an explicit zero — the canonical CSC
    * contract.
    */
  def hasCanonicalFormat: Boolean =
    t.hasCanonicalFormat

  /** Order the row indices ascending within each column, carrying values along.
    * Duplicates and explicit zeros are preserved.
    */
  def sortedIndices: CSC =
    t.sortedIndices.t

  /** Drop entries whose value is an exact zero, preserving order and structure. */
  def pruneZeros: CSC =
    t.pruneZeros.t

  /** Drop entries whose magnitude is `<= absBelow`, preserving survivor order. */
  def prune(absBelow: Double): CSC =
    t.prune(absBelow).t

  /** Fully canonicalize: sort rows within each column, sum duplicate rows, and
    * drop entries that are (or sum to) zero. Satisfies [[hasCanonicalFormat]].
    */
  def canonicalize: CSC =
    t.canonicalize.t

  // Elementwise arithmetic reuses the CSR implementations through the zero-copy
  // transpose: `(A op B) == (Aᵀ op Bᵀ)ᵀ`. Each CSR op builds fresh, canonical
  // arrays, so transposing the result back to CSC is safe and stays canonical.

  def +(that: CSC): CSC =
    (t + that.t).t

  def -(that: CSC): CSC =
    (t - that.t).t

  def *(alpha: Double): CSC =
    (t * alpha).t

  def mapValues(f: Double => Double): CSC =
    t.mapValues(f).t

  def zipValues(that: CSC)(f: (Double, Double) => Double): CSC =
    t.zipValues(that.t)(f).t

  /** Dense materialisation by column scatter, producing a contiguous row-major
    * matrix directly (rather than a transposed view of the CSR densification).
    */
  def toDense(maxEntries: Int = Int.MaxValue): DMat =
    val entries = rows.toLong * cols.toLong
    if entries > maxEntries.toLong then
      throw LinAlgError.UnsupportedOperation(s"dense conversion would allocate $entries entries")
    val out = DMat.zeros(rows, cols)
    val cPtr = colPtr
    val rIdx = rowIdx
    val vData = values
    val n = cols
    var col = 0
    while col < cols do
      var p = cPtr(col)
      val end = cPtr(col + 1)
      while p < end do
        out.data(rIdx(p) * n + col) = vData(p)
        p += 1
      col += 1
    out

  /** Main diagonal. `Aᵀ` shares `A`'s diagonal, so this reuses the CSR path. */
  def diagonal: DVec =
    t.diagonal

  def trace: Double =
    t.trace

final class Diagonal private[gale] (private val diagonal: DoubleArray)
    extends SparseMatrix[Double]
    with DoubleLinearOperator:
  def rows: Int =
    diagonal.length

  def cols: Int =
    diagonal.length

  def nnz: Int =
    var count = 0
    var i = 0
    while i < diagonal.length do
      if diagonal(i) != 0.0 then count += 1
      i += 1
    count

  def apply(row: Int, col: Int): Double =
    check(row, rows)
    check(col, cols)
    if row == col then diagonal(row) else 0.0

  def row(index: Int): DVec =
    check(index, rows)
    val out = DVec.zeros(cols)
    out.data(index) = diagonal(index)
    out

  def col(index: Int): DVec =
    row(index)

  def t: Diagonal =
    this

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    if x.length != cols then throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if into.length != rows then throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(into.length), Cols(1)))
    val intoData = into.data
    val intoOff = into.offset.value
    val intoStep = into.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val diag = diagonal
    var i = 0
    while i < diag.length do
      intoData(intoOff + i * intoStep) = diag(i) * xData(xOff + i * xStep)
      i += 1

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    applyTo(x, into)

  def toCSR: CSR =
    val builder = Sparse.coo(rows, cols)
    var i = 0
    while i < diagonal.length do
      if diagonal(i) != 0.0 then builder.add(i, i, diagonal(i))
      i += 1
    builder.toCSR()

final class Identity private[gale] (val rows: Int)
    extends SparseMatrix[Double]
    with DoubleLinearOperator:
  def cols: Int =
    rows

  def nnz: Int =
    rows

  def apply(row: Int, col: Int): Double =
    check(row, rows)
    check(col, cols)
    if row == col then 1.0 else 0.0

  def row(index: Int): DVec =
    check(index, rows)
    val out = DVec.zeros(cols)
    out.data(index) = 1.0
    out

  def col(index: Int): DVec =
    row(index)

  def t: Identity =
    this

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    if x.length != cols then throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if into.length != rows then throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(into.length), Cols(1)))
    val intoData = into.data
    val intoOff = into.offset.value
    val intoStep = into.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    var i = 0
    while i < rows do
      intoData(intoOff + i * intoStep) = xData(xOff + i * xStep)
      i += 1

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    applyTo(x, into)

final class Zero private[gale] (val rows: Int, val cols: Int)
    extends SparseMatrix[Double]
    with DoubleLinearOperator:
  require(rows >= 0 && cols >= 0, "zero matrix shape must be non-negative")

  def nnz: Int =
    0

  def apply(row: Int, col: Int): Double =
    check(row, rows)
    check(col, cols)
    0.0

  def row(index: Int): DVec =
    check(index, rows)
    DVec.zeros(cols)

  def col(index: Int): DVec =
    check(index, cols)
    DVec.zeros(rows)

  def t: Zero =
    new Zero(cols, rows)

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    if x.length != cols then throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if into.length != rows then throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(into.length), Cols(1)))
    into.clear()

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    t.applyTo(x, into)

final class Permutation private[gale] (private val columnsByRow: Array[Int])
    extends SparseMatrix[Double]
    with DoubleLinearOperator:
  val rows: Int =
    columnsByRow.length

  def cols: Int =
    rows

  def nnz: Int =
    rows

  def apply(row: Int, col: Int): Double =
    check(row, rows)
    check(col, cols)
    if columnsByRow(row) == col then 1.0 else 0.0

  def row(index: Int): DVec =
    check(index, rows)
    val out = DVec.zeros(cols)
    out.data(columnsByRow(index)) = 1.0
    out

  def col(index: Int): DVec =
    check(index, cols)
    val out = DVec.zeros(rows)
    var row = 0
    while row < rows do
      if columnsByRow(row) == index then out.data(row) = 1.0
      row += 1
    out

  def t: Permutation =
    val inverse = new Array[Int](columnsByRow.length)
    var row = 0
    while row < columnsByRow.length do
      inverse(columnsByRow(row)) = row
      row += 1
    new Permutation(inverse)

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    if x.length != cols then throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if into.length != rows then throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(into.length), Cols(1)))
    val intoData = into.data
    val intoOff = into.offset.value
    val intoStep = into.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    val colsByRow = columnsByRow
    var row = 0
    while row < rows do
      intoData(intoOff + row * intoStep) = xData(xOff + colsByRow(row) * xStep)
      row += 1

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    t.applyTo(x, into)

  /** The permutation as row → column: `toArray(row)` is the column carrying the 1
    * in that row. A fresh `Array[Int]` copy (safe to mutate); see [[toIndexSeq]]
    * for an immutable snapshot. Index arrays stay public under P4's storage rule,
    * which restricts only `Double` storage to the interop doorway.
    */
  def toArray: Array[Int] =
    columnsByRow.clone()

  /** An immutable row → column snapshot of the permutation. */
  def toIndexSeq: IndexedSeq[Int] =
    columnsByRow.toIndexedSeq

object CSR:
  private[gale] def fromCOO(coo: COO): CSR =
    val rowIndices = coo.rawRows
    val colIndices = coo.rawCols
    val entryValues = coo.rawValues
    val counts = new Array[Int](coo.rows)
    var i = 0
    while i < rowIndices.length do
      counts(rowIndices(i)) += 1
      i += 1
    val rowPtr = prefixCounts(counts)
    val next = rowPtr.clone()
    // Counting-sort straight into platform storage; the CSR adopts these
    // outputs, so the only copies are the defensive COO reads above.
    val colIdx = IndexArray.alloc(entryValues.length)
    val values = DoubleArray.alloc(entryValues.length)
    i = 0
    while i < entryValues.length do
      val row = rowIndices(i)
      val dest = next(row)
      colIdx(dest) = colIndices(i)
      values(dest) = entryValues(i)
      next(row) += 1
      i += 1
    new CSR(coo.rows, coo.cols, toIndexArray(rowPtr, coo.rows + 1), colIdx, values)

object CSC:
  private[gale] def fromCOO(coo: COO): CSC =
    coo.toCSR.toCSC

object Sparse:
  def coo(rows: Int, cols: Int): COOBuilder =
    new COOBuilder(rows, cols)

  def diagonal(values: Double*): Diagonal =
    new Diagonal(DoubleArray.fromArray(values.toArray))

  def identity(size: Int): Identity =
    require(size >= 0, "size must be non-negative")
    new Identity(size)

  def zero(rows: Int, cols: Int): Zero =
    new Zero(rows, cols)

  def permutation(columnsByRow: Int*): Permutation =
    val seen = Array.fill(columnsByRow.length)(false)
    var i = 0
    while i < columnsByRow.length do
      val col = columnsByRow(i)
      if col < 0 || col >= columnsByRow.length then
        throw LinAlgError.IndexOutOfBounds(col, columnsByRow.length)
      if seen(col) then
        throw LinAlgError.InvalidArgument(s"duplicate permutation target $col")
      seen(col) = true
      i += 1
    new Permutation(columnsByRow.toArray)

  def fromCOO(rows: Int, cols: Int, entries: Seq[(Int, Int, Double)]): COO =
    val builder = coo(rows, cols)
    entries.foreach { case (row, col, value) => builder.add(row, col, value) }
    builder.toCOO()

  /** A [[Banded]] matrix from diagonals keyed by offset (`0` main, `d>0` super,
    * `d<0` sub). See [[Banded.fromDiagonals]].
    */
  def banded(rows: Int, cols: Int, diagonals: Map[Int, Seq[Double]]): Banded =
    Banded.fromDiagonals(rows, cols, diagonals)

  /** A [[Banded]] matrix packed from a dense matrix, detecting its bandwidth. */
  def banded(dense: DMat): Banded =
    Banded.fromDense(dense)

private def check(index: Int, bound: Int): Unit =
  if index < 0 || index >= bound then
    throw LinAlgError.IndexOutOfBounds(index, bound)

/** Copy the first `length` entries of a plain builder array into an exact-size
  * platform [[IndexArray]] (typed-array-backed on JS, `Array[Int]` on the JVM).
  */
private def toIndexArray(src: Array[Int], length: Int): IndexArray =
  val out = IndexArray.alloc(length)
  var i = 0
  while i < length do
    out(i) = src(i)
    i += 1
  out

/** Copy the first `length` entries of a plain builder array into an exact-size
  * platform [[DoubleArray]].
  */
private def toDoubleArray(src: Array[Double], length: Int): DoubleArray =
  val out = DoubleArray.alloc(length)
  var i = 0
  while i < length do
    out(i) = src(i)
    i += 1
  out

/** Stable insertion sort of the parallel `(keys, vals)` arrays over `[start, end)`
  * by ascending `keys`. Insertion sort suits the short per-row/per-column slices
  * that canonicalization operates on and keeps equal keys in their original order,
  * so summed duplicates are deterministic.
  */
private def insertionSortRange(keys: Array[Int], vals: Array[Double], start: Int, end: Int): Unit =
  var i = start + 1
  while i < end do
    val k = keys(i)
    val v = vals(i)
    var j = i - 1
    while j >= start && keys(j) > k do
      keys(j + 1) = keys(j)
      vals(j + 1) = vals(j)
      j -= 1
    keys(j + 1) = k
    vals(j + 1) = v
    i += 1

private def prefixCounts(counts: Array[Int]): Array[Int] =
  val out = new Array[Int](counts.length + 1)
  var i = 0
  var sum = 0
  while i < counts.length do
    out(i) = sum
    sum += counts(i)
    i += 1
  out(counts.length) = sum
  out
