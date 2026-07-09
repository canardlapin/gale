package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import scala.annotation.targetName

trait Matrix[A] extends LinearOperator[A]:
  def rows: Int
  def cols: Int
  def shape: Shape = Shape(Rows(rows), Cols(cols))
  def apply(row: Int, col: Int): A
  def row(index: Int): Vec[A]
  def col(index: Int): Vec[A]
  def t: Matrix[A]

/** A dense `Double` matrix: a rows×cols window over platform storage described by
  * an offset and independent (strictly positive) row and column strides.
  *
  *   - '''Views vs copies.''' [[row]], [[col]], and [[t]] (transpose) are `O(1)`
  *     aliasing views sharing the backing storage — `t` just swaps the strides.
  *     [[updated]] and the `Array`/`Seq` exporters return independent data. `DMat`
  *     exposes no element mutators; the write path is `mulInto`-style APIs that
  *     take a caller-supplied [[MutableDVec]] destination.
  *   - '''Positive strides only,''' arbitrary otherwise: element `(i, j)` lives at
  *     `offset + i*rowStride + j*colStride`, so row-major, column-major, and
  *     strided submatrix layouts are all valid inputs to the kernels.
  *   - '''NaN / `beta == 0` policy.''' When a kernel forms `y := alpha*A*x` with an
  *     implicit `beta == 0` (as `mulInto`/`*` do), the destination is '''assigned''',
  *     never read-and-scaled — so a pre-existing `NaN`/`Inf` in the destination
  *     buffer cannot poison the result via `0*NaN`.
  */
final class DMat private[gale] (
    private[gale] val data: DoubleArray,
    private[gale] val offsetValue: Offset,
    private[gale] val rowsValue: Rows,
    private[gale] val colsValue: Cols,
    private[gale] val rowStrideValue: Stride,
    private[gale] val colStrideValue: Stride
) extends Matrix[Double]
    with DoubleLinearOperator:
  def rows: Int = rowsValue.value
  def cols: Int = colsValue.value
  def offset: Offset = offsetValue
  def rowStride: Stride = rowStrideValue
  def colStride: Stride = colStrideValue
  override def shape: Shape = Shape(rowsValue, colsValue)

  /** Rows are stored contiguously: walking a row touches adjacent storage. */
  def isContiguousRowMajor: Boolean =
    colStride.value == 1 && rowStride.value == cols

  /** Columns are stored contiguously: walking a column touches adjacent storage. */
  def isContiguousColMajor: Boolean =
    rowStride.value == 1 && colStride.value == rows

  def apply(row: Int, col: Int): Double =
    checkRow(row)
    checkCol(col)
    data(index(row, col))

  override def row(index: Int): DVec =
    checkRow(index)
    new DVec(
      data,
      Offset.unsafe(offset.value + index * rowStride.value),
      Length.unsafe(cols),
      colStride
    )

  override def col(index: Int): DVec =
    checkCol(index)
    new DVec(
      data,
      Offset.unsafe(offset.value + index * colStride.value),
      Length.unsafe(rows),
      rowStride
    )

  override def t: DMat =
    new DMat(
      data,
      offset,
      Rows.unsafe(cols),
      Cols.unsafe(rows),
      colStride,
      rowStride
    )

  def updated(row: Int, col: Int, value: Double): DMat =
    checkRow(row)
    checkCol(col)
    val out = toDoubleArrayCopyRowMajor
    out(row * cols + col) = value
    DMat.fromDoubleArrayOwned(rows, cols, out)

  private[gale] def toArrayRowMajor: Array[Double] =
    val out = new Array[Double](rows * cols)
    var i = 0
    while i < rows do
      var j = 0
      var aij = offset.value + i * rowStride.value
      while j < cols do
        out(i * cols + j) = data(aij)
        aij += colStride.value
        j += 1
      i += 1
    out

  /** Immutable copy of the elements in row-major logical order. */
  def valuesRowMajor: Seq[Double] =
    val builder = Vector.newBuilder[Double]
    builder.sizeHint(rows * cols)
    var i = 0
    while i < rows do
      var j = 0
      var aij = offset.value + i * rowStride.value
      while j < cols do
        builder += data(aij)
        aij += colStride.value
        j += 1
      i += 1
    builder.result()

  /** Contiguous row-major platform copy owned by the caller (offset 0,
    * colStride 1). One copy, even from a strided or transposed view.
    */
  private[gale] def toDoubleArrayCopyRowMajor: DoubleArray =
    val out = DoubleArray.alloc(rows * cols)
    var i = 0
    while i < rows do
      var j = 0
      var aij = offset.value + i * rowStride.value
      while j < cols do
        out(i * cols + j) = data(aij)
        aij += colStride.value
        j += 1
      i += 1
    out

  override def *(x: DVec): DVec =
    if cols != x.length then
      throw LinAlgError.VectorLengthMismatch(cols, x.length)
    val out = MutableDVec.zeros(rows)
    mulInto(x, out)
    out.asVec

  def mulInto(x: DVec, y: MutableDVec): Unit =
    if cols != x.length then
      throw LinAlgError.VectorLengthMismatch(cols, x.length)
    if rows != y.length then
      throw LinAlgError.VectorLengthMismatch(rows, y.length)
    if DoubleArray.sameStorage(x.data, y.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val rowStep = rowStride.value
    val colStep = colStride.value
    val xStep = x.stride.value
    val yStep = y.stride.value
    if colStep == 1 && xStep == 1 then
      DoubleKernels.dgemvRowMajor(
        rows,
        cols,
        1.0,
        data,
        offset.value,
        rowStep,
        x.data,
        x.offset.value,
        0.0,
        y.data,
        y.offset.value,
        yStep
      )

    else if rowStep == 1 then
      DoubleKernels.dgemvColMajor(
        rows,
        cols,
        1.0,
        data,
        offset.value,
        colStep,
        x.data,
        x.offset.value,
        xStep,
        0.0,
        y.data,
        y.offset.value,
        yStep
      )
    else
      DoubleKernels.dgemv(
        rows,
        cols,
        1.0,
        data,
        offset.value,
        rowStep,
        colStep,
        x.data,
        x.offset.value,
        xStep,
        0.0,
        y.data,
        y.offset.value,
        yStep
      )

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    mulInto(x, into)

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    t.mulInto(x, into)

  def asLinearOperator: DoubleLinearOperator =
    this

  def *(that: DMat): DMat =
    if cols != that.rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(that.cols)), Shape(Rows(that.rows), Cols(that.cols)))
    val out = DMat.zeros(rows, that.cols)
    DoubleKernels.dgemm(
      rows,
      that.cols,
      cols,
      1.0,
      data,
      offset.value,
      rowStride.value,
      colStride.value,
      that.data,
      that.offset.value,
      that.rowStride.value,
      that.colStride.value,
      0.0,
      out.data,
      out.offset.value,
      out.rowStride.value,
      out.colStride.value
    )
    out

  def +(that: DMat): DMat =
    requireSameShape(that)
    addSub(that, subtract = false)

  def -(that: DMat): DMat =
    requireSameShape(that)
    addSub(that, subtract = true)

  /** Elementwise add/subtract through the `dadd`/`dsub` kernels. When both
    * operands are contiguous row-major the whole block is a single kernel call;
    * otherwise each row is one strided call, honouring arbitrary layouts.
    */
  private def addSub(that: DMat, subtract: Boolean): DMat =
    val out = DMat.zeros(rows, cols)
    val outData = out.data
    if isContiguousRowMajor && that.isContiguousRowMajor then
      val n = rows * cols
      if subtract then
        DoubleKernels.dsub(n, data, offset.value, 1, that.data, that.offset.value, 1, outData, 0, 1)
      else
        DoubleKernels.dadd(n, data, offset.value, 1, that.data, that.offset.value, 1, outData, 0, 1)
    else
      val ncols = cols
      val aColStep = colStride.value
      val bColStep = that.colStride.value
      val aRowStep = rowStride.value
      val bRowStep = that.rowStride.value
      var i = 0
      var aRow = offset.value
      var bRow = that.offset.value
      var outRow = 0
      if subtract then
        while i < rows do
          DoubleKernels.dsub(ncols, data, aRow, aColStep, that.data, bRow, bColStep, outData, outRow, 1)
          aRow += aRowStep
          bRow += bRowStep
          outRow += ncols
          i += 1
      else
        while i < rows do
          DoubleKernels.dadd(ncols, data, aRow, aColStep, that.data, bRow, bColStep, outData, outRow, 1)
          aRow += aRowStep
          bRow += bRowStep
          outRow += ncols
          i += 1
    out

  def lu: Either[LinAlgError, LU] =
    DenseDecompositions.lu(this)

  def cholesky: Either[LinAlgError, Cholesky] =
    DenseDecompositions.cholesky(this)

  def qr: QR =
    DenseDecompositions.qr(this)

  def qrWith(workspace: DenseWorkspace): QR =
    DenseDecompositions.qr(this, workspace)

  def leastSquares(b: DVec): Either[LinAlgError, DVec] =
    qr.solveLeastSquares(b)

  def rankEstimate: Int =
    DenseDecompositions.rankEstimate(this)

  def conditionEstimate: Either[LinAlgError, Double] =
    DenseDecompositions.conditionEstimate(this)

  def det: Either[LinAlgError, Double] =
    lu.flatMap(_.det)

  def solve(b: DVec): Either[LinAlgError, DVec] =
    DenseDecompositions.solve(this, b)

  private inline def index(row: Int, col: Int): Int =
    offset.value + row * rowStride.value + col * colStride.value

  private def requireSameShape(that: DMat): Unit =
    if rows != that.rows || cols != that.cols then
      throw LinAlgError.DimensionMismatch(shape, that.shape)

  private def checkRow(row: Int): Unit =
    if row < 0 || row >= rows then
      throw LinAlgError.IndexOutOfBounds(row, rows)

  private def checkCol(col: Int): Unit =
    if col < 0 || col >= cols then
      throw LinAlgError.IndexOutOfBounds(col, cols)

object DMat:
  /** Reject shapes whose element count is negative or overflows an Int index,
    * so we never allocate a wrapped-around buffer while claiming the full shape.
    */
  private def requireStorable(rows: Int, cols: Int): Unit =
    require(rows >= 0 && cols >= 0, "rows and cols must be non-negative")
    if rows.toLong * cols.toLong > Int.MaxValue.toLong then
      throw LinAlgError.InvalidArgument(
        s"matrix size ${rows}x${cols} exceeds ${Int.MaxValue} storable elements"
      )

  def zeros(rows: Int, cols: Int): DMat =
    requireStorable(rows, cols)
    new DMat(
      DoubleArray.alloc(rows * cols),
      Offset.unsafe(0),
      Rows.unsafe(rows),
      Cols.unsafe(cols),
      Stride.unsafe(if cols == 0 then 1 else cols),
      Stride.unsafe(1)
    )

  def eye(size: Int): DMat =
    require(size >= 0, "size must be non-negative")
    val out = zeros(size, size)
    var i = 0
    while i < size do
      out.data(i * size + i) = 1.0
      i += 1
    out

  def dense(rows: Int, cols: Int, values: Seq[Double]): DMat =
    requireStorable(rows, cols)
    require(values.length == rows * cols, s"expected ${rows * cols} values, got ${values.length}")
    val out = zeros(rows, cols)
    var i = 0
    values.foreach { value =>
      out.data(i) = value
      i += 1
    }
    out

  private[gale] def fromArrayRowMajor(rows: Int, cols: Int, values: Array[Double]): DMat =
    requireStorable(rows, cols)
    require(values.length == rows * cols, s"expected ${rows * cols} values, got ${values.length}")
    fromDoubleArrayOwned(rows, cols, DoubleArray.fromArray(values))

  /** Wrap a contiguous row-major platform array as a matrix without copying;
    * the caller transfers ownership of `data` and must not mutate it afterwards.
    */
  private[gale] def fromDoubleArrayOwned(rows: Int, cols: Int, data: DoubleArray): DMat =
    requireStorable(rows, cols)
    require(data.length == rows * cols, s"expected ${rows * cols} values, got ${data.length}")
    new DMat(
      data,
      Offset.unsafe(0),
      Rows.unsafe(rows),
      Cols.unsafe(cols),
      Stride.unsafe(if cols == 0 then 1 else cols),
      Stride.unsafe(1)
    )

  def tabulate(rows: Int, cols: Int)(f: (Int, Int) => Double): DMat =
    val out = zeros(rows, cols)
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        out.data(i * cols + j) = f(i, j)
        j += 1
      i += 1
    out

object Matrix:
  def zeros(rows: Int, cols: Int): DMat =
    DMat.zeros(rows, cols)

  def eye(size: Int): DMat =
    DMat.eye(size)

  def dense(rows: Int, cols: Int, values: Seq[Double]): DMat =
    DMat.dense(rows, cols, values)

  @targetName("denseVarargs")
  def dense(rows: Int, cols: Int)(values: Double*): DMat =
    DMat.dense(rows, cols, values)

  def tabulate(rows: Int, cols: Int)(f: (Int, Int) => Double): DMat =
    DMat.tabulate(rows, cols)(f)
