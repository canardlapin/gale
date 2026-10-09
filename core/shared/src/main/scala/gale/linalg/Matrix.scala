package gale.linalg

import gale.backend.Backend
import gale.backend.PureBackend
import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.spectral.SVD
import gale.spectral.SingularSelection
import gale.spectral.SpectralBackend
import gale.spectral.Svds
import gale.spectral.SvdCutoff
import gale.spectral.TruncatedSvd
import gale.spectral.MinimumNormSolution
import gale.spectral.MinimumNormSolutions
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
  *   - '''Views vs copies.''' [[row]], [[col]], [[slice]], and [[t]] (transpose)
  *     are `O(1)` aliasing views sharing the backing storage — `t` just swaps the
  *     strides. [[gatherRows]], [[gatherColumns]], [[updated]], and [[toBuilder]]
  *     return independently owned data. `DMat` exposes no element mutators; the
  *     write path is `mulInto`-style APIs that take a caller-supplied
  *     [[MutableDVec]] destination.
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

  /** Contiguous half-open row/column window as an `O(1)` strided view.
    *
    * The result aliases this matrix storage and retains its row/column strides;
    * use [[gatherRows]], [[gatherColumns]], or [[toBuilder]] when an owned copy is
    * required. Empty windows are supported.
    */
  def slice(rowFrom: Int, rowUntil: Int, colFrom: Int, colUntil: Int): DMat =
    if rowFrom < 0 || rowUntil < rowFrom || rowUntil > rows ||
      colFrom < 0 || colUntil < colFrom || colUntil > cols
    then
      throw LinAlgError.InvalidArgument(
        s"invalid matrix slice rows [$rowFrom, $rowUntil), columns [$colFrom, $colUntil) for ${rows}x${cols}"
      )
    new DMat(
      data,
      Offset.unsafe(offset.value + rowFrom * rowStride.value + colFrom * colStride.value),
      Rows.unsafe(rowUntil - rowFrom),
      Cols.unsafe(colUntil - colFrom),
      rowStride,
      colStride
    )

  /** Gather rows in caller-specified order into one independently owned,
    * contiguous row-major matrix. Repeated indices repeat rows.
    */
  def gatherRows(indices: IndexedSeq[Int]): DMat =
    var i = 0
    while i < indices.length do
      checkRow(indices(i))
      i += 1
    val out = DMatBuilder.zeros(indices.length, cols)
    var row = 0
    while row < indices.length do
      val sourceRow = indices(row)
      var col = 0
      while col < cols do
        out(row, col) = apply(sourceRow, col)
        col += 1
      row += 1
    out.result()

  /** Gather columns in caller-specified order into one independently owned,
    * contiguous row-major matrix. Repeated indices repeat columns.
    */
  def gatherColumns(indices: IndexedSeq[Int]): DMat =
    var i = 0
    while i < indices.length do
      checkCol(indices(i))
      i += 1
    val out = DMatBuilder.zeros(rows, indices.length)
    var row = 0
    while row < rows do
      var col = 0
      while col < indices.length do
        out(row, col) = apply(row, indices(col))
        col += 1
      row += 1
    out.result()

  def updated(row: Int, col: Int, value: Double): DMat =
    checkRow(row)
    checkCol(col)
    val out = toDoubleArrayCopyRowMajor
    out(row * cols + col) = value
    DMat.fromDoubleArrayOwned(rows, cols, out)

  /** One owned logical row-major copy, ready for in-place construction edits.
    * Calling [[DMatBuilder.result]] transfers that storage without another copy.
    */
  def toBuilder: DMatBuilder =
    DMatBuilder.from(this)

  /** Add `amount` to every element on the main diagonal, returning one owned
    * copy. Rectangular matrices update `min(rows, cols)` diagonal elements.
    */
  def addToDiagonal(amount: Double): DMat =
    val out = toBuilder
    var i = 0
    while i < math.min(rows, cols) do
      out(i, i) = out(i, i) + amount
      i += 1
    out.result()

  /** Average a square matrix with its transpose:
    * `B(i,j) = 0.5*A(i,j) + 0.5*A(j,i)`.
    *
    * The explicit name distinguishes this from lower- or upper-triangle
    * symmetrization conventions used by factorization kernels.
    */
  def symmetrizedAverage: DMat =
    if rows != cols then throw LinAlgError.NonSquareMatrix(shape)
    val out = DMatBuilder.zeros(rows, cols)
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        out(i, j) = 0.5 * apply(i, j) + 0.5 * apply(j, i)
        j += 1
      i += 1
    out.result()

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

  /** Visit logical elements in row-major order without allocating an
    * intermediate collection. The callback receives primitive `Double` values.
    */
  def foreachRowMajor(f: Double => Unit): Unit =
    var i = 0
    while i < rows do
      var j = 0
      var aij = offset.value + i * rowStride.value
      while j < cols do
        f(data(aij))
        aij += colStride.value
        j += 1
      i += 1

  /** Copy logical elements into caller-owned primitive storage in row-major
    * order. No Gale backing storage is exposed or adopted.
    */
  def copyRowMajorTo(destination: Array[Double], destinationOffset: Int = 0): Unit =
    val required = rows.toLong * cols.toLong
    val end = destinationOffset.toLong + required
    if destinationOffset < 0 || end > destination.length.toLong then
      throw LinAlgError.InvalidArgument(
        s"destination range [$destinationOffset, $end) does not fit array length ${destination.length}"
      )
    var write = destinationOffset
    var i = 0
    while i < rows do
      var j = 0
      var aij = offset.value + i * rowStride.value
      while j < cols do
        destination(write) = data(aij)
        write += 1
        aij += colStride.value
        j += 1
      i += 1

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

  override def *(x: DVec)(using backend: Backend): DVec =
    if cols != x.length then
      throw LinAlgError.VectorLengthMismatch(cols, x.length)
    val out = MutableDVec.zeros(rows)
    mulInto(x, out)
    out.asVec

  def mulInto(x: DVec, y: MutableDVec)(using backend: Backend): Unit =
    mulIntoWithBackend(x, y, backend)

  private def mulIntoWithBackend(x: DVec, y: MutableDVec, backend: Backend): Unit =
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
    if backend.acceleratesGemv &&
      rows.toLong * cols.toLong >= backend.thresholds.nativeGemvMinWork
    then
      backend.denseDouble.gemv(
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
    else if colStep == 1 && xStep == 1 then
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
    // LinearOperator application is also used inside iterative hot loops; keep that
    // contract statically pure. The standalone `DMat * DVec` facade is the seam.
    mulIntoWithBackend(x, into, PureBackend)

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    t.mulIntoWithBackend(x, into, PureBackend)

  def asLinearOperator: DoubleLinearOperator =
    this

  /** Matrix product `this * that`. The coarse gemm dispatch seam (doc §A.2): with no
    * acceleration import the `given` resolves to [[gale.backend.Backend.pure]], whose
    * `acceleratesGemm` is `false`, so the pure `DoubleKernels` path runs — byte-identical
    * to before the seam. An imported accelerating backend routes the general product to
    * `backend.denseDouble.gemm` once the work clears its measured `nativeGemmMinFlops`.
    * The structural `AᵀA` fast-path routes to `backend.denseDouble.syrk` under the same
    * gate, and stays on the dedicated pure symmetric kernel below it.
    */
  def *(that: DMat)(using backend: Backend): DMat =
    requireProductShape(that)
    val out = DMat.zeros(rows, that.cols)
    gemmToStorage(
      that,
      out.data,
      out.offset.value,
      out.rowStride.value,
      out.colStride.value,
      alpha = 1.0,
      beta = 0.0,
      backend = backend
    )
    out

  /** General matrix product into caller-owned single-owner storage:
    * `destination := alpha * this * that + beta * destination`.
    *
    * `beta == 0.0` is replacement semantics: every destination element is
    * assigned without reading its old value, so pre-existing NaN or infinity
    * cannot poison the result. Any nonzero `beta` is accumulation semantics and
    * reads the previous destination. Inputs may be row-major, transposed, or
    * strided; the destination is the open contiguous row-major storage owned by
    * a [[DMatBuilder]]. No intermediate matrix is allocated.
    */
  def gemmInto(
      that: DMat,
      destination: DMatBuilder,
      alpha: Double = 1.0,
      beta: Double = 0.0
  )(using backend: Backend): Unit =
    val out = destination.writableData
    requireProductShape(that)
    requireDestinationShape(destination, rows, that.cols)
    if DoubleArray.sameStorage(data, out) || DoubleArray.sameStorage(that.data, out) then
      throw LinAlgError.UnsupportedOperation("aliased gemmInto destination")
    gemmToStorage(
      that,
      out,
      outOffset = 0,
      outRowStride = destination.cols,
      outColStride = 1,
      alpha = alpha,
      beta = beta,
      backend = backend
    )

  /** Fused elementwise replacement:
    * `destination := alpha * this + beta * that`.
    *
    * Unlike GEMM's `beta`, both coefficients scale immutable inputs; the old
    * destination is never read. A zero coefficient also suppresses reading its
    * corresponding input, so `0 * NaN` cannot contaminate the result. Structural
    * aliasing with either input is rejected even though public builders normally
    * make such an alias impossible.
    */
  def linearCombinationInto(
      that: DMat,
      destination: DMatBuilder,
      alpha: Double,
      beta: Double
  ): Unit =
    val out = destination.writableData
    requireSameShape(that)
    requireDestinationShape(destination, rows, cols)
    if DoubleArray.sameStorage(data, out) || DoubleArray.sameStorage(that.data, out) then
      throw LinAlgError.UnsupportedOperation("aliased linearCombinationInto destination")

    val aData = data
    val bData = that.data
    val aRowStride = rowStride.value
    val aColStride = colStride.value
    val bRowStride = that.rowStride.value
    val bColStride = that.colStride.value
    var row = 0
    var write = 0
    if alpha == 0.0 && beta == 0.0 then
      while write < rows * cols do
        out(write) = 0.0
        write += 1
    else if beta == 0.0 then
      while row < rows do
        var aIndex = offset.value + row * aRowStride
        var col = 0
        while col < cols do
          out(write) = alpha * aData(aIndex)
          write += 1
          aIndex += aColStride
          col += 1
        row += 1
    else if alpha == 0.0 then
      while row < rows do
        var bIndex = that.offset.value + row * bRowStride
        var col = 0
        while col < cols do
          out(write) = beta * bData(bIndex)
          write += 1
          bIndex += bColStride
          col += 1
        row += 1
    else
      while row < rows do
        var aIndex = offset.value + row * aRowStride
        var bIndex = that.offset.value + row * bRowStride
        var col = 0
        while col < cols do
          out(write) = alpha * aData(aIndex) + beta * bData(bIndex)
          write += 1
          aIndex += aColStride
          bIndex += bColStride
          col += 1
        row += 1

  private def gemmToStorage(
      that: DMat,
      out: DoubleArray,
      outOffset: Int,
      outRowStride: Int,
      outColStride: Int,
      alpha: Double,
      beta: Double,
      backend: Backend
  ): Unit =
    // `Aᵀ · A` with a row-major `A` (the common Gram/normal-equations product):
    // a general gemm would traverse the transposed left operand column-strided and
    // cache-hostile, so route it to the dedicated symmetric rank-k kernel (half the
    // flops, unit-stride throughout). Detected structurally: `this` is exactly the
    // transpose view of `that`, and `that` has unit column stride.
    if isTransposeOf(that) && that.colStride.value == 1 && alpha == 1.0 && beta == 0.0 then
      if backend.routesGemm(that.rows, that.cols, that.cols) then
        backend.denseDouble.syrk(
          that.rows,
          that.cols,
          that.data,
          that.offset.value,
          that.rowStride.value,
          out,
          outOffset,
          outRowStride
        )
      else
        DoubleKernels.dsyrkRowMajor(
          that.rows,
          that.cols,
          that.data,
          that.offset.value,
          that.rowStride.value,
          out,
          outOffset,
          outRowStride
        )
    else if backend.routesGemm(rows, that.cols, cols) then
      backend.denseDouble.gemm(
        rows,
        that.cols,
        cols,
        alpha,
        data,
        offset.value,
        rowStride.value,
        colStride.value,
        that.data,
        that.offset.value,
        that.rowStride.value,
        that.colStride.value,
        beta,
        out,
        outOffset,
        outRowStride,
        outColStride
      )
    else
      DoubleKernels.dgemm(
        rows,
        that.cols,
        cols,
        alpha,
        data,
        offset.value,
        rowStride.value,
        colStride.value,
        that.data,
        that.offset.value,
        that.rowStride.value,
        that.colStride.value,
        beta,
        out,
        outOffset,
        outRowStride,
        outColStride
      )

  private def requireProductShape(that: DMat): Unit =
    if cols != that.rows then
      throw LinAlgError.DimensionMismatch(
        Shape(Rows(cols), Cols(that.cols)),
        Shape(Rows(that.rows), Cols(that.cols))
      )

  private def requireDestinationShape(destination: DMatBuilder, expectedRows: Int, expectedCols: Int): Unit =
    if destination.rows != expectedRows || destination.cols != expectedCols then
      throw LinAlgError.DimensionMismatch(
        Shape(Rows(expectedRows), Cols(expectedCols)),
        Shape(Rows(destination.rows), Cols(destination.cols))
      )

  /** True when `this` is exactly the transpose view of `that` — same backing
    * storage and offset, swapped shape and strides — so `this * that` is `AᵀA`.
    */
  private def isTransposeOf(that: DMat): Boolean =
    DoubleArray.sameStorage(data, that.data) &&
      offset.value == that.offset.value &&
      rows == that.cols && cols == that.rows &&
      rowStride.value == that.colStride.value &&
      colStride.value == that.rowStride.value

  def +(that: DMat): DMat =
    requireSameShape(that)
    addSub(that, subtract = false)

  def -(that: DMat): DMat =
    requireSameShape(that)
    addSub(that, subtract = true)

  /** Scale every entry by `alpha`. The result is an owned contiguous row-major
    * matrix, including when `this` is a transpose or slice view. Same arithmetic
    * as [[DVec.*]]: `0 * Inf` is `NaN`, and `1.0` still copies.
    */
  def *(alpha: Double): DMat =
    val n = rows * cols
    val out = toDoubleArrayCopyRowMajor
    if n > 0 then DoubleKernels.dscal(n, alpha, out, 0, 1)
    DMat.fromDoubleArrayOwned(rows, cols, out)

  /** Elementwise (Hadamard) product of two same-shape matrices, as an owned
    * row-major result; the kernel behind `a.pointwise * b`.
    */
  private[gale] def hadamard(that: DMat): DMat =
    requireSameShape(that)
    zipKernel(that)(DoubleKernels.dmul)

  /** Elementwise quotient of two same-shape matrices; the kernel behind
    * `a.pointwise / b`.
    */
  private[gale] def elementwiseQuotient(that: DMat): DMat =
    requireSameShape(that)
    zipKernel(that)(DoubleKernels.ddiv)

  private def addSub(that: DMat, subtract: Boolean): DMat =
    if subtract then zipKernel(that)(DoubleKernels.dsub) else zipKernel(that)(DoubleKernels.dadd)

  /** Elementwise binary op through a `dadd`-shaped kernel. When both operands
    * are contiguous row-major the whole block is a single kernel call;
    * otherwise each row is one strided call, honouring arbitrary layouts.
    */
  private inline def zipKernel(that: DMat)(
      inline kernel: (Int, DoubleArray, Int, Int, DoubleArray, Int, Int, DoubleArray, Int, Int) => Unit
  ): DMat =
    val out = DMat.zeros(rows, cols)
    val outData = out.data
    if isContiguousRowMajor && that.isContiguousRowMajor then
      kernel(rows * cols, data, offset.value, 1, that.data, that.offset.value, 1, outData, 0, 1)
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
      while i < rows do
        kernel(ncols, data, aRow, aColStep, that.data, bRow, bColStep, outData, outRow, 1)
        aRow += aRowStep
        bRow += bRowStep
        outRow += ncols
        i += 1
    out

  /** The factorization dispatch gate in one place: the backend's provider, iff it
    * accelerates factorizations and the work clears the routine's size threshold.
    * (`acceleratesFactorizations` implies `denseFactorizations.isDefined`.) Structural
    * validation stays in the facade, BEFORE this gate, so a provider only ever sees
    * inputs the pure path would accept.
    */
  private def factorizationProvider(size: Int, minSize: Int)(using backend: Backend) =
    if backend.acceleratesFactorizations && size >= minSize then backend.denseFactorizations
    else None

  def lu(using backend: Backend): Either[LinAlgError, LU] =
    if rows != cols then Left(LinAlgError.NonSquareMatrix(shape))
    else
      factorizationProvider(rows, backend.thresholds.nativeLuMinSize) match
        case Some(provider) => provider.lu(this)
        case None           => DenseDecompositions.lu(this)

  def cholesky(using backend: Backend): Either[LinAlgError, Cholesky] =
    if rows != cols then Left(LinAlgError.NonSquareMatrix(shape))
    else
      factorizationProvider(rows, backend.thresholds.nativeCholeskyMinSize) match
        case Some(provider) => provider.cholesky(this)
        case None           => DenseDecompositions.cholesky(this)

  /** Cholesky with an explicit absolute pivot tolerance. Explicit numerical
    * policy is handled by the portable implementation because native provider
    * contracts do not currently expose a matching tolerance parameter.
    */
  def cholesky(options: CholeskyOptions)(using Backend): Either[LinAlgError, Cholesky] =
    DenseDecompositions.cholesky(this, options)

  /** QR is a total facade — it always returns a `QR`, exactly as the pure Householder
    * path always succeeds. A provider that declines the input (`Left`) is therefore a
    * fallback, not a failure: the pure path computes the answer instead.
    */
  def qr(using backend: Backend): QR =
    factorizationProvider(math.max(rows, cols), backend.thresholds.nativeQrMinSize) match
      case Some(provider) => provider.qr(this).getOrElse(DenseDecompositions.qr(this))
      case None           => DenseDecompositions.qr(this)

  /** QR with explicit pivoting and rank policy. Explicit policy uses the
    * portable factorization so every platform and backend observes the same
    * permutation and rank decision.
    */
  def qr(options: QROptions)(using Backend): QR =
    DenseDecompositions.qr(this, options)

  /** QR into a caller-supplied workspace. This is the '''allocation-controlled''' facade:
    * it always runs the pure kernels against `workspace` and never routes to a native
    * provider (whose factor storage gale cannot place in the workspace) — routing would
    * silently void the reuse contract the caller allocated for. Use [[qr]] for the
    * backend-routed path.
    */
  def qrWith(workspace: DenseWorkspace)(using Backend): QR =
    DenseDecompositions.qr(this, workspace)

  def qrWith(options: QROptions, workspace: DenseWorkspace)(using Backend): QR =
    DenseDecompositions.qr(this, options, workspace)

  /** Factor the algebraic row-scaled matrix `diag(scales) * this` without
    * materializing that matrix. Zero and negative finite scales are supported;
    * non-finite scales are rejected with a typed error.
    */
  def qrScaledRows(scales: DVec)(using Backend): Either[LinAlgError, QR] =
    DenseDecompositions.qrScaledRows(this, scales, QROptions.Default, DenseWorkspace.forQR(rows, cols))

  def qrScaledRows(scales: DVec, options: QROptions)(using Backend): Either[LinAlgError, QR] =
    DenseDecompositions.qrScaledRows(this, scales, options, DenseWorkspace.forQR(rows, cols, options))

  def qrScaledRows(scales: DVec, workspace: DenseWorkspace)(using Backend): Either[LinAlgError, QR] =
    DenseDecompositions.qrScaledRows(this, scales, QROptions.Default, workspace)

  def qrScaledRows(
      scales: DVec,
      options: QROptions,
      workspace: DenseWorkspace
  )(using Backend): Either[LinAlgError, QR] =
    DenseDecompositions.qrScaledRows(this, scales, options, workspace)

  def leastSquares(b: DVec)(using Backend): Either[LinAlgError, DVec] =
    qr.solveLeastSquares(b)

  def leastSquares(b: DVec, options: QROptions)(using Backend): Either[LinAlgError, DVec] =
    qr(options).solveLeastSquares(b)

  def leastSquares(b: DMat)(using Backend): Either[LinAlgError, DMat] =
    qr.solveLeastSquares(b)

  def leastSquares(b: DMat, options: QROptions)(using Backend): Either[LinAlgError, DMat] =
    qr(options).solveLeastSquares(b)

  /** Rank from a QR factorization on the SAME dispatch policy as [[qr]]'s pure path
    * (the ambient backend drives blocked QR's internal gemms), so `rankEstimate` and
    * `qr.diagnostics.rank` agree under any imported gemm-accelerating backend.
    */
  def rankEstimate(using Backend): Int =
    DenseDecompositions.rankEstimate(this)

  /** Deliberately pinned to the pure, deterministic LU path: the estimate is
    * gemm-free, so an accelerating backend could not change its cost profile —
    * only its reproducibility.
    */
  def conditionEstimate: Either[LinAlgError, Double] =
    DenseDecompositions.conditionEstimate(this)

  def det(using Backend): Either[LinAlgError, Double] =
    lu.flatMap(_.det)

  /** Solves `A x = b` for one vector right-hand side.
    *
    * The matrix is factored once for this call. Retain [[lu]] when solving
    * several independent systems with the same matrix.
    */
  def solve(b: DVec)(using Backend): Either[LinAlgError, DVec] =
    lu.flatMap(_.solve(b))

  /** Solves `A X = B` for all columns of a matrix right-hand side.
    *
    * The matrix is factored once, then the same factor is applied to every
    * column of `b`. Shape and singularity failures remain typed
    * [[LinAlgError]] values.
    */
  def solve(b: DMat)(using Backend): Either[LinAlgError, DMat] =
    lu.flatMap(_.solve(b))

  /** Full (economy) singular value decomposition `A = U Σ Vᵀ`: singular values
    * '''descending''', `U` `m×k`, `Vᵀ` `k×n` with `k = min(m, n)`. The spectral
    * dispatch gate follows `docs/spectral-backend-boundary.md` (seam S7): a
    * [[gale.spectral.SpectralCapability.DenseSvd]]-capable `given SpectralBackend`
    * computes the raw factors, and with no import the pure bidiagonal kernel runs —
    * unlike the kernel-`Backend` factorization gates there is no size threshold, and
    * canonical order, residuals, and rank are always the facade's. `Left` on an
    * empty dimension or (in practice unreachable) kernel non-convergence.
    */
  def svd(using SpectralBackend): Either[LinAlgError, SVD] =
    Svds.svd(this, SingularSelection.All)

  /** Moore–Penrose pseudo-inverse (`n×m` for an `m×n` matrix) via the economy
    * [[svd]], on the same spectral dispatch gate. Singular values at or below
    * the MATLAB/SciPy-convention cutoff `max(m, n) · ε · σ_max` are treated as zero (see
    * [[gale.spectral.Svds.pinv]] for the exact convention), so a rank-deficient —
    * even all-zero — matrix pseudo-inverts cleanly rather than failing. `Left`
    * on non-finite input entries or an underlying [[svd]] failure.
    */
  def pinv(using SpectralBackend): Either[LinAlgError, DMat] =
    Svds.pinv(this)

  def pinv(cutoff: SvdCutoff)(using SpectralBackend): Either[LinAlgError, DMat] =
    Svds.pinv(this, cutoff)

  /** Retain one policy-selected SVD for repeated solves and subspace operations. */
  def truncatedSvd(cutoff: SvdCutoff = SvdCutoff.Default)(using SpectralBackend): Either[LinAlgError, TruncatedSvd] =
    TruncatedSvd.factor(this, cutoff)

  def minimumNormLeastSquares(b: DVec, cutoff: SvdCutoff = SvdCutoff.Default)(using
      SpectralBackend
  ): Either[LinAlgError, MinimumNormSolution] =
    Svds.minimumNormLeastSquares(this, b, cutoff)

  def minimumNormLeastSquares(b: DMat)(using SpectralBackend): Either[LinAlgError, MinimumNormSolutions] =
    Svds.minimumNormLeastSquares(this, b, SvdCutoff.Default)

  def minimumNormLeastSquares(b: DMat, cutoff: SvdCutoff)(using
      SpectralBackend
  ): Either[LinAlgError, MinimumNormSolutions] =
    Svds.minimumNormLeastSquares(this, b, cutoff)

  /** Kronecker product `this ⊗ that`: the `(m·p)×(n·q)` block matrix whose
    * `(i, j)` block is `this(i, j) * that`. Total on every shape (including
    * empty operands) like the other structural products; the only throw is the
    * standard storable-size guard when the result would exceed `Int.MaxValue`
    * elements. Strided/transposed views are read through their strides — no
    * copy of either operand.
    */
  def kron(that: DMat): DMat =
    val outRowsL = rows.toLong * that.rows.toLong
    val outColsL = cols.toLong * that.cols.toLong
    if outRowsL > Int.MaxValue.toLong || outColsL > Int.MaxValue.toLong || outRowsL * outColsL > Int.MaxValue.toLong
    then
      throw LinAlgError.InvalidArgument(
        s"Kronecker product size ${outRowsL}x$outColsL exceeds ${Int.MaxValue} storable elements"
      )
    val bRows = that.rows
    val bCols = that.cols
    val outCols = outColsL.toInt
    val out = DMat.zeros(outRowsL.toInt, outCols)
    val outData = out.data
    val bData = that.data
    val bBase = that.offset.value
    val bRowStep = that.rowStride.value
    val bColStep = that.colStride.value
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        val aij = data(index(i, j))
        val blockRow = i * bRows
        val blockCol = j * bCols
        var p = 0
        while p < bRows do
          var outIdx = (blockRow + p) * outCols + blockCol
          var bIdx = bBase + p * bRowStep
          var q = 0
          while q < bCols do
            outData(outIdx) = aij * bData(bIdx)
            outIdx += 1
            bIdx += bColStep
            q += 1
          p += 1
        j += 1
      i += 1
    out

  // ---------------------------------------------------------------------------
  // Reductions and norms. Whole-matrix reductions read a contiguous block (row-
  // or column-major) in one kernel call; otherwise they walk lines along the
  // smaller stride. Per-axis reductions follow [[Axis]]: `Axis.Rows` gives one
  // value per row, `Axis.Cols` one value per column.
  // ---------------------------------------------------------------------------

  /** Fast sum of all entries; `0.0` for an empty matrix. May reassociate, so
    * the final bits are build- and platform-specific; see [[sumExact]].
    */
  def sum: Double =
    if isContiguousBlock then DoubleKernels.dsum(rows * cols, data, offset.value, 1)
    else
      var acc = 0.0
      var line = 0
      while line < traversalLines do
        acc += DoubleKernels.dsum(traversalLength, data, offset.value + line * traversalLineStep, traversalElementStep)
        line += 1
      acc

  /** Sum of all entries rounded once to the nearest `Double` via
    * [[gale.numeric.ExactSum]]: independent of layout, order, and platform.
    */
  def sumExact: Double =
    DVec.exactSum(data, offset.value, traversalLength, traversalElementStep, traversalLines, traversalLineStep)

  /** Arithmetic mean of all entries, `sum / (rows * cols)`, with the same
    * overflow repair as [[DVec.mean]]: a non-finite sum is recomputed as
    * `sum(a_ij / (rows * cols))`, so entries near `Double.MaxValue` give a finite
    * mean while infinite and NaN entries keep their IEEE results. Throws
    * [[LinAlgError.EmptyInput]] when the matrix has no entries.
    */
  def mean: Double =
    requireNonEmpty("mean")
    val count = rows.toDouble * cols.toDouble
    val total = sum
    if total.isFinite then total / count
    else
      var acc = 0.0
      var hasInfinity = false
      var line = 0
      while line < traversalLines do
        val start = offset.value + line * traversalLineStep
        acc += DoubleKernels.dsumDivided(traversalLength, data, start, traversalElementStep, count)
        hasInfinity ||= DoubleKernels.dcontainsInfinity(traversalLength, data, start, traversalElementStep)
        line += 1
      DoubleKernels.clampMeanOverflow(acc, hasInfinity)

  /** Largest entry; NaN if any entry is NaN. Throws [[LinAlgError.EmptyInput]]
    * when the matrix has no entries.
    */
  def max: Double =
    requireNonEmpty("max")
    val (row, col) = extremePosition(largest = true)
    data(index(row, col))

  /** Smallest entry; NaN if any entry is NaN. Throws [[LinAlgError.EmptyInput]]
    * when the matrix has no entries.
    */
  def min: Double =
    requireNonEmpty("min")
    val (row, col) = extremePosition(largest = false)
    data(index(row, col))

  /** `(row, col)` of the first largest entry in row-major order, or of the
    * first NaN in row-major order if any entry is NaN — independent of the
    * storage layout. Throws [[LinAlgError.EmptyInput]] when empty.
    */
  def argmax: (Int, Int) =
    requireNonEmpty("argmax")
    extremePosition(largest = true)

  /** `(row, col)` of the first smallest entry in row-major order, or of the
    * first NaN if any. Throws [[LinAlgError.EmptyInput]] when empty.
    */
  def argmin: (Int, Int) =
    requireNonEmpty("argmin")
    extremePosition(largest = false)

  /** Per-axis sums: `Axis.Rows` gives the length-`rows` vector of row sums,
    * `Axis.Cols` the length-`cols` vector of column sums. Empty lines sum to
    * `0.0`.
    */
  def sum(axis: Axis): DVec =
    val lines = axisLines(axis)
    val length = axisLength(axis)
    val lineStep = axisLineStep(axis)
    val elementStep = axisElementStep(axis)
    val out = DVec.zeros(lines)
    val outData = out.data
    if streamsAxis(axis) then
      // The kept axis is unit-stride: accumulate whole slices into `out` so the
      // reduction streams through storage instead of striding down each line.
      var k = 0
      while k < length do
        DoubleKernels.dadd(lines, outData, 0, 1, data, offset.value + k * elementStep, 1, outData, 0, 1)
        k += 1
    else
      var line = 0
      while line < lines do
        outData(line) = DoubleKernels.dsum(length, data, offset.value + line * lineStep, elementStep)
        line += 1
    out

  /** Per-axis means: the per-axis sum divided by the line length. A line whose
    * sum is not finite is recomputed as in [[DVec.mean]], so it overflows only
    * when its mean does. Throws [[LinAlgError.EmptyInput]] when the result is
    * non-empty but each reduced line is empty (for example `Axis.Rows` on an
    * `n×0` matrix, `n > 0`).
    */
  def mean(axis: Axis): DVec =
    requireNonEmptyLines(axis, "mean")
    val out = sum(axis)
    val count = axisLength(axis)
    val length = count.toDouble
    val outData = out.data
    val lineStep = axisLineStep(axis)
    val elementStep = axisElementStep(axis)
    var line = 0
    while line < out.length do
      val total = outData(line)
      outData(line) =
        if total.isFinite then total / length
        else DoubleKernels.dmeanOfNonFiniteSum(count, data, offset.value + line * lineStep, elementStep)
      line += 1
    out

  /** Per-axis maxima, each with the semantics of [[DVec.max]] (NaN propagates).
    * Throws [[LinAlgError.EmptyInput]] when the reduced lines are empty.
    */
  def max(axis: Axis): DVec =
    requireNonEmptyLines(axis, "max")
    axisExtremes(axis, largest = true)

  /** Per-axis minima, each with the semantics of [[DVec.min]]. Throws
    * [[LinAlgError.EmptyInput]] when the reduced lines are empty.
    */
  def min(axis: Axis): DVec =
    requireNonEmptyLines(axis, "min")
    axisExtremes(axis, largest = false)

  /** Matrix 1-norm: the maximum absolute column sum. `0.0` when empty; NaN if
    * any entry is NaN.
    */
  def norm1: Double =
    maxAbsLineSum(Axis.Cols)

  /** Matrix ∞-norm: the maximum absolute row sum. `0.0` when empty; NaN if any
    * entry is NaN.
    */
  def normInf: Double =
    maxAbsLineSum(Axis.Rows)

  /** Frobenius norm `sqrt(sum a_ij^2)`, computed with a max-scaled second pass
    * when needed so it neither overflows (entries near `1e300`) nor loses tiny
    * entries to underflow. `0.0` when empty; NaN if any entry is NaN, otherwise
    * `+Inf` if any entry is infinite.
    */
  def normFrobenius: Double =
    DoubleKernels.dnrmFrobenius(rows, cols, data, offset.value, rowStride.value, colStride.value)

  private def isContiguousBlock: Boolean =
    isContiguousRowMajor || isContiguousColMajor

  // Whole-matrix traversal: lines along the smaller stride, except that a
  // single row or column is always one line (one kernel call).
  private def traverseRows: Boolean =
    if rows == 1 then true else if cols == 1 then false else colStride.value <= rowStride.value
  private[gale] def traversalLines: Int = if traverseRows then rows else cols
  private[gale] def traversalLength: Int = if traverseRows then cols else rows
  private[gale] def traversalLineStep: Int = if traverseRows then rowStride.value else colStride.value
  private[gale] def traversalElementStep: Int = if traverseRows then colStride.value else rowStride.value

  // Per-axis geometry: `axisLines(axis)` results, each reducing `axisLength(axis)` entries.
  private[gale] def axisLines(axis: Axis): Int = if axis == Axis.Rows then rows else cols
  private[gale] def axisLength(axis: Axis): Int = if axis == Axis.Rows then cols else rows
  private[gale] def axisLineStep(axis: Axis): Int =
    if axis == Axis.Rows then rowStride.value else colStride.value
  private[gale] def axisElementStep(axis: Axis): Int =
    if axis == Axis.Rows then colStride.value else rowStride.value

  /** True when a per-axis reduction should stream whole slices into per-line
    * accumulators: the kept axis is unit-stride while each line is strided (for
    * example `Axis.Cols` on a row-major matrix), so reading slice by slice is
    * sequential in memory. Requires at least one entry per line.
    */
  private[gale] def streamsAxis(axis: Axis): Boolean =
    axisLines(axis) > 1 && axisLength(axis) > 0 && axisLineStep(axis) == 1 && axisElementStep(axis) > 1

  /** Streamed per-line extremes for [[streamsAxis]] geometry, written to
    * `out(0 until lines)`. Equivalent to the per-line `dmaxIndex`/`dminIndex`
    * value: the first strict improvement wins ties, and a NaN, once seen, sticks.
    */
  private[gale] def streamedExtremes(axis: Axis, largest: Boolean, out: DoubleArray): Unit =
    if largest then streamExtremes(axis, out)(math.max)(DoubleKernels.dmax)
    else streamExtremes(axis, out)(math.min)(DoubleKernels.dmin)

  // The streamed pass is a branch-free `math.max`/`math.min` per slice, which C2
  // vectorizes. It can differ from the first-occurrence value only in which
  // signed zero or NaN payload it keeps, so lines ending at zero or NaN are
  // recomputed by the per-line kernel.
  private inline def streamExtremes(axis: Axis, out: DoubleArray)(inline pick: (Double, Double) => Double)(
      inline exact: (Int, DoubleArray, Int, Int) => Double
  ): Unit =
    val lines = axisLines(axis)
    val length = axisLength(axis)
    val elementStep = axisElementStep(axis)
    val base = offset.value
    var line = 0
    while line < lines do
      out(line) = data(base + line)
      line += 1
    var k = 1
    while k < length do
      val start = base + k * elementStep
      line = 0
      while line < lines do
        out(line) = pick(out(line), data(start + line))
        line += 1
      k += 1
    line = 0
    while line < lines do
      val value = out(line)
      if value == 0.0 || value != value then out(line) = exact(length, data, base + line, elementStep)
      line += 1

  private def extremePosition(largest: Boolean): (Int, Int) =
    if isContiguousRowMajor then
      val n = rows * cols
      val k =
        if largest then DoubleKernels.dmaxIndex(n, data, offset.value, 1)
        else DoubleKernels.dminIndex(n, data, offset.value, 1)
      (k / cols, k % cols)
    else if cols == 1 then
      val k =
        if largest then DoubleKernels.dmaxIndex(rows, data, offset.value, rowStride.value)
        else DoubleKernels.dminIndex(rows, data, offset.value, rowStride.value)
      (k, 0)
    else
      val colStep = colStride.value
      var bestRow = -1
      var bestCol = -1
      var best = 0.0
      var row = 0
      while row < rows do
        val start = offset.value + row * rowStride.value
        val col =
          if largest then DoubleKernels.dmaxIndex(cols, data, start, colStep)
          else DoubleKernels.dminIndex(cols, data, start, colStep)
        val value = data(start + col * colStep)
        if value.isNaN then return (row, col)
        if bestRow < 0 || (largest && value > best) || (!largest && value < best) then
          bestRow = row
          bestCol = col
          best = value
        row += 1
      (bestRow, bestCol)

  private def axisExtremes(axis: Axis, largest: Boolean): DVec =
    val lines = axisLines(axis)
    val length = axisLength(axis)
    val lineStep = axisLineStep(axis)
    val elementStep = axisElementStep(axis)
    val out = DVec.zeros(lines)
    if streamsAxis(axis) then streamedExtremes(axis, largest, out.data)
    else
      var line = 0
      while line < lines do
        val start = offset.value + line * lineStep
        out.data(line) =
          if largest then DoubleKernels.dmax(length, data, start, elementStep)
          else DoubleKernels.dmin(length, data, start, elementStep)
        line += 1
    out

  /** Largest absolute line sum, the lines being given by `axis`. */
  private def maxAbsLineSum(axis: Axis): Double =
    val lines = axisLines(axis)
    val length = axisLength(axis)
    val lineStep = axisLineStep(axis)
    val elementStep = axisElementStep(axis)
    var out = 0.0
    if streamsAxis(axis) then
      // Accumulate |slice| into per-line sums so storage is read sequentially.
      val sums = DoubleArray.alloc(lines)
      var k = 0
      while k < length do
        val start = offset.value + k * elementStep
        var line = 0
        while line < lines do
          sums(line) = sums(line) + math.abs(data(start + line))
          line += 1
        k += 1
      var line = 0
      while line < lines do
        out = math.max(out, sums(line))
        line += 1
    else
      var line = 0
      while line < lines do
        out = math.max(out, DoubleKernels.dasum(length, data, offset.value + line * lineStep, elementStep))
        line += 1
    out

  private def requireNonEmpty(operation: String): Unit =
    if rows == 0 || cols == 0 then throw LinAlgError.EmptyInput(s"DMat.$operation")

  private[gale] def requireNonEmptyLines(axis: Axis, operation: String): Unit =
    if axisLines(axis) > 0 && axisLength(axis) == 0 then
      throw LinAlgError.EmptyInput(s"DMat.$operation(Axis.$axis)")

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
  private[gale] def requireStorable(rows: Int, cols: Int): Unit =
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

  /** Allocate a single-owner row-major builder whose [[DMatBuilder.result]]
    * transfers storage into an immutable matrix without copying.
    */
  def newBuilder(rows: Int, cols: Int): DMatBuilder =
    DMatBuilder.zeros(rows, cols)

  /** Create a builder containing one owned logical row-major copy of `matrix`. */
  def builderFrom(matrix: DMat): DMatBuilder =
    DMatBuilder.from(matrix)

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
  /** Construct an owned dense matrix from row-major literal values.
    *
    * This is the concise spelling of [[dense]] for handwritten matrices. The
    * number of values must equal `rows * cols`.
    */
  def apply(rows: Int, cols: Int)(values: Double*): DMat =
    DMat.dense(rows, cols, values)

  def zeros(rows: Int, cols: Int): DMat =
    DMat.zeros(rows, cols)

  def newBuilder(rows: Int, cols: Int): DMatBuilder =
    DMat.newBuilder(rows, cols)

  def builderFrom(matrix: DMat): DMatBuilder =
    DMat.builderFrom(matrix)

  def eye(size: Int): DMat =
    DMat.eye(size)

  def dense(rows: Int, cols: Int, values: Seq[Double]): DMat =
    DMat.dense(rows, cols, values)

  @targetName("denseVarargs")
  def dense(rows: Int, cols: Int)(values: Double*): DMat =
    DMat.dense(rows, cols, values)

  def tabulate(rows: Int, cols: Int)(f: (Int, Int) => Double): DMat =
    DMat.tabulate(rows, cols)(f)
