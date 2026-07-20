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
    if cols != that.rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(that.cols)), Shape(Rows(that.rows), Cols(that.cols)))
    val out = DMat.zeros(rows, that.cols)
    // `Aᵀ · A` with a row-major `A` (the common Gram/normal-equations product):
    // a general gemm would traverse the transposed left operand column-strided and
    // cache-hostile, so route it to the dedicated symmetric rank-k kernel (half the
    // flops, unit-stride throughout). Detected structurally: `this` is exactly the
    // transpose view of `that`, and `that` has unit column stride.
    if isTransposeOf(that) && that.colStride.value == 1 then
      if backend.routesGemm(that.rows, that.cols, that.cols) then
        backend.denseDouble.syrk(
          that.rows,
          that.cols,
          that.data,
          that.offset.value,
          that.rowStride.value,
          out.data,
          out.offset.value,
          out.rowStride.value
        )
      else
        DoubleKernels.dsyrkRowMajor(
          that.rows,
          that.cols,
          that.data,
          that.offset.value,
          that.rowStride.value,
          out.data,
          out.offset.value,
          out.rowStride.value
        )
    else if backend.routesGemm(rows, that.cols, cols) then
      backend.denseDouble.gemm(
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
    else
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

  def solve(b: DVec)(using Backend): Either[LinAlgError, DVec] =
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
    * exactly when [[svd]] is.
    */
  def pinv(using SpectralBackend): Either[LinAlgError, DMat] =
    Svds.pinv(this)

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

  def newBuilder(rows: Int, cols: Int): DMatBuilder =
    DMat.newBuilder(rows, cols)

  def eye(size: Int): DMat =
    DMat.eye(size)

  def dense(rows: Int, cols: Int, values: Seq[Double]): DMat =
    DMat.dense(rows, cols, values)

  @targetName("denseVarargs")
  def dense(rows: Int, cols: Int)(values: Double*): DMat =
    DMat.dense(rows, cols, values)

  def tabulate(rows: Int, cols: Int)(f: (Int, Int) => Double): DMat =
    DMat.tabulate(rows, cols)(f)
