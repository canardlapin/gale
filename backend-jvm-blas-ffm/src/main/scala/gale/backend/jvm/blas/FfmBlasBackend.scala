package gale.backend.jvm.blas

import gale.backend.*
import gale.backend.jvm.`native`.{Layout, NativeDMat}
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

import java.lang.foreign.{Arena, MemorySegment}
import scala.util.control.NonFatal

final case class BlasLoadError(message: String, cause: Throwable)
    extends RuntimeException(message, cause)

/** Conservative heap-copy-inclusive defaults pending platform JMH sweeps. */
final case class FfmBlasThresholds(
    nativeGemmMinFlops: Long = Long.MaxValue,
    nativeGemvMinWork: Long = Long.MaxValue
) extends BackendThresholds:
  val nativeFactorizationMinSize: Int = Int.MaxValue
  require(nativeGemmMinFlops >= 0L)
  require(nativeGemvMinWork >= 0L)

object FfmBlasThresholds:
  /** Only known optimized families dispatch by default. An unknown/reference
    * `libblas` remains direct-callable but cannot regress Gale's public route.
    */
  def forLibrary(name: String): FfmBlasThresholds =
    val normalized = name.toLowerCase(java.util.Locale.ROOT)
    val optimized =
      normalized.contains("accelerate") || normalized.contains("openblas") || normalized.contains("mkl")
    if optimized then FfmBlasThresholds(nativeGemmMinFlops = 256L * 256L * 256L)
    else FfmBlasThresholds()

final case class BlasLibraryInfo(
    name: String,
    hasThreadControl: Boolean,
    threadSetter: Option[String]
)

final class FfmBlasBackend private[blas] (
    private val bindings: CblasBindings,
    val config: BackendConfig,
    val thresholds: BackendThresholds
) extends Backend, AutoCloseable:
  val libraryInfo: BlasLibraryInfo = BlasLibraryInfo(
    bindings.libraryName(),
    bindings.hasThreadControl(),
    Option(bindings.threadSetterName())
  )
  val name: String = s"jvm-blas-ffm:${libraryInfo.name}"
  val capabilities: Set[Capability] =
    Set(Capability.NativeBlas) ++
      (if config.nativeThreads > 1 then Set(Capability.Multithreaded) else Set.empty)
  val denseDouble: DenseDoubleKernel = new FfmDenseDoubleKernel(bindings)

  /** Copy-free explicit native GEMM. The output layout selects CBLAS row/column
    * major; operands in the other layout are represented with a transpose flag,
    * so mixed layouts need no repacking.
    */
  def gemm(
      a: NativeDMat,
      b: NativeDMat,
      c: NativeDMat,
      alpha: Double = 1.0,
      beta: Double = 0.0
  ): Unit =
    require(a.cols == b.rows, s"native GEMM inner mismatch: ${a.cols} != ${b.rows}")
    require(c.rows == a.rows && c.cols == b.cols,
      s"native GEMM output mismatch: ${c.rows} x ${c.cols} != ${a.rows} x ${b.cols}")
    val cblasLayout = if c.layout == Layout.RowMajor then CblasBindings.ROW_MAJOR else CblasBindings.COL_MAJOR
    def operand(matrix: NativeDMat): (Int, Int) =
      if matrix.layout == c.layout then (CblasBindings.NO_TRANS, matrix.leadingDimension)
      else (CblasBindings.TRANS, matrix.leadingDimension)
    val (transA, lda) = operand(a)
    val (transB, ldb) = operand(b)
    bindings.dgemm(
      cblasLayout, transA, transB,
      a.rows, b.cols, a.cols, alpha,
      a.memory, lda, b.memory, ldb, beta, c.memory, c.leadingDimension
    )

  override def close(): Unit = bindings.close()

object FfmBlasBackend:
  def load(
      config: BackendConfig = BackendConfig.singleThreaded,
      thresholds: Option[BackendThresholds] = None
  ): Either[BlasLoadError, FfmBlasBackend] =
    try
      val bindings = CblasBindings.loadDefault()
      val threadCountConfigured =
        if bindings.hasThreadControl() then bindings.configureThreads(config.nativeThreads)
        else false
      if config.nativeThreads != 1 && !threadCountConfigured then
        bindings.close()
        Left(BlasLoadError(
          s"${bindings.libraryName()} exposes no supported thread-count setter; " +
            "use nativeThreads=1 or configure the vendor before JVM startup",
          IllegalStateException("native thread control unavailable")
        ))
      else
        val selectedThresholds = thresholds.getOrElse(FfmBlasThresholds.forLibrary(bindings.libraryName()))
        Right(Backend.requireValid(new FfmBlasBackend(bindings, config, selectedThresholds)))
    catch
      case NonFatal(error) => Left(BlasLoadError(error.getMessage, error))

  lazy val default: FfmBlasBackend =
    load().fold(throw _, identity)

private final class FfmDenseDoubleKernel(bindings: CblasBindings) extends DenseDoubleKernel:
  export PureDenseDoubleKernel.{gemm => _, gemv => _, syrk => _, *}

  private final case class MatrixLayout(transpose: Int, leadingDimension: Int)

  def gemm(
      rows: Int, cols: Int, shared: Int, alpha: Double,
      a: DoubleArray, aOffset: Int, aRowStride: Int, aColStride: Int,
      b: DoubleArray, bOffset: Int, bRowStride: Int, bColStride: Int,
      beta: Double, c: DoubleArray, cOffset: Int, cRowStride: Int, cColStride: Int
  ): Unit =
    val aLayout = matrixLayout(rows, shared, aRowStride, aColStride)
    val bLayout = matrixLayout(shared, cols, bRowStride, bColStride)
    val work = rows.toLong * cols.toLong * shared.toLong
    val copiedElements = a.length.toLong + b.length.toLong + c.length.toLong
    // A scalar product can exceed the facade's cubic threshold while doing only
    // O(n) arithmetic for O(n) copying. Require enough arithmetic intensity for
    // the heap/native boundary; otherwise the pure kernel is structurally better.
    val copyAmortized = work >= 8L * copiedElements
    if rows == 0 || cols == 0 || shared == 0 || !copyAmortized || aLayout.isEmpty || bLayout.isEmpty || cColStride != 1 then
      PureDenseDoubleKernel.gemm(
        rows, cols, shared, alpha,
        a, aOffset, aRowStride, aColStride,
        b, bOffset, bRowStride, bColStride,
        beta, c, cOffset, cRowStride, cColStride
      )
    else
      withNative3(a, b, c) { (aSegment, bSegment, cSegment) =>
        bindings.dgemm(
          CblasBindings.ROW_MAJOR, aLayout.get.transpose, bLayout.get.transpose,
          rows, cols, shared, alpha,
          atOffset(aSegment, aOffset), aLayout.get.leadingDimension,
          atOffset(bSegment, bOffset), bLayout.get.leadingDimension,
          beta, atOffset(cSegment, cOffset), cRowStride
        )
      }

  def gemv(
      rows: Int, cols: Int, alpha: Double,
      a: DoubleArray, aOffset: Int, rowStride: Int, colStride: Int,
      x: DoubleArray, xOffset: Int, xStride: Int,
      beta: Double, y: DoubleArray, yOffset: Int, yStride: Int
  ): Unit =
    val layout = matrixLayout(rows, cols, rowStride, colStride)
    if rows == 0 || cols == 0 || layout.isEmpty then
      PureDenseDoubleKernel.gemv(
        rows, cols, alpha, a, aOffset, rowStride, colStride,
        x, xOffset, xStride, beta, y, yOffset, yStride
      )
    else
      withNative3(a, x, y) { (aSegment, xSegment, ySegment) =>
        val transposed = layout.get.transpose == CblasBindings.TRANS
        val physicalRows = if transposed then cols else rows
        val physicalCols = if transposed then rows else cols
        bindings.dgemv(
          layout.get.transpose, physicalRows, physicalCols, alpha,
          atOffset(aSegment, aOffset), layout.get.leadingDimension,
          atOffset(xSegment, xOffset), xStride,
          beta, atOffset(ySegment, yOffset), yStride
        )
      }

  def syrk(
      m: Int, k: Int, a: DoubleArray, aOffset: Int, aRowStride: Int,
      c: DoubleArray, cOffset: Int, cRowStride: Int
  ): Unit =
    if m == 0 || k == 0 then
      PureDenseDoubleKernel.syrk(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)
    else
      withNative2(a, c) { (aSegment, cSegment) =>
        bindings.dsyrk(
          k, m, 1.0, atOffset(aSegment, aOffset), aRowStride,
          0.0, atOffset(cSegment, cOffset), cRowStride
        )
      }
      var i = 0
      while i < k do
        var j = i + 1
        while j < k do
          c(cOffset + j * cRowStride + i) = c(cOffset + i * cRowStride + j)
          j += 1
        i += 1

  private def matrixLayout(rows: Int, cols: Int, rowStride: Int, colStride: Int): Option[MatrixLayout] =
    if colStride == 1 && rowStride >= math.max(1, cols) then
      Some(MatrixLayout(CblasBindings.NO_TRANS, rowStride))
    else if rowStride == 1 && colStride >= math.max(1, rows) then
      Some(MatrixLayout(CblasBindings.TRANS, colStride))
    else None

  private def withNative2(a: DoubleArray, b: DoubleArray)(operation: (MemorySegment, MemorySegment) => Unit): Unit =
    val arena = Arena.ofConfined()
    val aArray = DoubleArray.asArray(a)
    val bArray = DoubleArray.asArray(b)
    try
      val aSegment = copyIn(aArray, arena)
      val bSegment = copyIn(bArray, arena)
      operation(aSegment, bSegment)
      copyOut(bSegment, bArray)
    finally arena.close()

  private def withNative3(a: DoubleArray, b: DoubleArray, c: DoubleArray)(
      operation: (MemorySegment, MemorySegment, MemorySegment) => Unit
  ): Unit =
    val arena = Arena.ofConfined()
    val aArray = DoubleArray.asArray(a)
    val bArray = DoubleArray.asArray(b)
    val cArray = DoubleArray.asArray(c)
    try
      val aSegment = copyIn(aArray, arena)
      val bSegment = copyIn(bArray, arena)
      val cSegment = copyIn(cArray, arena)
      operation(aSegment, bSegment, cSegment)
      copyOut(cSegment, cArray)
    finally arena.close()

  private def copyIn(array: Array[Double], arena: Arena): MemorySegment =
    val bytes = array.length.toLong * java.lang.Double.BYTES
    val segment = arena.allocate(math.max(1L, bytes), java.lang.Double.BYTES)
    if bytes > 0 then segment.asSlice(0L, bytes).copyFrom(MemorySegment.ofArray(array))
    segment

  private def copyOut(segment: MemorySegment, array: Array[Double]): Unit =
    val bytes = array.length.toLong * java.lang.Double.BYTES
    if bytes > 0 then MemorySegment.ofArray(array).copyFrom(segment.asSlice(0L, bytes))

  private inline def atOffset(segment: MemorySegment, offset: Int): MemorySegment =
    segment.asSlice(offset.toLong * java.lang.Double.BYTES)

/** Explicit import point. Evaluating it loads the first conforming CBLAS candidate. */
given ffmBlasBackend: Backend = FfmBlasBackend.default
