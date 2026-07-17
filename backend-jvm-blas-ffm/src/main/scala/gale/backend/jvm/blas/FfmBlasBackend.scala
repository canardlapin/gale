package gale.backend.jvm.blas

import gale.backend.*
import gale.backend.jvm.`native`.{Layout, NativeDMat}
import gale.linalg.{Cholesky, DMat, DVec, DenseDecompositions, FactorizationDiagnostics, LU, LinAlgError, Matrix, PivotVector, QR}
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.spectral.{RawSymmetricEigen, SpectralBackend, SpectralCapability}

import java.lang.foreign.{Arena, MemorySegment, ValueLayout}
import scala.util.control.NonFatal

final case class BlasLoadError(message: String, cause: Throwable)
    extends RuntimeException(message, cause)

/** Conservative heap-copy-inclusive defaults pending platform JMH sweeps. */
final case class FfmBlasThresholds(
    nativeGemmMinFlops: Long = Long.MaxValue,
    nativeGemvMinWork: Long = Long.MaxValue,
    override val nativeLuMinSize: Int = Int.MaxValue,
    override val nativeCholeskyMinSize: Int = Int.MaxValue,
    override val nativeQrMinSize: Int = Int.MaxValue
) extends BackendThresholds:
  val nativeFactorizationMinSize: Int = math.min(nativeLuMinSize, math.min(nativeCholeskyMinSize, nativeQrMinSize))
  require(nativeGemmMinFlops >= 0L)
  require(nativeGemvMinWork >= 0L)
  require(nativeLuMinSize >= 0)
  require(nativeCholeskyMinSize >= 0)
  require(nativeQrMinSize >= 0)

object FfmBlasThresholds:
  /** Only known optimized families dispatch by default. An unknown/reference
    * `libblas` remains direct-callable but cannot regress Gale's public route.
    */
  def forLibrary(name: String): FfmBlasThresholds =
    val normalized = name.toLowerCase(java.util.Locale.ROOT)
    if normalized.contains("accelerate") then
      // JDK 22 / Apple ARM64, two-fork copy-inclusive JMH (2026-07-17):
      // LU wins from n=64 onward, so n=128 leaves one measured size of margin.
      // Accelerate QR and Cholesky are non-monotone through this heap/FFM route;
      // keep both disabled unless a caller explicitly supplies thresholds.
      FfmBlasThresholds(
        nativeGemmMinFlops = 256L * 256L * 256L,
        nativeLuMinSize = 128
      )
    else if normalized.contains("openblas") || normalized.contains("mkl") then
      // GEMM has portable evidence; factorization defaults require a sweep on
      // the actual library family because copy cost and vendor kernels differ.
      FfmBlasThresholds(nativeGemmMinFlops = 256L * 256L * 256L)
    else FfmBlasThresholds()

final case class BlasLibraryInfo(
    name: String,
    hasLapack: Boolean,
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
    bindings.hasLapack(),
    bindings.hasThreadControl(),
    Option(bindings.threadSetterName())
  )
  val name: String = s"jvm-blas-ffm:${libraryInfo.name}"
  val capabilities: Set[Capability] =
    Set(Capability.NativeBlas) ++
      (if bindings.hasLapack() then Set(Capability.NativeLapack) else Set.empty) ++
      (if config.nativeThreads > 1 then Set(Capability.Multithreaded) else Set.empty)
  val denseDouble: DenseDoubleKernel = new FfmDenseDoubleKernel(bindings)
  override val denseFactorizations: Option[DenseDoubleFactorizations] =
    Option.when(bindings.hasLapack())(new FfmDenseDoubleFactorizations(bindings))
  override val spectral: Option[SpectralBackend] =
    Option.when(bindings.hasLapack())(new FfmLapackSpectralBackend(bindings, name))

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

private object FfmLapackMemory:
  def nativeDoubles(arena: Arena, length: Int): MemorySegment =
    arena.allocate(math.max(1L, length.toLong * java.lang.Double.BYTES), java.lang.Double.BYTES)

  def nativeInts(arena: Arena, length: Int): MemorySegment =
    arena.allocate(math.max(1L, length.toLong * java.lang.Integer.BYTES), java.lang.Integer.BYTES)

  def copyIn(values: Array[Double], arena: Arena): MemorySegment =
    val target = nativeDoubles(arena, values.length)
    if values.nonEmpty then target.asSlice(0L, values.length.toLong * java.lang.Double.BYTES).copyFrom(MemorySegment.ofArray(values))
    target

  def copyOut(source: MemorySegment, values: Array[Double]): Unit =
    if values.nonEmpty then
      MemorySegment.ofArray(values).copyFrom(source.asSlice(0L, values.length.toLong * java.lang.Double.BYTES))

  def toColumnMajor(a: DMat): Array[Double] =
    val out = new Array[Double](a.rows * a.cols)
    var j = 0
    while j < a.cols do
      var i = 0
      while i < a.rows do
        out(j * a.rows + i) = a(i, j)
        i += 1
      j += 1
    out

  def fromColumnMajor(rows: Int, cols: Int, values: Array[Double]): DMat =
    Matrix.tabulate(rows, cols)((i, j) => values(j * rows + i))

private final class FfmDenseDoubleFactorizations(bindings: CblasBindings) extends DenseDoubleFactorizations:
  import FfmLapackMemory.*

  def lu(a: DMat): Either[LinAlgError, LU] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else if a.rows == 0 then DenseDecompositions.lu(a)
    else
      val n = a.rows
      val packedColumnMajor = toColumnMajor(a)
      val arena = Arena.ofConfined()
      try
        val matrix = copyIn(packedColumnMajor, arena)
        val ipiv = nativeInts(arena, n)
        val info = bindings.dgetrf(n, n, matrix, n, ipiv)
        if info < 0 then Left(LinAlgError.InvalidArgument(s"native dgetrf rejected argument ${-info}"))
        else if info > 0 then Left(LinAlgError.SingularMatrix(info - 1))
        else
          copyOut(matrix, packedColumnMajor)
          val permutation = Array.tabulate(n)(i => i)
          var parity = 1
          var k = 0
          while k < n do
            val pivot = ipiv.getAtIndex(ValueLayout.JAVA_INT, k.toLong) - 1
            if pivot != k then
              val tmp = permutation(k)
              permutation(k) = permutation(pivot)
              permutation(pivot) = tmp
              parity = -parity
            k += 1
          Right(LU(
            packed = fromColumnMajor(n, n, packedColumnMajor),
            pivots = PivotVector.fromArray(permutation),
            parity = parity,
            diagnostics = FactorizationDiagnostics()
          ))
      finally arena.close()

  def cholesky(a: DMat): Either[LinAlgError, Cholesky] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else if a.rows == 0 then DenseDecompositions.cholesky(a)
    else
      val n = a.rows
      val packedColumnMajor = toColumnMajor(a)
      val arena = Arena.ofConfined()
      try
        val matrix = copyIn(packedColumnMajor, arena)
        val info = bindings.dpotrf('L'.toByte, n, matrix, n)
        if info < 0 then Left(LinAlgError.InvalidArgument(s"native dpotrf rejected argument ${-info}"))
        else if info > 0 then Left(LinAlgError.NotPositiveDefinite(info - 1))
        else
          copyOut(matrix, packedColumnMajor)
          val lower = Matrix.tabulate(n, n)((i, j) => if i >= j then packedColumnMajor(j * n + i) else 0.0)
          Right(Cholesky(lower, FactorizationDiagnostics()))
      finally arena.close()

  def qr(a: DMat): Either[LinAlgError, QR] =
    if a.rows == 0 || a.cols == 0 then Right(DenseDecompositions.qr(a)(using PureBackend))
    else
      val m = a.rows
      val n = a.cols
      val limit = math.min(m, n)
      val packedColumnMajor = toColumnMajor(a)
      val arena = Arena.ofConfined()
      try
        val matrix = copyIn(packedColumnMajor, arena)
        val tauSegment = nativeDoubles(arena, limit)
        val query = nativeDoubles(arena, 1)
        val queryInfo = bindings.dgeqrf(m, n, matrix, m, tauSegment, query, -1)
        if queryInfo != 0 then Left(LinAlgError.InvalidArgument(s"native dgeqrf workspace query failed with info=$queryInfo"))
        else
          val lwork = math.max(1, math.ceil(query.get(ValueLayout.JAVA_DOUBLE, 0L)).toInt)
          val work = nativeDoubles(arena, lwork)
          val info = bindings.dgeqrf(m, n, matrix, m, tauSegment, work, lwork)
          if info != 0 then Left(LinAlgError.InvalidArgument(s"native dgeqrf failed with info=$info"))
          else
            copyOut(matrix, packedColumnMajor)
            val tauArray = new Array[Double](limit)
            copyOut(tauSegment, tauArray)
            val reflectors = Matrix.tabulate(m, limit)((i, j) =>
              if i < j then 0.0 else if i == j then 1.0 else packedColumnMajor(j * m + i)
            )
            val r = Matrix.tabulate(m, n)((i, j) => if i <= j then packedColumnMajor(j * m + i) else 0.0)
            Right(QR(
              reflectors,
              DoubleArray.adopt(tauArray),
              r,
              FactorizationDiagnostics(rank = Some(DenseDecompositions.rankFromMatrix(r)))
            ))
      finally arena.close()

private final class FfmLapackSpectralBackend(bindings: CblasBindings, backendName: String) extends SpectralBackend:
  import FfmLapackMemory.*

  val name: String = s"$backendName:lapack"
  val capabilities: Set[SpectralCapability] = Set(SpectralCapability.DenseSymmetricEigen)

  override def denseSymmetricEigen(a: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
    if a.rows != a.cols then Left(LinAlgError.NonSquareMatrix(a.shape))
    else if a.rows == 0 then Right(RawSymmetricEigen(DVec.zeros(0), DMat.zeros(0, 0)))
    else
      val n = a.rows
      val packedColumnMajor = toColumnMajor(a)
      val arena = Arena.ofConfined()
      try
        val matrix = copyIn(packedColumnMajor, arena)
        val eigenvalues = nativeDoubles(arena, n)
        val query = nativeDoubles(arena, 1)
        val jobz = if wantVectors then 'V'.toByte else 'N'.toByte
        val queryInfo = bindings.dsyev(jobz, 'L'.toByte, n, matrix, n, eigenvalues, query, -1)
        if queryInfo != 0 then Left(LinAlgError.InvalidArgument(s"native dsyev workspace query failed with info=$queryInfo"))
        else
          val lwork = math.max(1, math.ceil(query.get(ValueLayout.JAVA_DOUBLE, 0L)).toInt)
          val work = nativeDoubles(arena, lwork)
          val info = bindings.dsyev(jobz, 'L'.toByte, n, matrix, n, eigenvalues, work, lwork)
          if info < 0 then Left(LinAlgError.InvalidArgument(s"native dsyev rejected argument ${-info}"))
          else if info > 0 then Left(LinAlgError.DidNotConverge(info, Double.NaN))
          else
            val values = new Array[Double](n)
            copyOut(eigenvalues, values)
            if wantVectors then copyOut(matrix, packedColumnMajor)
            val vectors = if wantVectors then fromColumnMajor(n, n, packedColumnMajor) else DMat.zeros(n, 0)
            Right(RawSymmetricEigen(DVec.fromArray(values), vectors))
      finally arena.close()

/** Explicit import point. Evaluating it loads the first conforming BLAS/LAPACK candidate. */
given ffmBlasBackend: Backend = FfmBlasBackend.default
