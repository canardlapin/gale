package gale.linalg

import gale.TestAccess
import gale.backend.*
import gale.platform.DoubleArray

class DenseDestinationSuite extends munit.FunSuite:
  private val a = Matrix.dense(2, 3)(
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0
  )
  private val b = Matrix.dense(3, 2)(
    2.0, -1.0,
    0.5, 3.0,
    -2.0, 4.0
  )

  test("gemmInto accepts row-major, transposed, and strided inputs") {
    val aTransposed = Matrix.dense(3, 2)(1.0, 4.0, 2.0, 5.0, 3.0, 6.0).t
    val aStrided = TestAccess.mat(
      TestAccess.doubleArray(1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0, 0.0, 6.0, 0.0),
      offset = 0,
      rows = 2,
      cols = 3,
      rowStride = 6,
      colStride = 2
    )
    val bTransposed = Matrix.dense(2, 3)(2.0, 0.5, -2.0, -1.0, 3.0, 4.0).t
    val bStrided = TestAccess.mat(
      TestAccess.doubleArray(2.0, 0.0, -1.0, 0.0, 0.5, 0.0, 3.0, 0.0, -2.0, 0.0, 4.0, 0.0),
      offset = 0,
      rows = 3,
      cols = 2,
      rowStride = 4,
      colStride = 2
    )
    val expected = a * b

    for
      (aLabel, left) <- Seq("row" -> a, "transpose" -> aTransposed, "strided" -> aStrided)
      (bLabel, right) <- Seq("row" -> b, "transpose" -> bTransposed, "strided" -> bStrided)
    do
      val destination = DMatBuilder.zeros(2, 2)
      left.gemmInto(right, destination)
      assertMatrixClose(destination.result(), expected, 1e-12, s"$aLabel/$bLabel")
  }

  test("gemmInto makes replacement and accumulation semantics explicit") {
    val identity = Matrix.eye(2)
    val source = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
    val replacement = DMatBuilder.zeros(2, 2)
    replacement.fill(Double.NaN)
    source.gemmInto(identity, replacement, alpha = 1.0, beta = 0.0)
    assertEquals(replacement.result().valuesRowMajor, source.valuesRowMajor)

    val accumulated = DMatBuilder.zeros(2, 2)
    accumulated.fill(1.0)
    source.gemmInto(identity, accumulated, alpha = 2.0, beta = 3.0)
    assertEquals(accumulated.result().valuesRowMajor, Seq(5.0, 7.0, 9.0, 11.0))

    val emptyProduct = DMatBuilder.zeros(2, 2)
    emptyProduct.fill(Double.NaN)
    Matrix.zeros(2, 0).gemmInto(Matrix.zeros(0, 2), emptyProduct, beta = 0.0)
    assertEquals(emptyProduct.result().valuesRowMajor, Seq.fill(4)(0.0))
  }

  test("linearCombinationInto fuses two strided inputs and ignores old destination") {
    val left = Matrix.dense(3, 2)(1.0, 4.0, 2.0, 5.0, 3.0, 6.0).t
    val right = TestAccess.mat(
      TestAccess.doubleArray(6.0, 0.0, 5.0, 0.0, 4.0, 0.0, 3.0, 0.0, 2.0, 0.0, 1.0, 0.0),
      offset = 0,
      rows = 2,
      cols = 3,
      rowStride = 6,
      colStride = 2
    )
    val destination = DMatBuilder.zeros(2, 3)
    destination.fill(Double.NaN)
    left.linearCombinationInto(right, destination, alpha = 0.5, beta = -2.0)

    val expected = Matrix.tabulate(2, 3)((row, col) => 0.5 * left(row, col) - 2.0 * right(row, col))
    assertMatrixClose(destination.result(), expected, 0.0, "linear combination")
  }

  test("zero linear-combination coefficients suppress reads of NaN inputs") {
    val finite = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
    val nan = Matrix.tabulate(2, 2)((_, _) => Double.NaN)

    val ignoreLeft = DMatBuilder.zeros(2, 2)
    nan.linearCombinationInto(finite, ignoreLeft, alpha = 0.0, beta = 2.0)
    assertEquals(ignoreLeft.result().valuesRowMajor, Seq(2.0, 4.0, 6.0, 8.0))

    val ignoreRight = DMatBuilder.zeros(2, 2)
    finite.linearCombinationInto(nan, ignoreRight, alpha = -1.0, beta = 0.0)
    assertEquals(ignoreRight.result().valuesRowMajor, Seq(-1.0, -2.0, -3.0, -4.0))

    val ignoreBoth = DMatBuilder.zeros(2, 2)
    nan.linearCombinationInto(nan, ignoreBoth, alpha = 0.0, beta = 0.0)
    assertEquals(ignoreBoth.result().valuesRowMajor, Seq.fill(4)(0.0))
  }

  test("destination operations reject shape mismatch and closed builders") {
    intercept[LinAlgError.DimensionMismatch] {
      a.gemmInto(Matrix.zeros(4, 2), DMatBuilder.zeros(2, 2))
    }
    intercept[LinAlgError.DimensionMismatch] {
      a.gemmInto(b, DMatBuilder.zeros(2, 3))
    }
    intercept[LinAlgError.DimensionMismatch] {
      a.linearCombinationInto(Matrix.zeros(2, 2), DMatBuilder.zeros(2, 3), 1.0, 1.0)
    }

    val closed = DMatBuilder.zeros(2, 2)
    closed.result()
    intercept[LinAlgError.UnsupportedOperation] {
      a.gemmInto(b, closed)
    }
  }

  test("destination operations reject internal aliases deterministically") {
    val destination = DMatBuilder.zeros(2, 2)
    val alias = TestAccess.mat(
      TestAccess.dmatBuilderStorage(destination),
      offset = 0,
      rows = 2,
      cols = 2,
      rowStride = 2,
      colStride = 1
    )
    val identity = Matrix.eye(2)

    intercept[LinAlgError.UnsupportedOperation](alias.gemmInto(identity, destination))
    intercept[LinAlgError.UnsupportedOperation](identity.gemmInto(alias, destination))
    intercept[LinAlgError.UnsupportedOperation](alias.linearCombinationInto(identity, destination, 1.0, 1.0))
  }

  test("gemmInto preserves backend gemm and syrk routing") {
    for capability <- Seq(Capability.Vectorized, Capability.NativeBlas) do
      RecordingKernel.reset()
      val destination = DMatBuilder.zeros(2, 2)
      a.gemmInto(b, destination)(using RecordingBackend(capability))
      assertEquals(RecordingKernel.gemmCalls, 1)
      assertMatrixClose(destination.result(), a.*(b)(using PureBackend), 1e-12, capability.toString)

    RecordingKernel.reset()
    val gramInput = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val gram = DMatBuilder.zeros(2, 2)
    gramInput.t.gemmInto(gramInput, gram)(using RecordingBackend(Capability.Vectorized))
    assertEquals(RecordingKernel.syrkCalls, 1)
    assertMatrixClose(gram.result(), gramInput.t * gramInput, 1e-12, "syrk")
  }

  private object RecordingKernel extends DenseDoubleKernel:
    export PureDenseDoubleKernel.{gemm => _, syrk => _, *}
    var gemmCalls = 0
    var syrkCalls = 0

    def reset(): Unit =
      gemmCalls = 0
      syrkCalls = 0

    def gemm(
        rows: Int,
        cols: Int,
        shared: Int,
        alpha: Double,
        a: DoubleArray,
        aOffset: Int,
        aRowStride: Int,
        aColStride: Int,
        b: DoubleArray,
        bOffset: Int,
        bRowStride: Int,
        bColStride: Int,
        beta: Double,
        c: DoubleArray,
        cOffset: Int,
        cRowStride: Int,
        cColStride: Int
    ): Unit =
      gemmCalls += 1
      PureDenseDoubleKernel.gemm(
        rows,
        cols,
        shared,
        alpha,
        a,
        aOffset,
        aRowStride,
        aColStride,
        b,
        bOffset,
        bRowStride,
        bColStride,
        beta,
        c,
        cOffset,
        cRowStride,
        cColStride
      )

    def syrk(
        m: Int,
        k: Int,
        a: DoubleArray,
        aOffset: Int,
        aRowStride: Int,
        c: DoubleArray,
        cOffset: Int,
        cRowStride: Int
    ): Unit =
      syrkCalls += 1
      PureDenseDoubleKernel.syrk(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)

  private final case class RecordingBackend(capability: Capability) extends Backend:
    val name: String = s"recording-$capability"
    val capabilities: Set[Capability] = Set(capability)
    val denseDouble: DenseDoubleKernel = RecordingKernel
    val thresholds: BackendThresholds = new BackendThresholds:
      def nativeGemmMinFlops: Long = 1L
      def nativeGemvMinWork: Long = 1L
      def nativeFactorizationMinSize: Int = Int.MaxValue
    val config: BackendConfig = BackendConfig.singleThreaded

  private def assertMatrixClose(actual: DMat, expected: DMat, tolerance: Double, clue: String): Unit =
    assertEquals(actual.rows, expected.rows, clue)
    assertEquals(actual.cols, expected.cols, clue)
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assert(math.abs(actual(row, col) - expected(row, col)) <= tolerance, s"$clue ($row,$col)")
        col += 1
      row += 1
