package gale.backend

import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.spectral.SpectralBackend
import gale.spectral.SpectralCapability

/** The coarse gemm dispatch seam on `DMat.*` (doc §A.2). With no acceleration import the
  * `given` resolves to `Backend.pure` and the pure kernel runs — the entire existing test
  * suite is the byte-identical witness for that; these tests add the ROUTING behaviour: an
  * imported accelerating backend takes over the general product above its threshold, a
  * below-threshold product stays pure, and every documented coarse path is observable.
  */
class BackendSeamSuite extends munit.FunSuite:

  /** A kernel that computes `2·(A·B)` for `gemm` — a deliberately WRONG product, so a
    * routed call is unmistakably distinguishable from the true pure result — and forwards
    * every other method to the pure kernels.
    */
  private object DoublingKernel extends DenseDoubleKernel:
    export PureDenseDoubleKernel.{gemm => _, gemv => _, syrk => _, *}
    var gemmCalls: Int = 0
    def reset(): Unit = gemmCalls = 0
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
        rows, cols, shared, 2.0 * alpha,
        a, aOffset, aRowStride, aColStride,
        b, bOffset, bRowStride, bColStride,
        beta, c, cOffset, cRowStride, cColStride
      )

    def gemv(
        rows: Int, cols: Int, alpha: Double,
        a: DoubleArray, aOffset: Int, rowStride: Int, colStride: Int,
        x: DoubleArray, xOffset: Int, xStride: Int,
        beta: Double, y: DoubleArray, yOffset: Int, yStride: Int
    ): Unit =
      PureDenseDoubleKernel.gemv(
        rows, cols, 2.0 * alpha,
        a, aOffset, rowStride, colStride,
        x, xOffset, xStride,
        beta, y, yOffset, yStride
      )

    def syrk(
        m: Int, k: Int, a: DoubleArray, aOffset: Int, aRowStride: Int,
        c: DoubleArray, cOffset: Int, cRowStride: Int
    ): Unit =
      PureDenseDoubleKernel.syrk(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)
      var i = 0
      while i < k do
        var j = 0
        while j < k do
          val index = cOffset + i * cRowStride + j
          c(index) = 2.0 * c(index)
          j += 1
        i += 1

  private object LowThreshold extends BackendThresholds:
    def nativeGemmMinFlops: Long = 1L
    def nativeGemvMinWork: Long = 1L
    def nativeFactorizationMinSize: Int = 1

  private object HighThreshold extends BackendThresholds:
    def nativeGemmMinFlops: Long = Long.MaxValue
    def nativeGemvMinWork: Long = Long.MaxValue
    def nativeFactorizationMinSize: Int = Int.MaxValue

  /** An accelerating backend (advertises `Vectorized`) whose gemm doubles, with a tunable
    * threshold so both branches of the seam are reachable.
    */
  private final class DoublingBackend(th: BackendThresholds) extends Backend:
    val name: String = "doubling"
    val capabilities: Set[Capability] = Set(Capability.Vectorized)
    val denseDouble: DenseDoubleKernel = DoublingKernel
    val thresholds: BackendThresholds = th
    val config: BackendConfig = BackendConfig.singleThreaded

  /** Same doubling gemm, but advertised via the OTHER `acceleratesGemm` disjunct. */
  private final class DoublingBlasBackend(th: BackendThresholds) extends Backend:
    val name: String = "doubling-blas"
    val capabilities: Set[Capability] = Set(Capability.NativeBlas)
    val denseDouble: DenseDoubleKernel = DoublingKernel
    val thresholds: BackendThresholds = th
    val config: BackendConfig = BackendConfig.singleThreaded

  private object RecordingFactorizations extends DenseDoubleFactorizations:
    var luCalls = 0
    var choleskyCalls = 0
    var qrCalls = 0
    def reset(): Unit =
      luCalls = 0; choleskyCalls = 0; qrCalls = 0
    def lu(a: DMat): Either[LinAlgError, LU] =
      luCalls += 1; DenseDecompositions.lu(a)
    def cholesky(a: DMat): Either[LinAlgError, Cholesky] =
      choleskyCalls += 1; DenseDecompositions.cholesky(a)
    def qr(a: DMat): Either[LinAlgError, QR] =
      qrCalls += 1; Right(DenseDecompositions.qr(a)(using PureBackend))

  private object FactorBackend extends Backend:
    val name = "recording-lapack"
    val capabilities = Set(Capability.NativeLapack)
    val denseDouble = PureDenseDoubleKernel
    val thresholds = LowThreshold
    val config = BackendConfig.singleThreaded
    override val denseFactorizations = Some(RecordingFactorizations)
    override val spectral = Some(new SpectralBackend:
      val name = "recording-spectral"
      val capabilities = Set(SpectralCapability.DenseSymmetricEigen)
    )

  /** A provider that declines every input, so the facades' fallback/validation
    * behaviour is observable in isolation from any real factorization.
    */
  private object DecliningFactorizations extends DenseDoubleFactorizations:
    def lu(a: DMat): Either[LinAlgError, LU] = Left(LinAlgError.UnsupportedOperation("declined lu"))
    def cholesky(a: DMat): Either[LinAlgError, Cholesky] = Left(LinAlgError.UnsupportedOperation("declined cholesky"))
    def qr(a: DMat): Either[LinAlgError, QR] = Left(LinAlgError.UnsupportedOperation("declined qr"))

  private object DecliningBackend extends Backend:
    val name = "declining-lapack"
    val capabilities = Set(Capability.NativeLapack)
    val denseDouble = PureDenseDoubleKernel
    val thresholds = LowThreshold
    val config = BackendConfig.singleThreaded
    override val denseFactorizations = Some(DecliningFactorizations)
    override val spectral = Some(new SpectralBackend:
      val name = "declining-spectral"
      val capabilities = Set(SpectralCapability.DenseSymmetricEigen)
    )

  private val a = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
  private val b = Matrix.dense(2, 2)(5.0, 6.0, 7.0, 8.0)

  test("no backend import → pure path computes the true product") {
    val c = a * b // resolves the companion `given pure`
    assertEqualsDouble(c(0, 0), 19.0, 1e-12) // 1·5 + 2·7
    assertEqualsDouble(c(1, 1), 50.0, 1e-12) // 3·6 + 4·8
  }

  test("an accelerating backend above threshold routes the general product to its gemm") {
    val pure = a * b
    val routed = a.*(b)(using DoublingBackend(LowThreshold))
    assertEqualsDouble(routed(0, 0), 2.0 * pure(0, 0), 1e-12) // doubled ⇒ it routed
    assertEqualsDouble(routed(1, 1), 2.0 * pure(1, 1), 1e-12)
  }

  test("an accelerating backend below threshold stays on the pure path") {
    val pure = a * b
    val notRouted = a.*(b)(using DoublingBackend(HighThreshold))
    assertEqualsDouble(notRouted(0, 0), pure(0, 0), 1e-12) // NOT doubled ⇒ pure
    assertEqualsDouble(notRouted(1, 1), pure(1, 1), 1e-12)
  }

  test("the AᵀA fast-path routes through backend syrk above threshold") {
    val m = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val expected = m.t * m // pure syrk
    val viaBackend = m.t.*(m)(using DoublingBackend(LowThreshold))
    assertEqualsDouble(viaBackend(0, 0), 2.0 * expected(0, 0), 1e-12)
    assertEqualsDouble(viaBackend(1, 0), 2.0 * expected(1, 0), 1e-12)
    assertEqualsDouble(viaBackend(1, 1), 2.0 * expected(1, 1), 1e-12)
  }

  test("the AᵀA fast-path below threshold stays on the pure symmetric kernel") {
    val m = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val expected = m.t * m // pure syrk
    val notRouted = m.t.*(m)(using DoublingBackend(HighThreshold))
    assertEqualsDouble(notRouted(0, 0), expected(0, 0), 1e-12) // NOT doubled ⇒ pure syrk ran
    assertEqualsDouble(notRouted(1, 0), expected(1, 0), 1e-12)
    assertEqualsDouble(notRouted(1, 1), expected(1, 1), 1e-12)
  }

  test("a routed GENERAL product with a strided (transposed) operand doubles correctly") {
    // m1.t is a transposed VIEW (non-unit column stride); m1.t * m2 is NOT AᵀA (different
    // storage), so it takes the general routed path — guarding argument order under strides.
    val m1 = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val m2 = Matrix.dense(3, 2)(1.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    val pure = m1.t * m2
    val routed = m1.t.*(m2)(using DoublingBackend(LowThreshold))
    var i = 0
    while i < 2 do
      var j = 0
      while j < 2 do
        assertEqualsDouble(routed(i, j), 2.0 * pure(i, j), 1e-12)
        j += 1
      i += 1
  }

  test("a `given` accelerating backend in scope routes a plain `a * b` (the ergonomic pattern)") {
    given Backend = DoublingBackend(LowThreshold) // shadows the companion `given pure`
    val pure = a.*(b)(using PureBackend)
    val routed = a * b // resolves the local given, not pure
    assertEqualsDouble(routed(0, 0), 2.0 * pure(0, 0), 1e-12)
    assertEqualsDouble(routed(1, 1), 2.0 * pure(1, 1), 1e-12)
  }

  test("the NativeBlas disjunct also routes (not just Vectorized)") {
    val pure = a * b
    val routed = a.*(b)(using DoublingBlasBackend(LowThreshold))
    assertEqualsDouble(routed(0, 0), 2.0 * pure(0, 0), 1e-12)
  }

  test("standalone DMat * DVec routes through backend gemv above threshold") {
    val x = Vec(2.0, -1.0)
    val pure = a.*(x)(using PureBackend)
    val routed = a.*(x)(using DoublingBackend(LowThreshold))
    assertEqualsDouble(routed(0), 2.0 * pure(0), 1e-12)
    assertEqualsDouble(routed(1), 2.0 * pure(1), 1e-12)
  }

  test("blocked QR routes its two trailing-update products through backend gemm") {
    DoublingKernel.reset()
    // The blocked path is selected by min(rows, cols) >= 96.
    val large = Matrix.tabulate(128, 96)((i, j) => if i == j then 3.0 else ((i * 7 + j * 5) % 11).toDouble / 11.0)
    large.qr(using DoublingBackend(LowThreshold))
    assert(DoublingKernel.gemmCalls > 0, "blocked QR never reached the backend gemm seam")
  }

  test("native factorization provider receives LU, Cholesky, QR, solve, and least-squares facades") {
    RecordingFactorizations.reset()
    val spd = Matrix.dense(2, 2)(4.0, 1.0, 1.0, 3.0)
    val rhs = Vec(1.0, 2.0)
    assert(spd.lu(using FactorBackend).isRight)
    assert(spd.cholesky(using FactorBackend).isRight)
    spd.qr(using FactorBackend)
    assert(spd.solve(rhs)(using FactorBackend).isRight)
    assert(spd.leastSquares(rhs)(using FactorBackend).isRight)
    assertEquals(RecordingFactorizations.luCalls, 2)
    assertEquals(RecordingFactorizations.choleskyCalls, 1)
    assertEquals(RecordingFactorizations.qrCalls, 2)
  }

  test("qr stays total: a provider Left falls back to the pure path instead of throwing") {
    val m = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val pure = m.qr(using PureBackend)
    val viaDeclining = m.qr(using DecliningBackend) // must not throw
    assertEquals(viaDeclining.diagnostics.rank, pure.diagnostics.rank)
    assert(m.leastSquares(Vec(1.0, 2.0, 3.0))(using DecliningBackend).isRight)
  }

  test("qrWith never routes to a provider — the workspace facade is allocation-controlled") {
    RecordingFactorizations.reset()
    val m = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    m.qrWith(DenseWorkspace.forQR(3, 2))(using FactorBackend)
    assertEquals(RecordingFactorizations.qrCalls, 0)
  }

  test("lu/cholesky on a non-square matrix return the typed Left before any provider is invoked") {
    RecordingFactorizations.reset()
    val rect = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    rect.lu(using FactorBackend) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected Left(NonSquareMatrix), got $other")
    rect.cholesky(using FactorBackend) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected Left(NonSquareMatrix), got $other")
    assertEquals(RecordingFactorizations.luCalls, 0)
    assertEquals(RecordingFactorizations.choleskyCalls, 0)
  }

  test("rankEstimate agrees with qr.diagnostics.rank under an accelerating backend") {
    given Backend = DoublingBackend(LowThreshold)
    // Blocked path (min(rows, cols) >= 96), so the ambient backend's gemm is in play.
    val large = Matrix.tabulate(128, 96)((i, j) => if i == j then 3.0 else ((i * 7 + j * 5) % 11).toDouble / 11.0)
    assertEquals(large.rankEstimate, large.qr.diagnostics.rank.get)
  }
