package gale.backend

import gale.linalg.*
import gale.spectral.SpectralBackend
import gale.spectral.SpectralCapability

/** The core backend contract (doc §A): the default `given` is the pure backend,
  * capability discovery and `compose` behave as specified, and `PureBackend`'s kernel
  * surface computes the same results as the public operators it forwards to. No facade
  * seam is wired into `DMat.*` yet — this is the contract foundation the acceleration
  * modules (vector, FFM) build on, so the no-import path is provably unchanged.
  */
class BackendContractSuite extends munit.FunSuite:

  /** Distinct, finite thresholds so `compose` routing is observable by value (pure's
    * are all `MaxValue`).
    */
  private object FakeThresholds extends BackendThresholds:
    def nativeGemmMinFlops: Long = 1000L
    def nativeGemvMinWork: Long = 1000L
    def nativeFactorizationMinSize: Int = 64

  /** A test-only backend that merely CLAIMS `Vectorized` (its kernels are still pure) —
    * enough to exercise capability/dispatch logic without a real SIMD implementation.
    */
  private object FakeVectorBackend extends Backend:
    val name: String = "fake-vector"
    val capabilities: Set[Capability] = Set(Capability.Vectorized)
    val denseDouble: DenseDoubleKernel = PureDenseDoubleKernel
    val thresholds: BackendThresholds = FakeThresholds
    val config: BackendConfig = BackendConfig.singleThreaded

  /** Claims `NativeBlas` (the OTHER `acceleratesGemm` disjunct). */
  private object FakeNativeBlasBackend extends Backend:
    val name: String = "fake-blas"
    val capabilities: Set[Capability] = Set(Capability.NativeBlas)
    val denseDouble: DenseDoubleKernel = PureDenseDoubleKernel
    val thresholds: BackendThresholds = FakeThresholds
    val config: BackendConfig = BackendConfig.singleThreaded

  /** A minimal spectral provider (every op defaults to Left(unsupported)). */
  private object FakeSpectral extends SpectralBackend:
    val name: String = "fake-spectral"
    val capabilities: Set[SpectralCapability] = Set(SpectralCapability.GeneralizedNonsymmetricEigen)

  /** A backend that also provides the spectral half (the A.5 bridge). */
  private object FakeLapackBackend extends Backend:
    val name: String = "fake-lapack"
    val capabilities: Set[Capability] = Set(Capability.NativeLapack)
    val denseDouble: DenseDoubleKernel = PureDenseDoubleKernel
    val thresholds: BackendThresholds = FakeThresholds
    val config: BackendConfig = BackendConfig.singleThreaded
    override val spectral: Option[SpectralBackend] = Some(FakeSpectral)

  test("the default given Backend is PureBackend") {
    assert(summon[Backend] eq PureBackend)
  }

  test("PureBackend: pure-path capabilities, does not accelerate gemm, infinite thresholds") {
    assertEquals(PureBackend.capabilities, Set(Capability.WasmFriendly, Capability.Deterministic))
    assert(!PureBackend.acceleratesGemm, "the pure path is not a gemm accelerator")
    assert(PureBackend.spectral.isEmpty, "pure provides no spectral half")
    assertEquals(PureBackend.config, BackendConfig.singleThreaded)
    assertEquals(PureBackend.thresholds.nativeGemmMinFlops, Long.MaxValue)
    assertEquals(PureBackend.thresholds.nativeFactorizationMinSize, Int.MaxValue)
  }

  test("acceleratesGemm triggers on either disjunct: Vectorized and NativeBlas") {
    assert(FakeVectorBackend.acceleratesGemm, "Vectorized must trigger")
    assert(FakeNativeBlasBackend.acceleratesGemm, "NativeBlas must trigger")
  }

  test("BackendConfig factories") {
    assertEquals(BackendConfig.singleThreaded, BackendConfig(1, 1))
    assertEquals(BackendConfig.dedicated(8), BackendConfig(8, 8))
    assertEquals(BackendConfig.singleThreaded.allowNestedParallelism, false)
  }

  test("compose: unions capabilities and routes the coarse surface to the accelerating part") {
    val composed = Backend.compose(PureBackend, FakeVectorBackend)
    assert(composed.capabilities.contains(Capability.Vectorized))
    assert(composed.capabilities.contains(Capability.WasmFriendly))
    assert(composed.acceleratesGemm, "composite advertises the accelerator")
    assert(composed.spectral.isEmpty, "no part provides spectral")
    // routed to the accelerating part's thresholds (finite), not pure's MaxValue:
    assertEquals(composed.thresholds.nativeGemmMinFlops, 1000L)
  }

  test("compose of a single backend is that backend") {
    assert(Backend.compose(PureBackend) eq PureBackend)
  }

  test("compose: multiple parts union all capabilities; all-pure compose accelerates nothing") {
    val composed = Backend.compose(PureBackend, FakeNativeBlasBackend, FakeLapackBackend)
    assert(composed.capabilities.contains(Capability.NativeBlas))
    assert(composed.capabilities.contains(Capability.NativeLapack))
    assert(composed.capabilities.contains(Capability.WasmFriendly))
    assert(composed.acceleratesGemm)
    // gemmProvider = first accelerating part (FakeNativeBlasBackend) → its finite threshold
    assertEquals(composed.thresholds.nativeGemmMinFlops, 1000L)

    // A compose where no part accelerates gemm falls back to the primary's kernels.
    val allPure = Backend.compose(PureBackend, PureBackend)
    assert(!allPure.acceleratesGemm)
    assertEquals(allPure.thresholds.nativeGemmMinFlops, Long.MaxValue)
  }

  test("compose surfaces the spectral bridge from a providing part (A.5)") {
    val composed = Backend.compose(PureBackend, FakeLapackBackend)
    assert(composed.spectral.isDefined, "the composite must expose the spectral half")
    assert(composed.spectral.get.capabilities.contains(SpectralCapability.GeneralizedNonsymmetricEigen))
  }

  test("PureBackend.denseDouble.gemm computes the same product as DMat.*") {
    val a = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val b = Matrix.dense(3, 2)(7.0, 8.0, 9.0, 10.0, 11.0, 12.0)
    val expected = a * b
    val c = Matrix.dense(2, 2)(0.0, 0.0, 0.0, 0.0)
    PureBackend.denseDouble.gemm(
      a.rows, b.cols, a.cols, 1.0,
      a.data, a.offset.value, a.rowStride.value, a.colStride.value,
      b.data, b.offset.value, b.rowStride.value, b.colStride.value,
      0.0,
      c.data, c.offset.value, c.rowStride.value, c.colStride.value
    )
    var i = 0
    while i < 2 do
      var j = 0
      while j < 2 do
        assertEqualsDouble(c(i, j), expected(i, j), 1e-12)
        j += 1
      i += 1
  }

  test("PureBackend.denseDouble.syrk computes the same as A.t * A (upper triangle)") {
    val a = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val expected = a.t * a // 2x2 symmetric
    val c = Matrix.dense(2, 2)(0.0, 0.0, 0.0, 0.0)
    // C(k×k) = AᵀA where A is m rows × k cols: here m = a.rows, k = a.cols.
    PureBackend.denseDouble.syrk(
      a.rows, a.cols,
      a.data, a.offset.value, a.rowStride.value,
      c.data, c.offset.value, c.rowStride.value
    )
    // dsyrkRowMajor computes the UPPER triangle then mirrors it into the lower, so
    // the result is the FULL symmetric AᵀA (both triangles populated).
    assertEqualsDouble(c(0, 0), expected(0, 0), 1e-12)
    assertEqualsDouble(c(0, 1), expected(0, 1), 1e-12)
    assertEqualsDouble(c(1, 0), expected(1, 0), 1e-12) // mirrored lower
    assertEqualsDouble(c(1, 1), expected(1, 1), 1e-12)
  }

  test("PureBackend.denseDouble.gemm honours a strided (transposed) operand") {
    val a = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val at = a.t // 2x3 VIEW with non-unit column stride — exercises dgemmStrided
    val b = Matrix.dense(3, 2)(1.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    val expected = at * b // 2x2
    val c = Matrix.dense(2, 2)(0.0, 0.0, 0.0, 0.0)
    PureBackend.denseDouble.gemm(
      at.rows, b.cols, at.cols, 1.0,
      at.data, at.offset.value, at.rowStride.value, at.colStride.value,
      b.data, b.offset.value, b.rowStride.value, b.colStride.value,
      0.0,
      c.data, c.offset.value, c.rowStride.value, c.colStride.value
    )
    var i = 0
    while i < 2 do
      var j = 0
      while j < 2 do
        assertEqualsDouble(c(i, j), expected(i, j), 1e-12)
        j += 1
      i += 1
  }

  test("PureBackend.denseDouble L1/L2 forwards match the public ops") {
    val x = Vec(1.0, 2.0, 3.0, 4.0)
    val y = Vec(5.0, 6.0, 7.0, 8.0)
    val n = x.length

    // dot / nrm2
    assertEqualsDouble(
      PureBackend.denseDouble.dot(n, x.data, x.offset.value, x.stride.value, y.data, y.offset.value, y.stride.value),
      x.dot(y),
      1e-12
    )
    assertEqualsDouble(PureBackend.denseDouble.nrm2(n, x.data, x.offset.value, x.stride.value), x.norm2, 1e-12)

    // gemv: A(2x4) * x
    val mat = Matrix.dense(2, 4)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
    val expectedMv = mat * x
    val yv = Vec.fill(2)(0.0)
    PureBackend.denseDouble.gemv(
      mat.rows, mat.cols, 1.0,
      mat.data, mat.offset.value, mat.rowStride.value, mat.colStride.value,
      x.data, x.offset.value, x.stride.value,
      0.0, yv.data, yv.offset.value, yv.stride.value
    )
    assertEqualsDouble(yv(0), expectedMv(0), 1e-12)
    assertEqualsDouble(yv(1), expectedMv(1), 1e-12)

    // axpy: yb := 2*x + yb (fresh buffer)
    val yb = Vec(5.0, 6.0, 7.0, 8.0)
    PureBackend.denseDouble.axpy(
      n, 2.0, x.data, x.offset.value, x.stride.value, yb.data, yb.offset.value, yb.stride.value
    )
    assertEqualsDouble(yb(0), 2.0 * 1.0 + 5.0, 1e-12)
    assertEqualsDouble(yb(3), 2.0 * 4.0 + 8.0, 1e-12)

    // scal: xb := 3*xb
    val xb = Vec(1.0, 2.0, 3.0, 4.0)
    PureBackend.denseDouble.scal(n, 3.0, xb.data, xb.offset.value, xb.stride.value)
    assertEqualsDouble(xb(2), 9.0, 1e-12)

    // copy: yc := x
    val yc = Vec.fill(n)(0.0)
    PureBackend.denseDouble.copy(n, x.data, x.offset.value, x.stride.value, yc.data, yc.offset.value, yc.stride.value)
    assertEqualsDouble(yc(1), 2.0, 1e-12)
  }
