package gale.backend

import gale.linalg.*
import gale.platform.DoubleArray

/** The coarse gemm dispatch seam on `DMat.*` (doc §A.2). With no acceleration import the
  * `given` resolves to `Backend.pure` and the pure kernel runs — the entire existing test
  * suite is the byte-identical witness for that; these tests add the ROUTING behaviour: an
  * imported accelerating backend takes over the general product above its threshold, a
  * below-threshold product stays pure, and the structural `AᵀA` fast-path is never routed.
  */
class BackendSeamSuite extends munit.FunSuite:

  /** A kernel that computes `2·(A·B)` for `gemm` — a deliberately WRONG product, so a
    * routed call is unmistakably distinguishable from the true pure result — and forwards
    * every other method to the pure kernels.
    */
  private object DoublingKernel extends DenseDoubleKernel:
    export PureDenseDoubleKernel.{gemm => _, *}
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
      PureDenseDoubleKernel.gemm(
        rows, cols, shared, 2.0 * alpha,
        a, aOffset, aRowStride, aColStride,
        b, bOffset, bRowStride, bColStride,
        beta, c, cOffset, cRowStride, cColStride
      )

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

  test("the AᵀA fast-path is never routed, even with an accelerating backend below-or-above threshold") {
    val m = Matrix.dense(3, 2)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val expected = m.t * m // pure syrk
    val viaBackend = m.t.*(m)(using DoublingBackend(LowThreshold))
    // structural AᵀA takes the dedicated pure symmetric kernel, so it is NOT doubled.
    assertEqualsDouble(viaBackend(0, 0), expected(0, 0), 1e-12)
    assertEqualsDouble(viaBackend(1, 0), expected(1, 0), 1e-12)
    assertEqualsDouble(viaBackend(1, 1), expected(1, 1), 1e-12)
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
