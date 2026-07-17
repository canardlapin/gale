package gale.laws

import gale.backend.*
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Reusable, capability-driven conformance suite for optional [[Backend]] modules.
  *
  * The kernel checks use independently accumulated `BigDecimal` references rather
  * than Gale's pure kernel, so a defect shared by two implementations cannot make
  * the suite self-confirming. Fixtures cover offsets, padded leading dimensions,
  * transposed storage, non-unit vector strides, `alpha`/`beta`, and both triangles
  * of `syrk`. A backend advertising `NativeLapack` additionally has to reconstruct
  * LU/Cholesky/QR factors and solve a system with a small residual.
  */
abstract class BackendConformanceSuite extends munit.FunSuite:
  def backend: Backend

  /** Backends may loosen this only with a documented numerical reason. */
  def relativeTolerance: Double = 2e-10
  def absoluteTolerance: Double = 2e-12

  private def close(actual: Double, expected: Double, clue: String): Unit =
    val tolerance = absoluteTolerance + relativeTolerance * math.max(math.abs(actual), math.abs(expected))
    assert(actual.isFinite == expected.isFinite, s"$clue: finiteness mismatch: $actual vs $expected")
    assert(math.abs(actual - expected) <= tolerance,
      s"$clue: $actual vs $expected (|delta|=${math.abs(actual - expected)}, tolerance=$tolerance)")

  private def decimalDot(n: Int)(left: Int => Double, right: Int => Double): Double =
    var total = BigDecimal(0)
    var i = 0
    while i < n do
      total += BigDecimal(left(i)) * BigDecimal(right(i))
      i += 1
    total.toDouble

  test("backend metadata and capability invariants conform"):
    assertEquals(Backend.validationErrors(backend), Nil)
    assert(Backend.requireValid(backend) eq backend)

  test("dense gemm conforms against an independent strided alpha/beta oracle"):
    val rows = 4
    val cols = 3
    val shared = 5
    val aOffset = 2
    val aRowStride = 8
    val aColStride = 1
    val bOffset = 1
    val bRowStride = 1       // logical B is a transpose view of row-major storage
    val bColStride = 7
    val cOffset = 2
    val cRowStride = 6
    val cColStride = 1
    val alpha = 0.75
    val beta = -0.5
    val a = DoubleArray.alloc(40)
    val b = DoubleArray.alloc(30)
    val c = DoubleArray.alloc(30)
    val before = DoubleArray.alloc(30)

    var i = 0
    while i < rows do
      var k = 0
      while k < shared do
        a(aOffset + i * aRowStride + k * aColStride) = (i * 7 - k * 3 + 2).toDouble / 5.0
        k += 1
      i += 1
    var k = 0
    while k < shared do
      var j = 0
      while j < cols do
        b(bOffset + k * bRowStride + j * bColStride) = (k * 5 + j * 11 - 4).toDouble / 7.0
        j += 1
      k += 1
    i = 0
    while i < rows do
      var j = 0
      while j < cols do
        val value = (i - 2 * j + 1).toDouble / 3.0
        c(cOffset + i * cRowStride + j) = value
        before(cOffset + i * cRowStride + j) = value
        j += 1
      i += 1

    backend.denseDouble.gemm(
      rows, cols, shared, alpha,
      a, aOffset, aRowStride, aColStride,
      b, bOffset, bRowStride, bColStride,
      beta, c, cOffset, cRowStride, cColStride
    )

    i = 0
    while i < rows do
      var j = 0
      while j < cols do
        val sum = decimalDot(shared)(
          kk => a(aOffset + i * aRowStride + kk * aColStride),
          kk => b(bOffset + kk * bRowStride + j * bColStride)
        )
        val expected = alpha * sum + beta * before(cOffset + i * cRowStride + j)
        close(c(cOffset + i * cRowStride + j), expected, s"gemm($i,$j)")
        j += 1
      i += 1

  test("dense gemv conforms with padded matrix and non-unit vector strides"):
    val rows = 5
    val cols = 7
    val aOffset = 3
    val rowStride = 10
    val xOffset = 1
    val xStride = 2
    val yOffset = 2
    val yStride = 3
    val alpha = -1.25
    val beta = 0.4
    val a = DoubleArray.alloc(60)
    val x = DoubleArray.alloc(20)
    val y = DoubleArray.alloc(20)
    val before = DoubleArray.alloc(20)
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        a(aOffset + i * rowStride + j) = (i * 13 + j * 3 - 8).toDouble / 9.0
        j += 1
      val yValue = (i - 2).toDouble / 6.0
      y(yOffset + i * yStride) = yValue
      before(yOffset + i * yStride) = yValue
      i += 1
    var j = 0
    while j < cols do
      x(xOffset + j * xStride) = (j * 5 - 4).toDouble / 8.0
      j += 1

    backend.denseDouble.gemv(
      rows, cols, alpha, a, aOffset, rowStride, 1,
      x, xOffset, xStride, beta, y, yOffset, yStride
    )

    i = 0
    while i < rows do
      val sum = decimalDot(cols)(
        jj => a(aOffset + i * rowStride + jj),
        jj => x(xOffset + jj * xStride)
      )
      val expected = alpha * sum + beta * before(yOffset + i * yStride)
      close(y(yOffset + i * yStride), expected, s"gemv($i)")
      i += 1

  test("dense syrk conforms against an independent full-symmetric oracle"):
    val m = 7
    val k = 4
    val aOffset = 1
    val aRowStride = 7
    val cOffset = 2
    val cRowStride = 6
    val a = DoubleArray.alloc(52)
    val c = DoubleArray.alloc(30)
    var row = 0
    while row < m do
      var col = 0
      while col < k do
        a(aOffset + row * aRowStride + col) = (row * 7 - col * 5 + 3).toDouble / 11.0
        col += 1
      row += 1

    backend.denseDouble.syrk(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)

    var i = 0
    while i < k do
      var j = 0
      while j < k do
        val expected = decimalDot(m)(
          r => a(aOffset + r * aRowStride + i),
          r => a(aOffset + r * aRowStride + j)
        )
        close(c(cOffset + i * cRowStride + j), expected, s"syrk($i,$j)")
        j += 1
      i += 1

  test("advertised native factorizations reconstruct and solve"):
    if backend.capabilities.contains(Capability.NativeLapack) then
      val a = Matrix.dense(3, 3)(4.0, 2.0, -1.0, 1.0, 5.0, 2.0, 2.0, 1.0, 6.0)
      val b = Vec(7.0, -1.0, 4.0)
      val lu = backend.denseFactorizations.get.lu(a).orThrow
      val l = Matrix.tabulate(3, 3)((i, j) => if i == j then 1.0 else if i > j then lu.packed(i, j) else 0.0)
      val u = Matrix.tabulate(3, 3)((i, j) => if i <= j then lu.packed(i, j) else 0.0)
      val pa = Matrix.tabulate(3, 3)((i, j) => a(lu.pivots(i), j))
      MatrixLaws.assertCloseRel(l.*(u)(using PureBackend), pa, relativeTolerance)
      val x = lu.solve(b).orThrow
      VecLaws.assertCloseRel(a.*(x)(using PureBackend), b, relativeTolerance)

      val spd = Matrix.dense(3, 3)(6.0, 2.0, 1.0, 2.0, 5.0, 2.0, 1.0, 2.0, 4.0)
      val chol = backend.denseFactorizations.get.cholesky(spd).orThrow
      MatrixLaws.assertCloseRel(chol.lower.*(chol.lower.t)(using PureBackend), spd, relativeTolerance)

      val rectangular = Matrix.dense(4, 3)(
        1.0, 2.0, -1.0,
        3.0, 0.5, 4.0,
        -2.0, 1.0, 3.0,
        0.25, -1.5, 2.0
      )
      val qr = backend.denseFactorizations.get.qr(rectangular).orThrow
      MatrixLaws.assertCloseRel(qr.q.*(qr.r)(using PureBackend), rectangular, 5.0 * relativeTolerance)
