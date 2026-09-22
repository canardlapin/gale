package gale.laws

import gale.backend.*
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Reusable, capability-driven conformance suite for optional `Backend` modules.
  *
  * The kernel checks use independently accumulated `BigDecimal` references rather
  * than Gale's pure kernel, so a defect shared by two implementations cannot make
  * the suite self-confirming. Fixtures cover offsets, padded leading dimensions,
  * transposed storage, non-unit vector strides, `alpha`/`beta`, and both triangles
  * of `syrk`. A backend advertising `NativeLapack` additionally has to reconstruct
  * LU/Cholesky/QR factors and solve a system with a small residual. Every backend
  * must also agree with `PureBackend` on facade residuals and on the `LinAlgError`
  * class for IEEE-exact singular / non-SPD plants; a `NativeLapack` provider is
  * checked on that same residual and error-class contract directly.
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

  test("facade products and solves agree with PureBackend within residual tolerance"):
    val a = Matrix.dense(4, 4)(
      4.0, 1.0, 0.0, 2.0,
      1.0, 5.0, 1.0, 0.0,
      0.0, 1.0, 6.0, 1.0,
      2.0, 0.0, 1.0, 7.0
    )
    val c = Matrix.dense(4, 3)(
      1.0, -1.0, 0.5,
      0.0, 2.0, 1.0,
      -0.5, 0.0, 3.0,
      1.5, 0.25, -1.0
    )
    val b = Vec(3.0, -1.0, 2.0, 4.0)
    MatrixLaws.assertCloseRel(a.*(c)(using backend), a.*(c)(using PureBackend), relativeTolerance)
    VecLaws.assertCloseRel(a.*(b)(using backend), a.*(b)(using PureBackend), relativeTolerance)

    val xBackend = a.solve(b)(using backend).orThrow
    val xPure = a.solve(b)(using PureBackend).orThrow
    VecLaws.assertCloseRel(xBackend, xPure, relativeTolerance)
    assert(
      solveResidual(a, xBackend, b) <= 1e-12,
      s"backend solve residual ${solveResidual(a, xBackend, b)}"
    )
    assert(
      solveResidual(a, xPure, b) <= 1e-12,
      s"pure solve residual ${solveResidual(a, xPure, b)}"
    )

    val spd = Matrix.dense(3, 3)(6.0, 2.0, 1.0, 2.0, 5.0, 2.0, 1.0, 2.0, 4.0)
    val rhs = Vec(1.0, 0.0, -1.0)
    val cholBackend = spd.cholesky(using backend).orThrow.solve(rhs).orThrow
    val cholPure = spd.cholesky(using PureBackend).orThrow.solve(rhs).orThrow
    VecLaws.assertCloseRel(cholBackend, cholPure, relativeTolerance)
    assert(solveResidual(spd, cholBackend, rhs) <= 1e-12)

    val tall = Matrix.dense(5, 3)(
      1.0, 0.0, 2.0,
      0.5, 3.0, -1.0,
      2.0, 1.0, 0.0,
      -1.0, 0.5, 1.5,
      0.0, 2.0, 1.0
    )
    val observations = Vec(1.0, -2.0, 0.5, 3.0, 0.25)
    val lsBackend = tall.qr(using backend).solveLeastSquares(observations).orThrow
    val lsPure = tall.qr(using PureBackend).solveLeastSquares(observations).orThrow
    VecLaws.assertCloseRel(lsBackend, lsPure, 5.0 * relativeTolerance)

  test("facade typed errors agree with PureBackend on IEEE-exact plants"):
    val singular = Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      2.0, 4.0, 6.0,
      3.0, 6.0, 9.0
    )
    assertSameErrorClass(
      singular.lu(using backend),
      singular.lu(using PureBackend),
      "LU exact rank-1"
    )
    val indefinite = Matrix.dense(2, 2)(1.0, 2.0, 2.0, 1.0)
    assertSameErrorClass(
      indefinite.cholesky(using backend),
      indefinite.cholesky(using PureBackend),
      "Cholesky indefinite"
    )
    val rectangular = Matrix.zeros(2, 3)
    assertSameErrorClass(
      rectangular.lu(using backend),
      rectangular.lu(using PureBackend),
      "LU non-square"
    )

  test("native factorizations match PureBackend residuals and error classes"):
    if backend.capabilities.contains(Capability.NativeLapack) then
      val provider = backend.denseFactorizations.get
      val a = Matrix.dense(4, 4)(
        4.0, 1.0, 0.0, 2.0,
        1.0, 5.0, 1.0, 0.0,
        0.0, 1.0, 6.0, 1.0,
        2.0, 0.0, 1.0, 7.0
      )
      val b = Vec(3.0, -1.0, 2.0, 4.0)
      val xNative = provider.lu(a).orThrow.solve(b).orThrow
      val xPure = a.solve(b)(using PureBackend).orThrow
      VecLaws.assertCloseRel(xNative, xPure, relativeTolerance)
      assert(solveResidual(a, xNative, b) <= 1e-12, s"native LU residual ${solveResidual(a, xNative, b)}")
      assert(solveResidual(a, xPure, b) <= 1e-12, s"pure LU residual ${solveResidual(a, xPure, b)}")

      val spd = Matrix.dense(3, 3)(6.0, 2.0, 1.0, 2.0, 5.0, 2.0, 1.0, 2.0, 4.0)
      val rhs = Vec(1.0, 0.0, -1.0)
      val xCholNative = provider.cholesky(spd).orThrow.solve(rhs).orThrow
      val xCholPure = spd.cholesky(using PureBackend).orThrow.solve(rhs).orThrow
      VecLaws.assertCloseRel(xCholNative, xCholPure, relativeTolerance)
      assert(solveResidual(spd, xCholNative, rhs) <= 1e-12)

      val rectangular = Matrix.dense(5, 3)(
        1.0, 2.0, -1.0,
        3.0, 0.5, 4.0,
        -2.0, 1.0, 3.0,
        0.25, -1.5, 2.0,
        1.5, 0.0, -0.5
      )
      val qrNative = provider.qr(rectangular).orThrow
      val qrPure = rectangular.qr(using PureBackend)
      MatrixLaws.assertCloseRel(qrNative.q.*(qrNative.r)(using PureBackend), rectangular, 5.0 * relativeTolerance)
      MatrixLaws.assertCloseRel(qrPure.q.*(qrPure.r)(using PureBackend), rectangular, 5.0 * relativeTolerance)

      val singular = Matrix.dense(3, 3)(
        1.0, 2.0, 3.0,
        2.0, 4.0, 6.0,
        3.0, 6.0, 9.0
      )
      assertSameErrorClass(provider.lu(singular), singular.lu(using PureBackend), "native LU exact rank-1")
      val indefinite = Matrix.dense(2, 2)(1.0, 2.0, 2.0, 1.0)
      assertSameErrorClass(
        provider.cholesky(indefinite),
        indefinite.cholesky(using PureBackend),
        "native Cholesky indefinite"
      )

  private def frobenius(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  private def solveResidual(a: DMat, x: DVec, b: DVec): Double =
    val residual = a.*(x)(using PureBackend) - b
    residual.norm2 / (1.0 + frobenius(a) * x.norm2 + b.norm2)

  private def assertSameErrorClass[A](
      actual: Either[LinAlgError, A],
      expected: Either[LinAlgError, A],
      clue: String
  ): Unit =
    (actual, expected) match
      case (Left(got), Left(want)) =>
        assertEquals(got.getClass, want.getClass, s"$clue: ${got.getClass.getName} vs ${want.getClass.getName}")
      case (Right(_), Right(_)) =>
        fail(s"$clue: both sides succeeded; expected a typed error")
      case (got, want) =>
        fail(s"$clue: error-class mismatch: $got vs $want")
