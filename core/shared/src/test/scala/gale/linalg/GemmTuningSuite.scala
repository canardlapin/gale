package gale.linalg

/** Correctness guards for the tuned dense kernels: the `Aᵀ·A` syrk fast path
  * (routed from `a.t * a`) and the unroll-and-jam gemm across shapes whose row
  * count is not a multiple of the 4-row unroll (exercising the scalar tail).
  */
class GemmTuningSuite extends munit.FunSuite:

  private def randomMat(m: Int, n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(m, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  private def naiveMul(a: DMat, b: DMat): DMat =
    Matrix.tabulate(a.rows, b.cols): (i, j) =>
      var s = 0.0
      var k = 0
      while k < a.cols do
        s += a(i, k) * b(k, j)
        k += 1
      s

  private def assertMatrixClose(actual: DMat, expected: DMat, tol: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        assert(math.abs(actual(i, j) - expected(i, j)) < tol, s"[$i,$j] ${actual(i, j)} != ${expected(i, j)}")
        j += 1
      i += 1

  test("syrk path: a.t * a matches the naive transpose product (tall/square/wide)") {
    for (m, n, seed) <- Seq((12, 4, 1L), (5, 5, 2L), (3, 7, 3L), (17, 6, 4L), (11, 9, 5L)) do
      val a = randomMat(m, n, seed)
      val gram = a.t * a // routes to dsyrk when a is row-major
      val reference = naiveMul(a.t, a)
      assertEquals(gram.rows, n)
      assertEquals(gram.cols, n)
      assertMatrixClose(gram, reference, 1e-12)
      // Symmetry is exact (mirrored, not recomputed).
      var i = 0
      while i < n do
        var j = 0
        while j < n do
          assertEquals(gram(i, j), gram(j, i), s"asymmetry at [$i,$j]")
          j += 1
        i += 1
  }

  test("syrk tiled path obeys quadratic scaling at extreme finite magnitudes") {
    val base = randomMat(9, 7, 20260710L)
    for scale <- Seq(1e100, 1e-100) do
      val scaled = Matrix.tabulate(base.rows, base.cols)((i, j) => scale * base(i, j))
      val actual = scaled.t * scaled
      val expected = naiveMul(base.t, base)
      val scale2 = scale * scale
      var i = 0
      while i < actual.rows do
        var j = 0
        while j < actual.cols do
          val e = scale2 * expected(i, j)
          val err = math.abs(actual(i, j) - e)
          val tolerance = 2e-14 * math.max(math.abs(e), scale2)
          assert(actual(i, j).isFinite, s"non-finite gram entry at [$i,$j], scale=$scale")
          assert(err <= tolerance, s"[$i,$j] $err > $tolerance at scale=$scale")
          j += 1
        i += 1
  }

  test("syrk path: strided (non-unit-column) left operand falls back correctly") {
    // A transposed view whose base is itself a transpose is not unit-column-major,
    // so `*` must fall back to the general gemm and still be correct.
    val base = randomMat(6, 4, 9L)
    val a = base.t // 4x6, rowStride 1, colStride 6 (not unit column stride)
    val gram = a.t * a // (a.t is 6x4 view of base); general gemm path
    assertMatrixClose(gram, naiveMul(a.t, a), 1e-12)
  }

  test("gemm all-tail partition: whole matrix smaller than the 4x4 tile") {
    // rows <= 3 AND cols <= 3: the register-tile main loops never run, so the
    // entire product goes through the leftover-column and scalar-row tails.
    for (m, k, n, seed) <- Seq((3, 3, 3, 10L), (2, 5, 3, 11L), (1, 4, 2, 12L), (3, 1, 1, 13L)) do
      val a = randomMat(m, k, seed)
      val b = randomMat(k, n, seed + 100)
      assertMatrixClose(a * b, naiveMul(a, b), 1e-12)
  }

  test("gemm unroll-and-jam: correct for row counts not divisible by 4") {
    // shapes covering row-remainders 0..3 (4, 5, 6, 7, 13 rows)
    for (m, k, n) <- Seq((4, 3, 5), (5, 4, 6), (6, 7, 2), (7, 5, 8), (13, 9, 11)) do
      val a = randomMat(m, k, m * 31L + k)
      val b = randomMat(k, n, n * 17L + k)
      assertMatrixClose(a * b, naiveMul(a, b), 1e-11)
  }
