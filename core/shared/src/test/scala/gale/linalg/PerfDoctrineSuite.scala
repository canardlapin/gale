package gale.linalg

import gale.TestAccess
import scala.util.Random

/** Equivalence tests for the P5 dense performance-doctrine rewrites: matrix
  * add/subtract kernels (K1), scaled `norm2` (K8), the gemm fast paths (K9), and
  * the hoisted factorization back/forward substitutions (K2).
  */
class PerfDoctrineSuite extends munit.FunSuite:
  // ---- K1: DMat + / - through dadd/dsub kernels ----------------------------

  private def refAddSub(a: DMat, b: DMat, subtract: Boolean): DMat =
    Matrix.tabulate(a.rows, a.cols)((i, j) => if subtract then a(i, j) - b(i, j) else a(i, j) + b(i, j))

  private def assertAddSubMatches(a: DMat, b: DMat): Unit =
    // Elementwise, so no summation reordering: results must be bit-exact.
    assertEquals((a + b).valuesRowMajor, refAddSub(a, b, subtract = false).valuesRowMajor)
    assertEquals((a - b).valuesRowMajor, refAddSub(a, b, subtract = true).valuesRowMajor)

  private def randStrided(
      rng: Random,
      rows: Int,
      cols: Int,
      rowStride: Int,
      colStride: Int,
      offset: Int
  ): DMat =
    val maxIdx = offset + math.max(rows - 1, 0) * rowStride + math.max(cols - 1, 0) * colStride
    val backing = TestAccess.doubleArray(Seq.fill(maxIdx + 1)(rng.nextDouble() * 2.0 - 1.0)*)
    TestAccess.mat(backing, offset, rows, cols, rowStride, colStride)

  test("DMat + and - match a dense recompute across storage layouts") {
    val rng = new Random(1L)

    // Both contiguous row-major -> single flat kernel call.
    val a = Matrix.dense(3, 4, Seq.fill(12)(rng.nextDouble() * 2.0 - 1.0))
    val b = Matrix.dense(3, 4, Seq.fill(12)(rng.nextDouble() * 2.0 - 1.0))
    assertAddSubMatches(a, b)

    // One transposed (column-major view) -> per-row strided kernel call.
    val c = Matrix.dense(4, 3, Seq.fill(12)(rng.nextDouble() * 2.0 - 1.0)).t
    assertAddSubMatches(a, c)
    assertAddSubMatches(c, a)

    // Both strided sub-grids of larger buffers.
    val s1 = randStrided(rng, 3, 4, rowStride = 9, colStride = 2, offset = 1)
    val s2 = randStrided(rng, 3, 4, rowStride = 7, colStride = 1, offset = 3)
    assertAddSubMatches(s1, s2)
    assertAddSubMatches(a, s1)
  }

  // ---- K8: scaled norm2 (dnrm2) --------------------------------------------

  test("norm2 does not overflow for very large elements") {
    val huge = DVec.fill(4)(1e155)
    val n = huge.norm2
    assert(n.isFinite, s"norm2 overflowed to $n")
    // sqrt(4) * 1e155 == 2e155.
    assert(math.abs(n - 2e155) <= 1e-13 * 2e155, s"norm2=$n")
  }

  test("norm2 does not underflow to zero for very small elements") {
    val tiny = DVec.fill(4)(1e-170)
    val n = tiny.norm2
    assert(n > 0.0, s"norm2 underflowed to $n")
    assert(math.abs(n - 2e-170) <= 1e-13 * 2e-170, s"norm2=$n")
  }

  test("norm2 matches sqrt(dot) to 1e-15 relative on ordinary values") {
    val rng = new Random(7L)
    var trial = 0
    while trial < 20 do
      val x = DVec.fromSeq(Seq.fill(37)(rng.nextDouble() * 200.0 - 100.0))
      val expected = math.sqrt(x.dot(x))
      val actual = x.norm2
      assert(math.abs(actual - expected) <= 1e-15 * expected, s"norm2=$actual sqrt(dot)=$expected")
      trial += 1
  }

  // ---- K9: dgemm fast paths vs the naive product ---------------------------

  private def naiveGemm(a: DMat, b: DMat): DMat =
    Matrix.tabulate(a.rows, b.cols) { (i, j) =>
      var acc = 0.0
      var k = 0
      while k < a.cols do
        acc += a(i, k) * b(k, j)
        k += 1
      acc
    }

  private def assertGemmClose(a: DMat, b: DMat, rel: Double): Unit =
    val actual = a * b
    val expected = naiveGemm(a, b)
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        val e = expected(i, j)
        val tol = rel * math.max(1.0, math.abs(e))
        assert(math.abs(actual(i, j) - e) <= tol, s"gemm mismatch at ($i, $j): ${actual(i, j)} != $e")
        j += 1
      i += 1

  private def randMat(rng: Random, rows: Int, cols: Int): DMat =
    Matrix.dense(rows, cols, Seq.fill(rows * cols)(rng.nextDouble() * 2.0 - 1.0))

  test("dgemm blocked and row-major fast paths match the naive product") {
    val rng = new Random(11L)
    // Large square -> blocked path (>= 64^3 elements).
    assertGemmClose(randMat(rng, 96, 96), randMat(rng, 96, 96), 1e-10)
    // Large rectangular -> blocked path.
    assertGemmClose(randMat(rng, 128, 96), randMat(rng, 96, 112), 1e-10)
    // Small -> unblocked i-k-j row-major path.
    assertGemmClose(randMat(rng, 32, 40), randMat(rng, 40, 24), 1e-10)
  }

  test("dgemm strided fallback (transposed operand) matches the naive product") {
    val rng = new Random(13L)
    // A.t has non-unit column stride, forcing the general strided kernel.
    val a = randMat(rng, 96, 100).t // 100 x 96
    val b = randMat(rng, 96, 70)
    assertGemmClose(a, b, 1e-10)
  }

  // ---- K2: hoisted factorization substitutions still solve correctly -------

  private def diagDominant(rng: Random, n: Int): DMat =
    Matrix.tabulate(n, n) { (i, j) =>
      if i == j then n.toDouble + 1.0 + rng.nextDouble()
      else rng.nextDouble() * 2.0 - 1.0
    }

  test("LU solve (hoisted substitution) recovers the right-hand side") {
    val rng = new Random(17L)
    val n = 9
    val A = diagDominant(rng, n)
    val truth = DVec.fromSeq(Seq.fill(n)(rng.nextDouble() * 4.0 - 2.0))
    val b = A * truth
    val x = A.solve(b).orThrow
    var i = 0
    while i < n do
      assert(math.abs(x(i) - truth(i)) < 1e-9, s"index $i: ${x(i)} != ${truth(i)}")
      i += 1
  }

  test("Cholesky solve (hoisted substitution) recovers the right-hand side") {
    val rng = new Random(19L)
    val n = 8
    // SPD: A = BᵀB + n*I.
    val B = Matrix.dense(n, n, Seq.fill(n * n)(rng.nextDouble() * 2.0 - 1.0))
    val shift = Matrix.tabulate(n, n)((i, j) => if i == j then n.toDouble else 0.0)
    val A = (B.t * B) + shift
    val truth = DVec.fromSeq(Seq.fill(n)(rng.nextDouble() * 4.0 - 2.0))
    val b = A * truth
    val x = A.cholesky.orThrow.solve(b).orThrow
    var i = 0
    while i < n do
      assert(math.abs(x(i) - truth(i)) < 1e-8, s"index $i: ${x(i)} != ${truth(i)}")
      i += 1
  }
