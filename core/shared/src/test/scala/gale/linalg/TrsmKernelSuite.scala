package gale.linalg

import gale.TestAccess
import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** The blocked multi-right-hand-side triangular kernel `dtrsmLeft` against the
  * single-column `dtrsv` it generalizes: sub-block offsets and leading
  * dimensions, both triangles, implicit unit diagonals, right-hand-side counts
  * on and off the four-column groups, and sizes on both sides of the block.
  */
class TrsmKernelSuite extends munit.FunSuite:
  private val Sentinel = -7.0

  /** A well-conditioned triangle at `aOffset` with leading dimension `lda`. The
    * opposite triangle is NaN, and so is the diagonal when `unit`, so a read of
    * either poisons the solution.
    */
  private def triangle(n: Int, lower: Boolean, unit: Boolean, aOffset: Int, lda: Int, seed: Int): DoubleArray =
    val rng = new scala.util.Random(seed)
    val a = TestAccess.filled(aOffset + n * lda, Double.NaN)
    var i = 0
    while i < n do
      var j = 0
      while j < n do
        val inTriangle = if lower then j < i else j > i
        if inTriangle then a(aOffset + i * lda + j) = (rng.nextDouble() * 2.0 - 1.0) / n
        else if i == j && !unit then a(aOffset + i * lda + j) = (if rng.nextBoolean() then 1.0 else -1.0) * (1.0 + rng.nextDouble())
        j += 1
      i += 1
    a

  /** Right-hand sides at `bOffset` with leading dimension `ldb`; every cell
    * outside the `n × nrhs` view holds the sentinel.
    */
  private def rhs(n: Int, nrhs: Int, bOffset: Int, ldb: Int, seed: Int): DoubleArray =
    val rng = new scala.util.Random(seed)
    val b = TestAccess.filled(bOffset + n * ldb + 3, Sentinel)
    var i = 0
    while i < n do
      var c = 0
      while c < nrhs do
        b(bOffset + i * ldb + c) = rng.nextDouble() * 2.0 - 1.0
        c += 1
      i += 1
    b

  private def copyOf(x: DoubleArray): DoubleArray =
    val out = DoubleArray.alloc(x.length)
    var i = 0
    while i < x.length do
      out(i) = x(i)
      i += 1
    out

  test("dtrsmLeft matches dtrsv per column across sizes, triangles, unit diagonals, and transposed storage") {
    for
      n <- Seq(1, 3, 8, 9, 37, 100)
      nrhs <- Seq(1, 3, 4, 7, 9)
      lower <- Seq(true, false)
      unit <- Seq(true, false)
      transposed <- Seq(false, true)
    do
      val clue = s"n=$n nrhs=$nrhs lower=$lower unit=$unit transposed=$transposed"
      val (aOffset, lda, bOffset, ldb) = (3, n + 2, 5, nrhs + 3)
      // Transposed: T is the opposite triangle stored row-major, read with
      // swapped strides (the Cholesky `Lᵀ` layout).
      val a = triangle(n, lower != transposed, unit, aOffset, lda, seed = 31 * n + nrhs)
      val (rowStep, colStep) = if transposed then (1, lda) else (lda, 1)
      val b = rhs(n, nrhs, bOffset, ldb, seed = 17 * n + nrhs)
      val expected = copyOf(b)
      var c = 0
      while c < nrhs do
        assertEquals(DoubleKernels.dtrsv(n, lower, unit, 0.0, a, aOffset, rowStep, colStep, expected, bOffset + c, ldb), -1)
        c += 1
      assertEquals(DoubleKernels.dtrsmLeft(lower, unit, n, nrhs, a, aOffset, rowStep, colStep, b, bOffset, ldb), -1, clue)
      var i = 0
      while i < b.length do
        val inView = i >= bOffset && (i - bOffset) / ldb < n && (i - bOffset) % ldb < nrhs
        if !inView then assertEquals(b(i), Sentinel, s"$clue: cell $i outside the view was written")
        else
          val (got, want) = (b(i), expected(i))
          assert(got.isFinite, s"$clue: non-finite at $i")
          // One block: dtrsv's exact sequence. Blocked: reassociated, so a
          // relative bound against the solution scale.
          if n <= 8 then assertEquals(got, want, s"$clue at $i")
          else assert(math.abs(got - want) <= 1e-13 * (1.0 + math.abs(want)), s"$clue at $i: $got vs $want")
        i += 1
  }

  test("dtrsmLeft reports the first zero diagonal in substitution order and leaves B untouched") {
    val n = 20
    for lower <- Seq(true, false) do
      val a = triangle(n, lower, unit = false, aOffset = 0, lda = n, seed = 5)
      a(3 * n + 3) = 0.0
      a(15 * n + 15) = -0.0
      val b = rhs(n, 6, bOffset = 0, ldb = 6, seed = 6)
      val before = copyOf(b)
      val info = DoubleKernels.dtrsmLeft(lower, unit = false, n, 6, a, 0, n, 1, b, 0, 6)
      assertEquals(info, if lower then 3 else 15)
      var i = 0
      while i < b.length do
        assertEquals(java.lang.Double.doubleToRawLongBits(b(i)), java.lang.Double.doubleToRawLongBits(before(i)))
        i += 1
      // An implicit unit diagonal never inspects the stored zeros.
      assertEquals(DoubleKernels.dtrsmLeft(lower, unit = true, n, 6, a, 0, n, 1, b, 0, 6), -1)
  }

  test("dtrsmLeft with no right-hand sides or an empty system is a no-op") {
    val a = TestAccess.doubleArray(2.0)
    val b = TestAccess.doubleArray(Sentinel)
    assertEquals(DoubleKernels.dtrsmLeft(lower = true, unit = false, 1, 0, a, 0, 1, 1, b, 0, 0), -1)
    assertEquals(DoubleKernels.dtrsmLeft(lower = false, unit = false, 0, 3, a, 0, 1, 1, b, 0, 3), -1)
    assertEquals(b(0), Sentinel)
  }
