package gale.laws

import gale.linalg.*
import munit.Assertions

/** Reusable algebraic laws for [[gale.linalg.DMat]], expressed against the public
  * API. See [[VecLaws]] for the throwing/`Unit` convention that lets these be
  * used from both munit tests and ScalaCheck property bodies.
  */
object MatrixLaws extends Assertions:
  /** Two matrices are elementwise equal within an absolute `tolerance`. */
  def assertClose(actual: DMat, expected: DMat, tolerance: Double = 1e-10): Unit =
    assertEquals(actual.rows, expected.rows, "row count mismatch")
    assertEquals(actual.cols, expected.cols, "col count mismatch")
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        assert(
          math.abs(actual(i, j) - expected(i, j)) <= tolerance,
          s"matrix mismatch at ($i, $j): ${actual(i, j)} != ${expected(i, j)} (tol $tolerance)"
        )
        j += 1
      i += 1

  /** Two matrices are elementwise equal within a `rel`-relative tolerance. */
  def assertCloseRel(actual: DMat, expected: DMat, rel: Double = 1e-9): Unit =
    assertEquals(actual.rows, expected.rows, "row count mismatch")
    assertEquals(actual.cols, expected.cols, "col count mismatch")
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        val e = expected(i, j)
        val tol = rel * math.max(1.0, math.abs(e))
        assert(math.abs(actual(i, j) - e) <= tol, s"matrix mismatch at ($i, $j): ${actual(i, j)} != $e (tol $tol)")
        j += 1
      i += 1

  /** Two matrices are bit-for-bit elementwise equal (no arithmetic performed). */
  def assertExact(actual: DMat, expected: DMat): Unit =
    assertEquals(actual.rows, expected.rows, "row count mismatch")
    assertEquals(actual.cols, expected.cols, "col count mismatch")
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        assert(actual(i, j) == expected(i, j), s"matrix mismatch at ($i, $j): ${actual(i, j)} != ${expected(i, j)}")
        j += 1
      i += 1

  /** Transpose is an involution: `A.t.t == A` elementwise. */
  def transposeIsInvolution(a: DMat): Unit =
    assertExact(a.t.t, a)

  /** Matrix-vector product is linear in the vector:
    * `A (x + y) == A x + A y`.
    */
  def matVecIsLinear(a: DMat, x: DVec, y: DVec): Unit =
    VecLaws.assertCloseRel(a * (x + y), (a * x) + (a * y))

  /** Matrix multiplication is associative: `A (B C) == (A B) C`.
    *
    * The tolerance scales with the '''intermediate''' magnitude
    * `maxAbs(A)·maxAbs(B)·maxAbs(C)·(inner dims)`, not with the result entry: a
    * result entry that nearly cancels to zero out of large intermediates carries
    * roundoff proportional to those intermediates, and the two association
    * orders round differently. A result-relative tolerance (`rel·max(1, |e|)`)
    * intermittently fails such draws — the source of a rare ScalaCheck flake
    * before this scaling.
    */
  def multiplicationAssociates(a: DMat, b: DMat, c: DMat): Unit =
    val left = a * (b * c)
    val right = (a * b) * c
    val scale =
      maxAbs(a) * maxAbs(b) * maxAbs(c) *
        math.max(1, a.cols).toDouble * math.max(1, b.cols).toDouble
    val tolerance = 1e-12 * math.max(1.0, scale)
    assertClose(left, right, tolerance)

  private def maxAbs(m: DMat): Double =
    var out = 0.0
    var i = 0
    while i < m.rows do
      var j = 0
      while j < m.cols do
        val v = math.abs(m(i, j))
        if v > out then out = v
        j += 1
      i += 1
    out
