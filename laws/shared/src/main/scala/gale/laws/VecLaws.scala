package gale.laws

import gale.linalg.*
import munit.Assertions

/** Reusable algebraic laws for [[gale.linalg.DVec]], expressed entirely against
  * the public API so they hold identically on every platform.
  *
  * Each method throws a munit failure (via the inherited [[munit.Assertions]])
  * when the law is violated and returns `Unit` otherwise, so it drops straight
  * into a plain munit `test` body or a ScalaCheck property body — munit-scalacheck
  * treats a non-throwing `Unit` as a passing property.
  */
object VecLaws extends Assertions:
  /** Two vectors are elementwise equal within an absolute `tolerance`. */
  def assertClose(actual: DVec, expected: DVec, tolerance: Double = 1e-10): Unit =
    assertEquals(actual.length, expected.length, "vector length mismatch")
    var i = 0
    while i < actual.length do
      assert(
        math.abs(actual(i) - expected(i)) <= tolerance,
        s"vector mismatch at $i: ${actual(i)} != ${expected(i)} (tol $tolerance)"
      )
      i += 1

  /** Two vectors are elementwise equal within a `rel`-relative tolerance,
    * floored so near-zero expected values still admit a small absolute slack.
    */
  def assertCloseRel(actual: DVec, expected: DVec, rel: Double = 1e-9): Unit =
    assertEquals(actual.length, expected.length, "vector length mismatch")
    var i = 0
    while i < actual.length do
      val e = expected(i)
      val tol = rel * math.max(1.0, math.abs(e))
      assert(math.abs(actual(i) - e) <= tol, s"vector mismatch at $i: ${actual(i)} != $e (tol $tol)")
      i += 1

  /** Vector addition commutes: `x + y == y + x`. */
  def additionCommutes(x: DVec, y: DVec): Unit =
    assertClose(x + y, y + x)

  /** Scalar multiplication distributes over vector addition:
    * `alpha * (x + y) == alpha * x + alpha * y`.
    */
  def scalarDistributesOverAddition(alpha: Double, x: DVec, y: DVec): Unit =
    assertCloseRel(alpha * (x + y), (alpha * x) + (alpha * y))
