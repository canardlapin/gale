package gale.linalg

/** Plain regression suite for dense shape/edge-case behaviour that stays in core.
  *
  * The algebraic laws (commutativity, distributivity, linearity, transpose
  * involution, associativity, strided/transposed-view equivalence) moved to the
  * gale-laws module (`gale.laws.DenseLawSuite`), where they are expressed through
  * the reusable [[gale.laws.VecLaws]] / [[gale.laws.MatrixLaws]] bundles. What
  * remains here are the dimension-guard and zero-dimension regressions, which are
  * edge-case checks rather than algebraic laws.
  */
class DenseLawSuite extends munit.FunSuite:
  test("shape checks reject invalid matvec and vector dot inputs") {
    intercept[LinAlgError.VectorLengthMismatch] {
      Matrix.eye(2) * Vec(1.0, 2.0, 3.0)
    }
    intercept[LinAlgError.VectorLengthMismatch] {
      Vec(1.0, 2.0).dot(Vec(1.0))
    }
  }

  test("zero-dimension shapes construct, multiply, and transpose without throwing") {
    val a0n = Matrix.zeros(0, 3)
    val an0 = Matrix.zeros(3, 0)
    val empty = DVec.zeros(0)

    assertEquals((a0n * Vec(1.0, 2.0, 3.0)).length, 0)
    assertEquals((an0 * empty).length, 3)
    assertEquals(a0n.t.rows, 3)
    assertEquals(a0n.t.cols, 0)
    assertEquals((an0.t * Vec(1.0, 2.0, 3.0)).length, 0)

    val z = Matrix.zeros(0, 0)
    assertEquals((z * empty).length, 0)
    assertEquals(z.t.rows, 0)
    assertEquals(z.t.cols, 0)
  }
