package gale.sized

import gale.linalg.DVec
import gale.linalg.Matrix
import gale.linalg.Vec

/** Tests for the optional compile-time-dimension layer: dimension-correct ops
  * compile and compute the same as their `DVec`/`DMat` equivalents; dimension-
  * mismatched ops '''do not compile''' (the layer's whole point, checked with
  * munit's `compileErrors`); runtime adoption goes through `Either`; lowering is
  * zero-copy; and the tiny fixed-size det/inverse kernels are exact.
  */
class SizedSuite extends munit.FunSuite:

  private def assertVec(actual: DVec, expected: Seq[Double])(using munit.Location): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < expected.length do
      assertEqualsDouble(actual(i), expected(i), 1e-12)
      i += 1

  // --- compile-time-correct ops compute right --------------------------------

  test("vector add / sub / scale / dot match DVec") {
    val a = SVec.vec3(1.0, 2.0, 3.0)
    val b = SVec.vec3(4.0, 5.0, 6.0)
    assertVec((a + b).toDVec, Seq(5.0, 7.0, 9.0))
    assertVec((b - a).toDVec, Seq(3.0, 3.0, 3.0))
    assertVec((a * 2.0).toDVec, Seq(2.0, 4.0, 6.0))
    assertEqualsDouble(a.dot(b), 32.0, 1e-12)
    assertEqualsDouble(a.norm2, math.sqrt(14.0), 1e-12)
  }

  test("matvec: SMat[2,3] * SVec[3] = SVec[2], matches DMat * DVec") {
    val a: SMat[2, 3] = SMat.sized[2, 3](1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val x: SVec[3] = SVec.vec3(1.0, 0.0, -1.0)
    val y: SVec[2] = a * x
    assertVec(y.toDVec, Seq(-2.0, -2.0))
    val dref = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0) * Vec(1.0, 0.0, -1.0)
    assertVec(y.toDVec, (0 until 2).map(dref(_)))
  }

  test("matmul: SMat[2,3] * SMat[3,2] = SMat[2,2], matches DMat * DMat; transpose") {
    val a: SMat[2, 3] = SMat.sized[2, 3](1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val b: SMat[3, 2] = SMat.sized[3, 2](7.0, 8.0, 9.0, 10.0, 11.0, 12.0)
    val c: SMat[2, 2] = a * b
    val ref = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0) * Matrix.dense(3, 2)(7.0, 8.0, 9.0, 10.0, 11.0, 12.0)
    assertEqualsDouble(c(0, 0), ref(0, 0), 1e-12)
    assertEqualsDouble(c(1, 1), ref(1, 1), 1e-12)
    // Transpose flips the static shape and the values.
    val at: SMat[3, 2] = a.t
    assertEquals(at.rows, 3)
    assertEquals(at.cols, 2)
    assertEqualsDouble(at(0, 1), 4.0, 1e-12) // a(1,0)
  }

  test("matrix add / sub / scale") {
    val a = SMat.sized[2, 2](1.0, 2.0, 3.0, 4.0)
    val b = SMat.sized[2, 2](4.0, 3.0, 2.0, 1.0)
    assertEqualsDouble((a + b)(0, 0), 5.0, 1e-12)
    assertEqualsDouble((a - b)(1, 1), 3.0, 1e-12)
    assertEqualsDouble(a.scale(3.0)(0, 1), 6.0, 1e-12)
  }

  // --- compile-time NEGATIVES (the whole point) ------------------------------

  test("dimension-mismatched ops do NOT compile") {
    val a23: SMat[2, 3] = SMat.sized[2, 3](1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val a24: SMat[2, 4] = SMat.sized[2, 4](1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
    val v2: SVec[2] = SVec.vec2(1.0, 2.0)
    val v3: SVec[3] = SVec.vec3(1.0, 2.0, 3.0)

    // matmul: the right factor's rows must equal the left's cols (3), but a24 is 2×4.
    assert(compileErrors("a23 * a24").nonEmpty, "mismatched matmul must not compile")
    // vector add: lengths must match.
    assert(compileErrors("v2 + v3").nonEmpty, "mismatched add must not compile")
    // matvec: the vector length must equal the matrix's cols (3), but v2 is length 2.
    assert(compileErrors("a23 * v2").nonEmpty, "mismatched matvec must not compile")
    // matrix add: shapes must match.
    assert(compileErrors("a23 + a24").nonEmpty, "mismatched matrix add must not compile")

    // Positive controls: the dimension-correct forms DO compile.
    assertEquals(compileErrors("val ok: SVec[2] = a23 * v3"), "")
    assertEquals(compileErrors("val ok: SMat[2, 4] = a23 * SMat.sized[3, 4](1.0)"), "")
  }

  // --- runtime adoption via Either -------------------------------------------

  test("fromDVec / fromDMat: Right on a shape match, Left on a mismatch") {
    assert(SVec.fromDVec[3](Vec(1.0, 2.0, 3.0)).isRight)
    assert(SVec.fromDVec[3](Vec(1.0, 2.0)).isLeft)
    assert(SMat.fromDMat[2, 2](Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)).isRight)
    assert(SMat.fromDMat[2, 2](Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)).isLeft)
  }

  test("sized throws on an arity mismatch (a programming error)") {
    intercept[IllegalArgumentException](SVec.sized[3](1.0, 2.0))
    intercept[IllegalArgumentException](SMat.sized[2, 2](1.0, 2.0, 3.0))
  }

  // --- zero-copy lowering -----------------------------------------------------

  test("toDVec / toDMat are zero-copy — the same instance that was adopted") {
    val d = Vec(1.0, 2.0, 3.0)
    val s = SVec.fromDVec[3](d).toOption.get
    assert(s.toDVec eq d, "SVec lowering must return the adopted DVec instance")
    val m = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
    val sm = SMat.fromDMat[2, 2](m).toOption.get
    assert(sm.toDMat eq m, "SMat lowering must return the adopted DMat instance")
  }

  // --- tiny fixed-size kernels: det / inverse --------------------------------

  test("2×2 det and inverse (m · m⁻¹ = I); singular ⇒ Left") {
    val m = SMat.mat2(4.0, 7.0, 2.0, 6.0)
    assertEqualsDouble(m.det, 10.0, 1e-12)
    val inv = m.inverse.toOption.get
    val id = m * inv
    assertEqualsDouble(id(0, 0), 1.0, 1e-12)
    assertEqualsDouble(id(0, 1), 0.0, 1e-12)
    assertEqualsDouble(id(1, 0), 0.0, 1e-12)
    assertEqualsDouble(id(1, 1), 1.0, 1e-12)
    assert(SMat.mat2(1.0, 2.0, 2.0, 4.0).inverse.isLeft) // det 0
  }

  test("3×3 det and inverse (m · m⁻¹ = I); singular ⇒ Left") {
    // Chosen so the full determinant (27) differs from the leading 2×2 minor
    // (1·5 − 2·4 = −3): a kernel that dropped the cofactor terms would still
    // agree on a degenerate matrix, so this discriminates the full expansion.
    val m = SMat.mat3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.0)
    assertEqualsDouble(m.det, 27.0, 1e-9)
    val inv = m.inverse.toOption.get
    val id = m * inv
    var i = 0
    while i < 3 do
      var j = 0
      while j < 3 do
        assertEqualsDouble(id(i, j), if i == j then 1.0 else 0.0, 1e-9)
        j += 1
      i += 1
    assert(SMat.mat3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0).inverse.isLeft) // det 0
  }
