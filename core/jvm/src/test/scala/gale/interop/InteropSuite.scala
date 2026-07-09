// NOTE: package `userland` (not gale.*) is intentional — this suite must run as
// a downstream user would, so private[gale] storage is genuinely out of scope
// and only the public interop extensions are reachable. The file lives under a
// gale/ directory only because the sandbox forbade relocating it; package, not
// path, is what the compiler and this test care about.
package userland

import gale.linalg.*

class InteropSuite extends munit.FunSuite:
  test("Vec.fromArrayUnsafe adopts the array without copying") {
    val raw = Array(1.0, 2.0, 3.0)
    val v = Vec.fromArrayUnsafe(raw)
    assertEquals(v.toSeq, Seq(1.0, 2.0, 3.0))
    // Adoption contract: mutating the source array changes the vector.
    raw(0) = 9.0
    assertEquals(v(0), 9.0)
  }

  test("DVec.toArray returns an independent copy") {
    val v = Vec(1.0, 2.0, 3.0)
    val out = v.toArray
    out(0) = 9.0
    assertEquals(v(0), 1.0)
    assertEquals(out.toSeq, Seq(9.0, 2.0, 3.0))
  }

  test("DVec.toArray materialises a strided slice in logical order") {
    val v = Vec(1.0, 2.0, 3.0, 4.0).slice(1, 4)
    assertEquals(v.toArray.toSeq, Seq(2.0, 3.0, 4.0))
  }

  test("Matrix.fromArrayUnsafe adopts row-major storage") {
    val raw = Array(1.0, 2.0, 3.0, 4.0)
    val m = Matrix.fromArrayUnsafe(2, 2, raw)
    assertEquals(m(0, 1), 2.0)
    assertEquals(m(1, 0), 3.0)
    raw(3) = 8.0
    assertEquals(m(1, 1), 8.0)
  }

  test("DMat.toArrayRowMajor returns an independent copy and materialises transposes") {
    val m = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val out = m.toArrayRowMajor
    out(0) = 9.0
    assertEquals(m(0, 0), 1.0)
    assertEquals(m.t.toArrayRowMajor.toSeq, Seq(1.0, 4.0, 2.0, 5.0, 3.0, 6.0))
  }
