// NOTE: package `userland` (not gale.*) is intentional — see the JVM twin. This
// suite runs as a downstream user would, so only the public JS interop
// extensions are reachable; private[gale] typed-array storage stays hidden.
package userland

import gale.linalg.*
import scala.scalajs.js.typedarray.Float64Array

class InteropSuite extends munit.FunSuite:
  test("Vec.fromFloat64ArrayUnsafe adopts the typed array without copying") {
    val raw = new Float64Array(3)
    raw(0) = 1.0
    raw(1) = 2.0
    raw(2) = 3.0
    val v = Vec.fromFloat64ArrayUnsafe(raw)
    assertEquals(v.toSeq, Seq(1.0, 2.0, 3.0))
    raw(0) = 9.0
    assertEquals(v(0), 9.0)
  }

  test("DVec.toFloat64Array returns an independent copy") {
    val v = Vec(1.0, 2.0, 3.0)
    val out = v.toFloat64Array
    out(0) = 9.0
    assertEquals(v(0), 1.0)
    assertEquals(out(0), 9.0)
  }

  test("DVec.toFloat64Array materialises a strided slice in logical order") {
    val out = Vec(1.0, 2.0, 3.0, 4.0).slice(1, 4).toFloat64Array
    assertEquals(out.length, 3)
    assertEquals(out(0), 2.0)
    assertEquals(out(2), 4.0)
  }

  test("Matrix.fromFloat64ArrayUnsafe adopts row-major storage") {
    val raw = new Float64Array(4)
    raw(0) = 1.0
    raw(1) = 2.0
    raw(2) = 3.0
    raw(3) = 4.0
    val m = Matrix.fromFloat64ArrayUnsafe(2, 2, raw)
    assertEquals(m(1, 0), 3.0)
    raw(3) = 8.0
    assertEquals(m(1, 1), 8.0)
  }

  test("DMat.toFloat64Array returns a row-major copy that materialises transposes") {
    val m = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val out = m.toFloat64Array
    out(0) = 9.0
    assertEquals(m(0, 0), 1.0)
    val transposed = m.t.toFloat64Array
    assertEquals(transposed.length, 6)
    assertEquals(transposed(1), 4.0)
  }
