package userland

import gale.linalg.*

/** Consumer-vantage checks for the mutable/immutable ownership boundary. */
class MutableViewOpacitySuite extends munit.FunSuite:
  private val mutable = MutableVec.from(Vec(1.0, 2.0))

  test("a mutable vector has no public aliasing conversion to DVec"):
    assert(compileErrors("val leaked: DVec = mutable.asVec").nonEmpty)
    assertEquals(compileErrors("val snapshot: DVec = mutable.toVec"), "")

  test("the public conversion is an independent snapshot"):
    val snapshot = mutable.toVec
    mutable(0) = 99.0
    assertEquals(snapshot.toSeq, Seq(1.0, 2.0))
