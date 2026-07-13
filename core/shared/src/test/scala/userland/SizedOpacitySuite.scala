package userland

import gale.linalg.*
import gale.sized.*

/** Consumer-vantage guarantees for the optional sized layer. This suite lives
  * '''outside''' `gale.sized` on purpose: an opaque type is only opaque beyond its
  * defining scope, so the anti-leak property (an `SVec`/`SMat` cannot be used as
  * the unsized `DVec`/`DMat`, and the trusted `unsafe`/`raw` constructors are
  * unreachable) can only be checked from a real downstream package like this one.
  * The in-package `gale.sized.SizedSuite` sees through the opacity and cannot.
  */
class SizedOpacitySuite extends munit.FunSuite:

  private val s: SVec[3] = SVec.vec3(1.0, 2.0, 3.0)
  private val m: SMat[2, 2] = SMat.mat2(1.0, 2.0, 3.0, 4.0)

  test("SVec/SMat do not leak into the unsized DVec/DMat API") {
    // Opaque: a sized value is NOT a subtype of the unsized runtime type, so it
    // cannot silently flow into unsized code without the sanctioned lowering.
    assert(compileErrors("val d: DVec = s").nonEmpty, "SVec must not be usable as a DVec")
    assert(compileErrors("val dm: DMat = m").nonEmpty, "SMat must not be usable as a DMat")

    // Nor can a runtime value be forged into a sized one without the checked adopter.
    assert(compileErrors("val bad: SVec[3] = Vec(1.0, 2.0, 3.0)").nonEmpty)
    assert(compileErrors("val bad: SMat[2, 2] = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)").nonEmpty)
  }

  test("trusted (unchecked) constructors are unreachable from userland") {
    // `unsafe` / `raw` are private[sized]: a consumer must go through the checked
    // `fromDVec`/`fromDMat` (Either) or the arity-checked factories.
    assert(compileErrors("SVec.unsafe[3](Vec(1.0, 2.0, 3.0))").nonEmpty)
    assert(compileErrors("SVec.raw[3](s)").nonEmpty)
    assert(compileErrors("SMat.unsafe[2, 2](Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0))").nonEmpty)
  }

  test("the sanctioned lowering DOES compile and is zero-copy") {
    // toDVec / toDMat are the public doorway back to the unsized world.
    assertEquals(compileErrors("val d: DVec = s.toDVec"), "")
    assertEquals(compileErrors("val dm: DMat = m.toDMat"), "")
    // And it is the identity on the wrapped value (opaque alias, no allocation).
    val d = Vec(4.0, 5.0, 6.0)
    val adopted = SVec.fromDVec[3](d).toOption.get
    assert(adopted.toDVec eq d)
  }
