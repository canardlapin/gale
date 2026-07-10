package gale.kernel

import gale.linalg.*
import gale.platform.PlatformMath

/** JVM-only path guards for the fused arithmetic used by dense kernels. */
class PlatformFmaSuite extends munit.FunSuite:
  private val delta = java.lang.Math.scalb(1.0, -27)
  private val a = 1.0 + delta
  private val b = 1.0 - delta
  private val fused = java.lang.Math.fma(a, b, -1.0)

  test("platform fma has one-rounding JVM semantics") {
    assertEquals(a * b - 1.0, 0.0)
    assert(fused != 0.0)
    assertEquals(PlatformMath.fma(a, b, -1.0), fused)
  }

  test("contiguous axpy uses the platform fused primitive") {
    val y = MutableVec.from(Vec(-1.0))
    y.axpyInPlace(a, Vec(b))
    assertEquals(y(0), fused)
  }

  test("row-major gemm accumulation uses the platform fused primitive") {
    val left = Matrix.dense(1, 2)(1.0, a)
    val right = Matrix.dense(2, 1)(-1.0, b)
    assertEquals((left * right)(0, 0), fused)
  }
