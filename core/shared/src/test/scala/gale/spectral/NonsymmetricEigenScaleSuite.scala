package gale.spectral

import gale.linalg.*

class NonsymmetricEigenScaleSuite extends munit.FunSuite:
  private val eps = math.ulp(1.0)
  private def power(k: Int): Double =
    if k >= -1022 then java.lang.Double.longBitsToDouble((1023L + k) << 52)
    else power(k + 600) * power(-600)

  for exponent <- Seq(-1030, -700, -550, -257, -256, 0, 256, 257, 550, 700) do
    test(s"nonsymmetric rotation at 2^$exponent preserves values, pairs and decoded vectors") {
      val scale = power(exponent)
      val a = Matrix.tabulate(2, 2)((i, j) => (if i == j then 1.0 else if i < j then -2.0 else 2.0) * scale)
      val before = Array.tabulate(2, 2)((i, j) => java.lang.Double.doubleToRawLongBits(a(i, j)))
      for mode <- Seq(EigenVectors.Right, EigenVectors.ValuesOnly) do
        val d = Eigen.eigNonsymmetric(a, EigenSelection.All, mode).orThrow
        assertEquals(d.size, 2)
        assert(d.diagnostics.allConverged)
        val positive = d.eigenvalue(0)
        val negative = d.eigenvalue(1)
        assert(positive.im > 0.0)
        assertEquals(negative.re, positive.re)
        assertEquals(negative.im, -positive.im)
        val quantum = Double.MinPositiveValue / scale
        assert(math.abs(positive.re / scale - 1.0) <= 64.0 * eps + 4.0 * quantum)
        assert(math.abs(positive.im / scale - 2.0) <= 64.0 * eps + 4.0 * quantum)
        if mode == EigenVectors.Right then
          for k <- 0 until 2 do
            val lambda = d.eigenvalue(k)
            val (vr, vi) = d.eigenvector(k)
            var norm = 0.0
            var error = 0.0
            for i <- 0 until 2 do
              val ar = (a(i, 0) / scale) * vr(0) + (a(i, 1) / scale) * vr(1)
              val ai = (a(i, 0) / scale) * vi(0) + (a(i, 1) / scale) * vi(1)
              error = math.hypot(error, ar - (lambda.re / scale) * vr(i) + (lambda.im / scale) * vi(i))
              error = math.hypot(error, ai - (lambda.re / scale) * vi(i) - (lambda.im / scale) * vr(i))
              norm = math.hypot(norm, math.hypot(vr(i), vi(i)))
            assert(norm > 0.0 && norm.isFinite)
            assert(error / norm <= 128.0 * eps + 8.0 * quantum)
      for i <- 0 until 2; j <- 0 until 2 do
        assertEquals(java.lang.Double.doubleToRawLongBits(a(i, j)), before(i)(j))
    }

  test("unrepresentable eigenvalues and collapsed imaginary parts fail explicitly") {
    val huge = Matrix.tabulate(2, 2)((_, _) => Double.MaxValue)
    // Characteristic polynomial x^3-16x^2-1 has a pair with |Im| about0.25.
    // Scaling by the smallest subnormal rounds that pair's imaginary parts to0.
    val companion = Array(Array(0.0, 0.0, 1.0), Array(1.0, 0.0, 0.0), Array(0.0, 1.0, 16.0))
    val normal = Matrix.tabulate(3, 3)((i, j) => companion(i)(j))
    val check = Eigen.eigNonsymmetric(normal, EigenSelection.All, EigenVectors.Right).orThrow
    assert((0 until check.size).exists(k => check.eigenvalue(k).im != 0.0))
    val tiny = Matrix.tabulate(3, 3)((i, j) => companion(i)(j) * Double.MinPositiveValue)
    for a <- Seq(huge, tiny); mode <- Seq(EigenVectors.Right, EigenVectors.ValuesOnly) do
      Eigen.eigNonsymmetric(a, EigenSelection.All, mode) match
        case Left(_: LinAlgError.DidNotConverge) => ()
        case other => fail(s"expected a typed unrepresentable-spectrum failure, got $other")
  }
