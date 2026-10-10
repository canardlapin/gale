package gale.spectral

import gale.golden.GoldenTolerance.{Eps, svdFactorRatios}
import gale.linalg.*

class FullSvdScaleSuite extends munit.FunSuite:
  private def power(k: Int): Double =
    if k >= -1022 then java.lang.Double.longBitsToDouble((1023L + k) << 52)
    else power(k + 600) * power(-600)

  // Orthogonal columns of length 5 and 2, then a right rotation.
  private val base = Array(Array(1.8, 2.4), Array(-1.6, 1.2), Array(2.4, 3.2))

  for exponent <- Seq(-1030, -700, -550, -257, -256, -255, 0, 255, 256, 257, 550, 700); wide <- Seq(false, true) do
    test(s"full SVD scale 2^$exponent wide=$wide, vectors and values-only") {
      val scale = power(exponent)
      val a = Matrix.tabulate(if wide then 2 else 3, if wide then 3 else 2)((i, j) => (if wide then base(j)(i) else base(i)(j)) * scale)
      val before = Array.tabulate(a.rows, a.cols)((i, j) => java.lang.Double.doubleToRawLongBits(a(i, j)))
      val d = a.svd.orThrow
      val values = Svds.svd(a, SingularSelection.All, EigenVectors.ValuesOnly).orThrow
      assert(d.diagnostics.allConverged)
      // Restored subnormal singular values have unavoidable quantization error.
      val quantum = Double.MinPositiveValue / scale
      for (expected, i) <- Seq(5.0, 2.0).zipWithIndex do
        assert(math.abs(d.singularValues(i) / scale - expected) <= 128.0 * Eps + 8.0 * quantum)
        assert(math.abs(values.singularValues(i) / scale - expected) <= 128.0 * Eps + 8.0 * quantum)
      val bound = 32.0 + 8.0 * quantum / (3.0 * Eps)
      for (name, units) <- svdFactorRatios(a.apply, a.rows, a.cols, d.singularValues.toArray, d.u.apply, d.vt.apply) do
        assert(units <= bound, s"$name: $units units, bound=$bound")
      for i <- 0 until a.rows; j <- 0 until a.cols do
        assertEquals(java.lang.Double.doubleToRawLongBits(a(i, j)), before(i)(j), "input is caller-owned")
    }

  test("zero/rank-deficient inputs retain complete finite factors") {
    for scale <- Seq(0.0, power(-700), power(700)) do
      val a = Matrix.tabulate(4, 3)((i, j) => if j == 0 then (i + 1) * scale else 0.0)
      val d = a.svd.orThrow
      assertEquals(d.size, 3)
      assert(d.singularValues.toArray.forall(_.isFinite))
      for (name, units) <- svdFactorRatios(a.apply, a.rows, a.cols, d.singularValues.toArray, d.u.apply, d.vt.apply) do
        assert(units <= 32.0, s"$name: $units units")
  }

  test("an unrepresentable spectrum is a typed failure, never successful infinities") {
    val a = Matrix.tabulate(2, 2)((_, _) => Double.MaxValue)
    for mode <- Seq(EigenVectors.Right, EigenVectors.ValuesOnly) do
      Svds.svd(a, SingularSelection.All, mode) match
        case Left(_: LinAlgError.DidNotConverge) => ()
        case other => fail(s"expected typed non-finite-arithmetic failure, got $other")
  }
