package gale.laws

import gale.linalg.*

/** Reversal, positive scaling and superposition do not depend on a second
  * Cholesky implementation. Dense products supply independent residuals.
  */
class BandedCholeskyLawSuite extends munit.FunSuite:
  private def packed(a: DMat, b: Int): DMat =
    Matrix.tabulate(a.rows, b + 1)((i, d) => if d <= i then a(i, i - d) else 0.0)

  private def close(x: DVec, y: DVec): Unit =
    for i <- 0 until x.length do assertEqualsDouble(x(i), y(i), 3e-11)

  for n <- List(3, 9, 27); b <- List(0, 1, 2) do
    test(s"banded SPD solve laws n=$n bandwidth=$b"):
      val a = Matrix.tabulate(n, n): (i, j) =>
        if i == j then 2.0 * b + 1.0 + i / n.toDouble
        else if math.abs(i - j) <= b then math.sin(i + j * 1.0)
        else 0.0
      val f = BandedCholesky.factorLower(packed(a, b)).orThrow
      val u = Vec.tabulate(n)(i => math.sin(i * 2.0))
      val v = Vec.tabulate(n)(i => math.cos(i / 3.0))
      val x = f.solve(u).orThrow
      close(a * x, u)
      close(f.solve(u * 2.0 + v * -3.0).orThrow, x * 2.0 + f.solve(v).orThrow * -3.0)
      close(f.solveTranspose(u).orThrow, x)
      val reversed = Matrix.tabulate(n, n)((i, j) => a(n - 1 - i, n - 1 - j))
      val reversedFactor = BandedCholesky.factorLower(packed(reversed, b)).orThrow
      close(reversedFactor.solve(Vec.tabulate(n)(i => u(n - 1 - i))).orThrow, Vec.tabulate(n)(i => x(n - 1 - i)))
      assertEqualsDouble(reversedFactor.logDet, f.logDet, 3e-12)
      for scale <- List(0.001, 7.0, 1000.0) do
        val scaled = BandedCholesky.factorLower(packed(a * scale, b)).orThrow
        close(scaled.solve(u * scale).orThrow, x)
        assertEqualsDouble(scaled.logDet, f.logDet + n * math.log(scale), 3e-11)
      val inverse = f.solve(Matrix.eye(n)).orThrow
      val inverseNorm = (0 until n).map(j => (0 until n).map(i => math.abs(inverse(i, j))).sum).max
      assert(f.conditioning.inverseOneNormUpperBound >= inverseNorm * (1.0 - 1e-13))
      assert(f.conditioning.conditionNumberUpperBound >= 1.0)
