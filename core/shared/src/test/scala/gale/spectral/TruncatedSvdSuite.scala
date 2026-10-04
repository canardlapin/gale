package gale.spectral

import gale.linalg.*

class TruncatedSvdSuite extends munit.FunSuite:
  private def close(actual: DVec, expected: DVec, tolerance: Double = 1e-12): Unit =
    assertEquals(actual.length, expected.length)
    assert((actual - expected).norm2 <= tolerance, s"$actual != $expected")

  test("relative and absolute cutoffs exclude equality and share the kept rank"):
    val a = Matrix.dense(3, 3)(4.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.25)
    for cutoff <- Seq(SvdCutoff.Relative(0.25), SvdCutoff.Absolute(1.0)) do
      val factor = a.truncatedSvd(cutoff).orThrow
      assertEquals(factor.rank, 1)
      assertEquals(factor.cutoff, 1.0)
      val fit = factor.solve(Vec(8.0, 3.0, 5.0)).orThrow
      close(fit.solution, Vec(2.0, 0.0, 0.0))
      close(fit.residual, Vec(0.0, 3.0, 5.0))
      close(a.pinv(cutoff).orThrow * Vec(8.0, 3.0, 5.0), fit.solution)

  test("caller policy does not reuse SVD.rank and default preserves pinv"):
    val a = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1e-12)
    assertEquals(a.svd.orThrow.rank, 1)
    assertEquals(a.truncatedSvd().orThrow.rank, 2)
    assertEquals(a.truncatedSvd(SvdCutoff.Relative(1e-8)).orThrow.rank, 1)
    assertEquals(a.pinv.orThrow(1, 1), 1e12)
    assertEquals(a.pinv(SvdCutoff.Default).orThrow(1, 1), 1e12)

  test("rotated truncation has the analytical inverse and projected residual"):
    // Eigenvectors (1,1)/sqrt(2), (1,-1)/sqrt(2); singular values 4 and 0.1.
    val a = Matrix.dense(2, 2)(2.05, 1.95, 1.95, 2.05)
    val cutoff = SvdCutoff.Relative(0.1)
    val factor = a.truncatedSvd(cutoff).orThrow
    val expectedInverse = Matrix.dense(2, 2)(0.125, 0.125, 0.125, 0.125)
    for j <- 0 until 2 do close(factor.pinv.col(j), expectedInverse.col(j))
    val fit = factor.solve(Vec(2.0, 0.0)).orThrow
    close(fit.solution, Vec(0.25, 0.25))
    close(fit.residual, Vec(1.0, -1.0))
    val retained = Matrix.dense(2, 2)(2.0, 2.0, 2.0, 2.0)
    for j <- 0 until 2 do
      close((retained * factor.pinv * retained).col(j), retained.col(j))
      close((factor.pinv * retained * factor.pinv).col(j), factor.pinv.col(j))

  test("rank-deficient wide solve has the analytical minimum norm"):
    val a = Matrix.dense(2, 3)(1.0, 1.0, 0.0, 0.0, 0.0, 2.0)
    val fit = a.minimumNormLeastSquares(Vec(2.0, 4.0)).orThrow
    close(fit.solution, Vec(1.0, 1.0, 2.0))
    close(fit.residual, DVec.zeros(2))
    assertEquals(fit.rank, 2)
    val nullDirection = Vec(1.0, -1.0, 0.0)
    assertEqualsDouble(fit.solution.dot(nullDirection), 0.0, 1e-12)
    assert((fit.solution + nullDirection).norm2 > fit.solution.norm2)

  test("reusable matrix solve agrees with vector solves including strided RHS"):
    val a = Matrix.dense(3, 3)(1.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0)
    val b = Matrix.dense(2, 3)(2.0, 4.0, 3.0, -2.0, 6.0, 7.0).t
    val factor = a.truncatedSvd().orThrow
    val result = factor.solve(b).orThrow
    for j <- 0 until b.cols do
      val vector = factor.solve(b.col(j)).orThrow
      close(result.solution.col(j), vector.solution)
      close(result.residual.col(j), vector.residual)
      close(result.residual.col(j), b.col(j) - a * result.solution.col(j))
      close(a.t * result.residual.col(j), DVec.zeros(3))
    assertEquals(factor.solve(DMat.zeros(3, 0)).orThrow.solution.cols, 0)

  test("relative cutoff is invariant under common scaling of matrix and RHS"):
    val a = Matrix.dense(2, 2)(2.0, 0.0, 0.0, 1e-7)
    val b = Vec(4.0, 3.0)
    for scale <- Seq(1e-90, 1.0, 1e90) do
      val fit = (a * scale).minimumNormLeastSquares(b * scale, SvdCutoff.Relative(1e-6)).orThrow
      assertEquals(fit.rank, 1)
      close(fit.solution, Vec(2.0, 0.0))

  test("zero rank and typed validation"):
    val zero = DMat.zeros(2, 3).truncatedSvd().orThrow
    val fit = zero.solve(Vec(1.0, 2.0)).orThrow
    assertEquals(fit.rank, 0)
    close(fit.solution, DVec.zeros(3))
    close(fit.residual, Vec(1.0, 2.0))
    assert(zero.solve(DVec.zeros(1)).isLeft)
    assert(zero.solve(DMat.zeros(1, 2)).isLeft)
    assert(zero.solve(Vec(Double.NaN, 0.0)).isLeft)
    assert(zero.solve(Matrix.dense(2, 1)(0.0, Double.PositiveInfinity)).isLeft)
    for value <- Seq(-1.0, Double.NaN, Double.PositiveInfinity) do
      assert(Matrix.eye(2).truncatedSvd(SvdCutoff.Relative(value)).isLeft)
      assert(Matrix.eye(2).truncatedSvd(SvdCutoff.Absolute(value)).isLeft)
    assert(Matrix.dense(1, 1)(Double.NaN).truncatedSvd().isLeft)
    assert(DMat.zeros(0, 2).truncatedSvd().isLeft)
    assertEquals(Matrix.eye(2).truncatedSvd(SvdCutoff.Relative(0.0)).orThrow.rank, 2)
