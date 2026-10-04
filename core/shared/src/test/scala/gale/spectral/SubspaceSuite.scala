package gale.spectral

import gale.linalg.*

class SubspaceSuite extends munit.FunSuite:
  private def close(actual: DVec, expected: DVec): Unit =
    assert((actual - expected).norm2 < 1e-12)

  test("row space projects onto equal first coordinates without a square projector"):
    val a = Matrix.dense(2, 3)(1.0, 1.0, 0.0, 0.0, 0.0, 2.0)
    val space = Subspace.rowSpace(a).orThrow
    assertEquals(space.dimension, 3)
    assertEquals(space.rank, 2)
    val b = Vec(3.0, 1.0, 7.0)
    val projected = space.project(b).orThrow
    close(projected, Vec(2.0, 2.0, 7.0))
    close(space.project(projected).orThrow, projected)
    close(space.residual(b).orThrow, Vec(1.0, -1.0, 0.0))
    assertEqualsDouble(space.distance(b).orThrow, math.sqrt(2.0), 1e-12)
    close(space.basis.t * space.residual(b).orThrow, DVec.zeros(space.rank))

  test("column space residualizes each column of a strided RHS"):
    val a = Matrix.dense(3, 2)(1.0, 0.0, 1.0, 0.0, 0.0, 2.0)
    val space = Subspace.columnSpace(a).orThrow
    val b = Matrix.dense(2, 3)(3.0, 1.0, 7.0, -1.0, 5.0, 0.0).t
    val projected = space.project(b).orThrow
    val residual = space.residual(b).orThrow
    for j <- 0 until b.cols do
      close(projected.col(j), space.project(b.col(j)).orThrow)
      close(projected.col(j) + residual.col(j), b.col(j))
    close(projected.col(0), Vec(2.0, 2.0, 7.0))
    assertEquals(space.project(DMat.zeros(3, 0)).orThrow.cols, 0)

  test("retained factors expose singular metadata and share cutoff across spaces"):
    val a = Matrix.dense(2, 2)(4.0, 0.0, 0.0, 1.0)
    val factor = a.truncatedSvd(SvdCutoff.Relative(0.25)).orThrow
    for space <- Seq(factor.rowSpace, factor.columnSpace) do
      assertEquals(space.rank, 1)
      assertEquals(space.sigmaMax, 4.0)
      assertEquals(space.cutoff, 1.0)
      assertEquals(space.singularValues(0), 4.0)
      close(space.project(Vec(2.0, 3.0)).orThrow, Vec(2.0, 0.0))

  test("zero and full spaces, dimension errors, and non-finite inputs"):
    val zero = Subspace.rowSpace(DMat.zeros(2, 3)).orThrow
    close(zero.project(Vec(1.0, 2.0, 3.0)).orThrow, DVec.zeros(3))
    close(zero.residual(Vec(1.0, 2.0, 3.0)).orThrow, Vec(1.0, 2.0, 3.0))
    val full = Subspace.columnSpace(Matrix.eye(3)).orThrow
    close(full.project(Vec(1.0, 2.0, 3.0)).orThrow, Vec(1.0, 2.0, 3.0))
    assert(full.project(DVec.zeros(2)).isLeft)
    assert(full.residual(DMat.zeros(2, 3)).isLeft)
    assert(full.distance(Vec(1.0, Double.NaN, 0.0)).isLeft)
