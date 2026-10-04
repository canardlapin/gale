package gale.linalg

import gale.spectral.SvdCutoff

class RowReductionSuite extends munit.FunSuite:
  private def maxAbs(a: DMat): Double =
    (for i <- 0 until a.rows; j <- 0 until a.cols yield math.abs(a(i, j))).maxOption.getOrElse(0.0)

  test("ordinary RREF has ordered pivots, unit pivot columns and known rows"):
    val a = Matrix.dense(3, 3)(1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 0.0, 1.0, 1.0)
    val reduced = RowReduction.reducedRowEchelon(a, 0.0).orThrow
    assertEquals(reduced.pivotColumns, Vector(0, 1))
    assert(maxAbs(reduced.matrix - Matrix.dense(3, 3)(1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0)) < 1e-12)
    assertEquals(a(0, 1), 2.0)
    val twice = RowReduction.reducedRowEchelon(reduced.matrix, 0.0).orThrow
    assertEquals(twice.pivotColumns, reduced.pivotColumns)
    assertEquals(maxAbs(twice.matrix - reduced.matrix), 0.0)
    val scaled = RowReduction.reducedRowEchelon(Matrix.dense(1, 2)(2e4, 4e4), 1e3).orThrow
    assertEquals(scaled.matrix.valuesRowMajor, IndexedSeq(1.0, 2.0))

  test("spectral cutoff and elimination tolerance deliberately have different contracts"):
    val a = Matrix.dense(2, 2)(1.0, 1.0, 0.0, 0.09)
    val factor = a.truncatedSvd(SvdCutoff.Relative(0.05)).orThrow
    assertEquals(factor.rank, 1)
    assertEquals(RowReduction.reducedRowEchelon(a, factor.cutoff).orThrow.rank, 2)
    val reduced = RowReduction.reducedRowEchelon(a, SvdCutoff.Relative(0.05)).orThrow
    assertEquals(reduced.rank, 1)
    assertEquals(reduced.matrix.rows, 2)
    assertEquals(maxAbs(reduced.matrix.slice(1, 2, 0, 2)), 0.0)

  test("sparse null equations use the null-projector RREF convention"):
    val a = Matrix.dense(1, 3)(1.0, 1.0, 1.0)
    val result = RowReduction.nullSpaceBasis(a).orThrow
    val expected = Matrix.dense(3, 2)(1.0, 0.0, 0.0, 1.0, -1.0, -1.0)
    assertEquals(result.nullity, 2)
    assert(maxAbs(result.basis - expected) < 1e-12)
    assert(maxAbs(a * result.basis) < 1e-12)

  test("wide matrix includes every missing null direction in both forms"):
    val a = Matrix.dense(2, 5)(1.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 3.0)
    for form <- NullSpaceForm.values do
      val result = RowReduction.nullSpaceBasis(a, form = form).orThrow
      assertEquals(result.rank, 2)
      assertEquals(result.nullity, 3)
      assert(maxAbs(a * result.basis) < 1e-12)
      assertEquals(result.basis.rankEstimate, 3)
      if form == NullSpaceForm.Orthonormal then assert(maxAbs(result.basis.t * result.basis - Matrix.eye(3)) < 1e-12)

  test("relative policy scales consistently and numerical null vectors obey the residual bound"):
    val a = Matrix.dense(2, 3)(1.0, 1.0, 0.0, 0.0, 1e-8, 0.0)
    for scale <- Seq(1e-40, 1.0, 1e40) do
      val result = RowReduction.nullSpaceBasis(a * scale, SvdCutoff.Relative(1e-6), NullSpaceForm.Orthonormal).orThrow
      assertEquals(result.rank, 1)
      assertEquals(result.nullity, 2)
      assert(maxAbs((a * scale) * result.basis) <= 1e-6 * scale)
      assert(maxAbs(result.basis.t * result.basis - Matrix.eye(2)) < 1e-12)

  test("seeded integer low-rank plants preserve nullity and complete span"):
    val rng = new scala.util.Random(804L)
    for (m, n, r) <- Seq((3, 7, 2), (8, 5, 3), (4, 4, 1)) do
      val left = Matrix.tabulate(m, r)((i, j) => if i == j then 1.0 else (rng.nextInt(7) - 3).toDouble)
      val right =
        Matrix.tabulate(r, n)((i, j) => if j < r then (if i == j then 1.0 else 0.0) else (rng.nextInt(7) - 3).toDouble)
      val a = left * right
      val factor = a.truncatedSvd(SvdCutoff.Relative(1e-10)).orThrow
      assertEquals(factor.rank, r)
      for form <- NullSpaceForm.values do
        val result = RowReduction.nullSpaceBasis(factor, form).orThrow
        assertEquals(result.nullity, n - r)
        for j <- 0 until result.nullity do
          assert((a * result.basis.col(j)).norm2 <= 1e-11 * math.max(1.0, maxAbs(a)) * result.basis.col(j).norm2)
        assertEquals(result.basis.rankEstimate, n - r)
        if form == NullSpaceForm.Orthonormal then
          assert(maxAbs(result.basis.t * result.basis - Matrix.eye(n - r)) < 1e-12)
          val reconstructed = factor.rowSpace.basis * factor.rowSpace.basis.t + result.basis * result.basis.t
          assert(maxAbs(reconstructed - Matrix.eye(n)) < 1e-12)

  test("zero and full rank null spaces, empty algebraic RREF and typed errors"):
    val zero = RowReduction.nullSpaceBasis(DMat.zeros(2, 3)).orThrow
    assertEquals(zero.nullity, 3)
    assertEquals(maxAbs(zero.basis - Matrix.eye(3)), 0.0)
    for form <- NullSpaceForm.values do
      assertEquals(RowReduction.nullSpaceBasis(Matrix.eye(3), form = form).orThrow.nullity, 0)
    assertEquals(RowReduction.reducedRowEchelon(DMat.zeros(0, 3), 0.0).orThrow.rank, 0)
    assertEquals(RowReduction.reducedRowEchelon(DMat.zeros(3, 0), 0.0).orThrow.rank, 0)
    assert(RowReduction.reducedRowEchelon(Matrix.eye(2), -1.0).isLeft)
    assert(RowReduction.reducedRowEchelon(Matrix.dense(1, 1)(Double.NaN), 0.0).isLeft)
    assert(RowReduction.nullSpaceBasis(Matrix.dense(1, 1)(Double.PositiveInfinity)).isLeft)
