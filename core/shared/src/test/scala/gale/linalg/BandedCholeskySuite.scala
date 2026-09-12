package gale.linalg

import gale.sparse.Banded

class BandedCholeskySuite extends munit.FunSuite:
  private def factor(bands: DMat, options: CholeskyOptions = CholeskyOptions.Default): BandedCholesky =
    BandedCholesky.factorLower(bands, options).fold(throw _, identity)

  private def close(actual: DMat, expected: DMat, tolerance: Double = 1e-11): Unit =
    assertEquals(actual.shape, expected.shape)
    for i <- 0 until actual.rows; j <- 0 until actual.cols do
      assertEqualsDouble(actual(i, j), expected(i, j), tolerance)

  private def tridiagonal(n: Int): DMat =
    Matrix.tabulate(n, math.min(2, math.max(1, n)))((_, d) => if d == 0 then 2.0 else -1.0)

  test("analytic discrete Laplacian factors, determinant, inverse and solution"):
    for n <- List(1, 2, 7, 31, 300) do
      val f = factor(tridiagonal(n))
      assertEqualsDouble(f.logDet, math.log(n + 1.0), 2e-12)
      val x = f.solve(Vec.fill(n)(1.0)).fold(throw _, identity)
      for i <- 0 until n do
        assertEqualsDouble(f.lowerBands(i, 0), math.sqrt((i + 2.0) / (i + 1.0)), 2e-14)
        if i > 0 then assertEqualsDouble(f.lowerBands(i, 1), -math.sqrt(i.toDouble / (i + 1.0)), 2e-14)
        assertEqualsDouble(x(i), (i + 1.0) * (n - i) / 2.0, 3e-8)
      if n <= 31 then
        val inverse = f.solve(Matrix.eye(n)).fold(throw _, identity)
        for i <- 0 until n; j <- 0 until n do
          val expected = (math.min(i, j) + 1.0) * (n - math.max(i, j)) / (n + 1.0)
          assertEqualsDouble(inverse(i, j), expected, 2e-12)

  test("known lower factors reconstruct and solve every supported bandwidth"):
    val n = 13
    for b <- List(0, 1, 2, 6, n - 1) do
      val lower = Matrix.tabulate(n, n): (i, j) =>
        if i == j then 2.0 + i / 10.0
        else if j < i && i - j <= b then math.sin(i * 3.0 + j) / 8.0
        else 0.0
      val a = lower * lower.t
      val bands = Matrix.tabulate(n, b + 1)((i, d) => if d <= i then a(i, i - d) else Double.NaN)
      val f = factor(bands)
      for i <- 0 until n; d <- 0 to b do
        assertEqualsDouble(f.lowerBands(i, d), if d <= i then lower(i, i - d) else 0.0, 2e-14)
      val expected = Matrix.tabulate(n, 4)((i, j) => math.cos(i + j * 2.0))
      val rhs = a * expected
      close(f.solve(rhs).fold(throw _, identity), expected)
      close(f.solveTranspose(rhs).fold(throw _, identity), expected)
      close(f.solveLowerTranspose(f.solveLower(rhs).fold(throw _, identity)).fold(throw _, identity), expected)
      assertEqualsDouble(f.logDet, 2.0 * (0 until n).map(i => math.log(lower(i, i))).sum, 2e-13)

  test("strided bands and RHS views are read logically and preserved"):
    val bands = DMatBuilder.from(tridiagonal(9).t).result().t
    val f = factor(bands)
    val source = Matrix.tabulate(12, 5)((i, j) => i + j * 0.25)
    val rhs = source.slice(1, 10, 1, 4)
    val snapshot = DMatBuilder.from(rhs).result()
    val x = f.solve(rhs).fold(throw _, identity)
    val dense = Matrix.tabulate(9, 9)((i, j) => if i == j then 2.0 else if math.abs(i - j) == 1 then -1.0 else 0.0)
    close(dense * x, rhs)
    close(rhs, snapshot, 0.0)
    val v = rhs.col(1)
    val solved = f.solve(v).fold(throw _, identity)
    for i <- 0 until 9 do assertEqualsDouble(solved(i), x(i, 1), 1e-12)
    close(bands, tridiagonal(9), 0.0)

  test("general band adapter uses the lower triangle and rejects nonsquare input"):
    val a = Matrix.tabulate(8, 8): (i, j) =>
      if i == j then 2.0 else if i == j + 1 then -1.0 else if j == i + 1 then 123.0 else 0.0
    val f = BandedCholesky.factorLower(Banded.fromDense(a)).fold(throw _, identity)
    close(f.lowerBands, factor(tridiagonal(8)).lowerBands)
    assert(BandedCholesky.factorLower(Banded.fromDense(Matrix.zeros(2, 3))).isLeft)

  test("consuming a builder transfers ownership on success and failure"):
    val builder = DMatBuilder.from(tridiagonal(5))
    val f = builder.consumeBandedCholesky().fold(throw _, identity)
    intercept[LinAlgError.UnsupportedOperation](builder(0, 0))
    intercept[LinAlgError.UnsupportedOperation](builder.update(0, 0, 7.0))
    intercept[LinAlgError.UnsupportedOperation](builder.result())
    assertEqualsDouble(f.logDet, math.log(6.0), 1e-14)
    val bad = DMatBuilder.zeros(3, 2)
    assert(bad.consumeBandedCholesky().isLeft)
    intercept[LinAlgError.UnsupportedOperation](bad.fill(1.0))
    val badShape = DMatBuilder.zeros(2, 3)
    assert(badShape.consumeBandedCholesky().isLeft)
    intercept[LinAlgError.UnsupportedOperation](badShape.result())

  test("in-place strided vector and multiple RHS solves preserve other destinations"):
    val f = factor(tridiagonal(8))
    val rhs = Matrix.tabulate(8, 3)((i, j) => i + j + 1.0)
    val expected = f.solve(rhs).fold(throw _, identity)
    val builder = DMatBuilder.from(rhs)
    assert(f.solveInPlace(builder.mutableColumn(1)).isRight)
    for i <- 0 until 8; j <- 0 until 3 do
      assertEqualsDouble(builder(i, j), if j == 1 then expected(i, j) else rhs(i, j), 1e-12)
    val split = DMatBuilder.from(rhs)
    assert(f.solveLowerInPlace(split).isRight)
    assert(f.solveLowerTransposeInPlace(split).isRight)
    close(split.result(), expected)
    val full = DMatBuilder.from(rhs)
    assert(f.solveTransposeInPlace(full).isRight)
    close(full.result(), expected)
    val mutable = MutableDVec.from(rhs.col(2))
    assert(f.solveLowerInPlace(mutable).isRight)
    assert(f.solveLowerTransposeInPlace(mutable).isRight)
    for i <- 0 until 8 do assertEqualsDouble(mutable(i), expected(i, 2), 1e-12)

  test("invalid RHS dimensions and nonfinite values are rejected before writes"):
    val f = factor(tridiagonal(4))
    assert(f.solve(Vec.zeros(3)).isLeft)
    assert(f.solve(Matrix.zeros(3, 2)).isLeft)
    val wrong = DMatBuilder.from(Matrix.tabulate(3, 2)((_, _) => 7.0))
    assert(f.solveInPlace(wrong).isLeft)
    for i <- 0 until 3; j <- 0 until 2 do assertEqualsDouble(wrong(i, j), 7.0, 0.0)
    val bad = DMatBuilder.from(Matrix.tabulate(4, 2)((i, j) => if i == 3 && j == 1 then Double.PositiveInfinity else 2.0))
    assert(f.solveInPlace(bad).isLeft)
    for i <- 0 until 4; j <- 0 until 2 if i != 3 || j != 1 do assertEqualsDouble(bad(i, j), 2.0, 0.0)
    val closed = DMatBuilder.zeros(4, 1)
    closed.result()
    intercept[LinAlgError.UnsupportedOperation](f.solveInPlace(closed))

  test("empty, singleton, zero-column RHS, padding, and invalid storage shapes"):
    val empty = factor(Matrix.zeros(0, 1))
    assertEquals(empty.size, 0)
    assertEqualsDouble(empty.logDet, 0.0, 0.0)
    assertEqualsDouble(empty.conditioning.conditionNumberUpperBound, 1.0, 0.0)
    assertEquals(empty.solve(Vec.zeros(0)).fold(throw _, identity).length, 0)
    assertEquals(empty.solve(Matrix.zeros(0, 3)).fold(throw _, identity).cols, 3)
    val f = factor(Matrix.dense(1, 1)(4.0))
    assertEqualsDouble(f.solve(Vec(6.0)).fold(throw _, identity)(0), 1.5, 0.0)
    assertEquals(f.solve(Matrix.zeros(1, 0)).fold(throw _, identity).cols, 0)
    val padded = Matrix.tabulate(3, 3)((i, d) => if d > i then Double.NaN else if d == 0 then 2.0 else 0.0)
    assertEqualsDouble(factor(padded).lowerBands(0, 2), 0.0, 0.0)
    for bands <- List(Matrix.zeros(2, 0), Matrix.zeros(2, 3), Matrix.zeros(0, 2)) do
      assert(BandedCholesky.factorLower(bands).isLeft)

  test("nonpositive pivots, finite input and absolute pivot tolerance are explicit"):
    val singular = Matrix.dense(2, 2)(1.0, 0.0, 1.0, 1.0)
    BandedCholesky.factorLower(singular) match
      case Left(LinAlgError.NotPositiveDefinite(index)) => assertEquals(index, 1)
      case other => fail(s"expected pivot 1 failure, got $other")
    for value <- List(Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity) do
      assert(BandedCholesky.factorLower(Matrix.dense(1, 1)(value)).isLeft)
    val tiny = Matrix.dense(1, 1)(1e-200)
    assert(BandedCholesky.factorLower(tiny).isRight)
    assert(BandedCholesky.factorLower(tiny, CholeskyOptions(pivotTolerance = 1e-199)).isLeft)
    assert(factor(tiny).solve(Vec(1e200)).isLeft)

  test("conditioning bounds agree with analytic diagonal and Laplacian inverses"):
    val diagonal = factor(Matrix.dense(3, 1)(2.0, 5.0, 8.0)).conditioning
    assertEqualsDouble(diagonal.matrixOneNorm, 8.0, 0.0)
    assertEqualsDouble(diagonal.inverseOneNormUpperBound, 0.5, 1e-15)
    assertEqualsDouble(diagonal.conditionNumberUpperBound, 4.0, 2e-15)
    for n <- List(2, 7, 30); scale <- List(1e-150, 1.0, 1e150) do
      val f = factor(Matrix.tabulate(n, 2)((_, d) => scale * (if d == 0 then 2.0 else -1.0)))
      val inverseNorm = (0 until n).map(i => (i + 1.0) * (n - i) / 2.0).max
      val norm = if n == 2 then 3.0 else 4.0
      assertEqualsDouble(f.conditioning.conditionNumberUpperBound, norm * inverseNorm, 2e-10)
      assertEqualsDouble(f.conditioning.reciprocalConditionLowerBound, 1.0 / (norm * inverseNorm), 2e-14)
      assert(f.conditioning.pivotRatio > f.conditioning.reciprocalConditionLowerBound)

  test("conditioning overflow is explicit and zero bands do not propagate NaN"):
    val f = factor(Matrix.dense(2, 2)(1e-320, 0.0, 1e300, 0.0))
    assertEquals(f.conditioning.inverseOneNormUpperBound, Double.PositiveInfinity)
    assertEquals(f.conditioning.conditionNumberUpperBound, Double.PositiveInfinity)
    assertEqualsDouble(f.conditioning.reciprocalConditionLowerBound, 0.0, 0.0)
    assertEqualsDouble(f.conditioning.pivotRatio, 0.0, 0.0)

  test("100000-row narrow-band system uses packed storage and solves without dense expansion"):
    val n = 100000
    val bands = Matrix.tabulate(n, 3)((_, d) => if d == 0 then 5.0 else -0.25)
    val f = factor(bands)
    assertEquals(f.lowerBands.rows * f.lowerBands.cols, 300000)
    val rhs = Vec.tabulate(n)(i => 5.0 - 0.25 * (math.min(i, 2) + math.min(n - i - 1, 2)))
    val x = f.solve(rhs).fold(throw _, identity)
    for i <- 0 until n do assertEqualsDouble(x(i), 1.0, 2e-15)
