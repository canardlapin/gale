package gale.linalg

class DenseCholeskyWorkspaceSuite extends munit.FunSuite:
  private val a = Array(4.0, 2.0, -2.0, 2.0, 10.0, 2.0, -2.0, 2.0, 6.0)
  // Independent exact factor L = [[2,0,0],[1,3,0],[-1,1,2]].
  private val lower = Array(2.0, 0.0, 0.0, 1.0, 3.0, 0.0, -1.0, 1.0, 2.0)

  test("factor and solve match an exact SPD oracle across units") {
    val workspace = new DenseCholeskyWorkspace(3)
    for scale <- Vector(1e-150, 1.0, 1e150) do
      val matrix = a.map(_ * scale)
      assert(workspace.factorLowerInPlace(3, matrix).isRight)
      for r <- 0 until 3; c <- 0 to r do
        assertEqualsDouble(matrix(r*3+c)/math.sqrt(scale), lower(r*3+c), 1e-14)
      // A [1,-2,3]^T = [-6,-12,12]^T.
      val rhs = Array(-6.0, -12.0, 12.0).map(_ * scale)
      assert(workspace.solveLowerInPlace(3, matrix, rhs).isRight)
      rhs.zip(Array(1.0, -2.0, 3.0)).foreach((x, y) => assertEqualsDouble(x, y, 1e-13))
  }

  test("offset spans, upper garbage and factor reuse preserve caller storage") {
    val workspace = new DenseCholeskyWorkspace(3)
    val matrix = Array(99.0) ++ a ++ Array(88.0)
    matrix(2) = Double.NaN
    assert(workspace.factorLowerInPlace(3, matrix, 1).isRight)
    assert(matrix(2).isNaN)
    assertEqualsDouble(matrix.head, 99.0, 0.0)
    assertEqualsDouble(matrix.last, 88.0, 0.0)
    val snapshot = matrix.clone()
    for _ <- 0 until 3 do
      val rhs = Array(77.0, -6.0, -12.0, 12.0, 66.0)
      assert(workspace.solveLowerInPlace(3, matrix, rhs, 1, 1).isRight)
      assertEqualsDouble(rhs(1), 1.0, 1e-14)
      assertEqualsDouble(rhs.head, 77.0, 0.0)
      assertEqualsDouble(rhs.last, 66.0, 0.0)
    assertEquals(matrix.toVector.map(java.lang.Double.doubleToLongBits), snapshot.toVector.map(java.lang.Double.doubleToLongBits))
  }

  test("failure is transactional and workspace can be reused after failure") {
    val workspace = new DenseCholeskyWorkspace(3)
    for bad <- Vector(Array(0.0), Array(-1.0), Array(Double.NaN), Array(Double.PositiveInfinity), Array(1.0, 2.0, 2.0, 1.0)) do
      val before = bad.clone()
      val size = if bad.length == 1 then 1 else 2
      assert(workspace.factorLowerInPlace(size, bad).isLeft)
      assertEquals(bad.toVector.map(java.lang.Double.doubleToLongBits), before.toVector.map(java.lang.Double.doubleToLongBits))
    val matrix = a.clone()
    assert(workspace.factorLowerInPlace(3, matrix).isRight)
    val rhs = Array(1.0, Double.NaN, 2.0)
    assert(workspace.solveLowerInPlace(3, matrix, rhs).isLeft)
    assertEqualsDouble(rhs(0), 1.0, 0.0)
    assert(rhs(1).isNaN)
    assert(workspace.checkPositiveDefinite(3, a).isRight)
  }

  test("capacity, spans, tolerance and empty shape have explicit behavior") {
    val workspace = new DenseCholeskyWorkspace(3, CholeskyOptions(1e-8))
    assert(workspace.factorLowerInPlace(1, Array(1e-10)).isLeft)
    assert(workspace.factorLowerInPlace(4, new Array[Double](16)).isLeft)
    assert(workspace.factorLowerInPlace(3, a, Int.MaxValue).isLeft)
    assert(workspace.solveLowerInPlace(3, lower, Array(1.0)).isLeft)
    assert(workspace.solveLowerInPlace(3, lower, Array(1.0, 2.0, 3.0), offset = -1).isLeft)
    assert(workspace.factorLowerInPlace(0, Array.emptyDoubleArray).isRight)
    assert(workspace.solveLowerInPlace(0, Array.emptyDoubleArray, Array.emptyDoubleArray).isRight)
  }

  test("overlapping factor and rhs spans use a complete input snapshot") {
    val workspace = new DenseCholeskyWorkspace(3)
    val shared = lower.clone()
    // A one-dimensional overlapping span exercises mutation of its own factor.
    val one = Array(2.0)
    assert(workspace.solveLowerInPlace(1, one, one).isRight)
    assertEqualsDouble(one(0), 0.5, 0.0)
    assert(workspace.solveLowerInPlace(3, shared, shared).isRight)
    val expected = aToMatrix.cholesky.toOption.get.solve(DVec.fromSeq(lower.take(3).toSeq)).toOption.get
    for i <- 0 until 3 do assertEqualsDouble(shared(i), expected(i), 1e-14)
  }

  private def aToMatrix: DMat = Matrix.tabulate(3, 3)((r, c) => a(r*3+c))

  test("allocation-free classification distinguishes failure kinds without mutating input") {
    val workspace = new DenseCholeskyWorkspace(3)
    assertEquals(workspace.testPositiveDefinite(3, a), DensePositiveDefiniteness.PositiveDefinite)
    assertEquals(workspace.testPositiveDefinite(1, Array(-1.0)), DensePositiveDefiniteness.NotPositiveDefinite)
    assertEquals(workspace.testPositiveDefinite(1, Array(Double.NaN)), DensePositiveDefiniteness.NonFiniteInput)
    assertEquals(workspace.testPositiveDefinite(4, a), DensePositiveDefiniteness.InvalidSpan)
    assertEquals(a.toVector, Vector(4.0, 2.0, -2.0, 2.0, 10.0, 2.0, -2.0, 2.0, 6.0))
  }

  test("non-finite computed solution is refused before writing the right hand side") {
    val workspace = new DenseCholeskyWorkspace(1)
    val rhs = Array(Double.MaxValue)
    assert(workspace.solveLowerInPlace(1, Array(1e-200), rhs).isLeft)
    assertEqualsDouble(rhs(0), Double.MaxValue, 0.0)
    intercept[IllegalArgumentException](new DenseCholeskyWorkspace(46341))
  }

  test("immutable Cholesky shares the finite-pivot refusal policy") {
    for value <- Vector(Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity) do
      assert(Matrix.tabulate(1, 1)((_, _) => value).cholesky.isLeft)
    // Finite input whose Schur update overflows must also fail closed.
    val overflow = Matrix.tabulate(2, 2)((r, c) => if r == c then 1.0 else Double.MaxValue)
    assert(overflow.cholesky.isLeft)
  }
