package gale.linalg

class RankConditionSuite extends munit.FunSuite:
  test("rankEstimate distinguishes full rank and rank-deficient matrices") {
    val full = Matrix.dense(3, 2)(
      1.0, 0.0,
      1.0, 1.0,
      1.0, 2.0
    )
    val deficient = Matrix.dense(3, 2)(
      1.0, 2.0,
      2.0, 4.0,
      3.0, 6.0
    )

    assertEquals(full.rankEstimate, 2)
    assertEquals(deficient.rankEstimate, 1)
  }

  test("conditionEstimate is exact for simple diagonal matrices") {
    val A = Matrix.dense(2, 2)(
      1.0, 0.0,
      0.0, 0.001
    )

    assert(math.abs(Matrix.eye(3).conditionEstimate.orThrow - 1.0) < 1e-12)
    assert(math.abs(A.conditionEstimate.orThrow - 1000.0) < 1e-9)
  }

  test("conditionEstimate uses the 1-norm (column sums), not the inf-norm") {
    // A and A^{-1} whose 1-norm condition number (4) differs from the inf-norm
    // condition number (9), so the assertion pins the documented 1-norm
    // semantics. A = [[1,0,0],[0,1,0],[1,1,1]], A^{-1} = [[1,0,0],[0,1,0],[-1,-1,1]].
    // ||A||_1 = 2, ||A^{-1}||_1 = 2  ->  kappa_1 = 4.
    val A = Matrix.dense(3, 3)(
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      1.0, 1.0, 1.0
    )
    val exact1Norm = 4.0
    val cond = A.conditionEstimate.orThrow
    assert(math.abs(cond - exact1Norm) < 1e-9, s"cond=$cond")
    // Hager's estimator is a lower bound on the true 1-norm; it must never
    // overshoot the exact value.
    assert(cond <= exact1Norm * (1.0 + 1e-8), s"cond=$cond")
  }

  test("conditionEstimate reports non-square and singular cases") {
    val rectangular = Matrix.zeros(2, 3)
    val singular = Matrix.dense(2, 2)(
      1.0, 2.0,
      2.0, 4.0
    )

    assert(rectangular.conditionEstimate.left.exists(_.isInstanceOf[LinAlgError.NonSquareMatrix]))
    assert(singular.conditionEstimate.orThrow.isPosInfinity)
  }
