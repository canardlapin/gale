package gale.linalg

class IllConditionedSuite extends munit.FunSuite:
  // Item 5: the rank tolerance must scale with the matrix magnitude, not with a
  // hard floor at 1.0. A uniformly tiny but well-conditioned matrix is still
  // full rank.
  test("rankEstimate is scale-invariant for a uniformly tiny full-rank matrix") {
    val tiny = Matrix.dense(3, 3)(
      1e-20, 0.0, 0.0,
      0.0, 1e-20, 0.0,
      0.0, 0.0, 1e-20
    )
    assertEquals(tiny.rankEstimate, 3)
  }

  test("rankEstimate still detects genuine rank deficiency and zero matrices") {
    val mixed = Matrix.dense(3, 3)(
      1.0, 0.0, 0.0,
      0.0, 1e-25, 0.0,
      0.0, 0.0, 1e-25
    )
    assertEquals(mixed.rankEstimate, 1)
    assertEquals(Matrix.zeros(3, 3).rankEstimate, 0)
  }

  // Item 12: dense direct solve / conditioning on the notoriously ill-conditioned
  // Hilbert matrix, plus a genuinely singular LU.
  private def hilbert(n: Int): DMat =
    Matrix.tabulate(n, n)((i, j) => 1.0 / (i + j + 1))

  test("Hilbert 8x8 solve recovers the exact right-hand side") {
    val h = hilbert(8)
    val ones = Vec.fill(8)(1.0)
    val b = h * ones
    val x = h.solve(b).orThrow
    var i = 0
    var maxErr = 0.0
    while i < 8 do
      maxErr = math.max(maxErr, math.abs(x(i) - 1.0))
      i += 1
    assert(maxErr < 1e-4, s"max error $maxErr")
  }

  test("conditionEstimate of the 5x5 Hilbert matrix is near the known value") {
    val cond = hilbert(5).conditionEstimate.orThrow
    val known = 4.77e5
    assert(cond > known / 20.0 && cond < known * 20.0, s"cond=$cond")
  }

  test("LU of a singular matrix returns Left(SingularMatrix)") {
    val singular = Matrix.dense(2, 2)(
      1.0, 2.0,
      2.0, 4.0
    )
    assert(singular.lu.left.exists(_.isInstanceOf[LinAlgError.SingularMatrix]))
  }
