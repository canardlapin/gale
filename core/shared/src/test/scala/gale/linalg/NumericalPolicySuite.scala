package gale.linalg

/** Cross-platform numerical-policy locks for LU singularity, QR rank, and
  * `pinv` cutoff.
  *
  * These fixtures must stay IEEE-exact (or far from a cutoff) so JVM, Scala.js,
  * and experimental Wasm observe the same decision. Scala.js `PlatformMath.fma`
  * is `a*b+c`, not IEEE FMA; a QR-reconstructed `UΣVᵀ` with a planted `σ=0` is
  * only approximately singular there. LU reports `SingularMatrix` only on an
  * exact-zero or NaN pivot (`Factorizations.lu`), and `conditionEstimate` maps
  * that failure to `Right(+∞)`.
  *
  * Do not plant singularity by reconstructing a reduced SVD. Use an exact
  * outer product or an exact zero row/column.
  */
class NumericalPolicySuite extends munit.FunSuite:

  private val machineEps = 2.220446049250313e-16

  /** Exact rank-1 outer product `v vᵀ` for `v = (1,2,3)`. Every entry is an
    * integer, so the dependence is bit-identical on every platform.
    */
  private def exactRank1: DMat =
    Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      2.0, 4.0, 6.0,
      3.0, 6.0, 9.0
    )

  private def frob(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  private def asymmetry(a: DMat): Double =
    frob(a - a.t)

  private def assertMoorePenrose(a: DMat, p: DMat, tol: Double, clue: String): Unit =
    val scale = math.max(1.0, frob(a))
    assert(frob(a * p * a - a) < tol * scale, s"$clue: A A+ A != A (${frob(a * p * a - a)})")
    assert(frob(p * a * p - p) < tol * math.max(1.0, frob(p)), s"$clue: A+ A A+ != A+ (${frob(p * a * p - p)})")
    assert(asymmetry(a * p) < tol * scale, s"$clue: A A+ not symmetric (${asymmetry(a * p)})")
    assert(asymmetry(p * a) < tol * scale, s"$clue: A+ A not symmetric (${asymmetry(p * a)})")

  test("LU reports SingularMatrix only for an IEEE-exact zero pivot") {
    exactRank1.lu match
      case Left(_: LinAlgError.SingularMatrix) => ()
      case other => fail(s"exact rank-1 outer product: expected SingularMatrix, got $other")

    Matrix.zeros(2, 2).lu match
      case Left(LinAlgError.SingularMatrix(index)) =>
        assertEquals(index, 0, "all-zero matrix is singular at the first pivot")
      case other =>
        fail(s"all-zero: expected SingularMatrix(0), got $other")

    val zeroRow = Matrix.dense(2, 2)(1.0, 2.0, 0.0, 0.0)
    zeroRow.lu match
      case Left(_: LinAlgError.SingularMatrix) => ()
      case other => fail(s"exact zero row: expected SingularMatrix, got $other")
  }

  test("near-singularity with a nonzero pivot is not SingularMatrix or +∞") {
    // 1e-12 above the exact rank-1 line: LU must succeed. This is the class of
    // matrix a reconstructed UΣVᵀ plant becomes on JS, and it must not be used
    // as a +∞ / SingularMatrix fixture.
    val near = Matrix.dense(2, 2)(1.0, 2.0, 2.0, 4.0 + 1e-12)
    near.lu match
      case Right(_) => ()
      case Left(error) => fail(s"near-singular with nonzero pivot must factor, got $error")
    val cond = near.conditionEstimate.orThrow
    assert(cond.isFinite, s"near-singular condition must stay finite, got $cond")
    assert(cond > 1e10, s"near-singular condition should be large, got $cond")
  }

  test("conditionEstimate is +∞ exactly when LU reports SingularMatrix") {
    assert(exactRank1.conditionEstimate.orThrow.isPosInfinity)
    assert(math.abs(Matrix.eye(3).conditionEstimate.orThrow - 1.0) < 1e-14)
    for (m, n) <- Seq((2, 3), (3, 1)) do
      Matrix.zeros(m, n).conditionEstimate match
        case Left(LinAlgError.NonSquareMatrix(shape)) =>
          assertEquals(shape.rows.value, m)
          assertEquals(shape.cols.value, n)
        case other =>
          fail(s"${m}x$n: expected NonSquareMatrix, got $other")
  }

  test("rankEstimate uses the QR cutoff on exact diagonals and exact outer products") {
    // QR rank tolerance is 2 · max(m, n) · ε · max|R_ii|. Gaps sit far from
    // that ~1e-15 scale so Householder rounding cannot flip the decision.
    val n = 4
    val kept = Matrix.tabulate(n, n)((i, j) => if i == j then (if i < n - 1 then 1.0 else 1e-8) else 0.0)
    val dropped = Matrix.tabulate(n, n)((i, j) => if i == j then (if i < n - 1 then 1.0 else 1e-20) else 0.0)
    assertEquals(kept.rankEstimate, n)
    assertEquals(dropped.rankEstimate, n - 1)
    assertEquals(exactRank1.rankEstimate, 1)
    assertEquals(Matrix.zeros(5, 3).rankEstimate, 0)
    val tinyFull = Matrix.tabulate(3, 3)((i, j) => if i == j then 1e-20 else 0.0)
    assertEquals(tinyFull.rankEstimate, 3)
  }

  test("pinv cutoff on a rectangular diagonal is the MATLAB/SciPy rule") {
    val m = 7
    val n = 5
    val sigmaMax = 2.0
    val cutoff = math.max(m, n).toDouble * machineEps * sigmaMax
    val above = 1e-8
    val below = cutoff * 0.25
    val sigmas = IndexedSeq(sigmaMax, 1.0, 0.5, above, below)
    val a = Matrix.tabulate(m, n)((i, j) => if i == j then sigmas(i) else 0.0)
    val p = a.pinv.orThrow
    assertEquals(p.rows, n)
    assertEquals(p.cols, m)
    assertMoorePenrose(a, p, 1e-8, "diagonal pinv")
    assert(above > cutoff, s"fixture above=$above must exceed cutoff=$cutoff")
    assert(below <= cutoff, s"fixture below=$below must sit at or below cutoff=$cutoff")

    val eTrunc = DVec.tabulate(m)(i => if i == 4 then 1.0 else 0.0)
    val leaked = p * eTrunc
    assert(leaked.norm2 < 1e-8, s"truncated left singular vector leaked into pinv: ${leaked.norm2}")

    val eKept = DVec.tabulate(m)(i => if i == 3 then 1.0 else 0.0)
    val recovered = p * (eKept * above)
    val vKept = DVec.tabulate(n)(i => if i == 3 then 1.0 else 0.0)
    val align = math.abs(recovered.dot(vKept)) / (recovered.norm2 * vKept.norm2)
    assert(math.abs(align - 1.0) < 1e-8, s"kept singular vector not inverted: align=$align")
  }
