package gale.linalg

/** Gale's own numerical contracts for `rankEstimate`, `pinv`, and
  * `conditionEstimate`. These are not Breeze differential checks: near-cutoff
  * rank / pseudo-inverse / condition behaviour is a deliberate non-goal of
  * `parity/` because the libraries use different norms and cutoffs. This suite
  * pins Gale's documented definitions.
  *
  *   - `rankEstimate` is QR numerical rank at
  *     `2 · max(m, n) · ε · max|R_ii|`.
  *   - `pinv` zeros singular values at or below the MATLAB/SciPy cutoff
  *     `max(m, n) · ε · σ_max` and otherwise satisfies Moore–Penrose.
  *   - `conditionEstimate` is a Hager/Higham lower bound on `‖A‖₁‖A⁻¹‖₁`;
  *     singular squares are `Right(+∞)` and rectangular inputs are
  *     `Left(NonSquareMatrix)`.
  */
class GaleNumericalContractSuite extends munit.FunSuite:

  private val machineEps = 2.220446049250313e-16

  private def randomMat(m: Int, n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(m, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  private def orthonormal(n: Int, seed: Long): DMat =
    randomMat(n, n, seed).qr.q

  private def fromSvd(u: DMat, sigmas: IndexedSeq[Double], v: DMat): DMat =
    val p = sigmas.length
    Matrix.tabulate(u.rows, v.rows): (i, j) =>
      var sum = 0.0
      var k = 0
      while k < p do
        sum += u(i, k) * sigmas(k) * v(j, k)
        k += 1
      sum

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

  private def oneNorm(a: DMat): Double =
    var worst = 0.0
    var j = 0
    while j < a.cols do
      var sum = 0.0
      var i = 0
      while i < a.rows do
        sum += math.abs(a(i, j))
        i += 1
      worst = math.max(worst, sum)
      j += 1
    worst

  private def exactCond1(a: DMat): Double =
    val inv = a.solve(Matrix.eye(a.rows)).orThrow
    oneNorm(a) * oneNorm(inv)

  private def asymmetry(a: DMat): Double =
    frob(a - a.t)

  private def assertMoorePenrose(a: DMat, p: DMat, tol: Double, clue: String): Unit =
    val scale = math.max(1.0, frob(a))
    assert(frob(a * p * a - a) < tol * scale, s"$clue: A A+ A != A (${frob(a * p * a - a)})")
    assert(frob(p * a * p - p) < tol * math.max(1.0, frob(p)), s"$clue: A+ A A+ != A+ (${frob(p * a * p - p)})")
    assert(asymmetry(a * p) < tol * scale, s"$clue: A A+ not symmetric (${asymmetry(a * p)})")
    assert(asymmetry(p * a) < tol * scale, s"$clue: A+ A not symmetric (${asymmetry(p * a)})")

  // --- rankEstimate ----------------------------------------------------------

  test("rankEstimate recovers a planted numerical rank when the gap is clear of the QR cutoff") {
    val u = orthonormal(8, 11L)
    val v = orthonormal(5, 12L)
    val full = fromSvd(u, IndexedSeq(4.0, 3.0, 2.0, 1.0, 0.5), v)
    val rank3 = fromSvd(u, IndexedSeq(4.0, 3.0, 2.0, 0.0, 0.0), v)
    val rank1 = fromSvd(u, IndexedSeq(5.0, 0.0, 0.0, 0.0, 0.0), v)
    assertEquals(full.rankEstimate, 5)
    assertEquals(rank3.rankEstimate, 3)
    assertEquals(rank1.rankEstimate, 1)
    assertEquals(DMat.zeros(6, 4).rankEstimate, 0)
  }

  test("rankEstimate drops a diagonal entry below the QR cutoff and keeps one above it") {
    // QR rank tolerance is 2 · max(m, n) · ε · max|R_ii|. Gaps here sit far
    // from that ~1e-15 scale so Householder rounding cannot flip the decision.
    val n = 4
    val kept = Matrix.tabulate(n, n)((i, j) => if i == j then (if i < n - 1 then 1.0 else 1e-8) else 0.0)
    val dropped = Matrix.tabulate(n, n)((i, j) => if i == j then (if i < n - 1 then 1.0 else 1e-20) else 0.0)
    assertEquals(kept.rankEstimate, n)
    assertEquals(dropped.rankEstimate, n - 1)
  }

  test("rankEstimate is scale-invariant for a uniformly tiny full-rank matrix") {
    val tiny = Matrix.tabulate(3, 3)((i, j) => if i == j then 1e-20 else 0.0)
    assertEquals(tiny.rankEstimate, 3)
  }

  // --- pinv cutoff -----------------------------------------------------------

  test("pinv inverts singular values above the MATLAB/SciPy cutoff and zeros those at or below it") {
    // Rectangular diagonal so the SVD is essentially exact and the cutoff
    // comparison is not lost in a random-basis reconstruction.
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
    assertMoorePenrose(a, p, 1e-8, "cutoff pinv")
    assert(above > cutoff, s"fixture above=$above must exceed cutoff=$cutoff")
    assert(below <= cutoff, s"fixture below=$below must sit at or below cutoff=$cutoff")

    // Standard-basis e₄ is the truncated left singular vector: A⁺ e₄ ≈ 0.
    val eTrunc = DVec.tabulate(m)(i => if i == 4 then 1.0 else 0.0)
    val image = p * eTrunc
    assert(image.norm2 < 1e-8, s"truncated left singular vector leaked into pinv: ${image.norm2}")

    // e₃ is kept: A⁺ (σ₃ e₃) ≈ e₃ in the n-space (the corresponding right vector).
    val eKept = DVec.tabulate(m)(i => if i == 3 then 1.0 else 0.0)
    val recovered = p * (eKept * above)
    val vKept = DVec.tabulate(n)(i => if i == 3 then 1.0 else 0.0)
    val align = math.abs(recovered.dot(vKept)) / (recovered.norm2 * vKept.norm2)
    assert(math.abs(align - 1.0) < 1e-8, s"kept singular vector not inverted: align=$align")
  }

  test("pinv of a planted rank-deficient matrix still satisfies Moore–Penrose") {
    val a = fromSvd(orthonormal(6, 31L), IndexedSeq(3.0, 2.0, 0.0, 0.0), orthonormal(4, 32L))
    assertMoorePenrose(a, a.pinv.orThrow, 1e-9, "planted rank-2")
  }

  // --- conditionEstimate -----------------------------------------------------

  test("conditionEstimate is a Hager lower bound on the exact 1-norm condition number") {
    val a = Matrix.dense(3, 3)(
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      1.0, 1.0, 1.0
    )
    val exact = exactCond1(a)
    val est = a.conditionEstimate.orThrow
    assert(est <= exact * (1.0 + 1e-8), s"Hager overshot: est=$est exact=$exact")
    assert(est > 0.0)
    assert(math.abs(est - 4.0) < 1e-9, s"this fixture is exactly κ₁=4, got $est")
  }

  test("conditionEstimate of a diagonally dominant random matrix stays at or below the exact 1-norm") {
    val n = 6
    val base = randomMat(n, n, 41L)
    val a = Matrix.tabulate(n, n)((i, j) => if i == j then base(i, j) + 8.0 else base(i, j))
    val exact = exactCond1(a)
    val est = a.conditionEstimate.orThrow
    assert(est <= exact * (1.0 + 1e-8), s"Hager overshot: est=$est exact=$exact")
    assert(est >= 1.0)
  }

  test("conditionEstimate is 1 for the identity and +∞ for a singular square") {
    assert(math.abs(Matrix.eye(4).conditionEstimate.orThrow - 1.0) < 1e-14)
    // LU reports SingularMatrix only on an exact-zero pivot. A QR-reconstructed
    // U Σ Vᵀ with a planted σ=0 is only approximately singular on Scala.js
    // (fma is a*b+c), so plant an IEEE-exact rank-1 outer product instead.
    val singular = Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      2.0, 4.0, 6.0,
      3.0, 6.0, 9.0
    )
    assert(singular.conditionEstimate.orThrow.isPosInfinity)
  }

  test("conditionEstimate is Left(NonSquareMatrix) for every rectangular shape") {
    for (m, n) <- Seq((2, 3), (5, 2), (1, 4), (4, 1)) do
      randomMat(m, n, m * 17L + n).conditionEstimate match
        case Left(LinAlgError.NonSquareMatrix(shape)) =>
          assertEquals(shape.rows.value, m)
          assertEquals(shape.cols.value, n)
        case other => fail(s"${m}x$n: expected NonSquareMatrix, got $other")
  }
