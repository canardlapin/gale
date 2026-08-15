package gale.spectral

import gale.linalg.*

/** Live-solver trust checks for the generalized symmetric-definite pencil
  * `A x = λ B x`: residuals, `B`-orthonormality, SPD-`B` failure, and the
  * `ResidualConverged` vs `ExtremeCertified` distinction on actual solver
  * results (not constructed diagnostics).
  */
class GeneralizedEigenTrustSuite extends munit.FunSuite:

  private def randomSymmetric(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(n, n): (i, j) =>
      val v = rng.nextDouble() * 2.0 - 1.0
      if i >= j then v else 0.0

  private def randomSpd(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    val raw = Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)
    val gram = raw * raw.t
    Matrix.tabulate(n, n)((i, j) => gram(i, j) + (if i == j then n.toDouble else 0.0))

  private def diagonal(values: IndexedSeq[Double]): DMat =
    Matrix.tabulate(values.length, values.length)((i, j) => if i == j then values(i) else 0.0)

  private def frobenius(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  private def generalizedResidual(a: DMat, b: DMat, d: EigenDecomposition, c: Int): Double =
    val x = d.eigenvectors.col(c)
    ((a * x) - (b * x) * d.eigenvalues(c)).norm2

  private def bOrthoError(d: EigenDecomposition, b: DMat): Double =
    val g = d.eigenvectors.t * (b * d.eigenvectors)
    frobenius(g - Matrix.eye(g.rows))

  private def diagonalOperator(values: IndexedSeq[Double]): DoubleLinearOperator =
    LinearOperator.fromFunction(values.length, values.length): (x, into) =>
      var i = 0
      while i < values.length do
        into(i) = values(i) * x(i)
        i += 1

  // --- dense live results ----------------------------------------------------

  test("dense generalized All is ExtremeCertified with tiny residuals and B-orthonormal vectors") {
    val n = 8
    val a = randomSymmetric(n, 101L)
    val b = randomSpd(n, 202L)
    val d = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).orThrow
    val aSym = Matrix.tabulate(n, n)((i, j) => if i >= j then a(i, j) else a(j, i))
    val scale = math.max(frobenius(aSym), frobenius(b))
    assertEquals(d.diagnostics.convergenceStatus, SpectralConvergenceStatus.ExtremeCertified)
    assert(d.requireExtremeCertified.isRight)
    assert(bOrthoError(d, b) < 1e-10, s"XᵀBX≠I: ${bOrthoError(d, b)}")
    var c = 0
    while c < n do
      val res = generalizedResidual(aSym, b, d, c)
      assert(res < 1e-8 * scale, s"residual $c = $res")
      assert(math.abs(d.diagnostics.residuals(c) - res) < 1e-12, s"diagnostics residual $c")
      c += 1
  }

  test("dense generalized with an ill-conditioned SPD B still reports honest residuals") {
    // κ₁(B) is large; the Cholesky reduction amplifies error, but the call is
    // a dense one-shot so membership is ExtremeCertified and residuals stay
    // the accuracy signal.
    val a = Matrix.tabulate(4, 4)((i, j) => if i == j then (i + 1).toDouble else 0.0)
    val b = Matrix.tabulate(4, 4)((i, j) => if i == j then math.pow(10.0, i.toDouble - 1.0) else 0.0)
    val d = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).orThrow
    assertEquals(d.diagnostics.convergenceStatus, SpectralConvergenceStatus.ExtremeCertified)
    assert(d.requireExtremeCertified.isRight)
    assert(bOrthoError(d, b) < 1e-8, s"B-ortho ${bOrthoError(d, b)}")
    var c = 0
    while c < 4 do
      assert(generalizedResidual(a, b, d, c) < 1e-8, s"residual $c")
      c += 1
    // Analytic λ_i = a_ii / b_ii = {10, 2, 0.3, 0.04} sorted ascending.
    assertEqualsDouble(d.eigenvalues(0), 0.04, 1e-10, "smallest eigenvalue")
    assertEqualsDouble(d.eigenvalues(3), 10.0, 1e-10, "largest eigenvalue")
  }

  test("dense generalized rejects a non-SPD B with NotPositiveDefinite") {
    val a = Matrix.eye(3)
    val indefinite = Matrix.tabulate(3, 3)((i, j) => if i == j then Array(1.0, -2.0, 1.0)(i) else 0.0)
    Eigen.eigSymmetricGeneralized(a, indefinite, EigenSelection.All) match
      case Left(_: LinAlgError.NotPositiveDefinite) => ()
      case other                                    => fail(s"expected NotPositiveDefinite, got $other")
  }

  // --- operator live results -------------------------------------------------

  test("LOBPCG residuals can pass without certifying the requested extreme") {
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
    val initial = DMat.tabulate(6, 2)((row, col) => if row == col then 1.0 else 0.0)
    val d = Eigen
      .eigSymmetricGeneralized(
        diagonalOperator(values).assumeSymmetricOperator,
        diagonalOperator(IndexedSeq.fill(6)(1.0)).assumePositiveDefiniteOperator,
        6,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(tolerance = 1e-12, initialSubspace = Some(initial))
      )
      .orThrow
    assert(d.requireConverged.isRight)
    assertEquals(d.diagnostics.convergenceStatus, SpectralConvergenceStatus.ResidualConverged)
    d.requireExtremeCertified match
      case Left(_: LinAlgError.SpectralExtremeNotCertified) => ()
      case other => fail(s"expected SpectralExtremeNotCertified, got $other")
    var c = 0
    while c < d.size do
      assert(d.diagnostics.residuals(c) < 1e-10, s"residual $c")
      c += 1
    assert(d.diagnostics.orthogonalityError < 1e-10, s"B-ortho ${d.diagnostics.orthogonalityError}")
  }

  test("LOBPCG certifies the extreme when the trial space is the whole ambient") {
    // n = 4, k = 2: one LOBPCG step concatenates the k-block with a complement
    // of width up to n-k, so the trial reaches n and extremalityCertified flips.
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0)
    val d = Eigen
      .eigSymmetricGeneralized(
        diagonalOperator(values).assumeSymmetricOperator,
        diagonalOperator(IndexedSeq.fill(4)(1.0)).assumePositiveDefiniteOperator,
        4,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(tolerance = 1e-12, maxIterations = 8)
      )
      .orThrow
    assert(d.requireConverged.isRight)
    assertEqualsDouble(d.eigenvalues(0), 1.0, 1e-10, "λ0")
    assertEqualsDouble(d.eigenvalues(1), 2.0, 1e-10, "λ1")
    assertEquals(d.diagnostics.convergenceStatus, SpectralConvergenceStatus.ExtremeCertified)
    assert(d.requireExtremeCertified.isRight)
  }

  test("operator path rejects encountered non-positive B geometry") {
    val initial = DMat.tabulate(4, 2)((row, col) => if row == col then 1.0 else 0.0)
    val result = Eigen.eigSymmetricGeneralized(
      diagonalOperator(IndexedSeq(1.0, 2.0, 3.0, 4.0)).assumeSymmetricOperator,
      diagonalOperator(IndexedSeq(1.0, -1.0, 2.0, 3.0)).assumePositiveDefiniteOperator,
      4,
      EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
      GeneralizedSpectralOptions(initialSubspace = Some(initial))
    )
    result match
      case Left(_: LinAlgError.NotPositiveDefinite) => ()
      case other                                    => fail(s"expected NotPositiveDefinite, got $other")
  }

  test("generalized Lanczos full-space projection is ExtremeCertified; a cap is not") {
    val n = 8
    val values = IndexedSeq.tabulate(n)(i => (i + 1).toDouble)
    val a = diagonal(values)
    val b = Matrix.eye(n)
    val solve =
      MetricSolveOperator
        .bind(b.assumePositiveDefinite, LinearSolveOperator.direct(b.cholesky.orThrow))
        .orThrow
    val full = Eigen
      .eigSymmetricGeneralizedLanczos(
        a.assumeSymmetric,
        solve,
        n,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedLanczosOptions(tolerance = 1e-10, maxIterations = 4, subspaceDimension = Some(n))
      )
      .orThrow
    assertEquals(full.diagnostics.convergenceStatus, SpectralConvergenceStatus.ExtremeCertified)
    assert(full.requireExtremeCertified.isRight)

    val capped = Eigen
      .eigSymmetricGeneralizedLanczos(
        a.assumeSymmetric,
        solve,
        n,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedLanczosOptions(tolerance = 1e-8, maxIterations = 20, subspaceDimension = Some(5))
      )
      .orThrow
    assert(capped.requireConverged.isRight, capped.diagnostics.toString)
    assertEquals(capped.diagnostics.convergenceStatus, SpectralConvergenceStatus.ResidualConverged)
    capped.requireExtremeCertified match
      case Left(_: LinAlgError.SpectralExtremeNotCertified) => ()
      case other => fail(s"expected SpectralExtremeNotCertified, got $other")
  }
