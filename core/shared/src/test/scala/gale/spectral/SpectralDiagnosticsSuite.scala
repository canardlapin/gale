package gale.spectral

import gale.linalg.*

class SpectralDiagnosticsSuite extends munit.FunSuite:
  private def diagnostics(requested: Int, converged: Int, residuals: Seq[Double], iterations: Int): SpectralDiagnostics =
    SpectralDiagnostics(
      requested = requested,
      converged = converged,
      residuals = Vec(residuals*),
      orthogonalityError = 0.0,
      iterations = iterations
    )

  test("convergence status separates residual convergence from extreme certification") {
    val full = diagnostics(3, 3, Seq(1e-12, 2e-12, 5e-13), 4)
    assert(full.allConverged)
    assertEquals(full.convergenceStatus, SpectralConvergenceStatus.ResidualConverged)
    assertEqualsDouble(full.worstResidual, 2e-12, 0.0)

    val certified = full.copy(extremalityCertified = true)
    assert(certified.allConverged)
    assertEquals(certified.convergenceStatus, SpectralConvergenceStatus.ExtremeCertified)

    val partial = diagnostics(3, 2, Seq(1e-12, 9e-3), 7)
    assert(!partial.allConverged)
    assertEquals(partial.convergenceStatus, SpectralConvergenceStatus.NotConverged)
    assertEqualsDouble(partial.worstResidual, 9e-3, 0.0)

    assertEqualsDouble(diagnostics(0, 0, Nil, 0).worstResidual, 0.0, 0.0)
  }

  test("symmetric strict enforcement distinguishes all three convergence states") {
    val values = Vec(1.0, 2.0)
    val vectors = Matrix.eye(2)

    val ok = EigenDecomposition(values, vectors, diagnostics(2, 2, Seq(0.0, 0.0), 3))
    assert(ok.requireConverged.isRight)
    assertEquals(ok.requireConverged.orThrow.size, 2)
    ok.requireExtremeCertified match
      case Left(LinAlgError.SpectralExtremeNotCertified(iters, residual)) =>
        assertEquals(iters, 3)
        assertEqualsDouble(residual, 0.0, 0.0)
      case other => fail(s"expected SpectralExtremeNotCertified, got $other")

    val certified = ok.copy(diagnostics = ok.diagnostics.copy(extremalityCertified = true))
    assert(certified.requireExtremeCertified.isRight)
    assertEquals(certified.requireExtremeCertified.orThrow.size, 2)

    val bad = EigenDecomposition(values, vectors, diagnostics(2, 1, Seq(0.0, 0.4), 9))
    bad.requireConverged match
      case Left(LinAlgError.DidNotConverge(iters, residual)) =>
        assertEquals(iters, 9)
        assertEqualsDouble(residual, 0.4, 0.0)
      case other => fail(s"expected DidNotConverge, got $other")
    bad.requireExtremeCertified match
      case Left(LinAlgError.DidNotConverge(iters, residual)) =>
        assertEquals(iters, 9)
        assertEqualsDouble(residual, 0.4, 0.0)
      case other => fail(s"expected DidNotConverge, got $other")
  }

  test("dense symmetric decomposition is extreme-certified") {
    val dense = Eigen
      .eigSymmetric(
        Matrix.tabulate(3, 3)((row, col) => if row == col then (row + 1).toDouble else 0.0),
        EigenSelection.All
      )
      .toOption
      .get
    assertEquals(dense.diagnostics.convergenceStatus, SpectralConvergenceStatus.ExtremeCertified)
    assert(dense.requireExtremeCertified.isRight)
  }

  test("SVD.requireConverged gates on convergence") {
    val ok = SVD(Vec(3.0, 1.0), Matrix.eye(2), Matrix.eye(2), rank = 2, diagnostics(2, 2, Seq(0.0, 0.0), 1))
    assert(ok.requireConverged.isRight)

    val bad = SVD(Vec(3.0), Matrix.eye(1), Matrix.eye(1), rank = 1, diagnostics(2, 1, Seq(0.5), 5))
    assert(bad.requireConverged.isLeft)
  }

  test("nonsymmetric requireConverged gates on convergence") {
    val re = Vec(2.0, -1.0)
    val im = Vec(0.0, 0.0)
    val ok = new NonsymmetricEigenDecomposition(re, im, DMat.zeros(2, 0), None, diagnostics(2, 2, Seq(0.0, 0.0), 0))
    assert(ok.requireConverged.isRight)

    val bad = new NonsymmetricEigenDecomposition(re, im, DMat.zeros(2, 0), None, diagnostics(2, 1, Seq(0.3), 6))
    assert(bad.requireConverged.isLeft)
  }

  test("strict enforcement is consistent across diagnostics-carrying result types") {
    val residualOnly = diagnostics(1, 1, Seq(0.0), 2)
    val certified = residualOnly.copy(extremalityCertified = true)

    val svd = SVD(Vec(1.0), Matrix.eye(1), Matrix.eye(1), rank = 1, residualOnly)
    assert(svd.requireExtremeCertified.left.exists(_.isInstanceOf[LinAlgError.SpectralExtremeNotCertified]))
    assert(svd.copy(diagnostics = certified).requireExtremeCertified.isRight)

    val nonsymmetric = new NonsymmetricEigenDecomposition(
      Vec(1.0),
      Vec(0.0),
      Matrix.eye(1),
      None,
      residualOnly
    )
    assert(nonsymmetric.requireExtremeCertified.left.exists(_.isInstanceOf[LinAlgError.SpectralExtremeNotCertified]))

    val generalizedEigen = new GeneralizedEigenDecomposition(
      Vec(1.0),
      Vec(0.0),
      Vec(1.0),
      Matrix.eye(1),
      None,
      certified
    )
    assert(generalizedEigen.requireExtremeCertified.isRight)

    val cs = 1.0 / math.sqrt(2.0)
    val generalizedSvd = GeneralizedSVD(
      Matrix.eye(1),
      Matrix.eye(1),
      Matrix.eye(1),
      Vec(cs),
      Vec(cs),
      IndexedSeq(GeneralizedSingularValue.Finite(1.0)),
      certified
    )
    assert(generalizedSvd.requireExtremeCertified.isRight)
  }

  test("deferred measurements run once, print as deferred, and match eager diagnostics") {
    var residualRuns = 0
    var orthogonalityRuns = 0
    val deferred = SpectralDiagnostics.deferred(
      requested = 2,
      converged = 2,
      residuals = () =>
        residualRuns += 1
        Vec(1e-14, 3e-14),
      orthogonalityError = () =>
        orthogonalityRuns += 1
        2e-15,
      iterations = 0,
      rank = None,
      extremalityCertified = true
    )
    assert(deferred.allConverged)
    assert(deferred.requireExtremeCertified(()).isRight)
    assertEquals(deferred.toString, "SpectralDiagnostics(2,2,<deferred>,<deferred>,0,None,true,None)")
    assertEquals((residualRuns, orthogonalityRuns), (0, 0))

    assertEqualsDouble(deferred.worstResidual, 3e-14, 0.0)
    assertEqualsDouble(deferred.worstResidual, 3e-14, 0.0)
    assertEquals((residualRuns, orthogonalityRuns), (1, 0))
    assert(deferred.toString.endsWith(",<deferred>,0,None,true,None)"), deferred.toString)

    // copy is strict: it measures what it carries over and shares the values.
    val copied = deferred.copy(iterations = 4)
    assertEquals((residualRuns, orthogonalityRuns), (1, 1))
    assertEquals(copied.iterations, 4)
    assert(copied.residuals eq deferred.residuals)
    assertEqualsDouble(copied.orthogonalityError, 2e-15, 0.0)
    assertEquals((residualRuns, orthogonalityRuns), (1, 1))

    // DVec compares by reference (as the former case class did), so share it.
    val eager = SpectralDiagnostics(2, 2, deferred.residuals, 2e-15, 0, None, extremalityCertified = true)
    assert(!eager.toString.contains("<deferred>"), eager.toString)
    assertEquals(deferred, eager)
    assertEquals(deferred.hashCode, eager.hashCode)
    assertEquals(deferred.toString, eager.toString)
    assertNotEquals(deferred, eager.copy(orthogonalityError = 1.0))
    assertEquals((residualRuns, orthogonalityRuns), (1, 1))
  }

  test("equality takes a deferred measurement") {
    var runs = 0
    def make() = SpectralDiagnostics.deferred(1, 1, () => { runs += 1; Vec(0.0) }, () => 0.0, 0, None, true)
    val first = make()
    assertNotEquals(first, make()) // distinct DVec instances compare by reference
    assertEquals(runs, 2)
  }

  test("equality is reflexive for NaN measurements without forcing deferred work") {
    val residuals = Vec(Double.NaN)
    val eager = SpectralDiagnostics(1, 1, residuals, Double.NaN, 0)
    assert(eager.equals(eager))
    // Preserve the former case class's field comparison for distinct objects.
    assert(!eager.equals(SpectralDiagnostics(1, 1, residuals, Double.NaN, 0)))
    assert(!eager.equals(null))
    var runs = 0
    val deferred = SpectralDiagnostics.deferred(
      1, 1,
      () => { runs += 1; residuals },
      () => { runs += 1; Double.NaN },
      0, None, true
    )
    assert(deferred.equals(deferred))
    assertEquals(runs, 0)
    assert(deferred.orthogonalityError.isNaN)
    assert(deferred.equals(deferred))
    assertEquals(runs, 1)
  }

  test("dense facades measure diagnostics on access to the same values as explicit products") {
    val n = 7
    val a = Matrix.tabulate(n, n)((r, c) => 1.0 / (1 + r + c) + (if r == c then r.toDouble else 0.0))
    val eig = Eigen.eigSymmetric(a, EigenSelection.All).toOption.get
    val v = eig.eigenvectors
    val av = a * v
    var worst = 0.0
    var c = 0
    while c < n do
      val r = (av.col(c) - v.col(c) * eig.eigenvalues(c)).norm2
      assertEqualsDouble(eig.diagnostics.residuals(c), r, 0.0)
      worst = math.max(worst, r)
      c += 1
    assertEqualsDouble(eig.diagnostics.worstResidual, worst, 0.0)
    assert(eig.diagnostics.orthogonalityError < 1e-13)

    val workspace = Eigen.eigSymmetricWith(a, EigenSelection.All, gale.linalg.DenseWorkspace.empty).toOption.get
    assertEquals(workspace.diagnostics.residuals.toSeq, eig.diagnostics.residuals.toSeq)
    assertEqualsDouble(workspace.diagnostics.orthogonalityError, eig.diagnostics.orthogonalityError, 0.0)

    val tall = Matrix.tabulate(9, 5)((r, c) => math.sin(1.0 + r * 5 + c))
    val svd = Svds.svd(tall, SingularSelection.All).toOption.get
    assertEquals(svd.diagnostics.residuals.length, 5)
    assert(svd.diagnostics.worstResidual < 1e-12)
    assert(svd.diagnostics.orthogonalityError < 1e-13)
  }

  /** A strongly non-normal matrix with two complex pairs (rotation blocks) and
    * real eigenvalues (diagonal ramp), scaled by `scale`.
    */
  private def nonNormal(n: Int, scale: Double): DMat =
    Matrix.tabulate(n, n)((r, c) =>
      val v =
        if r == c then 0.5 * r
        else if r < 4 && c < 4 && r / 2 == c / 2 then (if r < c then 2.0 else -2.0)
        else if c > r then 3.0 * math.cos(1.0 + r * n + c)
        else 0.01 * math.sin(1.0 + r * n + c)
      v * scale
    )

  /** The residual of eigenvalue `j` from the explicit `A V` columns, by the
    * definition: `‖A v − λ v‖` in real arithmetic, a pair sharing one value.
    */
  private def explicitResidual(av: DMat, v: DMat, re: DVec, im: DVec, j: Int): Double =
    if im(j) == 0.0 then (av.col(j) - v.col(j) * re(j)).norm2
    else
      val k = if im(j) > 0.0 then j else j - 1
      val vr = v.col(k)
      val vi = v.col(k + 1)
      val realPart = av.col(k) - (vr * re(k) - vi * im(k))
      val imagPart = av.col(k + 1) - (vi * re(k) + vr * im(k))
      Vec((realPart.toSeq ++ imagPart.toSeq)*).norm2

  test("dense nonsymmetric eigen measures residuals on access by their definition") {
    val n = 9
    val a = nonNormal(n, 1e6)
    val eig = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Right).toOption.get
    assert(eig.diagnostics.residuals ne null)
    val v = eig.rightVectorsPacked
    val im = eig.im
    assert((0 until n).count(j => im(j) != 0.0) == 4 && (0 until n).exists(j => im(j) == 0.0))
    // The same A·V product the facade forms, so both sides share its rounding
    // and the comparison sees only the residual formula and its reduction.
    val av = a * v
    var j = 0
    while j < n do
      val expected = explicitResidual(av, v, eig.re, im, j)
      val got = eig.diagnostics.residuals(j)
      assert(expected > 0.0, s"residual $j vanished")
      assert(math.abs(got - expected) <= 1e-12 * expected, s"residual $j: $got vs $expected")
      j += 1
    assert(eig.diagnostics.worstResidual < 1e-9 * 1e6)
    assertEqualsDouble(eig.diagnostics.orthogonalityError, 0.0, 0.0)

    val valuesOnly = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    assertEquals(valuesOnly.diagnostics.residuals.toSeq, Seq.fill(n)(0.0))
  }

  test("dense nonsymmetric eigen defers only the residuals") {
    val eig = Eigen.eigNonsymmetric(nonNormal(6, 1.0), EigenSelection.All, EigenVectors.Right).toOption.get
    assert(
      eig.diagnostics.toString.endsWith(",<deferred>,0.0,0,None,true,None)"),
      eig.diagnostics.toString
    )
    assert(eig.diagnostics.worstResidual < 1e-12)
    assert(!eig.diagnostics.toString.contains("<deferred>"), eig.diagnostics.toString)
  }

  test("dense nonsymmetric pair residuals stay finite near the overflow threshold") {
    val n = 9
    val scale = 1e300
    val a = nonNormal(n, scale)
    val eig = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Right).toOption.get
    val im = eig.im
    assert((0 until n).exists(j => im(j) != 0.0))
    var j = 0
    while j < n do
      val r = eig.diagnostics.residuals(j)
      assert(r.isFinite && r >= 0.0, s"residual $j = $r")
      assert(r < 1e-9 * scale, s"residual $j = $r")
      j += 1
  }
