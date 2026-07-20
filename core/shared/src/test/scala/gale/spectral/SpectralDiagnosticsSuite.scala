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
