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

  test("symmetric EigenDecomposition.requireConverged gates on convergence") {
    val values = Vec(1.0, 2.0)
    val vectors = Matrix.eye(2)

    val ok = EigenDecomposition(values, vectors, diagnostics(2, 2, Seq(0.0, 0.0), 3))
    assert(ok.requireConverged.isRight)
    assertEquals(ok.requireConverged.orThrow.size, 2)

    val bad = EigenDecomposition(values, vectors, diagnostics(2, 1, Seq(0.0, 0.4), 9))
    bad.requireConverged match
      case Left(LinAlgError.DidNotConverge(iters, residual)) =>
        assertEquals(iters, 9)
        assertEqualsDouble(residual, 0.4, 0.0)
      case other => fail(s"expected DidNotConverge, got $other")
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
