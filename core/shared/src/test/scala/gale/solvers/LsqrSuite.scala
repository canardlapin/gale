package gale.solvers

import gale.linalg.*

class LsqrSuite extends munit.FunSuite:
  private def assertVecClose(actual: DVec, expected: DVec, tol: Double): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assert(math.abs(actual(i) - expected(i)) <= tol, s"at $i: ${actual(i)} != ${expected(i)}")
      i += 1

  // Overdetermined 4x2 with a consistent right-hand side b = A x_true.
  private val a = Matrix.dense(4, 2)(
    1.0, 1.0,
    1.0, 2.0,
    1.0, 3.0,
    1.0, 4.0
  )

  test("consistent overdetermined system matches the QR least-squares solution") {
    val xTrue = Vec(1.0, 2.0)
    val b = a * xTrue
    val result = lsqr(a, b, SolverConfig(tolerance = 1e-12, maxIterations = 100))
    assert(result.converged, s"lsqr did not converge (residual ${result.residual})")
    assertVecClose(result.x, a.leastSquares(b).orThrow, 1e-8)
    assertVecClose(result.x, xTrue, 1e-8)
  }

  test("inconsistent system matches the normal-equations solution") {
    val b = Vec(3.0, 5.0, 7.0, 10.0) // not in the range of A
    // Normal equations: (AᵀA) x = Aᵀ b.
    val normalSolution = (a.t * a).solve(a.t * b).orThrow
    val result = lsqr(a, b, SolverConfig(tolerance = 1e-12, maxIterations = 100))
    assert(result.converged, s"lsqr did not converge (residual ${result.residual})")
    assertVecClose(result.x, normalSolution, 1e-7)
  }

  test("ill-conditioned tall matrix: lsqr needs no more iterations than cgnr") {
    // A 6x4 Vandermonde with column scaling makes AᵀA badly conditioned, which
    // hurts cgnr (it works on AᵀA, squaring the condition number) more than lsqr.
    val ts = Seq(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    val scales = Seq(1.0, 1e-1, 1e-2, 1e-3)
    val values =
      for
        t <- ts
        j <- 0 until 4
      yield math.pow(t, j.toDouble) * scales(j)
    val m = Matrix.dense(6, 4, values)
    val bTrue = Vec(1.0, -1.0, 2.0, -2.0, 3.0, -3.0)

    val config = SolverConfig(tolerance = 1e-8, maxIterations = 2000)
    val lsqrResult = lsqr(m, bTrue, config)
    val cgnrResult = cgnr(m, bTrue, config)

    assert(lsqrResult.converged, s"lsqr did not converge (residual ${lsqrResult.residual})")
    assert(cgnrResult.converged, s"cgnr did not converge (residual ${cgnrResult.residual})")

    // Both reach the same least-squares solution.
    assertVecClose(lsqrResult.x, m.leastSquares(bTrue).orThrow, 1e-5)

    assert(
      lsqrResult.iterations <= cgnrResult.iterations,
      s"expected lsqr iterations (${lsqrResult.iterations}) <= cgnr iterations (${cgnrResult.iterations})"
    )
  }
