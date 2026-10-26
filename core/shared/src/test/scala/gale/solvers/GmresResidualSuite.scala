package gale.solvers

import gale.linalg.*

class GmresResidualSuite extends munit.FunSuite:
  private def check(a: DMat, b: DVec, result: SolverResult, tolerance: Double): Unit =
    val actual = (b - a * result.x).norm2
    assertEqualsDouble(result.residual, actual, 1e-14 * math.max(1.0, actual))
    if result.converged then assert(actual <= tolerance, s"false convergence: $actual > $tolerance")

  for mode <- ToleranceMode.values do
    test(s"$mode: scaled Jacobi cannot accept an unsolved initial point") {
      for scale <- Seq(1e-12, 1.0, 1e12) do
        val a = Matrix(1, 1)(scale)
        val b = Vec(1.0)
        val result = gmres(a, b, SolverConfig(tolerance = 1e-10), Preconditioner.Jacobi(a), toleranceMode = mode)
        check(a, b, result, 1e-10)
        assert(result.converged)
        assertEquals(result.iterations, 1)
        assertEqualsDouble(result.x(0) * scale, 1.0, 1e-14)
    }

  test("restart cycles and exhausted budgets report the original-system residual") {
    val a = Matrix(2, 2)(4.0, 1.0, 1.0, 3.0)
    val b = Vec(1.0, 2.0)
    for budget <- Seq(0, 1, 2, 100) do
      val result = gmres(a, b, SolverConfig(tolerance = 1e-10, maxIterations = budget, restart = 1), Preconditioner.Jacobi(a))
      check(a, b, result, 1e-10)
      assert(result.iterations <= budget)
      if budget == 100 then assert(result.converged)
      if budget == 0 then assert(!result.converged)
  }

  test("annihilating and nonfinite preconditioners cannot claim success") {
    val a = Matrix.eye(2)
    val b = Vec(1.0, 2.0)
    for value <- Seq(0.0, Double.NaN, Double.PositiveInfinity) do
      val bad = new Preconditioner:
        def solve(r: DVec, into: MutableVec[Double]): Unit =
          for i <- 0 until r.length do into(i) = value
      val result = gmres(a, b, preconditioner = bad)
      assert(!result.converged)
      check(a, b, result, 1e-10)
  }

  test("exact initial guesses need only the initial verification matvec") {
    var calls = 0
    val op = LinearOperator.fromFunction(2, 2) { (x, into) =>
      calls += 1
      into(0) = x(0)
      into(1) = x(1)
    }
    val result = gmres(op, Vec(1.0, 2.0), initial = Some(Vec(1.0, 2.0)))
    assert(result.converged)
    assertEquals(result.iterations, 0)
    assertEquals(calls, 1)
  }

  test("one-step solve counts initial, Arnoldi and terminal verification matvecs") {
    var calls = 0
    val op = LinearOperator.fromFunction(1, 1) { (x, into) =>
      calls += 1
      into(0) = 1e12 * x(0)
    }
    val result = gmres(op, Vec(1.0), preconditioner = Preconditioner.Jacobi(Matrix(1, 1)(1e12)))
    assert(result.converged)
    assertEquals(result.iterations, 1)
    assertEquals(calls, 3)
  }

  test("happy breakdown is invariant to a large right-hand-side scale") {
    val a = Matrix(2, 2)(4.0, 1.0, 1.0, 3.0)
    val b = Vec(1e18, 2e18)
    val result = gmres(a, b, SolverConfig(tolerance = 1e-12, maxIterations = 2, restart = 2),
      toleranceMode = ToleranceMode.RelativeToRhs)
    check(a, b, result, 1e-12 * b.norm2)
    assert(result.converged)
  }

  test("nonfinite preconditioning during Arnoldi retains honest failure diagnostics") {
    var calls = 0
    val bad = new Preconditioner:
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        calls += 1
        for i <- 0 until r.length do into(i) = if calls == 1 then r(i) else Double.NaN
    val a = Matrix.eye(2)
    val b = Vec(1.0, 2.0)
    val result = gmres(a, b, preconditioner = bad)
    assert(!result.converged)
    check(a, b, result, 1e-10)
  }
