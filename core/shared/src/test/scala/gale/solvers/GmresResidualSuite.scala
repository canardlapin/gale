package gale.solvers

import gale.linalg.*

class GmresResidualSuite extends munit.FunSuite:
  private def scaledIdentity(scale: Double): Preconditioner =
    new Preconditioner:
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        var i = 0
        while i < r.length do
          into(i) = scale * r(i)
          i += 1

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

  for mode <- ToleranceMode.values do
    for scale <- Seq(1e-12, 1.0, 1e12) do
      test(s"preconditioner scaling preserves full Krylov convergence: $mode, $scale") {
        val A = Matrix.dense(2, 2)(2.0, 0.0, 0.0, 3.0)
        val b = Vec(2.0, 3.0)
        val result = gmres(
          A,
          b,
          SolverConfig(tolerance = 1e-10, maxIterations = 2, restart = 2),
          scaledIdentity(scale),
          toleranceMode = mode
        )
        val threshold = if mode == ToleranceMode.Absolute then 1e-10 else 1e-10 * b.norm2
        check(A, b, result, threshold)
        assert(result.converged, s"residual=${result.residual}")
        assertEquals(result.iterations, 2)
      }

      test(s"restart and iteration-limit results report the true residual: $mode, $scale") {
        val A = Matrix.dense(2, 2)(2.0, 0.0, 0.0, 3.0)
        val b = Vec(2.0, 3.0)
        val config = SolverConfig(tolerance = 1e-10, maxIterations = 2, restart = 1)
        val limited = gmres(A, b, config, scaledIdentity(scale), toleranceMode = mode)
        val threshold = if mode == ToleranceMode.Absolute then 1e-10 else 1e-10 * b.norm2
        check(A, b, limited, threshold)
        assert(!limited.converged)
        assertEquals(limited.iterations, 2)
        val solved = gmres(A, b, config.copy(maxIterations = 100), scaledIdentity(scale), toleranceMode = mode)
        check(A, b, solved, threshold)
        assert(solved.converged)
        assert(solved.iterations > 2 && solved.iterations <= 100)
      }

  test("a rejected inner convergence estimate retains the Krylov basis") {
    val A = Matrix.eye(2)
    val b = Vec(1.0, 1.0)
    val preconditioner = new Preconditioner:
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        into(0) = r(0)
        into(1) = 1e-8 * r(1)
    val result = gmres(A, b, SolverConfig(tolerance = 1e-6, maxIterations = 2, restart = 2), preconditioner)
    check(A, b, result, 1e-6)
    assert(result.converged)
    assertEquals(result.iterations, 2)
  }

  test("RHS scaling does not cause premature Arnoldi breakdown") {
    val A = Matrix.dense(2, 2)(2.0, 0.0, 0.0, 3.0)
    for scale <- Seq(1e-16, 1.0, 1e16) do
      val b = Vec(2.0 * scale, 3.0 * scale)
      val result = gmres(
        A,
        b,
        SolverConfig(tolerance = 1e-10, maxIterations = 2, restart = 2),
        toleranceMode = ToleranceMode.RelativeToRhs
      )
      check(A, b, result, 1e-10 * b.norm2)
      assert(result.converged, s"rhs scale=$scale, residual=${result.residual}")
      assertEquals(result.iterations, 2)
  }

  test("an acceptable initial guess converges even when preconditioning amplifies its residual") {
    val A = Matrix.eye(1)
    val b = Vec(1.0)
    val result = gmres(
      A,
      b,
      SolverConfig(tolerance = 1e-6),
      scaledIdentity(1e12),
      initial = Some(Vec(1.0 - 1e-8))
    )
    check(A, b, result, 1e-6)
    assert(result.converged)
    assertEquals(result.iterations, 0)
  }

  test("zero iteration budget reports the initial true residual") {
    val A = Matrix.eye(1)
    val b = Vec(1.0)
    val result = gmres(A, b, SolverConfig(maxIterations = 0), scaledIdentity(1e-12))
    check(A, b, result, 1e-10)
    assert(!result.converged)
    assertEquals(result.iterations, 0)
  }

  test("relative tolerance uses the original RHS norm at the initial guess") {
    val A = Matrix.eye(1)
    val b = Vec(1000.0)
    val config = SolverConfig(tolerance = 1e-8, maxIterations = 0)
    val guess = Some(Vec(1000.0 - 1e-6))
    val absolute = gmres(A, b, config, scaledIdentity(1e-12), initial = guess)
    val relative = gmres(
      A,
      b,
      config,
      scaledIdentity(1e12),
      initial = guess,
      toleranceMode = ToleranceMode.RelativeToRhs
    )
    check(A, b, absolute, 1e-8)
    check(A, b, relative, 1e-5)
    assert(!absolute.converged)
    assert(relative.converged)
    assertEquals(relative.iterations, 0)
  }

  test("zero RHS with an exact initial guess converges in both tolerance modes") {
    val A = Matrix.eye(2)
    val b = Vec(0.0, 0.0)
    for mode <- ToleranceMode.values do
      val result = gmres(A, b, preconditioner = scaledIdentity(1e12), toleranceMode = mode)
      check(A, b, result, 0.0)
      assert(result.converged)
      assertEquals(result.iterations, 0)
  }

  test("an inconsistent system cannot report convergence at numerical breakdown") {
    val A = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 0.0)
    val b = Vec(1.0, 1.0)
    val result = gmres(A, b, SolverConfig(maxIterations = 2), scaledIdentity(1e-12))
    check(A, b, result, 1e-10)
    assert(!result.converged)
    assert(result.iterations > 0 && result.iterations <= 2)
  }

  test("a zero preconditioned residual cannot certify a nonzero true residual") {
    val A = Matrix.eye(1)
    val b = Vec(1.0)
    val result = gmres(A, b, preconditioner = scaledIdentity(0.0))
    check(A, b, result, 1e-10)
    assert(!result.converged)
    assertEquals(result.iterations, 0)
  }

  test("hard breakdown counts attempted Arnoldi steps and residual checks are bounded") {
    var applications = 0
    val A = LinearOperator.fromFunction(2, 2): (x, into) =>
      applications += 1
      into(0) = x(1)
      into(1) = 0.0
    val result = gmres(A, Vec(1.0, 0.0), preconditioner = scaledIdentity(1e-12))
    assert(!result.converged)
    assertEquals(result.residual, 1.0)
    assertEquals(result.iterations, 1)
    assert(applications <= 1 + 2 * result.iterations, s"applications=$applications")
  }

  test("operator work is bounded by one initial check and two applications per iteration") {
    var applications = 0
    val A = LinearOperator.fromFunction(2, 2): (x, into) =>
      applications += 1
      into(0) = 2.0 * x(0)
      into(1) = 3.0 * x(1)
    val result = gmres(
      A,
      Vec(2.0, 3.0),
      SolverConfig(tolerance = 1e-10, maxIterations = 100, restart = 1),
      scaledIdentity(1e-12)
    )
    assert(result.converged)
    assert(result.iterations > 2 && result.iterations <= 100)
    assertEquals(applications, 1 + 2 * result.iterations)
  }
