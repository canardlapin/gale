package gale.optim

import gale.linalg.DMat

class AcceleratedProximalSuite extends munit.FunSuite:
  private val tolerance = FirstOrderTolerance.from(1e-7, 0.0).toOption.get
  private def quadratic(diagonal: Vector[Double], center: Vector[Double]): DifferentiableObjective =
    DifferentiableObjective(diagonal.size)(at =>
      val gradient = DMat.tabulate(at.rows, at.cols)((r, c) => diagonal(r) * (at(r, c) - center(r)))
      var f = 0.0
      var r = 0
      while r < at.rows do
        var c = 0
        while c < at.cols do
          f += .5 * diagonal(r) * Math.pow(at(r, c) - center(r), 2)
          c += 1
        r += 1
      Right(ObjectiveEvaluation(f, gradient))
    ).toOption.get

  test("backtracking accelerated L1 solves a heterogeneous diagonal KKT system"):
    val d = Vector(1.0, 10.0, 100.0)
    val center = Vector(2.0, -.01, -1.0)
    val term = ProximalTerms.l1(3, .25).toOption.get
    val solution = AcceleratedProximal
      .minimize(
        quadratic(d, center),
        term,
        DMat.zeros(3, 2),
        AcceleratedConfig(tolerance = tolerance, initialLipschitz = .1)
      )
      .toOption
      .get
    assertEquals(solution.status, FirstOrderStoppingStatus.Converged)
    for r <- 0 until 3; c <- 0 until 2 do
      val expected = Math.signum(center(r)) * Math.max(0.0, Math.abs(center(r)) - .25 / d(r))
      assertEqualsDouble(solution.primal(r, c), expected, 2e-7)
    assert(solution.certificate.binds(solution.primal, None))
    assert(solution.primalStep <= .01)
    assert(solution.evaluations.gradients > solution.certificate.iterations)

  test("convex projection uses feasible initial points and returns a feasible optimum"):
    val box = ProjectionSets.box(2, -.5, .5).toOption.get
    val objective = quadratic(Vector(1.0, 100.0), Vector(2.0, -1.0))
    val solution = AcceleratedProximal
      .minimize(objective, box.indicator, DMat.zeros(2, 1), AcceleratedConfig(tolerance = tolerance))
      .toOption
      .get
    assertEquals(solution.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solution.primal(0, 0), .5, 1e-8)
    assertEqualsDouble(solution.primal(1, 0), -.5, 1e-8)
    assert(AcceleratedProximal.minimize(objective, box.indicator, DMat.tabulate(2, 1)((_, _) => 2.0)).isLeft)

  test("budget and cancellation preserve the last checked point and count exact callbacks"):
    val objective = quadratic(Vector(1.0, 10.0), Vector(2.0, -1.0))
    val term = ProximalTerms.zero(2).toOption.get
    val limited = AcceleratedProximal
      .minimize(
        objective,
        term,
        DMat.zeros(2, 1),
        AcceleratedConfig(initialLipschitz = 10.0, control = SolverControl(maxEvaluations = 5))
      )
      .toOption
      .get
    assertEquals(limited.status, FirstOrderStoppingStatus.EvaluationLimit)
    assertEquals(limited.evaluations.callbacks, 5)
    assert(limited.certificate.binds(limited.primal, None))
    val stopped = AcceleratedProximal
      .minimize(
        objective,
        term,
        DMat.zeros(2, 1),
        AcceleratedConfig(control = SolverControl(progress = Some(_ => false), traceCapacity = 2))
      )
      .toOption
      .get
    assertEquals(stopped.status, FirstOrderStoppingStatus.Cancelled)
    assertEquals(stopped.trace.size, 1)

  test("gradient checks expose an incorrect analytic callback without solving"):
    val objective = quadratic(Vector(1.0, 3.0), Vector(2.0, 1.0))
    val at = DMat.zeros(2, 1)
    val direction = DMat.tabulate(2, 1)((_, _) => 1.0)
    assert(GradientCheck.directional(objective, at, direction).toOption.get.passed)
    val wrong =
      DifferentiableObjective(2)(x => objective.evaluate(x).map(_.copy(gradient = DMat.zeros(2, 1)))).toOption.get
    assert(!GradientCheck.directional(wrong, at, direction).toOption.get.passed)

  test("an underestimated initial Lipschitz estimate cannot certify a tiny constrained mapping"):
    val objective = quadratic(Vector(1.0), Vector(2.0))
    val box = ProjectionSets.box(1, -1.0, 1.0).toOption.get
    val solution = AcceleratedProximal
      .minimize(
        objective,
        box.indicator,
        DMat.zeros(1, 1),
        AcceleratedConfig(initialLipschitz = 1e-12, maxBacktracks = 60, tolerance = tolerance)
      )
      .toOption
      .get
    assertEquals(solution.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solution.primal(0, 0), 1.0, 1e-9)

  test("throwing progress hooks stay inside the typed failure boundary"):
    val solution = AcceleratedProximal.minimize(
      quadratic(Vector(1.0), Vector(2.0)),
      ProximalTerms.zero(1).toOption.get,
      DMat.zeros(1, 1),
      AcceleratedConfig(control = SolverControl(progress = Some(_ => throw new IllegalStateException("hook"))))
    )
    assert(solution.left.toOption.exists(_.isInstanceOf[FirstOrderError.OracleFailure]))

  test("objective offsets and affine objectives cannot shrink constrained stopping residuals through huge steps"):
    val box = ProjectionSets.box(1, -1.0, 1.0).toOption.get
    for curvature <- Vector(0.0, 1.0) do
      val objective = DifferentiableObjective(1)(x =>
        Right(
          ObjectiveEvaluation(
            1e16 + .5 * curvature * x(0, 0) * x(0, 0) - 2.0 * x(0, 0),
            DMat.tabulate(1, 1)((_, _) => curvature * x(0, 0) - 2.0)
          )
        )
      ).toOption.get
      val solved = AcceleratedProximal
        .minimize(
          objective,
          box.indicator,
          DMat.zeros(1, 1),
          AcceleratedConfig(initialLipschitz = 1e-12, maxBacktracks = 60, tolerance = tolerance)
        )
        .toOption
        .get
      assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
      assertEqualsDouble(solved.primal(0, 0), 1.0, 1e-8)
      assert(solved.primalResidualStep <= 1.0)

  test("initial checks require no trial objective evaluations"):
    val objective = quadratic(Vector(2.0), Vector(2.0))
    val term = ProximalTerms.zero(1).toOption.get
    for cancel <- Vector(false, true) do
      val result = AcceleratedProximal
        .minimize(
          objective,
          term,
          DMat.zeros(1, 1),
          AcceleratedConfig(
            maxIterations = if cancel then 10 else 0,
            initialLipschitz = .01,
            control = SolverControl(progress = if cancel then Some(_ => false) else None)
          )
        )
        .toOption
        .get
      assertEquals(
        result.status,
        if cancel then FirstOrderStoppingStatus.Cancelled else FirstOrderStoppingStatus.IterationLimit
      )
      assertEquals(result.evaluations.callbacks, 3)
      assertEquals(result.evaluations.gradients, 1)
      assertEqualsDouble(result.primal(0, 0), 0.0, 0.0)
      assertEqualsDouble(result.certificate.primalResidual, 4.0, 0.0)
      assertEqualsDouble(result.certificate.settings.primalResidualScale, 4.0, 0.0)

  test("a rounded capped check cannot prevent a representable large-step update"):
    val solved = AcceleratedProximal
      .minimize(
        quadratic(Vector(1e-20), Vector(0.0)),
        ProximalTerms.zero(1).toOption.get,
        DMat.tabulate(1, 1)((_, _) => 1.0),
        AcceleratedConfig(initialLipschitz = 1e-20, tolerance = FirstOrderTolerance.from(1e-22, 0.0).toOption.get)
      )
      .toOption
      .get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.primal(0, 0), 0.0, 1e-15)
    assertEquals(solved.certificate.iterations, 1)

  test("backtracking freezes the initial capped reference and reports the final-point map"):
    val solved = AcceleratedProximal
      .minimize(
        quadratic(Vector(10.0), Vector(2.0)),
        ProximalTerms.zero(1).toOption.get,
        DMat.zeros(1, 1),
        AcceleratedConfig(initialLipschitz = .1, tolerance = tolerance)
      )
      .toOption
      .get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.certificate.settings.primalResidualScale, 20.0, 0.0)
    val step = solved.primalResidualStep
    val x = solved.primal(0, 0)
    val mapped = x - step * 10.0 * (x - 2.0)
    assertEqualsDouble(solved.certificate.primalResidual, Math.abs(x - mapped) / step, 0.0)
    assert(Math.abs(10.0 * (x - 2.0)) <= 1e-7)

  test("both restart policies preserve box KKT conditions on free and active coordinates"):
    val d = Vector(1.0, 25.0, 1000.0)
    val center = Vector(.4, 2.0, .36)
    val box = ProjectionSets.box(3, -.5, .5).toOption.get
    for restart <- Vector(false, true) do
      val solved = AcceleratedProximal
        .minimize(
          quadratic(d, center),
          box.indicator,
          DMat.zeros(3, 1),
          AcceleratedConfig(tolerance = tolerance, restart = restart)
        )
        .toOption
        .get
      assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
      for i <- 0 until 3 do
        val x = solved.primal(i, 0)
        val projected = Math.max(-.5, Math.min(.5, x - d(i) * (x - center(i))))
        assert(Math.abs(x - projected) <= 2e-7)
        assert(x >= -.5 && x <= .5)
