package gale.optim

import gale.linalg.DMat

final class LBFGSBSuite extends munit.FunSuite:
  private def column(x: Double*): DMat = DMat.fromArrayRowMajor(x.size, 1, x.toArray)
  private val tolerance = FirstOrderTolerance.from(1e-8, 0.0).toOption.get
  private val config = LBFGSBConfig(tolerance = tolerance)
  private def quadratic(h: Vector[Vector[Double]], center: Vector[Double]): DifferentiableObjective =
    DifferentiableObjective(center.size): x =>
      val gradient = Array.tabulate(center.size)(i => center.indices.map(j => h(i)(j) * (x(j, 0) - center(j))).sum)
      val value = center.indices.map(i => 0.5 * (x(i, 0) - center(i)) * gradient(i)).sum
      Right(ObjectiveEvaluation(value, column(gradient.toSeq*)))
    .toOption.get

  test("rotated SPD box problem recovers analytic active and free coordinates"):
    val objective = quadratic(Vector(Vector(4.0, 1.0), Vector(1.0, 2.0)), Vector(2.0, 0.2))
    val bounds = MatrixBoxBounds.uniform(2, 0.0, 1.0).toOption.get
    for start <- Vector(column(0.0, 0.9), column(1.0, 0.0), column(0.1, 1.0)) do
      val solved = LBFGSB.minimize(objective, bounds, start, config).toOption.get
      assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
      assertEqualsDouble(solved.primal(0, 0), 1.0, 1e-8)
      assertEqualsDouble(solved.primal(1, 0), 0.7, 1e-8)
      assert(bounds.contains(solved.primal))
      assertEqualsDouble(solved.primalResidualStep, 1.0, 0.0)
      assert(solved.certificate.binds(solved.primal, None))

  test("linear objective accepts an Armijo step at the actual feasible endpoint"):
    var calls = 0
    val objective = DifferentiableObjective(1): x =>
      calls += 1
      assert(x(0, 0) >= 0.0 && x(0, 0) <= 1.0)
      Right(ObjectiveEvaluation(-x(0, 0), column(-1.0)))
    .toOption.get
    val solved =
      LBFGSB.minimize(objective, MatrixBoxBounds.uniform(1, 0.0, 1.0).toOption.get, column(0.25), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.primal(0, 0), 1.0, 0.0)
    assertEquals(solved.evaluations.callbacks, calls)

  test("fixed, infinite, and simultaneous bounds remain feasible at every callback"):
    val bounds =
      MatrixBoxBounds
        .from(column(0.5, Double.NegativeInfinity, 0.0), column(0.5, Double.PositiveInfinity, 1.0))
        .toOption
        .get
    val base =
      quadratic(Vector(Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0)), Vector(10.0, -2.0, 3.0))
    val objective = DifferentiableObjective(3)(x =>
      assert(bounds.contains(x))
      base.evaluate(x)
    ).toOption.get
    val solved = LBFGSB.minimize(objective, bounds, column(0.5, 0.0, 0.0), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.primal(0, 0), 0.5, 0.0)
    assertEqualsDouble(solved.primal(1, 0), -2.0, 1e-8)
    assertEqualsDouble(solved.primal(2, 0), 1.0, 0.0)
    val together = LBFGSB.cauchy(
      Array(0.0, 0.0),
      Array(-2.0, -2.0),
      Array(0.0, 0.0),
      Array(1.0, 1.0),
      LBFGSB.Model(1.0, Vector.empty)
    )
    assertEquals(together.toVector, Vector(1.0, 1.0))

  test("projected residual avoids false stationarity at translated points"):
    val x = Array(1e16)
    val g = Array(1.0)
    val residual = LBFGSB.projected(x, g, Array(Double.NegativeInfinity), Array(Double.PositiveInfinity))
    assertEqualsDouble(residual(0), 1.0, 0.0)
    val objective = DifferentiableObjective(1)(x => Right(ObjectiveEvaluation(x(0, 0), column(1.0)))).toOption.get
    val solved = LBFGSB
      .minimize(
        objective,
        MatrixBoxBounds.uniform(1, Double.NegativeInfinity, Double.PositiveInfinity).toOption.get,
        column(1e16),
        config.copy(maxIterations = 2)
      )
      .toOption
      .get
    assert(solved.status != FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.certificate.primalResidual, 1.0, 0.0)

  test("limited-memory Hessian agrees with independent dense BFGS after history eviction"):
    val pairs = Vector(
      LBFGSB.Pair(Array(1.0, 0.2, -0.1), Array(2.0, 0.3, -0.2)),
      LBFGSB.Pair(Array(0.1, 1.0, 0.2), Array(0.2, 3.0, 0.3)),
      LBFGSB.Pair(Array(-0.2, 0.1, 1.0), Array(-0.1, 0.2, 4.0))
    )
    for history <- Vector(pairs, pairs.takeRight(2), pairs.takeRight(1)) do
      val last = history.last
      val theta = last.y.map(v => v * v).sum / last.s.zip(last.y).map(_ * _).sum
      val dense = Array.tabulate(3, 3)((i, j) => if i == j then theta else 0.0)
      history.foreach: pair =>
        val bs = Array.tabulate(3)(i => (0 until 3).map(j => dense(i)(j) * pair.s(j)).sum)
        val denom = pair.s.zip(bs).map(_ * _).sum
        val curvature = pair.s.zip(pair.y).map(_ * _).sum
        for i <- 0 until 3; j <- 0 until 3 do dense(i)(j) += pair.y(i) * pair.y(j) / curvature - bs(i) * bs(j) / denom
      val v = Array(0.4, -0.7, 1.3)
      val actual = LBFGSB.model(history)(v)
      for i <- 0 until 3 do
        val expected = (0 until 3).map(j => dense(i)(j) * v(j)).sum
        assertEqualsDouble(actual(i), expected, 1e-12)

  test("generalized Cauchy point finds the first piecewise quadratic minimum"):
    // B = [[2,1],[1,2]], gradient [-4,-1], first coordinate hits its upper bound at t=1/8.
    // Along that segment the model still decreases; after fixing x0=.5, x1 minimizes at .25.
    val model = LBFGSB.Model(1.0, Vector(Array(1.0, 1.0) -> 1.0))
    val point = LBFGSB.cauchy(Array(0.0, 0.0), Array(-4.0, -1.0), Array(0.0, 0.0), Array(0.5, 2.0), model)
    assertEqualsDouble(point(0), 0.5, 1e-14)
    assertEqualsDouble(point(1), 0.25, 1e-14)

  test("restricted Hessian products match the full operator on free coordinates and reuse scratch"):
    val b = LBFGSB.model(
      Vector(
        LBFGSB.Pair(Array(1.0, 0.2, -0.1, 0.3), Array(2.0, 0.3, -0.2, 0.7)),
        LBFGSB.Pair(Array(0.1, 1.0, 0.2, -0.3), Array(0.2, 3.0, 0.3, -0.8))
      )
    )
    for free <- Vector(Array.emptyIntArray, Array(1), Array(0, 2), Array(0, 1, 2, 3)) do
      val scratch = new Array[Double](4)
      for factor <- Vector(1.0, -2.0) do
        val direction = Array.tabulate(4)(i => if free.contains(i) then factor * (i + 0.25) else 0.0)
        val full = b(direction)
        b.restricted(direction, free, scratch)
        for i <- 0 until 4 do assertEqualsDouble(scratch(i), if free.contains(i) then full(i) else 0.0, 0.0)

  test("breakpoint cancellation preserves smaller still-free directions"):
    for flip <- Vector(false, true) do
      val g = if flip then Array(-1.0, -1e8) else Array(-1e8, -1.0)
      val u = if flip then Array(2.0, 1.0) else Array(1.0, 2.0)
      val point = LBFGSB.cauchy(Array(0.0, 0.0), g, Array(0.0, 0.0), u, LBFGSB.Model(1.0, Vector.empty))
      assertEqualsDouble(point(0), 1.0, 1e-14)
      assertEqualsDouble(point(1), 1.0, 1e-14)

  test("history stays bounded and unconstrained Rosenbrock remains a valid box problem"):
    val objective = DifferentiableObjective(2): x =>
      val a = x(0, 0)
      val b = x(1, 0)
      val r = b - a * a
      Right(
        ObjectiveEvaluation(100.0 * r * r + (1.0 - a) * (1.0 - a), column(-400.0 * a * r - 2.0 * (1.0 - a), 200.0 * r))
      )
    .toOption.get
    val solved = LBFGSB
      .minimize(
        objective,
        MatrixBoxBounds.uniform(2, -2.0, 2.0).toOption.get,
        column(-1.2, 1.0),
        config.copy(historySize = 3)
      )
      .toOption
      .get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.primal(0, 0), 1.0, 1e-6)
    assertEqualsDouble(solved.primal(1, 0), 1.0, 1e-6)

  test("line search uses derivatives to resolve exact objective ties"):
    val objective = DifferentiableObjective(1)(x =>
      Right(ObjectiveEvaluation(1000.0 + 2.0 * x(0, 0) * x(0, 0), column(4.0 * x(0, 0))))
    ).toOption.get
    val solved = LBFGSB
      .minimize(
        objective,
        MatrixBoxBounds.uniform(1, -1.0, 1.0).toOption.get,
        column(1e-8),
        config.copy(tolerance = FirstOrderTolerance.from(1e-12, 0.0).toOption.get)
      )
      .toOption
      .get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.primal(0, 0), 0.0, 1e-12)

  test("budgets, cancellation, invalid inputs and callback failures are explicit"):
    val objective = quadratic(Vector(Vector(2.0)), Vector(0.2))
    val box = MatrixBoxBounds.uniform(1, 0.0, 1.0).toOption.get
    val limited = LBFGSB
      .minimize(objective, box, column(0.8), config.copy(control = SolverControl(maxEvaluations = 1)))
      .toOption
      .get
    assertEquals(limited.status, FirstOrderStoppingStatus.EvaluationLimit)
    assertEqualsDouble(limited.primal(0, 0), 0.8, 0.0)
    val cancelled = LBFGSB
      .minimize(objective, box, column(0.8), config.copy(control = SolverControl(progress = Some(_ => false))))
      .toOption
      .get
    assertEquals(cancelled.status, FirstOrderStoppingStatus.Cancelled)
    assert(LBFGSB.minimize(objective, box, column(2.0), config).isLeft)
    assert(MatrixBoxBounds.uniform(1, Double.NaN, 1.0).isLeft)
    assert(MatrixBoxBounds.uniform(1, 2.0, 1.0).isLeft)
