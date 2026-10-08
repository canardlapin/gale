package gale.optim

import gale.linalg.{DMat, DVec}

final class AugmentedLagrangianSuite extends munit.FunSuite:
  private def point(x: Double*): DMat = DMat.dense(x.size, 1, x)
  private def vec(x: Double*): DVec = DVec.fromSeq(x)
  private def get[A](value: Either[FirstOrderError, A]): A = value.fold(e => fail(e.message), identity)
  private def quadratic(center: DMat): DifferentiableObjective = get(DifferentiableObjective(center.rows): x =>
    val delta = DMat.tabulate(x.rows, 1)((i, _) => x(i, 0) - center(i, 0))
    Right(ObjectiveEvaluation(0.5 * (0 until x.rows).map(i => delta(i, 0) * delta(i, 0)).sum, delta)))
  private def equality(scale: Double = 1.0): NonlinearConstraints = get(NonlinearConstraints(2, 1, 0): x =>
    Right(ConstraintEvaluation(vec(scale * (x(0, 0) + x(1, 0) - 1.0)), DMat.dense(1, 2, Seq(scale, scale)))))
  private def solve(
      f: DifferentiableObjective,
      c: NonlinearConstraints,
      x: DMat,
      config: AugmentedLagrangianConfig = AugmentedLagrangianConfig(),
      bounds: Option[MatrixBoxBounds] = None,
      warm: Option[DVec] = None
  ): AugmentedLagrangianResult =
    get(AugmentedLagrangian.minimize(f, c, x, bounds, config, warm))
  private def converged(r: AugmentedLagrangianResult): Unit =
    assertEquals(r.status, AugmentedLagrangianStatus.Converged, r.toString)
    assert(r.diagnostics.satisfies(r.settings))

  test("cancellation among large multipliers cannot erase the original objective gradient"):
    val f = get(DifferentiableObjective(1)(x => Right(ObjectiveEvaluation(x(0, 0), point(1)))))
    val c = get(NonlinearConstraints(1, 2, 0)(x => Right(ConstraintEvaluation(vec(x(0, 0), x(0, 0)), point(1, 1)))))
    val r = solve(f, c, point(0), AugmentedLagrangianConfig(maxIterations = 0), warm = Some(vec(1e16, -1e16)))
    assertEquals(r.diagnostics.stationarity, 1.0)
    assertEquals(r.status, AugmentedLagrangianStatus.IterationLimit)

  test("overflowing multiplier-Jacobian products fail even at a fixed coordinate"):
    val f = quadratic(point(0))
    val c = get(NonlinearConstraints(1, 1, 0)(x => Right(ConstraintEvaluation(vec(x(0, 0)), point(1e200)))))
    val bounds = get(MatrixBoxBounds.uniform(1, 0, 0))
    val r = AugmentedLagrangian.minimize(f, c, point(0), Some(bounds), initialMultipliers = Some(vec(1e200)))
    assert(r.isLeft)

  test("equality constrained projection has the analytic primal and multiplier"):
    val r = solve(quadratic(point(3, 2)), equality(), point(0, 0))
    converged(r)
    assertEqualsDouble(r.primal(0, 0), 1.0, 1e-7)
    assertEqualsDouble(r.primal(1, 0), 0.0, 1e-7)
    assertEqualsDouble(r.multipliers(0), 2.0, 1e-7)
    assertEqualsDouble(r.objective, 4.0, 3e-8)
    assert(math.abs(r.primal(0, 0) + r.primal(1, 0) - 1.0) <= 1e-8)

  test("one subproblem agrees with its analytic solution; only the next multiplier estimate is clipped"):
    val c = get(NonlinearConstraints(1, 1, 0)(x => Right(ConstraintEvaluation(vec(x(0, 0)), point(1)))))
    val r = solve(
      quadratic(point(3)),
      c,
      point(0),
      AugmentedLagrangianConfig(maxIterations = 1, initialPenalty = 1, multiplierBound = 1),
      warm = Some(vec(10))
    )
    // The safeguarded input lambda=1 gives min .5*(x-3)^2 + x + .5*x^2 at x=1.
    // The update is lambda+=1+x=2, retained without clipping in the reported diagnostics.
    assertEqualsDouble(r.primal(0, 0), 1.0, 1e-10)
    assertEqualsDouble(r.multipliers(0), 2.0, 1e-10)
    assertEqualsDouble(r.diagnostics.stationarity, 0.0, 1e-10)
    assertEquals(r.status, AugmentedLagrangianStatus.IterationLimit)

  test("SPD diagonal quadratic agrees with the closed form equality projection"):
    val n = 12
    val q = Array.tabulate(n)(i => math.pow(20.0, i.toDouble / (n - 1)))
    val a = Array.tabulate(n)(i => (i + 1.0) / n)
    val center = Array.tabulate(n)(i => math.sin(i + 0.2))
    val lambda = (a.indices.map(i => a(i) * center(i)).sum - 1.0) / a.indices.map(i => a(i) * a(i) / q(i)).sum
    val f = get(DifferentiableObjective(n): x =>
      val g = DMat.tabulate(n, 1)((i, _) => q(i) * (x(i, 0) - center(i)))
      Right(ObjectiveEvaluation(0.5 * (0 until n).map(i => (x(i, 0) - center(i)) * g(i, 0)).sum, g)))
    val c = get(NonlinearConstraints(n, 1, 0): x =>
      Right(ConstraintEvaluation(vec((0 until n).map(i => a(i) * x(i, 0)).sum - 1.0), DMat.dense(1, n, a.toSeq))))
    val r = solve(f, c, DMat.zeros(n, 1))
    converged(r)
    for i <- 0 until n do assertEqualsDouble(r.primal(i, 0), center(i) - lambda * a(i) / q(i), 2e-7)
    assertEqualsDouble(r.multipliers(0), lambda, 2e-7)

  test("mixed active/inactive constraints and active bounds have correct signs"):
    val c = get(NonlinearConstraints(3, 1, 2): x =>
      Right(
        ConstraintEvaluation(
          vec(x(0, 0) + x(1, 0) - 1.0, -x(1, 0), x(2, 0) - 0.5),
          DMat.dense(3, 3, Seq(1, 1, 0, 0, -1, 0, 0, 0, 1))
        )
      ))
    val bounds = get(
      MatrixBoxBounds.from(
        point(Double.NegativeInfinity, Double.NegativeInfinity, 0),
        point(Double.PositiveInfinity, Double.PositiveInfinity, Double.PositiveInfinity)
      )
    )
    val r = solve(quadratic(point(3, 0, -1)), c, point(0, 0, 0.2), bounds = Some(bounds))
    converged(r)
    assertEqualsDouble(r.primal(0, 0), 1.0, 1e-7)
    assertEqualsDouble(r.primal(1, 0), 0.0, 1e-8)
    assertEquals(r.primal(2, 0), 0.0)
    assertEqualsDouble(r.multipliers(0), 2.0, 1e-7)
    assertEqualsDouble(r.multipliers(1), 2.0, 1e-7)
    assertEquals(r.multipliers(2), 0.0)
    assertEqualsDouble(r.objective, 2.5, 3e-8)

  test("nonlinear ball projection agrees with a geometric solution"):
    val n = 8
    val a = point((1 to n).map(_.toDouble / n)*)
    val length = math.sqrt((0 until n).map(i => a(i, 0) * a(i, 0)).sum)
    val c = get(NonlinearConstraints(n, 0, 1): x =>
      Right(
        ConstraintEvaluation(
          vec((0 until n).map(i => x(i, 0) * x(i, 0)).sum - 1.0),
          DMat.tabulate(1, n)((_, j) => 2 * x(j, 0))
        )
      ))
    val r = solve(quadratic(a), c, DMat.zeros(n, 1))
    converged(r)
    for i <- 0 until n do assertEqualsDouble(r.primal(i, 0), a(i, 0) / length, 1e-7)
    assertEqualsDouble(r.multipliers(0), (length - 1.0) / 2, 1e-7)

  test("constraint scaling preserves original-unit multipliers and raw feasibility remains a gate"):
    val f = quadratic(point(3, 2))
    val r = solve(
      f,
      equality(1000),
      point(0, 0),
      AugmentedLagrangianConfig(constraintScales = Some(vec(1000)), feasibilityTolerance = 1e-5)
    )
    converged(r)
    assertEqualsDouble(r.multipliers(0), 0.002, 1e-10)
    val unscaled = solve(f, equality(), point(0, 0))
    assertEqualsDouble(r.primal(0, 0), unscaled.primal(0, 0), 1e-12)
    val blocked = solve(
      quadratic(point(0.5, 0.5 + 1e-6)),
      equality(1e6),
      point(0.5, 0.5 + 1e-6),
      AugmentedLagrangianConfig(maxIterations = 0, constraintScales = Some(vec(1e12)))
    )
    assert(blocked.diagnostics.scaledFeasibility < 1e-8)
    assert(blocked.diagnostics.feasibility > 0.9)
    assertEquals(blocked.status, AugmentedLagrangianStatus.IterationLimit)

  test("warm multipliers use original units; inactive inequalities release their multipliers"):
    val f = quadratic(point(3, 2))
    val r = solve(
      f,
      equality(1000),
      point(1, 0),
      AugmentedLagrangianConfig(constraintScales = Some(vec(1000))),
      warm = Some(vec(0.002))
    )
    converged(r)
    assertEquals(r.iterations, 0)
    val c = get(NonlinearConstraints(1, 0, 1)(x => Right(ConstraintEvaluation(vec(x(0, 0) - 10), point(1)))))
    val inactive = solve(quadratic(point(0)), c, point(3), warm = Some(vec(5)))
    converged(inactive)
    assertEqualsDouble(inactive.primal(0, 0), 0.0, 1e-7)
    assertEquals(inactive.multipliers(0), 0.0)

  test("inner stationarity at an infeasible degenerate point never means convergence"):
    val c = get(NonlinearConstraints(1, 1, 0): x =>
      Right(ConstraintEvaluation(vec(x(0, 0) * x(0, 0) - 1.0), point(2 * x(0, 0)))))
    val r = solve(quadratic(point(0)), c, point(0), AugmentedLagrangianConfig(maximumPenalty = 100))
    assertEquals(r.status, AugmentedLagrangianStatus.PenaltyLimit)
    assertEquals(r.diagnostics.feasibility, 1.0)
    assertEquals(r.diagnostics.stationarity, 0.0)

  test("global work cap counts objective and constraints separately and retains checked diagnostics"):
    for budget <- 1 to 12 do
      var objectives = 0
      var constraintCalls = 0
      val base = quadratic(point(3, 2))
      val f = get(DifferentiableObjective(2)(x => { objectives += 1; base.evaluate(x) }))
      val eq = equality()
      val c = get(NonlinearConstraints(2, 1, 0)(x => { constraintCalls += 1; eq.evaluate(x) }))
      val r = AugmentedLagrangian.minimize(
        f,
        c,
        point(0, 0),
        config = AugmentedLagrangianConfig(control = SolverControl(maxEvaluations = budget))
      )
      assertEquals(objectives + constraintCalls, budget)
      if budget == 1 then
        assertEquals(r, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.EvaluationLimit)))
      else
        val out = get(r)
        assertEquals(out.status, AugmentedLagrangianStatus.EvaluationLimit)
        assertEquals(out.evaluations.callbacks, budget)
        assertEquals(out.evaluations.values, objectives)
        assertEquals(out.evaluations.jacobians, constraintCalls)
        assertEqualsDouble(out.constraints(0), out.primal(0, 0) + out.primal(1, 0) - 1, 0.0)

  test("progress cancellation and exceptions override convergence; trace is bounded"):
    val cancelled = solve(
      quadratic(point(0.5, 0.5)),
      equality(),
      point(0.5, 0.5),
      AugmentedLagrangianConfig(control = SolverControl(progress = Some(_ => false), traceCapacity = 2))
    )
    assertEquals(cancelled.status, AugmentedLagrangianStatus.Cancelled)
    assertEquals(cancelled.trace.size, 1)
    val failed = AugmentedLagrangian.minimize(
      quadratic(point(3, 2)),
      equality(),
      point(0, 0),
      config = AugmentedLagrangianConfig(control =
        SolverControl(progress = Some(_ => throw new IllegalStateException("hook")))
      )
    )
    assert(failed.left.toOption.exists(_.message.contains("hook")))
    val r = solve(
      quadratic(point(3, 2)),
      equality(),
      point(0, 0),
      AugmentedLagrangianConfig(control = SolverControl(traceCapacity = 2))
    )
    converged(r)
    assertEquals(r.trace.size, 2)
    assertEquals(r.trace.last.iteration, r.iterations)

  test("malformed input is rejected before callbacks; invalid callbacks remain typed failures"):
    var calls = 0
    val f = get(DifferentiableObjective(2)(_ => { calls += 1; throw new IllegalArgumentException("objective") }))
    for config <- Seq(
        AugmentedLagrangianConfig(initialPenalty = 0),
        AugmentedLagrangianConfig(constraintScales = Some(vec(0))),
        AugmentedLagrangianConfig(multiplierBound = Double.NaN)
      )
    do assert(AugmentedLagrangian.minimize(f, equality(), point(0, 0), config = config).isLeft)
    assertEquals(calls, 0)
    assert(
      AugmentedLagrangian.minimize(f, equality(), point(0, 0)).left.toOption.exists(_.message.contains("objective"))
    )
    assertEquals(calls, 1)
    val malformed = get(NonlinearConstraints(2, 1, 0)(_ => Right(ConstraintEvaluation(vec(0), point(1)))))
    assert(AugmentedLagrangian.minimize(quadratic(point(3, 2)), malformed, point(0, 0)).isLeft)
    val nonfinite = get(NonlinearConstraints(2, 1, 0): _ =>
      Right(ConstraintEvaluation(vec(Double.NaN), DMat.zeros(1, 2))))
    assert(AugmentedLagrangian.minimize(quadratic(point(3, 2)), nonfinite, point(0, 0)).isLeft)

  test("no nonlinear constraints uses the smooth kernel without calling the empty oracle"):
    val empty = get(NonlinearConstraints(2, 0, 0)(_ => fail("empty oracle must not be called")))
    val r = solve(quadratic(point(3, 2)), empty, point(0, 0))
    converged(r)
    assertEquals(r.evaluations.jacobians, 0)
    assertEqualsDouble(r.primal(0, 0), 3.0, 1e-8)

  test("fixed bounds are preserved; box feasibility is checked before callbacks"):
    val b = get(MatrixBoxBounds.from(point(1, Double.NegativeInfinity), point(1, Double.PositiveInfinity)))
    val base = quadratic(point(3, 2))
    val f = get(DifferentiableObjective(2)(x => { assertEquals(x(0, 0), 1.0); base.evaluate(x) }))
    val r = solve(f, equality(), point(1, 3), bounds = Some(b))
    converged(r)
    assertEqualsDouble(r.primal(1, 0), 0.0, 1e-8)
    assert(AugmentedLagrangian.minimize(f, equality(), point(0, 0), bounds = Some(b)).isLeft)

  test("adaptive tolerances tighten at feasible nonstationary points, including empty constraints"):
    val f = quadratic(point(0, 0))
    val equalityY =
      get(NonlinearConstraints(2, 1, 0)(x => Right(ConstraintEvaluation(vec(x(1, 0)), DMat.dense(1, 2, Seq(0, 1))))))
    val empty = get(NonlinearConstraints(2, 0, 0)(_ => fail("empty oracle")))
    for c <- Seq(equalityY, empty) do
      val r = solve(
        f,
        c,
        point(1e-3, 0),
        AugmentedLagrangianConfig(adaptiveInnerTolerance = true, control = SolverControl(traceCapacity = 100))
      )
      converged(r)
      assertEqualsDouble(r.primal(0, 0), 0.0, 1e-7)
      assert(r.innerTrace.head.iterations == 0)
      val eps = r.innerTrace.map(_.tolerance)
      assert(eps.forall(e => e.isFinite && e >= r.settings.stationarityTolerance))
      assert(eps.zip(eps.drop(1)).forall((a, b) => b <= a))
      assert(eps.last < eps.head)

  test("inner trace preserves numerics and accounts for interrupted attempts and endpoint checks"):
    val f = quadratic(point(3, 2))
    for budget <- Seq(3, 4, 5, 6, 9, 12, 500) do
      val config = AugmentedLagrangianConfig(control = SolverControl(maxEvaluations = budget))
      val plain = solve(f, equality(), point(0, 0), config)
      val traced = solve(f, equality(), point(0, 0), config.copy(control = config.control.copy(traceCapacity = 100)))
      assertEquals(traced.primal.col(0).toSeq, plain.primal.col(0).toSeq)
      assertEquals(traced.diagnostics, plain.diagnostics)
      assertEquals(traced.status, plain.status)
      assertEquals(traced.evaluations, plain.evaluations)
      assertEquals(traced.innerTrace.map(_.objectiveEvaluations).sum + 1L, traced.evaluations.values.toLong)
      assertEquals(traced.innerTrace.map(_.constraintEvaluations).sum + 1L, traced.evaluations.jacobians.toLong)
      assert(traced.innerTrace.forall(t => t.augmentedEvaluations >= t.iterations + 1L))
      if budget < 10 then assertEquals(traced.status, AugmentedLagrangianStatus.EvaluationLimit)
      if budget == 4 then
        assertEquals(traced.innerTrace.last.endpointChecked, false)
        assertEquals(traced.iterations, 0)
        assertEquals(traced.primal.col(0).toSeq, Seq(0.0, 0.0))
    val bounded =
      solve(f, equality(), point(0, 0), AugmentedLagrangianConfig(control = SolverControl(traceCapacity = 2)))
    assertEquals(bounded.innerTrace.size, 2)
    assertEquals(bounded.innerTrace.last.outerIteration, bounded.iterations)

  test("curvature scale resets after penalty changes and remains local to a solve"):
    val config = AugmentedLagrangianConfig(feasibilityContraction = 1e-4, control = SolverControl(traceCapacity = 100))
    val f = quadratic(point(3, 2))
    val first = solve(f, equality(), point(0, 0), config)
    converged(first)
    val transitions = first.innerTrace.zip(first.innerTrace.drop(1)).filter((a, b) => a.penalty != b.penalty)
    assert(transitions.nonEmpty)
    transitions.foreach((_, b) => assertEquals(b.initialCurvatureScale, 1.0))
    assert(
      first.innerTrace.forall(t =>
        t.finalCurvatureScale.isFinite && t.finalCurvatureScale > 0 && (1.0 / t.finalCurvatureScale).isFinite
      )
    )
    val second = solve(f, equality(), point(0, 0), config)
    assertEquals(second.innerTrace, first.innerTrace)
    assertEquals(second.evaluations, first.evaluations)

  test("near dependent equalities retain honest diagnostics and the default reaches the analytic projection"):
    val delta = 1e-3
    val c = get(
      NonlinearConstraints(2, 2, 0)(x =>
        Right(
          ConstraintEvaluation(
            vec(x(0, 0) + delta * x(1, 0) - 1, x(0, 0) - delta * x(1, 0) - 1),
            DMat.dense(2, 2, Seq(1, delta, 1, -delta))
          )
        )
      )
    )
    for adaptive <- Seq(false, true); reuse <- Seq(false, true) do
      val r = solve(
        quadratic(point(3, 2)),
        c,
        point(0, 0),
        AugmentedLagrangianConfig(adaptiveInnerTolerance = adaptive, reuseCurvatureScale = reuse)
      )
      val stationarity = math.max(
        math.abs(r.primal(0, 0) - 3 + r.multipliers(0) + r.multipliers(1)),
        math.abs(r.primal(1, 0) - 2 + delta * (r.multipliers(0) - r.multipliers(1)))
      )
      assertEqualsDouble(r.diagnostics.stationarity, stationarity, 1e-10)
      if adaptive && reuse then converged(r)
      if r.status == AugmentedLagrangianStatus.Converged then
        assert(stationarity <= 1e-7)
        assertEqualsDouble(r.primal(0, 0), 1.0, 1e-7)
        assertEqualsDouble(r.primal(1, 0), 0.0, 1e-5)
        assertEqualsDouble(r.multipliers(0) + r.multipliers(1), 2.0, 1e-6)
        assertEqualsDouble(delta * (r.multipliers(0) - r.multipliers(1)), 2.0, 1e-5)
      else assert(!r.diagnostics.satisfies(r.settings))

  test("scaled inequalities can change activity without exhausting the zoom budget"):
    for scale <- Seq(1e4, 1e5, 1e6, 1e7); start <- Seq(1.0, 1e-4); reuse <- Seq(false, true) do
      val c = get(NonlinearConstraints(1, 0, 1)(x => Right(ConstraintEvaluation(vec(scale * x(0, 0)), point(scale)))))
      val r = solve(
        quadratic(point(-1e-4)),
        c,
        point(start),
        AugmentedLagrangianConfig(adaptiveInnerTolerance = true, reuseCurvatureScale = reuse)
      )
      converged(r)
      assertEqualsDouble(r.primal(0, 0), -1e-4, 1e-7)
      assertEquals(r.multipliers(0), 0.0)

  test("roundoff-limited subproblem acceptance is visible and does not relax original stationarity"):
    val f = get(DifferentiableObjective(2)(x =>
      val delta = x(0, 0) - 1
      // Objective evaluation loses the small quadratic beneath a large constant, with one-ulp trial noise.
      val value = 1e16 + 0.5 * delta * delta + (if x(0, 0) == 0 then 0.0 else math.ulp(1e16))
      Right(ObjectiveEvaluation(value, point(delta, 0)))
    ))
    val c =
      get(NonlinearConstraints(2, 1, 0)(x => Right(ConstraintEvaluation(vec(x(1, 0)), DMat.dense(1, 2, Seq(0, 1))))))
    val config = AugmentedLagrangianConfig(control = SolverControl(traceCapacity = 4))
    val strict = solve(f, c, point(0, 0), config.copy(roundoffAwareLineSearch = false))
    assertEquals(strict.status, AugmentedLagrangianStatus.LineSearchFailed)
    assertEquals(strict.diagnostics.stationarity, 1.0)
    assertEquals(strict.innerTrace.map(_.approximateSteps).sum, 0)
    val aware = solve(f, c, point(0, 0), config)
    converged(aware)
    assertEquals(aware.primal(0, 0), 1.0)
    assertEquals(aware.diagnostics.stationarity, 0.0)
    assertEquals(aware.evaluations.callbacks, 4)
    assertEquals(aware.innerTrace.map(_.approximateSteps).sum, 1)
    val plain = solve(f, c, point(0, 0), config.copy(control = SolverControl()))
    assertEquals(plain.evaluations, aware.evaluations)
    assertEquals(plain.diagnostics, aware.diagnostics)
    val limited = solve(f, c, point(0, 0), config.copy(control = SolverControl(maxEvaluations = 3, traceCapacity = 4)))
    assertEquals(limited.status, AugmentedLagrangianStatus.EvaluationLimit)
    assertEquals(limited.evaluations.callbacks, 3)
    assertEquals(limited.innerTrace.map(_.approximateSteps).sum, 0)
