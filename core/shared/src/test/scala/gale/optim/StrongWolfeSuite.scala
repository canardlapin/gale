package gale.optim

import gale.linalg.DMat

final class StrongWolfeSuite extends munit.FunSuite:
  private def point(x: Double): DMat = DMat.dense(1, 1, Seq(x))
  private def quadratic(curvature: Double): DifferentiableObjective =
    DifferentiableObjective(1)(x =>
      Right(
        ObjectiveEvaluation(
          0.5 * curvature * x(0, 0) * x(0, 0) - x(0, 0),
          point(curvature * x(0, 0) - 1)
        )
      )
    ).toOption.get

  private def search(
      f: DifferentiableObjective,
      control: SolverControl = SolverControl(),
      trials: Int = 40,
      c2: Double = 0.1,
      maximum: Double = 1.0,
      interpolate: Boolean = true,
      approximateWolfe: Boolean = false,
      onApproximateAcceptance: () => Unit = StrongWolfe.ignoreApproximateAcceptance
  ): (Either[FirstOrderError, StrongWolfe.Accepted], OptimizationExecution) =
    val e = new OptimizationExecution(control)
    val initial = e.evaluate(f, point(0)).toOption.get
    (
      StrongWolfe.search(
        f,
        point(0),
        initial,
        point(1),
        e,
        1e-4,
        c2,
        trials,
        maximum,
        interpolate = interpolate,
        approximateWolfe = approximateWolfe,
        onApproximateAcceptance = onApproximateAcceptance
      ),
      e
    )

  private def checkWolfe(a: StrongWolfe.Accepted, curvature: Double, c2: Double): Unit =
    val x = a.point(0, 0)
    val value = 0.5 * curvature * x * x - x
    val derivative = curvature * x - 1
    assert(value <= -1e-4 * a.step)
    assert(math.abs(derivative) <= c2)

  test("stiff quadratic interpolation reaches analytic steps with two trial callbacks"):
    for curvature <- Seq(11.0, 1001.0, 1e6 + 1, 1e9 + 1) do
      val (r, e) = search(quadratic(curvature))
      val a = r.toOption.get
      assertEqualsDouble(a.step, 1.0 / curvature, 1e-15 / curvature)
      checkWolfe(a, curvature, 0.1)
      assertEquals(e.counts.callbacks, 3) // initial evaluation plus two line-search trials

  test("reversed zoom orientation keeps derivatives with the correct endpoints"):
    val (r, _) = search(quadratic(1.5))
    val a = r.toOption.get
    assertEqualsDouble(a.step, 2.0 / 3, 1e-15)
    checkWolfe(a, 1.5, 0.1)

  test("nonfinite high endpoints use bisection and never contaminate interpolation"):
    val base = quadratic(4)
    val f = DifferentiableObjective(1)(x =>
      if x(0, 0) > 0.4 then Right(ObjectiveEvaluation(Double.PositiveInfinity, point(Double.NaN)))
      else base.evaluate(x)
    ).toOption.get
    val (r, e) = search(f)
    checkWolfe(r.toOption.get, 4, 0.1)
    assertEquals(e.counts.callbacks, 4)

  test("zoom respects exact callback budget and cancellation during rejected trials"):
    val (limited, e) = search(quadratic(1e6), SolverControl(maxEvaluations = 2))
    assertEquals(limited, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.EvaluationLimit)))
    assertEquals(e.counts.callbacks, 2)
    var calls = 0
    val base = quadratic(1e6)
    val f = DifferentiableObjective(1)(x => { calls += 1; base.evaluate(x) }).toOption.get
    val (cancelled, work) = search(f, SolverControl(cancelled = () => calls >= 2))
    assertEquals(cancelled, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled)))
    assertEquals(work.counts.callbacks, 2)

  test("degenerate interpolation and exhausted trials return explicit failure"):
    // Inconsistent flat values and fixed gradient cannot pass sufficient decrease.
    val f = DifferentiableObjective(1)(_ => Right(ObjectiveEvaluation(0.0, point(-1)))).toOption.get
    val (r, e) = search(f, trials = 7)
    assertEquals(r, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed)))
    assertEquals(e.counts.callbacks, 8)

  test("public smooth solvers retain the bisection trial sequence"):
    val f = quadratic(11)
    val bounds = MatrixBoxBounds.from(point(Double.NegativeInfinity), point(Double.PositiveInfinity)).toOption.get
    val results = Seq(
      LBFGS.minimize(f, point(0), LBFGSConfig(maxIterations = 1)),
      LBFGSB.minimize(f, bounds, point(0), LBFGSBConfig(maxIterations = 1))
    )
    results.foreach: r =>
      val fit = r.toOption.get
      assertEquals(fit.primal(0, 0), 0.125)
      assertEquals(fit.evaluations.callbacks, 5)

  test("both zoom policies preserve nonfinite fallback, cancellation, and exact budgets"):
    for interpolate <- Seq(false, true) do
      val base = quadratic(4)
      val nonfinite = DifferentiableObjective(1)(x =>
        if x(0, 0) > 0.4 then Right(ObjectiveEvaluation(Double.PositiveInfinity, point(Double.NaN)))
        else base.evaluate(x)
      ).toOption.get
      val (finite, finiteWork) = search(nonfinite, interpolate = interpolate)
      checkWolfe(finite.toOption.get, 4, 0.1)
      assertEquals(finiteWork.counts.callbacks, 4)
      val (limited, work) = search(quadratic(1e6), SolverControl(maxEvaluations = 2), interpolate = interpolate)
      assertEquals(limited, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.EvaluationLimit)))
      assertEquals(work.counts.callbacks, 2)
      var calls = 0
      val f = DifferentiableObjective(1)(x => { calls += 1; base.evaluate(x) }).toOption.get
      val (cancelled, cancelledWork) = search(f, SolverControl(cancelled = () => calls >= 2), interpolate = interpolate)
      assertEquals(cancelled, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled)))
      assertEquals(cancelledWork.counts.callbacks, 2)

  test("roundoff acceptance recovers analytic minima across objective scales and offsets"):
    for offset <- Seq(1.0, -1.0, 1e-200, -1e-200, 1e200, -1e200) do
      val ulp = math.ulp(offset)
      val scale = 0.25 * ulp
      // One-ulp evaluation noise at trials, with the exact gradient of the underlying quadratic.
      val f = DifferentiableObjective(1)(x =>
        val delta = x(0, 0) - 1
        val noise = if x(0, 0) == 0 then 0.0 else ulp
        Right(ObjectiveEvaluation(offset + 0.5 * scale * delta * delta + noise, point(scale * delta)))
      ).toOption.get
      val (strict, _) = search(f, trials = 8)
      assertEquals(strict, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed)))
      var accepted = 0
      val (aware, counts) = search(f, approximateWolfe = true, onApproximateAcceptance = () => accepted += 1)
      assertEquals(aware.toOption.get.point(0, 0), 1.0)
      assertEquals(aware.toOption.get.directionalDerivative, 0.0)
      assertEquals(accepted, 1)
      assertEquals(counts.counts.callbacks, 2)

  test("roundoff acceptance rejects resolved increases and insufficient curvature"):
    val ulp = math.ulp(1.0)
    for (increase, gradient) <- Seq(
        (32 * ulp, (x: Double) => 0.25 * ulp * (x - 1)),
        (ulp, (_: Double) => -0.25 * ulp),
        (ulp, (x: Double) => -0.25 * ulp * (x + 1)),
        (ulp, (x: Double) => x - 1)
      )
    do
      val f = DifferentiableObjective(1)(x =>
        Right(ObjectiveEvaluation(if x(0, 0) == 0 then 1.0 else 1.0 + increase, point(gradient(x(0, 0)))))
      ).toOption.get
      var accepted = 0
      val (r, _) = search(f, trials = 8, approximateWolfe = true, onApproximateAcceptance = () => accepted += 1)
      assertEquals(r, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed)))
      assertEquals(accepted, 0)

  test("ordinary Wolfe remains preferred and approximate steps obey interruption budgets"):
    var accepted = 0
    val (ordinary, _) = search(quadratic(11), approximateWolfe = true, onApproximateAcceptance = () => accepted += 1)
    checkWolfe(ordinary.toOption.get, 11, 0.1)
    assertEquals(accepted, 0)
    val f = DifferentiableObjective(1)(x =>
      Right(ObjectiveEvaluation(if x(0, 0) == 0 then 1.0 else Math.nextUp(1.0), point(math.ulp(1.0) * (x(0, 0) - 1))))
    ).toOption.get
    val (limited, work) = search(
      f,
      SolverControl(maxEvaluations = 1),
      approximateWolfe = true,
      onApproximateAcceptance = () => accepted += 1
    )
    assertEquals(limited, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.EvaluationLimit)))
    assertEquals(work.counts.callbacks, 1)
    assertEquals(accepted, 0)
    var calls = 0
    val cancelledF = DifferentiableObjective(1)(x => { calls += 1; f.evaluate(x) }).toOption.get
    val (cancelled, _) = search(cancelledF, SolverControl(cancelled = () => calls >= 1), approximateWolfe = true)
    assertEquals(cancelled, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled)))
    assertEquals(calls, 1)

  test("approximate acceptance cannot use the feasible-boundary curvature exception"):
    val scale = math.ulp(1.0)
    val f = DifferentiableObjective(1)(x =>
      Right(ObjectiveEvaluation(if x(0, 0) == 0 then 1.0 else Math.nextUp(1.0), point(scale * (x(0, 0) - 1))))
    ).toOption.get
    val execution = new OptimizationExecution(SolverControl())
    val initial = execution.evaluate(f, point(0)).toOption.get
    val r = StrongWolfe.search(
      f,
      point(0),
      initial,
      point(1),
      execution,
      1e-4,
      0.1,
      8,
      0.25,
      feasibleBoundary = true,
      interpolate = true,
      approximateWolfe = true
    )
    assertEquals(r, Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.LineSearchFailed)))
