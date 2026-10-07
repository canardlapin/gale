package gale.optim

import gale.linalg.DMat
import munit.FunSuite

final class SmoothOptimizationSuite extends FunSuite:
  private val tolerance = FirstOrderTolerance.from(1e-6, 0.0).toOption.get

  test("L-BFGS solves a rotated SPD quadratic to an independently checked gradient"):
    val h = DMat.tabulate(2, 2): (r, c) =>
      (r, c) match
        case (0, 0) | (1, 1) => 500.5
        case _               => -499.5
    val center = column(0.4, -0.7)
    val objective = quadratic(h, center)
    val solved =
      LBFGS.minimize(objective, DMat.zeros(2, 1), LBFGSConfig(maxIterations = 2000, tolerance = tolerance)).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assert(normInf(gradient(objective, solved.primal)) <= 1e-6)
    assertEqualsDouble(solved.primal(0, 0), center(0, 0), 1e-5)
    assertEqualsDouble(solved.primal(1, 0), center(1, 0), 1e-5)

  test("L-BFGS reaches stationary Rosenbrock endpoints in dimensions 2 and 32"):
    for n <- Vector(2, 32) do
      val start = DMat.tabulate(n, 1)((r, _) => if r % 2 == 0 then -1.2 else 1.0)
      val objective = rosenbrock(n)
      val solved = LBFGS
        .minimize(objective, start, LBFGSConfig(maxIterations = 20000, tolerance = tolerance, historySize = 12))
        .toOption
        .get
      assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
      assert(normInf(gradient(objective, solved.primal)) <= 1e-6)

  test("strong-Wolfe accepted step independently satisfies sufficient decrease and curvature"):
    val objective = quadratic(DMat.eye(1), column(3.0))
    val initial = DMat.zeros(1, 1)
    val execution = new OptimizationExecution(SolverControl())
    val evaluation = execution.evaluate(objective, initial).toOption.get
    val direction = scale(evaluation.gradient, -1.0)
    val accepted =
      StrongWolfe.search(objective, initial, evaluation, direction, execution, 1e-4, .9, 40, 100.0).toOption.get
    val initialDerivative = dot(evaluation.gradient, direction)
    assert(accepted.evaluation.value <= evaluation.value + 1e-4 * accepted.step * initialDerivative)
    assert(Math.abs(dot(accepted.evaluation.gradient, direction)) <= -.9 * initialDerivative)

  test("flat objective reports line-search failure rather than a false convergence"):
    val flat = DifferentiableObjective(1)(at => Right(ObjectiveEvaluation(at(0, 0), column(1.0))))
    val result = LBFGS.minimize(flat.toOption.get, DMat.zeros(1, 1), LBFGSConfig(maxIterations = 10, maxLineSearch = 8))
    assertEquals(result.toOption.map(_.status), Some(FirstOrderStoppingStatus.LineSearchFailed))

  test("negative-curvature objective does not manufacture a BFGS curvature pair"):
    val concave = DifferentiableObjective(1): at =>
      val x = at(0, 0)
      Right(ObjectiveEvaluation(-0.5 * x * x, column(-x)))
    val result = LBFGS.minimize(concave.toOption.get, column(1.0), LBFGSConfig(maxIterations = 8, maxLineSearch = 8))
    assertEquals(result.toOption.map(_.status), Some(FirstOrderStoppingStatus.LineSearchFailed))

  test("evaluation budget and cancellation retain the last accepted finite point"):
    val objective = quadratic(DMat.eye(1), column(2.0))
    val budgeted = LBFGS
      .minimize(
        objective,
        DMat.zeros(1, 1),
        LBFGSConfig(maxIterations = 10, control = SolverControl(maxEvaluations = 1))
      )
      .toOption
      .get
    assertEquals(budgeted.status, FirstOrderStoppingStatus.EvaluationLimit)
    assertEqualsDouble(budgeted.primal(0, 0), 0.0, 0.0)
    var checks = 0
    val cancelled = LBFGS
      .minimize(
        objective,
        DMat.zeros(1, 1),
        LBFGSConfig(maxIterations = 10, control = SolverControl(cancelled = () => { checks += 1; checks >= 2 }))
      )
      .toOption
      .get
    assertEquals(cancelled.status, FirstOrderStoppingStatus.Cancelled)
    assert(cancelled.primal(0, 0).isFinite)

  test("results remain stable across later solves and directional gradient check is independent"):
    val objective = quadratic(DMat.eye(2), column(1.0, -2.0))
    val first = LBFGS.minimize(objective, DMat.zeros(2, 1), LBFGSConfig(tolerance = tolerance)).toOption.get
    val snapshot = DMat.tabulate(2, 1)((r, c) => first.primal(r, c))
    LBFGS.minimize(objective, column(20.0, 30.0), LBFGSConfig(tolerance = tolerance)).toOption.get
    for r <- 0 until snapshot.rows; c <- 0 until snapshot.cols do
      assertEqualsDouble(first.primal(r, c), snapshot(r, c), 0.0)
    val check = GradientCheck.directional(objective, column(0.2, -0.4), column(0.3, 0.7)).toOption.get
    assert(check.passed)

  private def rosenbrock(n: Int): DifferentiableObjective =
    DifferentiableObjective(n) { at =>
      val x = (0 until n).map(i => at(i, 0)).toArray
      val g = Array.fill(n)(0.0)
      var value = 0.0
      var i = 0
      while i < n - 1 do
        val residual = x(i + 1) - x(i) * x(i)
        value += 100.0 * residual * residual + (1.0 - x(i)) * (1.0 - x(i))
        g(i) += -400.0 * x(i) * residual - 2.0 * (1.0 - x(i))
        g(i + 1) += 200.0 * residual
        i += 1
      Right(ObjectiveEvaluation(value, DMat.fromArrayRowMajor(n, 1, g)))
    }.toOption.get

  private def quadratic(h: DMat, center: DMat): DifferentiableObjective =
    DifferentiableObjective(center.rows) { at =>
      val d = subtract(at, center)
      val g = multiply(h, d)
      Right(ObjectiveEvaluation(0.5 * dot(d, g), g))
    }.toOption.get

  private def gradient(objective: DifferentiableObjective, at: DMat): DMat =
    objective.evaluate(at).toOption.get.gradient
  private def column(values: Double*): DMat = DMat.fromArrayRowMajor(values.size, 1, values.toArray)
  private def scale(x: DMat, factor: Double): DMat = DMat.tabulate(x.rows, x.cols)((r, c) => factor * x(r, c))
  private def subtract(left: DMat, right: DMat): DMat =
    DMat.tabulate(left.rows, left.cols)((r, c) => left(r, c) - right(r, c))
  private def dot(left: DMat, right: DMat): Double = (0 until left.rows).map(r => left(r, 0) * right(r, 0)).sum
  private def normInf(value: DMat): Double = (0 until value.rows).map(r => Math.abs(value(r, 0))).max
  private def multiply(left: DMat, right: DMat): DMat =
    DMat.tabulate(left.rows, right.cols): (r, c) =>
      (0 until left.cols).map(k => left(r, k) * right(k, c)).sum

  test("L-BFGS reuses a curvature pair to solve a scalar quadratic in two updates"):
    val objective = DifferentiableObjective(1)(x =>
      Right(ObjectiveEvaluation(1.5 * x(0, 0) * x(0, 0), column(3 * x(0, 0))))
    ).toOption.get
    val solved = LBFGS
      .minimize(
        objective,
        column(1.0),
        LBFGSConfig(maxIterations = 2, tolerance = FirstOrderTolerance.from(1e-12, 0.0).toOption.get)
      )
      .toOption
      .get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.primal(0, 0), 0.0, 1e-12)

  test("a stationary initial point triggers exactly one cancellable progress event"):
    var calls = 0
    val objective = DifferentiableObjective(1)(x => Right(ObjectiveEvaluation(0.0, column(0.0)))).toOption.get
    val solved = LBFGS
      .minimize(
        objective,
        column(0.0),
        LBFGSConfig(control = SolverControl(progress = Some(_ => { calls += 1; false })))
      )
      .toOption
      .get
    assertEquals(calls, 1)
    assertEquals(solved.status, FirstOrderStoppingStatus.Cancelled)
