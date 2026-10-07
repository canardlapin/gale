package gale.optim

import gale.linalg.{DMat, DVec}
import gale.TestAccess

final class LevenbergMarquardtSuite extends munit.FunSuite:
  private def vector(x: Double*): DVec = DVec.fromSeq(x)
  private val tolerance = FirstOrderTolerance.from(1e-9, 0.0).toOption.get
  private val config = LevenbergMarquardtConfig(tolerance = tolerance)
  private def location = LeastSquaresObjective(1, 3)(
    x => Right(vector(x(0) - 1.0, x(0) - 2.0, x(0) - 6.0)),
    _ => Right(DMat.tabulate(3, 1)((_, _) => 1.0))
  ).toOption.get

  test("nonzero residual least squares returns the analytic mean and final raw stationarity"):
    val solved = LevenbergMarquardt.minimize(location, vector(0.0), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.parameters(0), 3.0, 1e-8)
    assertEqualsDouble(solved.objective, 7.0, 1e-12)
    assertEqualsDouble(
      solved.solution.certificate.primalResidual,
      Math.abs((0 until 3).map(solved.residuals(_)).sum),
      1e-14
    )
    assert(solved.solution.certificate.binds(solved.solution.primal, None))

  test("weighted rows recover the independently computed weighted mean"):
    val weighted = location.weighted(vector(1.0, 0.0, 4.0)).toOption.get
    val solved = LevenbergMarquardt.minimize(weighted, vector(0.0), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.parameters(0), 5.0, 1e-8)
    assertEqualsDouble(solved.objective, 10.0, 1e-12)
    assertEqualsDouble(solved.residuals(1), 0.0, 0.0)
    assert(location.weighted(vector(1.0, -1.0, 1.0)).isLeft)

  test("nonlinear exponential fitting recovers synthetic parameters"):
    val objective = LeastSquaresObjective(2, 12)(
      x => Right(DVec.tabulate(12)(i => x(0) * Math.exp(-x(1) * i / 3.0) - 2.0 * Math.exp(-0.4 * i / 3.0))),
      x =>
        Right(
          DMat.tabulate(12, 2)((r, c) =>
            if c == 0 then Math.exp(-x(1) * r / 3.0) else -x(0) * r / 3.0 * Math.exp(-x(1) * r / 3.0)
          )
        )
    ).toOption.get
    val solved = LevenbergMarquardt.minimize(objective, vector(1.0, 0.1), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.parameters(0), 2.0, 1e-7)
    assertEqualsDouble(solved.parameters(1), 0.4, 1e-7)

  test("rank deficient and underdetermined Jacobians are stabilized without a uniqueness claim"):
    val objective = LeastSquaresObjective(2, 1)(
      x => Right(vector(x(0) + x(1) - 3.0)),
      _ => Right(DMat.dense(1, 2, Seq(1.0, 1.0)))
    ).toOption.get
    val solved = LevenbergMarquardt.minimize(objective, vector(0.0, 0.0), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.parameters(0) + solved.parameters(1), 3.0, 1e-8)
    val zeroColumn = LeastSquaresObjective(2, 1)(
      x => Right(vector(x(0) - 2.0)),
      _ => Right(DMat.dense(1, 2, Seq(1.0, 0.0)))
    ).toOption.get
    val rank = LevenbergMarquardt.minimize(zeroColumn, vector(0.0, 7.0), config).toOption.get
    assertEquals(rank.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(rank.parameters(1), 7.0, 0.0)

  test("column scaling handles differently scaled parameters with raw gradient checks"):
    val objective = LeastSquaresObjective(2, 2)(
      x => Right(vector(1e-4 * x(0) - 1.0, 1e4 * x(1) - 1.0)),
      _ => Right(DMat.dense(2, 2, Seq(1e-4, 0.0, 0.0, 1e4)))
    ).toOption.get
    for scaling <- Vector(None, Some(vector(1e4, 1e-4))) do
      val solved =
        LevenbergMarquardt.minimize(objective, vector(0.0, 0.0), config.copy(parameterScale = scaling)).toOption.get
      assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
      assertEqualsDouble(solved.parameters(0), 1e4, 1e-5)
      assertEqualsDouble(solved.parameters(1), 1e-4, 1e-12)
      assert(solved.solution.certificate.primalResidual <= 1e-9)

  test("finite differences are explicit and every perturbation consumes budget"):
    val objective = LeastSquaresObjective.finiteDifferences(1, 1)(x => Right(vector(x(0) - 2.0))).toOption.get
    val initial = LevenbergMarquardt
      .minimize(objective, vector(0.0), config.copy(maxIterations = 0, control = SolverControl(maxEvaluations = 3)))
      .toOption
      .get
    assertEquals(initial.solution.evaluations.residuals, 3)
    assertEquals(initial.solution.evaluations.jacobians, 0)
    assertEquals(initial.status, FirstOrderStoppingStatus.IterationLimit)
    assert(
      LevenbergMarquardt
        .minimize(objective, vector(0.0), config.copy(control = SolverControl(maxEvaluations = 2)))
        .isLeft
    )
    val solved = LevenbergMarquardt.minimize(objective, vector(0.0), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(solved.parameters(0), 2.0, 1e-8)

  test("rejected trials are counted and final residuals describe the accepted point"):
    val objective = LeastSquaresObjective(1, 1)(
      x => Right(vector(x(0) * x(0) - 1.0)),
      x => Right(DMat.dense(1, 1, Seq(2.0 * x(0))))
    ).toOption.get
    val solved = LevenbergMarquardt.minimize(objective, vector(0.1), config).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.Converged)
    assert(solved.rejectedSteps > 0)
    assertEqualsDouble(solved.residuals(0), solved.parameters(0) * solved.parameters(0) - 1.0, 0.0)
    assert(solved.solution.evaluations.residuals > solved.solution.evaluations.jacobians)

  test("budgets and cancellation preserve the last completely evaluated accepted point"):
    val stopped = LevenbergMarquardt
      .minimize(location, vector(0.0), config.copy(control = SolverControl(maxEvaluations = 3)))
      .toOption
      .get
    assertEquals(stopped.status, FirstOrderStoppingStatus.EvaluationLimit)
    assertEqualsDouble(stopped.parameters(0), 0.0, 0.0)
    assertEqualsDouble(stopped.residuals(0), -1.0, 0.0)
    assertEquals(stopped.solution.evaluations.callbacks, 3)
    assertEquals(stopped.solution.certificate.iterations, 1)
    val cancelled = LevenbergMarquardt
      .minimize(location, vector(0.0), config.copy(control = SolverControl(progress = Some(_ => false))))
      .toOption
      .get
    assertEquals(cancelled.status, FirstOrderStoppingStatus.Cancelled)

  test("huge damping and an unrepresentable step do not certify convergence"):
    val solved = LevenbergMarquardt.minimize(location, vector(1e8), config.copy(initialDamping = 1e300)).toOption.get
    assertEquals(solved.status, FirstOrderStoppingStatus.NumericalStagnation)
    assert(solved.solution.certificate.primalResidual > 1.0)

  test("cancelling Jacobian products cannot manufacture stationarity"):
    for slopes <- Vector(Vector(1e16, 1.0, -1e16), Vector(1.0, -1e16, 1e16), Vector(-1e16, 1e16, 1.0)) do
      val objective = LeastSquaresObjective(1, 3)(
        x => Right(DVec.tabulate(3)(i => 1.0 + slopes(i) * x(0))),
        _ => Right(DMat.tabulate(3, 1)((r, _) => slopes(r)))
      ).toOption.get
      val solved = LevenbergMarquardt.minimize(objective, vector(0.0), config.copy(maxIterations = 0)).toOption.get
      assertEquals(solved.status, FirstOrderStoppingStatus.IterationLimit)
      assertEqualsDouble(solved.solution.certificate.primalResidual, 1.0, 0.0)

  test("invalid callbacks, dimensions and overflow fail through Either"):
    val bad = LeastSquaresObjective(1, 1)(
      _ => Right(vector(Double.MaxValue)),
      _ => Right(DMat.dense(1, 1, Seq(1.0)))
    ).toOption.get
    assert(LevenbergMarquardt.minimize(bad, vector(0.0)).isLeft)
    val shape = LeastSquaresObjective(1, 1)(_ => Right(vector(1.0, 2.0)), _ => Right(DMat.zeros(1, 1))).toOption.get
    assert(LevenbergMarquardt.minimize(shape, vector(0.0)).isLeft)
    val thrown = LeastSquaresObjective(1, 1)(
      _ => throw new IllegalStateException("fixture"),
      _ => Right(DMat.zeros(1, 1))
    ).toOption.get
    assert(
      LevenbergMarquardt.minimize(thrown, vector(0.0)).left.toOption.get.isInstanceOf[FirstOrderError.OracleFailure]
    )
    assert(LevenbergMarquardt.minimize(location, vector(0.0), config.copy(initialDamping = 0.0)).isLeft)

  test("large finite least-squares cost avoids overflowing the squared norm intermediate"):
    val objective = LeastSquaresObjective(1, 2)(
      _ => Right(vector(1e154, 1e154)),
      _ => Right(DMat.dense(2, 1, Seq(1e-154, 1e-154)))
    ).toOption.get
    val solved = LevenbergMarquardt.minimize(objective, vector(0.0), config.copy(maxIterations = 0)).toOption.get
    assert(Math.abs(solved.objective / 1e308 - 1.0) < 1e-15)
    assertEqualsDouble(solved.solution.certificate.primalResidual, 2.0, 1e-15)

  test("Jacobian snapshots preserve strided values and do not share callback storage"):
    val storage = DMat.dense(2, 2, Seq(1.0, 2.0, 3.0, 4.0))
    val objective = LeastSquaresObjective(2, 2)(_ => Right(vector(1.0, 2.0)), _ => Right(storage.t)).toOption.get
    val copied = LeastSquaresObjective
      .jacobian(objective, vector(0.0, 0.0), new OptimizationExecution(SolverControl()))
      .toOption
      .get
    TestAccess.dmatStorage(storage)(1) = 99.0
    assertEqualsDouble(copied(0, 0), 1.0, 0.0)
    assertEqualsDouble(copied(0, 1), 3.0, 0.0)
    assertEqualsDouble(copied(1, 0), 2.0, 0.0)
    assertEqualsDouble(copied(1, 1), 4.0, 0.0)
