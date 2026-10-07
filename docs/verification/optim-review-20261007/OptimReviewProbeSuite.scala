package gale.optim

import gale.linalg.DMat

/** Characterization probes for the planning review; not acceptance tests. */
class OptimReviewProbeSuite extends munit.FunSuite:
  private def smooth(center: Double = 0.0, bound: Double = 1.0): SmoothObjective =
    new SmoothObjective:
      val variableRows = 1
      val lipschitz = bound
      def value(at: DMat): Either[FirstOrderError, Double] =
        val delta = at(0, 0) - center
        Right(0.5 * delta * delta)
      def gradient(at: DMat): Either[FirstOrderError, DMat] =
        Right(DMat.tabulate(1, 1)((_, _) => at(0, 0) - center))

  private val zero = new ProximalTerm:
    val variableRows = 1
    def value(at: DMat): Either[FirstOrderError, Double] = Right(0.0)
    def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] = Right(at)

  test("summary binding accepts permuted distinct matrices"):
    val a = DMat.tabulate(2, 1)((r, _) => if r == 0 then 1.0 else 2.0)
    val b = DMat.tabulate(2, 1)((r, _) => if r == 0 then 2.0 else 1.0)
    assertEquals(ValueSummary.from(a), ValueSummary.from(b))
    val cert = FirstOrderCertificate(
      ValueSummary.from(a), None, 0.0, 0.0, 0.0, 0.0, 0,
      FirstOrderSettings(FirstOrderMethod.ProximalGradient, 1, FirstOrderTolerance.strict, 0.99, 1.0)
    )
    assert(cert.binds(b, None))

  test("single-step certificate residual is measured at the preceding iterate"):
    val config = FirstOrderConfig.from(1, FirstOrderTolerance.strict).toOption.get
    val result = FirstOrderSolvers.proximalGradient(smooth(), zero, DMat.eye(1), config).toOption.get
    assertEqualsDouble(result.primal(0, 0), 0.01, 1e-14)
    assertEqualsDouble(result.certificate.primalResidual, 1.0, 1e-14)
    assertEqualsDouble(smooth().gradient(result.primal).toOption.get(0, 0), 0.01, 1e-14)

  test("rounded-away step reports convergence despite a non-small gradient"):
    val objective = smooth(center = 1e8, bound = 1e12)
    val initial = DMat.tabulate(1, 1)((_, _) => 1e8 + 10.0)
    val result = FirstOrderSolvers.proximalGradient(objective, zero, initial).toOption.get
    assertEquals(result.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(result.certificate.primalResidual, 0.0, 0.0)
    val actualGradient = objective.gradient(result.primal).toOption.get(0, 0)
    assert(actualGradient > FirstOrderTolerance.strict.threshold(result.primal(0, 0)))

  test("valid finite configuration can underflow step and throw outside Either"):
    val config = FirstOrderConfig.from(1, FirstOrderTolerance.strict, stepSafety = 1e-320).toOption.get
    intercept[IllegalArgumentException]:
      FirstOrderSolvers.proximalGradient(smooth(bound = 1e308), zero, DMat.eye(1), config)

  test("null-space image check accepts an incomplete zero basis"):
    val basis = DMat.zeros(2, 1)
    val constraint = DMat.tabulate(1, 2)((_, c) => if c == 0 then 1.0 else -1.0)
    assert(ExactLinearReduction.verify(basis, constraint, FirstOrderTolerance.strict).isRight)

  test("Rayleigh accepts an asymmetric numerator and reports a nonstationary point as converged"):
    val a = DMat.tabulate(2, 2)((r, c) => if r == 0 && c == 0 then 1.0 else if r == 0 && c == 1 then 2.0 else 0.0)
    val result = ProjectedRayleigh.solve(
      a, DMat.eye(2), RayleighCone.NonnegativeOrthant,
      gale.linalg.Vec(1.0, 0.0)
    ).toOption.get
    assert(result.converged)
    assertEqualsDouble(result.certificate.stationarityResidual, 0.0, 0.0)
    // x(t) = (cos(t), sin(t)); q'(0) = a(0,1) + a(1,0) = 2 > 0.
    val angle = 1e-5
    val c = Math.cos(angle)
    val s = Math.sin(angle)
    val nearby = c * c + 2.0 * c * s
    assert(nearby > result.root)
