package gale.optim

import gale.linalg.{DMat, DVec, DoubleLinearOperator, LinAlgError, MutableDVec}

final class OptimizationContractsSuite extends munit.FunSuite:
  private val strict = FirstOrderTolerance.from(1e-8, 1e-8).toOption.get
  private def config(iterations: Int = 1000, theta: Double = 1.0) =
    FirstOrderConfig.from(iterations, strict, extrapolation = theta).toOption.get
  private def col(values: Double*) = DMat.fromArrayRowMajor(values.size, 1, values.toArray)

  test("both primal-dual APIs validate and honor warm starts"):
    val op = BoundedLinearOperator.from(zeroOperator, 0.0).toOption.get
    val badDual = Some(DMat.zeros(2, 1))
    assert(
      FirstOrderSolvers
        .linearCompositePrimalDual(proximalQuadratic(2.0), zeroFunctional, op, DMat.zeros(1, 1), initialDual = badDual)
        .isLeft
    )
    assert(
      FirstOrderSolvers
        .smoothCompositePrimalDual(
          quadratic(2.0),
          zeroTerm,
          zeroFunctional,
          op,
          DMat.zeros(1, 1),
          initialDual = badDual
        )
        .isLeft
    )
    val warm = Some(col(7.0))
    val linear = FirstOrderSolvers
      .linearCompositePrimalDual(proximalQuadratic(2.0), zeroFunctional, op, DMat.zeros(1, 1), config(), warm)
      .toOption
      .get
    val smooth = FirstOrderSolvers
      .smoothCompositePrimalDual(quadratic(2.0), zeroTerm, zeroFunctional, op, DMat.zeros(1, 1), config(), warm)
      .toOption
      .get
    assertEqualsDouble(linear.primal(0, 0), 2.0, 1e-6)
    assertEqualsDouble(smooth.primal(0, 0), 2.0, 1e-6)
    assertEquals(linear.dual.map(_(0, 0)), Some(0.0))
    assertEquals(smooth.dual.map(_(0, 0)), Some(0.0))

  test("primal-dual forward cache and exposed counts agree with instrumented operator calls"):
    val counted = new CountingIdentity
    val op = BoundedLinearOperator.from(counted, 1.0).toOption.get
    val solved = FirstOrderSolvers
      .linearCompositePrimalDual(proximalQuadratic(2.0), zeroFunctional, op, DMat.zeros(1, 1), config(4))
      .toOption
      .get
    val iterations = solved.certificate.iterations
    assertEquals(solved.evaluations.forwards, 1 + 2 * iterations)
    assertEquals(solved.evaluations.adjoints, 2 * iterations)
    assertEquals(counted.forwardCalls, solved.evaluations.forwards)
    assertEquals(counted.adjointCalls, solved.evaluations.adjoints)

  test("zero-norm operator is admitted while non-unit primal-dual extrapolation is rejected"):
    val op = BoundedLinearOperator.from(zeroOperator, 0.0).toOption.get
    assert(
      FirstOrderSolvers
        .smoothCompositePrimalDual(quadratic(1.0), zeroTerm, zeroFunctional, op, DMat.zeros(1, 1), config(10))
        .isRight
    )
    val noTheta = FirstOrderConfig.from(10, strict, extrapolation = 0.0).toOption.get
    assert(
      FirstOrderSolvers
        .smoothCompositePrimalDual(quadratic(1.0), zeroTerm, zeroFunctional, op, DMat.zeros(1, 1), noTheta)
        .isLeft
    )

  test("proximal gradient retains coupled matrix columns as one objective variable"):
    val coupled = new SmoothObjective:
      val variableRows = 1
      val lipschitz = 2.0
      def value(at: DMat) =
        val d = at(0, 0) - at(0, 1) - 2.0
        Right(0.5 * d * d)
      def gradient(at: DMat) =
        val d = at(0, 0) - at(0, 1) - 2.0
        Right(DMat.fromArrayRowMajor(1, 2, Array(d, -d)))
    val result = FirstOrderSolvers.proximalGradient(coupled, zeroTerm, DMat.zeros(1, 2), config()).toOption.get
    assertEqualsDouble(result.primal(0, 0) - result.primal(0, 1), 2.0, 1e-6)

  test("late infinite objective value becomes a typed error"):
    var calls = 0
    val lateInfinity = new SmoothObjective:
      val variableRows = 1
      val lipschitz = 1.0
      def value(at: DMat) =
        calls += 1
        Right(if calls == 1 then 0.5 * at(0, 0) * at(0, 0) else Double.PositiveInfinity)
      def gradient(at: DMat) = Right(at)
    val result = FirstOrderSolvers.proximalGradient(lateInfinity, zeroTerm, col(1.0), config(3))
    assert(result.left.toOption.exists(_.isInstanceOf[FirstOrderError.NonFiniteValue]))

  test("exact reduction rejects malformed and non-finite basis images"):
    val constraint = zeroOperator
    val malformed = new DoubleLinearOperator:
      val rows = 1; val cols = 1
      def applyTo(x: DVec, into: MutableDVec): Unit = into(0) = x(0)
      override def applyTo(input: DMat): Either[LinAlgError, DMat] = Right(DMat.zeros(1, input.cols + 1))
    val nonfinite = new DoubleLinearOperator:
      val rows = 1; val cols = 1
      def applyTo(x: DVec, into: MutableDVec): Unit = into(0) = Double.PositiveInfinity
    assert(ExactLinearReduction.verify(malformed, constraint, strict).isLeft)
    assert(
      ExactLinearReduction
        .verify(nonfinite, constraint, strict)
        .left
        .toOption
        .exists(_.isInstanceOf[FirstOrderError.NonFiniteValue])
    )

  test("DVec objective adapter gives the same scalar solution as a matrix objective"):
    val vector = DifferentiableObjective.vector(1)(x =>
      Right((0.5 * (x(0) - 3.0) * (x(0) - 3.0), DVec.tabulate(1)(i => x(i) - 3.0)))
    )
    val matrix = DifferentiableObjective(1)(x =>
      Right(ObjectiveEvaluation(0.5 * (x(0, 0) - 3.0) * (x(0, 0) - 3.0), col(x(0, 0) - 3.0)))
    )
    val a =
      LBFGS.minimize(vector.toOption.get, DVec.tabulate(1)(_ => 0.0), LBFGSConfig(tolerance = strict)).toOption.get
    val b = LBFGS.minimize(matrix.toOption.get, DMat.zeros(1, 1), LBFGSConfig(tolerance = strict)).toOption.get
    assertEqualsDouble(a.primal(0, 0), b.primal(0, 0), 1e-10)

  private val zeroOperator = new DoubleLinearOperator:
    val rows = 1; val cols = 1
    def applyTo(x: DVec, into: MutableDVec): Unit = into(0) = 0.0
    override def transposeApplyTo(x: DVec, into: MutableDVec): Unit = into(0) = 0.0

  private final class CountingIdentity extends DoubleLinearOperator:
    val rows = 1; val cols = 1
    var forwardCalls = 0; var adjointCalls = 0
    def applyTo(x: DVec, into: MutableDVec): Unit = { forwardCalls += 1; into(0) = x(0) }
    override def transposeApplyTo(x: DVec, into: MutableDVec): Unit = { adjointCalls += 1; into(0) = x(0) }

  private val zeroTerm = ProximalTerms.zero(1).toOption.get
  private val zeroFunctional = ProximalTerms.conjugate(zeroTerm)
  private def quadratic(center: Double) = new SmoothObjective:
    val variableRows = 1; val lipschitz = 1.0
    def value(at: DMat) = Right(0.5 * (at(0, 0) - center) * (at(0, 0) - center))
    def gradient(at: DMat) = Right(col(at(0, 0) - center))
  private def proximalQuadratic(center: Double) = new ProximalObjective:
    val variableRows = 1
    def value(at: DMat) = Right(0.5 * (at(0, 0) - center) * (at(0, 0) - center))
    def proximal(at: DMat, step: Double) = Right(col((at(0, 0) + step * center) / (1.0 + step)))
