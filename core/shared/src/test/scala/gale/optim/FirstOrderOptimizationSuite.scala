package gale.optim

import gale.linalg.DMat
import gale.linalg.LinAlgError
import gale.linalg.LinearOperator

class FirstOrderOptimizationSuite extends munit.FunSuite:

  test("capability selection is deterministic and fails before numerical execution"):
    val capabilities = FirstOrderCapabilities
      .from(Set(FirstOrderMethod.ProjectedGradient, FirstOrderMethod.ProximalGradient))
      .toOption
      .get
    val selected = capabilities.select(
      Vector(FirstOrderMethod.ProjectedGradient, FirstOrderMethod.ProximalGradient),
      SolverMethodRequest.Automatic
    )
    val missing = capabilities.select(
      Vector(FirstOrderMethod.LinearCompositePrimalDual),
      SolverMethodRequest.Require(FirstOrderMethod.LinearCompositePrimalDual)
    )

    assertEquals(selected, Right(FirstOrderMethod.ProximalGradient))
    assert(missing.left.toOption.exists(_.isInstanceOf[FirstOrderError.MissingCapability]))

  test("proximal gradient matches the separable quadratic plus l1 oracle"):
    val center = matrix(Vector(Vector(1.0), Vector(-2.0)))
    val smooth = quadratic(center)
    val l1 = l1Term(2, 0.25)
    val solution = FirstOrderSolvers
      .proximalGradient(smooth, l1, DMat.zeros(2, 1))
      .toOption
      .get

    assertMatrixClose(solution.primal, matrix(Vector(Vector(0.75), Vector(-1.75))), 1e-7)
    assertEquals(solution.status, FirstOrderStoppingStatus.Converged)
    assertEquals(solution.certificate.settings.method, FirstOrderMethod.ProximalGradient)
    assert(solution.certificate.binds(solution.primal, solution.dual))
    assert(!solution.certificate.binds(solution.primal.updated(0, 0, 0.0), solution.dual))

  test("projected gradient matches a box-constrained quadratic oracle"):
    val center = matrix(Vector(Vector(2.0), Vector(-3.0)))
    val box = new ProjectionSet:
      val variableRows = 2
      def project(at: DMat): Either[FirstOrderError, DMat] =
        Right(map(at)(value => Math.max(-1.0, Math.min(1.0, value))))
    val solution = FirstOrderSolvers
      .projectedGradient(quadratic(center), box, DMat.zeros(2, 1))
      .toOption
      .get

    assertMatrixClose(solution.primal, matrix(Vector(Vector(1.0), Vector(-1.0))), 1e-8)
    assertEquals(solution.status, FirstOrderStoppingStatus.Converged)
    assert(solution.certificate.primalResidual <= 1e-8)

  test("linear-composite primal-dual matches an independent diagonal lasso oracle"):
    val center = matrix(Vector(Vector(1.0), Vector(-2.0)))
    val diagonal = matrix(Vector(Vector(2.0, 0.0), Vector(0.0, 1.0)))
    val operator = BoundedLinearOperator.from(diagonal, Math.sqrt(5.0)).toOption.get
    val solution = FirstOrderSolvers
      .linearCompositePrimalDual(
        proximalQuadratic(center),
        l1Functional(2, 0.25),
        operator,
        DMat.zeros(2, 1)
      )
      .toOption
      .get

    assertMatrixClose(solution.primal, matrix(Vector(Vector(0.5), Vector(-1.75))), 1e-6)
    assertEquals(solution.status, FirstOrderStoppingStatus.Converged)
    assert(solution.certificate.primalResidual <= 1e-7)
    assert(solution.certificate.dualResidual <= 1e-7)
    assertEquals(solution.certificate.settings.method, FirstOrderMethod.LinearCompositePrimalDual)

  test("smooth composite primal-dual matches an independent sparse smooth TV oracle"):
    val center = matrix(Vector(Vector(2.0), Vector(0.0)))
    val difference = matrix(Vector(Vector(1.0, -1.0)))
    val operator = BoundedLinearOperator.from(difference, Math.sqrt(2.0)).toOption.get
    val solution = FirstOrderSolvers
      .smoothCompositePrimalDual(
        coupledQuadratic(center, smoothnessWeight = 1.0),
        l1Term(2, weight = 0.25),
        l1Functional(1, weight = 0.25),
        operator,
        DMat.zeros(2, 1)
      )
      .toOption
      .get

    assertMatrixClose(solution.primal, matrix(Vector(Vector(1.0), Vector(0.5))), 2e-6)
    assertEquals(solution.status, FirstOrderStoppingStatus.Converged)
    assert(solution.certificate.primalResidual <= 2e-7)
    assert(solution.certificate.dualResidual <= 2e-7)
    assertEquals(solution.certificate.settings.method, FirstOrderMethod.SmoothCompositePrimalDual)

  test("exact linear reduction verifies the declared null-space image"):
    val basis = matrix(Vector(Vector(1.0), Vector(1.0)))
    val constraint = matrix(Vector(Vector(1.0, -1.0)))
    val badConstraint = matrix(Vector(Vector(1.0, 1.0)))
    val proof = ExactLinearReduction.verify(basis, constraint, FirstOrderTolerance.strict).toOption.get

    assertEquals(proof.freeRows, 1)
    assertEquals(proof.semanticRows, 2)
    assertEquals(proof.constraintRows, 1)
    assertEqualsDouble(proof.residual, 0.0, 0.0)
    assert(ExactLinearReduction.verify(basis, badConstraint, FirstOrderTolerance.strict).isLeft)

  test("oracle shape and non-finite failures are typed"):
    val center = matrix(Vector(Vector(1.0), Vector(-2.0)))
    val bad = new ProximalTerm:
      val variableRows = 2
      def value(at: DMat): Either[FirstOrderError, Double] = Right(Double.NaN)
      def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] = Right(at)

    FirstOrderSolvers.proximalGradient(quadratic(center), bad, DMat.zeros(2, 1)) match
      case Left(FirstOrderError.NonFiniteValue("proximal term", 0, value)) => assert(value.isNaN)
      case other => fail(s"expected typed non-finite oracle failure, got $other")

  test("operator failures retain the typed Gale cause"):
    val failing = LinearOperator.fromFunctions(1, 1)(
      (_, _) => throw LinAlgError.InvalidArgument("deliberate forward failure"),
      (_, _) => throw LinAlgError.InvalidArgument("deliberate adjoint failure")
    )

    ExactLinearReduction.verify(failing, DMat.eye(1), FirstOrderTolerance.strict) match
      case Left(FirstOrderError.OperatorFailure("null-space basis", cause)) =>
        assertEquals(cause, LinAlgError.InvalidArgument("deliberate forward failure"))
      case other =>
        fail(s"expected a typed operator failure, got $other")

  test("single-step residual describes the returned point"):
    val config = FirstOrderConfig.from(1, FirstOrderTolerance.strict).toOption.get
    val result =
      FirstOrderSolvers.proximalGradient(quadratic(DMat.zeros(1, 1)), l1Term(1, 0.0), DMat.eye(1), config).toOption.get
    assertEquals(result.status, FirstOrderStoppingStatus.IterationLimit)
    assertEqualsDouble(result.primal(0, 0), 0.01, 1e-14)
    assertEqualsDouble(result.certificate.primalResidual, 0.01, 1e-14)

  test("rounded and partially rounded steps never certify stationarity"):
    for columns <- Vector(1, 2) do
      val objective = new SmoothObjective:
        val variableRows = 1
        val lipschitz = 1e12
        def value(at: DMat): Either[FirstOrderError, Double] =
          Right((0 until columns).map(c => 0.5 * math.pow(at(0, c) - (if c == 0 then 1e8 else 0.0), 2)).sum)
        def gradient(at: DMat): Either[FirstOrderError, DMat] =
          Right(DMat.tabulate(1, columns)((_, c) => at(0, c) - (if c == 0 then 1e8 else 0.0)))
      val initial = DMat.tabulate(1, columns)((_, c) => if c == 0 then 1e8 + 10.0 else 1.0)
      val config = FirstOrderConfig.from(2, FirstOrderTolerance.strict).toOption.get
      val result = FirstOrderSolvers.proximalGradient(objective, l1Term(1, 0.0), initial, config).toOption.get
      assert(result.status != FirstOrderStoppingStatus.Converged)
      assert(
        result.certificate.primalResolution > result.certificate.settings.tolerance
          .threshold(result.certificate.settings.primalResidualScale)
      )
      assertEqualsDouble(objective.gradient(result.primal).toOption.get(0, 0), 10.0, 0.0)

  test("derived zero step is a typed numerical failure"):
    val objective = new SmoothObjective:
      val variableRows = 1
      val lipschitz = 1e308
      def value(at: DMat): Either[FirstOrderError, Double] = Right(0.5 * at(0, 0) * at(0, 0))
      def gradient(at: DMat): Either[FirstOrderError, DMat] = Right(at)
    val config = FirstOrderConfig.from(1, FirstOrderTolerance.strict, stepSafety = 1e-320).toOption.get
    assert(
      FirstOrderSolvers
        .proximalGradient(objective, l1Term(1, 0.0), DMat.eye(1), config)
        .left
        .toOption
        .exists(_.isInstanceOf[FirstOrderError.NumericalFailure])
    )

  test("certificate binding rejects a permutation preserving all summaries"):
    val solution = FirstOrderSolvers
      .proximalGradient(quadratic(matrix(Vector(Vector(1.0), Vector(2.0)))), l1Term(2, 0.0), DMat.zeros(2, 1))
      .toOption
      .get
    val permuted = DMat.tabulate(2, 1)((r, _) => solution.primal(1 - r, 0))
    assertEquals(ValueSummary.from(permuted), ValueSummary.from(solution.primal))
    assert(!solution.certificate.binds(permuted, None))
    assert(solution.certificate.binds(DMat.tabulate(2, 1)((r, c) => solution.primal(r, c)), None))
    intercept[IllegalArgumentException](solution.copy(objective = solution.objective + 1.0))

  test("primal-dual rejects unqualified extrapolation before oracle execution"):
    val config = FirstOrderConfig.from(10, FirstOrderTolerance.strict, extrapolation = 0.0).toOption.get
    val operator = BoundedLinearOperator.from(DMat.eye(1), 1.0).toOption.get
    val center = DMat.eye(1)
    assert(
      FirstOrderSolvers
        .linearCompositePrimalDual(proximalQuadratic(center), l1Functional(1, 10.0), operator, DMat.zeros(1, 1), config)
        .isLeft
    )
    assert(
      FirstOrderSolvers
        .smoothCompositePrimalDual(
          quadratic(center),
          l1Term(1, 0.0),
          l1Functional(1, 10.0),
          operator,
          DMat.zeros(1, 1),
          config
        )
        .isLeft
    )

  test("nonconvex projections retain a local fixed-point path"):
    val center = matrix(Vector(Vector(2.0), Vector(1.0)))
    val circle = new ProjectionSet:
      val variableRows = 2
      def project(at: DMat): Either[FirstOrderError, DMat] =
        val norm = math.hypot(at(0, 0), at(1, 0))
        Right(at * (1.0 / norm))
    val result = FirstOrderSolvers
      .projectedGradient(quadratic(center), circle, matrix(Vector(Vector(1.0), Vector(0.0))))
      .toOption
      .get
    assertEquals(result.status, FirstOrderStoppingStatus.Converged)
    assertEqualsDouble(result.primal(0, 0), 2.0 / math.sqrt(5.0), 1e-7)
    assertEqualsDouble(result.primal(1, 0), 1.0 / math.sqrt(5.0), 1e-7)

  test("null-space verification certifies containment without completeness"):
    val constraint = matrix(Vector(Vector(1.0, -1.0)))
    val zeroBasis = DMat.zeros(2, 1)
    val certificate = ExactLinearReduction.verify(zeroBasis, constraint, FirstOrderTolerance.strict).toOption.get
    assertEqualsDouble(certificate.residual, 0.0, 0.0)
    assertEqualsDouble(certificate.basisImage.squaredNorm, 0.0, 0.0)

  private def quadratic(center: DMat): SmoothObjective =
    new SmoothObjective:
      val variableRows: Int = center.rows
      val lipschitz: Double = 1.0
      def value(at: DMat): Either[FirstOrderError, Double] =
        Right(0.5 * squaredNorm(subtract(at, center)))
      def gradient(at: DMat): Either[FirstOrderError, DMat] =
        Right(subtract(at, center))

  private def proximalQuadratic(center: DMat): ProximalObjective =
    new ProximalObjective:
      val variableRows: Int = center.rows
      def value(at: DMat): Either[FirstOrderError, Double] =
        Right(0.5 * squaredNorm(subtract(at, center)))
      def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] =
        Right(scale(add(at, scale(center, step)), 1.0 / (1.0 + step)))

  private def coupledQuadratic(center: DMat, smoothnessWeight: Double): SmoothObjective =
    new SmoothObjective:
      val variableRows: Int = 2
      val lipschitz: Double = 1.0 + 2.0 * smoothnessWeight
      def value(at: DMat): Either[FirstOrderError, Double] =
        val difference = at(0, 0) - at(1, 0)
        Right(0.5 * squaredNorm(subtract(at, center)) + 0.5 * smoothnessWeight * difference * difference)
      def gradient(at: DMat): Either[FirstOrderError, DMat] =
        val difference = at(0, 0) - at(1, 0)
        Right(
          matrix(
            Vector(
              Vector(at(0, 0) - center(0, 0) + smoothnessWeight * difference),
              Vector(at(1, 0) - center(1, 0) - smoothnessWeight * difference)
            )
          )
        )

  private def l1Term(rows: Int, weight: Double): ProximalTerm =
    new ProximalTerm:
      val variableRows: Int = rows
      def value(at: DMat): Either[FirstOrderError, Double] = Right(weight * l1(at))
      def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] =
        Right(map(at): value =>
          Math.signum(value) * Math.max(0.0, Math.abs(value) - weight * step))

  private def l1Functional(rows: Int, weight: Double): LinearCompositeFunctional =
    new LinearCompositeFunctional:
      val targetRows: Int = rows
      def value(at: DMat): Either[FirstOrderError, Double] = Right(weight * l1(at))
      def proximalConjugate(at: DMat, step: Double): Either[FirstOrderError, DMat] =
        Right(map(at)(value => Math.max(-weight, Math.min(weight, value))))

  private def matrix(rows: Vector[Vector[Double]]): DMat =
    val columnCount = rows.headOption.fold(0)(_.size)
    require(rows.forall(_.size == columnCount), "matrix rows must have equal length")
    DMat.tabulate(rows.size, columnCount): (row, column) =>
      rows(row)(column)

  private def map(value: DMat)(function: Double => Double): DMat =
    DMat.tabulate(value.rows, value.cols): (row, column) =>
      function(value(row, column))

  private def add(left: DMat, right: DMat): DMat =
    DMat.tabulate(left.rows, left.cols): (row, column) =>
      left(row, column) + right(row, column)

  private def subtract(left: DMat, right: DMat): DMat =
    add(left, scale(right, -1.0))

  private def scale(value: DMat, factor: Double): DMat =
    value * factor

  private def squaredNorm(value: DMat): Double =
    var result = 0.0
    var row = 0
    while row < value.rows do
      var column = 0
      while column < value.cols do
        val current = value(row, column)
        result += current * current
        column += 1
      row += 1
    result

  private def l1(value: DMat): Double =
    var result = 0.0
    var row = 0
    while row < value.rows do
      var column = 0
      while column < value.cols do
        result += Math.abs(value(row, column))
        column += 1
      row += 1
    result

  private def assertMatrixClose(actual: DMat, expected: DMat, tolerance: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var row = 0
    while row < actual.rows do
      var column = 0
      while column < actual.cols do
        assertEqualsDouble(actual(row, column), expected(row, column), tolerance)
        column += 1
      row += 1

  test("a valid tiny Lipschitz bound for an affine objective cannot dilute a box residual"):
    val smooth = new SmoothObjective:
      val variableRows = 1
      val lipschitz = 1e-12
      def value(x: DMat) = Right(1e16 - 2.0 * x(0, 0))
      def gradient(x: DMat) = Right(DMat.tabulate(1, 1)((_, _) => -2.0))
    val box = ProjectionSets.box(1, -1.0, 1.0).toOption.get
    val result = FirstOrderSolvers.projectedGradient(smooth, box, DMat.zeros(1, 1)).toOption.get
    assertEqualsDouble(result.primal(0, 0), 1.0, 1e-8)
    assertEquals(result.status, FirstOrderStoppingStatus.Converged)
    assertEquals(result.evaluations.projections, 4)
    assertEquals(result.evaluations.gradients, 4)
    assertEqualsDouble(result.primalResidualStep, 1.0, 0.0)
