package gale.optim

class BoxQuasiNewtonSuite extends munit.FunSuite:
  private def config(iterations: Int = 300, evaluations: Int = 1200, lines: Int = 30) =
    BoxQuasiNewtonConfig.from(iterations, evaluations, lines, 1e-9).toOption.get
  private def box(lower: Double, upper: Double, n: Int = 1) =
    BoxBounds.from(Vector.fill(n)(lower), Vector.fill(n)(upper)).toOption.get
  private def oracle(n: Int)(f: (Array[Double], Array[Double]) => Double) = new BoxDifferentiableObjective:
    val dimension = n
    def evaluate(x: Array[Double], g: Array[Double]): Either[BoxOptimizationError, Double] = Right(f(x, g))

  test("coupled quadratic reaches its interior solution with coherent value and gradient"):
    val objective = oracle(2): (x, g) =>
      val a = x(0) - 0.2
      val b = x(1) - 0.7
      g(0) = 4 * a + b
      g(1) = a + 2 * b
      2 * a * a + a * b + b * b
    val result = BoxQuasiNewton.minimize(objective, box(0, 1, 2), Vector(0.9, 0.1), config()).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.Stationary)
    val point = result.point.get
    assertEqualsDouble(point.coordinates(0), 0.2, 1e-8)
    assertEqualsDouble(point.coordinates(1), 0.7, 1e-8)
    assertEqualsDouble(point.value, 0.0, 1e-16)
    assert(point.projectedGradientNorm <= 1e-9)

  test("nonconvex negative-curvature start can reach a local minimum"):
    val objective = oracle(1): (x, g) =>
      g(0) = 4 * x(0) * (x(0) * x(0) - 1)
      math.pow(x(0) * x(0) - 1, 2)
    val result = BoxQuasiNewton.minimize(objective, box(-2, 2), Vector(0.2), config()).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.Stationary)
    assertEqualsDouble(result.point.get.coordinates.head, 1.0, 1e-8)

  test("coupled boundary optimum changes active set without leaving the box"):
    val objective = oracle(2): (x, g) =>
      assert(x.forall(v => v >= 0 && v <= 1))
      val a = x(0) + 0.2
      val b = x(1) - 0.4
      g(0) = 4 * a + b
      g(1) = a + 2 * b
      2 * a * a + a * b + b * b
    val result = BoxQuasiNewton.minimize(objective, box(0, 1, 2), Vector(0.9, 0.9), config()).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.Stationary)
    assertEqualsDouble(result.point.get.coordinates(0), 0.0, 1e-15)
    assertEqualsDouble(result.point.get.coordinates(1), 0.3, 1e-8)
    assert(result.point.get.gradient(0) > 0)
    assert(result.work.metricResets > 0)

  test("fixed coordinates are excluded from projected stationarity"):
    val objective = oracle(1): (x, g) =>
      g(0) = 2 * x(0)
      x(0) * x(0)
    val result = BoxQuasiNewton.minimize(objective, box(1, 1), Vector(1.0), config()).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.Stationary)
    assertEquals(result.work.evaluations, 1)

  test("evaluation and iteration limits retain exact counts and never claim convergence"):
    var calls = 0
    val objective = oracle(1): (x, g) =>
      calls += 1
      g(0) = 2 * x(0)
      x(0) * x(0)
    val limited = BoxQuasiNewton.minimize(objective, box(-2, 2), Vector(1.0), config(evaluations = 1)).toOption.get
    assertEquals(limited.status, BoxOptimizationStatus.EvaluationLimit)
    assertEquals(limited.work.evaluations, calls)
    assertEquals(calls, 1)
    val noSteps = BoxQuasiNewton.minimize(objective, box(-2, 2), Vector(1.0), config(iterations = 0)).toOption.get
    assertEquals(noSteps.status, BoxOptimizationStatus.IterationLimit)
    assertEquals(noSteps.work.evaluations, 1)

  test("rejected line search preserves the previous point and charges rejected calls"):
    val objective = oracle(1): (x, g) =>
      g(0) = 1.0 // deliberately inconsistent oracle forces rejection
      -x(0)
    val result = BoxQuasiNewton.minimize(objective, box(-2, 2), Vector(0.0), config(lines = 3)).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.LineSearchLimit)
    assertEquals(result.work.evaluations, 4)
    assertEquals(result.work.rejectedSteps, 3)
    assertEqualsDouble(result.point.get.coordinates.head, 0.0, 1e-15)

  test("refused and incomplete oracle calls are charged without stale gradients"):
    var calls = 0
    val objective = oracle(1): (x, g) =>
      calls += 1
      if calls == 1 then g(0) = 2 * x(0)
      x(0) * x(0)
    val result = BoxQuasiNewton.minimize(objective, box(-2, 2), Vector(1.0), config()).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.OracleRefused(BoxOptimizationError.NonFiniteOracle))
    assertEquals(result.work.evaluations, 2)
    assertEqualsDouble(result.point.get.coordinates.head, 1.0, 1e-15)
    val refused = new BoxDifferentiableObjective:
      val dimension = 1
      def evaluate(x: Array[Double], g: Array[Double]): Either[BoxOptimizationError, Double] =
        Left(BoxOptimizationError.OracleFailure("unavailable"))
    val empty = BoxQuasiNewton.minimize(refused, box(-2, 2), Vector(1.0), config()).toOption.get
    assertEquals(empty.point, None)
    assertEquals(empty.work.evaluations, 1)

  test("invalid input refuses before any oracle call"):
    val objective = oracle(1)((_, _) => fail("unexpected oracle call"))
    assert(BoxQuasiNewton.minimize(objective, box(0, 1), Vector(2.0), config()).isLeft)
    assert(BoxBounds.from(Vector(0.0), Vector(Double.PositiveInfinity)).isLeft)
    assert(BoxBounds.from(Vector(1.0), Vector(0.0)).isLeft)
    assert(BoxQuasiNewtonConfig.from(1, 0, 1, 1e-9).isLeft)
    assert(BoxQuasiNewtonConfig.from(1, 1, 1, Double.NaN).isLeft)

  test("Rosenbrock curved valley converges from an indefinite-Hessian region"):
    val objective = oracle(2): (x, g) =>
      val residual = x(1) - x(0) * x(0)
      g(0) = -400 * x(0) * residual + 2 * (x(0) - 1)
      g(1) = 200 * residual
      100 * residual * residual + math.pow(1 - x(0), 2)
    val result = BoxQuasiNewton.minimize(objective, box(-2, 2, 2), Vector(0.0, 1.0), config()).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.Stationary)
    result.point.get.coordinates.foreach(v => assertEqualsDouble(v, 1.0, 1e-7))

  test("an unrepresentable nonstationary step is not relabeled convergence"):
    val objective = oracle(1): (x, g) =>
      g(0) = 1e-300
      x(0) * 1e-300
    val settings = BoxQuasiNewtonConfig.from(10, 10, 3, 1e-310).toOption.get
    val result = BoxQuasiNewton.minimize(objective, box(0, 2), Vector(1.0), settings).toOption.get
    assertEquals(result.status, BoxOptimizationStatus.NoRepresentableStep)
    assertEquals(result.work.evaluations, 1)
