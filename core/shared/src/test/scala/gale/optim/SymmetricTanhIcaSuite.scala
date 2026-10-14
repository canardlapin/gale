package gale.optim

import gale.linalg.{DMat, Matrix}
import munit.FunSuite

class SymmetricTanhIcaSuite extends FunSuite:
  private val budget = BigInt(1024 * 1024)
  private val initial = Matrix.tabulate(2, 2)((i, j) => Vector(1.0, 0.1, 0.2, 1.3)(i * 2 + j))
  private val observations = Matrix.tabulate(2, 192): (i, t) =>
    val a = (t % 4 * 2 - 3) / math.sqrt(5.0)
    val b = ((t / 4) % 3 - 1) * math.sqrt(1.5)
    if i == 0 then math.cos(0.37) * a + math.sin(0.37) * b
    else -math.sin(0.37) * a + math.cos(0.37) * b
  private def ok[A](value: Either[?, A]): A = value.fold(e => fail(e.toString), identity)
  private def values(matrix: DMat): Vector[Double] =
    Vector.tabulate(matrix.rows * matrix.cols)(i => matrix(i / matrix.cols, i % matrix.cols))

  test("frozen NumPy 2.3.5 oracle, convergence and orthogonality on a Cartesian mixture"):
    val result = ok(SymmetricTanhIca.fitWhitened(observations, initial, budget))
    val oracle = Vector(0.9323273456060331, -0.36161543196496504, 0.3616154319649651, 0.9323273456060333)
    values(result.unmixing).zip(oracle).foreach((a, b) => assertEqualsDouble(a, b, 2e-10))
    assertEquals(result.iterations, 3)
    assert(result.change <= 1e-7)
    assert(result.orthogonalityError <= 1e-12)
    assertEquals(result.observationPasses, 4L)
    assertEquals(result.spectralBackend, "portable-pure")

  test("tile partition does not change arithmetic order or owned outputs"):
    val first = ok(SymmetricTanhIca.fitWhitened(observations, initial, budget, TanhIcaPolicy(tileSamples = 1)))
    val snapshot = values(first.unmixing)
    for tile <- Vector(7, 64, 192, 1024) do
      val next = ok(SymmetricTanhIca.fitWhitened(observations, initial, budget, TanhIcaPolicy(tileSamples = tile)))
      assertEquals(values(next.unmixing), snapshot)
      assertEquals(next.change, first.change)
      assertEquals(next.iterations, first.iterations)
    assertEquals(values(first.unmixing), snapshot)
    assertEquals(values(initial), Vector(1.0, 0.1, 0.2, 1.3))

  test("exhaustion returns nonconvergence rather than a usable rotation"):
    SymmetricTanhIca.fitWhitened(observations, initial, budget, TanhIcaPolicy(maxIterations = 1)) match
      case Left(TanhIcaFailure.NonConverged(iterations, change)) =>
        assertEquals(iterations, 1)
        assertEqualsDouble(change, 0.049883944726710316, 2e-12)
      case other => fail(other.toString)

  test("singular initialization is refused without a restart or rank reduction"):
    val zero = Matrix.zeros(2, 2)
    val duplicate = Matrix.tabulate(2, 2)((_, j) => j.toDouble + 1)
    for init <- Vector(zero, duplicate) do
      assert(SymmetricTanhIca.fitWhitened(observations, init, budget) match
        case Left(TanhIcaFailure.SingularRotation(0, _, _)) => true
        case _                                              => false)

  test("unwhitened, noncentered and nonfinite observations refuse"):
    val invalid = Vector(
      Matrix.tabulate(2, 192)((i, t) => observations(i, t) * 2),
      Matrix.tabulate(2, 192)((i, t) => observations(i, t) + 1),
      Matrix.tabulate(2, 192)((i, t) => if i == 0 && t == 3 then Double.NaN else observations(i, t)),
      Matrix.zeros(2, 192)
    )
    invalid.foreach(x => assert(SymmetricTanhIca.fitWhitened(x, initial, budget).isLeft))
    val nanInit = Matrix.tabulate(2, 2)((i, j) => if i == 0 then Double.NaN else initial(i, j))
    assert(SymmetricTanhIca.fitWhitened(observations, nanInit, budget).isLeft)

  test("budget is admitted before reading observations and covers the exact boundary"):
    val bytes = SymmetricTanhIca.numericalBytes(2, 192, 1024)
    val nan = Matrix.tabulate(2, 192)((_, _) => Double.NaN)
    assertEquals(SymmetricTanhIca.fitWhitened(nan, initial, bytes - 1), Left(TanhIcaFailure.Budget(bytes, bytes - 1)))
    assert(SymmetricTanhIca.fitWhitened(observations, initial, bytes).isRight)
    assert(SymmetricTanhIca.numericalBytes(Int.MaxValue, Int.MaxValue, Int.MaxValue) > BigInt(Long.MaxValue))

  test("cancellation before validation, during validation and during iteration publishes no model"):
    assertEquals(
      SymmetricTanhIca.fitWhitened(observations, initial, budget, cancelled = () => true),
      Left(TanhIcaFailure.Cancelled(0))
    )
    for limit <- Vector(3, 6, 8) do
      var calls = 0
      val result = SymmetricTanhIca.fitWhitened(
        observations,
        initial,
        budget,
        TanhIcaPolicy(tileSamples = 64),
        () => { calls += 1; calls >= limit }
      )
      assert(result match
        case Left(TanhIcaFailure.Cancelled(_)) => true
        case _                                 => false)

  test("shape and numerical policy failures precede allocations"):
    assert(SymmetricTanhIca.fitWhitened(Matrix.zeros(1, 4), Matrix.eye(1), budget).isLeft)
    assert(SymmetricTanhIca.fitWhitened(Matrix.eye(2), initial, budget).isLeft)
    for policy <- Vector(
        TanhIcaPolicy(maxIterations = 0),
        TanhIcaPolicy(tileSamples = 0),
        TanhIcaPolicy(tolerance = Double.NaN),
        TanhIcaPolicy(whiteningTolerance = -1),
        TanhIcaPolicy(orthogonalityTolerance = 1)
      )
    do assert(SymmetricTanhIca.fitWhitened(observations, initial, budget, policy).isLeft)
