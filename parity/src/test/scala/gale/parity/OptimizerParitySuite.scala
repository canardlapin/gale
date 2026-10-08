package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.eigSym
import breeze.optimize.DiffFunction
import breeze.optimize.FirstOrderMinimizer
import breeze.optimize.FirstOrderMinimizer.ConvergenceCheck
import breeze.optimize.FirstOrderMinimizer.ConvergenceReason
import breeze.optimize.FirstOrderMinimizer.State
import breeze.optimize.LBFGS as BreezeLBFGS
import breeze.optimize.LBFGSB as BreezeLBFGSB
import gale.linalg.DMat
import gale.optim.*

import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Paths
import scala.collection.mutable.ArrayBuffer

/** Solution-level parity of Gale's `LBFGS` / `LBFGSB` against Breeze 2.1
  * `breeze.optimize.LBFGS` / `LBFGSB`.
  *
  * Both libraries minimize the same objective callback, from the same start,
  * with the same history size and an equivalent stopping rule:
  *
  *   - Gale stops when `‖r(x)‖∞ ≤ abs + rel · max(1, ‖r(x₀)‖∞)`, where `r` is
  *     the gradient (L-BFGS) or the unit projected-gradient step
  *     `x − clamp(x − ∇f(x))` (L-BFGS-B).
  *   - Breeze's default checks are different (relative function improvement,
  *     an L2 gradient norm scaled by `|f|`, and for L-BFGS-B a hard-coded
  *     `‖r‖∞ < 1e-5`), so this suite gives Breeze L-BFGS a custom
  *     `ConvergenceCheck` computing Gale's threshold, and drives Breeze
  *     L-BFGS-B's iterator with the same predicate.
  *
  * Iterates and iteration counts are '''not''' compared. Assertions are on the
  * returned point and value, each bounded by what the shared stopping rule
  * implies:
  *
  *   - unconstrained, μ-strongly convex: `‖x − x*‖₂ ≤ √n τ / μ` per solver, so
  *     the two solvers agree within `2 √n τ / μ`;
  *   - box-constrained, μ-strongly convex, L-smooth: with `G` the unit
  *     projected-gradient step, the variational inequality gives
  *     `μ e² ≤ (1 + L) e ‖G‖ − ‖G‖²`, hence `e ≤ (1 + L) / μ · ‖G‖₂`. That
  *     bound is loose by about `κ`, so it also sets the active-set distance;
  *     once the known active set matches, the free coordinates are held to
  *     `(√|F| τ + L ‖d_A‖₂) / λmin(H_FF)` from the restricted Hessian;
  *   - values: `|f(a) − f(b)| ≤ ‖∇f(b)‖₂ d + L d² / 2` with `d = ‖a − b‖₂`,
  *     plus a rounding allowance;
  *   - Rosenbrock (nonconvex): the known optimum `1`, with the local strong
  *     convexity of the Hessian there and a factor 10 for nonquadratic terms.
  *
  * Wall time and callback counts are written to
  * `parity/target/optimizer-parity.md` for information only.
  */
class OptimizerParitySuite extends munit.FunSuite:

  // ---------------------------------------------------------------------------
  // Problems (library-agnostic, deterministic)
  // ---------------------------------------------------------------------------

  private type Oracle = Array[Double] => (Double, Array[Double])

  /** A minimization problem shared by both libraries.
    *
    * @param mu
    *   a lower bound on the strong-convexity constant, or for Rosenbrock the
    *   smallest Hessian eigenvalue at the optimum.
    * @param lipschitz
    *   an upper bound on the gradient Lipschitz constant (local for Rosenbrock).
    * @param bounds
    *   finite or infinite box bounds; `None` runs the unconstrained solvers.
    * @param known
    *   an analytic minimizer, when one exists.
    * @param hessian
    *   the constant Hessian of a quadratic, used for the free-coordinate bound.
    */
  private final case class Problem(
      name: String,
      n: Int,
      oracle: Oracle,
      start: Array[Double],
      mu: Double,
      lipschitz: Double,
      convex: Boolean = true,
      bounds: Option[(Array[Double], Array[Double])] = None,
      known: Option[Array[Double]] = None,
      knownActive: Option[Set[Int]] = None,
      hessian: Option[Array[Array[Double]]] = None,
      tolerance: FirstOrderTolerance = sharedTolerance
  )

  /** Product of three Householder reflections: an exactly specified orthogonal matrix. */
  private def orthogonal(n: Int, seed: Long): Array[Array[Double]] =
    val q = Array.tabulate(n, n)((i, j) => if i == j then 1.0 else 0.0)
    val rng = new scala.util.Random(seed)
    for _ <- 0 until 3 do
      val v = Array.fill(n)(rng.nextGaussian())
      val vv = v.map(x => x * x).sum
      // q ← q (I − 2 v vᵀ / vᵀv)
      for i <- 0 until n do
        var qv = 0.0
        var k = 0
        while k < n do
          qv += q(i)(k) * v(k)
          k += 1
        val scale = 2.0 * qv / vv
        k = 0
        while k < n do
          q(i)(k) -= scale * v(k)
          k += 1
    q

  /** `Q diag(λ) Qᵀ` with eigenvalues log-spaced on `[1, κ]`, so `μ = 1` and `L = κ`. */
  private def spdWithCondition(n: Int, kappa: Double, seed: Long): Array[Array[Double]] =
    val q = orthogonal(n, seed)
    val lambda = Array.tabulate(n)(i => if n == 1 then 1.0 else math.pow(kappa, i.toDouble / (n - 1)))
    val a = Array.ofDim[Double](n, n)
    for i <- 0 until n; j <- 0 to i do
      var s = 0.0
      var k = 0
      while k < n do
        s += q(i)(k) * lambda(k) * q(j)(k)
        k += 1
      a(i)(j) = s
      a(j)(i) = s
    a

  private def matVec(a: Array[Array[Double]], x: Array[Double]): Array[Double] =
    Array.tabulate(a.length): i =>
      val row = a(i)
      var s = 0.0
      var k = 0
      while k < x.length do
        s += row(k) * x(k)
        k += 1
      s

  /** `f(x) = ½ dᵀAd + g*ᵀd` with `d = x − x*`: `½ xᵀAx − bᵀx` for `b = Ax* − g*`, shifted so that `f(x*) = 0`.
    *
    * The shift matters. Strict sufficient decrease must resolve `f(x) − f(x + αp) ≈ ‖g‖² / L`, which falls below
    * `eps · |f|` before `‖g‖∞` reaches 1e-8 when `|f*|` is O(100); see the roundoff-floor test.
    */
  private def quadraticOracle(a: Array[Array[Double]], xStar: Array[Double], gStar: Array[Double]): Oracle = x =>
    val d = Array.tabulate(x.length)(i => x(i) - xStar(i))
    val ad = matVec(a, d)
    var value = 0.0
    var i = 0
    while i < x.length do
      value += d(i) * (0.5 * ad(i) + gStar(i))
      i += 1
    (value, Array.tabulate(x.length)(i => ad(i) + gStar(i)))

  /** The same quadratic in unshifted form `½ xᵀAx − bᵀx`. */
  private def unshiftedQuadraticOracle(a: Array[Array[Double]], b: Array[Double]): Oracle = x =>
    val ax = matVec(a, x)
    var value = 0.0
    var i = 0
    while i < x.length do
      value += x(i) * (0.5 * ax(i) - b(i))
      i += 1
    (value, Array.tabulate(x.length)(i => ax(i) - b(i)))

  private def unconstrainedQuadratic(n: Int, kappa: Double, seed: Long): Problem =
    val a = spdWithCondition(n, kappa, seed)
    val xStar = ParitySupport.vectorData(n, seed + 1).map(_ * 3.0)
    Problem(
      s"quadratic n=$n κ=${fmtKappa(kappa)}",
      n,
      quadraticOracle(a, xStar, Array.fill(n)(0.0)),
      Array.fill(n)(0.0),
      mu = 1.0,
      lipschitz = kappa,
      known = Some(xStar),
      hessian = Some(a)
    )

  /** Bound-constrained quadratic whose minimizer is fixed by construction through KKT.
    *
    * Coordinates are assigned to lower / upper / free; `x*` is placed on the bound or strictly inside, and the
    * gradient at `x*` is chosen strictly positive / negative / zero accordingly. Then `b = A x* − g*`, so `x*`
    * satisfies strict complementarity and is the unique minimizer of a strongly convex problem.
    */
  private def boxQuadratic(
      n: Int,
      kappa: Double,
      seed: Long,
      role: Int => Int, // 0 free, -1 lower active, +1 upper active
      start: Array[Double] => Array[Double] = _ => Array.empty
  ): Problem =
    val a = spdWithCondition(n, kappa, seed)
    val lower = Array.fill(n)(-1.0)
    val upper = Array.fill(n)(1.0)
    val rng = new scala.util.Random(seed + 7)
    val xStar = Array.tabulate(n): i =>
      role(i) match
        case -1 => lower(i)
        case 1  => upper(i)
        case _  => rng.nextDouble() * 1.6 - 0.8
    val gStar = Array.tabulate(n): i =>
      role(i) match
        case -1 => 0.5 + rng.nextDouble()
        case 1  => -(0.5 + rng.nextDouble())
        case _  => 0.0
    val chosen = start(lower)
    Problem(
      s"box quadratic n=$n κ=${fmtKappa(kappa)}",
      n,
      quadraticOracle(a, xStar, gStar),
      if chosen.isEmpty then Array.fill(n)(0.0) else chosen,
      mu = 1.0,
      lipschitz = kappa,
      bounds = Some((lower, upper)),
      known = Some(xStar),
      knownActive = Some((0 until n).filter(role(_) != 0).toSet),
      hessian = Some(a)
    )

  /** Chained Rosenbrock `Σ 100 (x_{i+1} − x_i²)² + (1 − x_i)²`. */
  private val rosenbrock: Oracle = x =>
    val n = x.length
    val g = new Array[Double](n)
    var value = 0.0
    var i = 0
    while i < n - 1 do
      val t = x(i + 1) - x(i) * x(i)
      val s = 1.0 - x(i)
      value += 100.0 * t * t + s * s
      g(i) += -400.0 * x(i) * t - 2.0 * s
      g(i + 1) += 200.0 * t
      i += 1
    (value, g)

  /** Smallest eigenvalue and spectral norm of the chained-Rosenbrock Hessian at `1`. */
  private def rosenbrockCurvatureAtOne(n: Int): (Double, Double) =
    val h = BDM.zeros[Double](n, n)
    for i <- 0 until n - 1 do
      // at x = 1: d²/dx_i² = 1200 x_i² − 400 x_{i+1} + 2 = 802, d²/dx_{i+1}² += 200, cross = −400 x_i = −400
      h(i, i) += 802.0
      h(i + 1, i + 1) += 200.0
      h(i, i + 1) += -400.0
      h(i + 1, i) += -400.0
    val values = eigSym.justEigenvalues(h)
    (values.toArray.min, values.toArray.max)

  private def rosenbrockProblem(n: Int): Problem =
    val (mu, l) = rosenbrockCurvatureAtOne(n)
    Problem(
      s"rosenbrock n=$n",
      n,
      rosenbrock,
      Array.tabulate(n)(i => if i % 2 == 0 then -1.2 else 1.0),
      mu = mu,
      lipschitz = l,
      convex = false,
      known = Some(Array.fill(n)(1.0))
    )

  /** Mean logistic loss plus `λ/2 ‖w‖²` on fixed synthetic data with labels `±1`. */
  private def logisticProblem(
      p: Int,
      m: Int,
      lambda: Double,
      seed: Long,
      box: Option[Double] = None
  ): Problem =
    val rng = new scala.util.Random(seed)
    val features = Array.fill(m, p)(rng.nextGaussian())
    val truth = Array.fill(p)(rng.nextGaussian())
    val labels = Array.tabulate(m): i =>
      val margin = features(i).indices.map(j => features(i)(j) * truth(j)).sum
      val probability = 1.0 / (1.0 + math.exp(-margin))
      if rng.nextDouble() < probability then 1.0 else -1.0
    val oracle: Oracle = w =>
      val g = new Array[Double](p)
      var loss = 0.0
      var i = 0
      while i < m do
        val row = features(i)
        var z = 0.0
        var j = 0
        while j < p do
          z += row(j) * w(j)
          j += 1
        z *= labels(i)
        loss += (if z > 0.0 then math.log1p(math.exp(-z)) else -z + math.log1p(math.exp(z)))
        // d/dz log(1 + e^{−z}) = −σ(−z)
        val weight = -labels(i) * (if z > 0.0 then math.exp(-z) / (1.0 + math.exp(-z)) else 1.0 / (1.0 + math.exp(z)))
        j = 0
        while j < p do
          g(j) += weight * row(j)
          j += 1
        i += 1
      var penalty = 0.0
      var j = 0
      while j < p do
        penalty += w(j) * w(j)
        g(j) = g(j) / m + lambda * w(j)
        j += 1
      (loss / m + 0.5 * lambda * penalty, g)
    val frobenius = features.map(_.map(x => x * x).sum).sum
    Problem(
      s"${if box.isDefined then "box " else ""}logistic p=$p m=$m",
      p,
      oracle,
      Array.fill(p)(0.0),
      mu = lambda,
      lipschitz = frobenius / (4.0 * m) + lambda,
      bounds = box.map(r => (Array.fill(p)(-r), Array.fill(p)(r)))
    )

  private def fmtKappa(kappa: Double): String =
    if kappa < 10 then f"$kappa%.0f" else s"1e${math.round(math.log10(kappa))}"

  // ---------------------------------------------------------------------------
  // Stopping rule shared by both libraries
  // ---------------------------------------------------------------------------

  private val absoluteTolerance = 1e-8
  private val relativeTolerance = 0.0
  private val history = 10
  private val iterationLimit = 20000
  private val sharedTolerance = FirstOrderTolerance.from(absoluteTolerance, relativeTolerance).toOption.get

  /** `‖∇f‖∞` unconstrained, `‖x − clamp(x − ∇f)‖∞` boxed. */
  private def residual(p: Problem, x: Array[Double], g: Array[Double]): Double =
    p.bounds match
      case None                 => g.map(math.abs).max
      case Some((lower, upper)) =>
        x.indices.map(i => math.abs(x(i) - math.max(lower(i), math.min(upper(i), x(i) - g(i))))).max

  /** Gale's threshold, computed from the start exactly as `LBFGS` / `LBFGSB` do. */
  private def threshold(p: Problem): Double =
    val (_, g0) = p.oracle(p.start)
    p.tolerance.threshold(math.max(1.0, residual(p, p.start, g0)))

  // ---------------------------------------------------------------------------
  // Drivers
  // ---------------------------------------------------------------------------

  private final case class Run(
      x: Array[Double],
      value: Double,
      converged: Boolean,
      status: String,
      callbacks: Int,
      millis: Double
  )

  private final class Counted(oracle: Oracle):
    var calls = 0
    def apply(x: Array[Double]): (Double, Array[Double]) =
      calls += 1
      oracle(x)

  private def column(x: Array[Double]): DMat = DMat.fromArrayRowMajor(x.length, 1, x.clone())

  private def runGale(p: Problem): Run =
    val counted = Counted(p.oracle)
    val objective = DifferentiableObjective(p.n): x =>
      val (value, g) = counted(Array.tabulate(p.n)(x(_, 0)))
      Right(ObjectiveEvaluation(value, column(g)))
    .toOption.get
    val started = System.nanoTime()
    val solved = p.bounds match
      case None =>
        LBFGS.minimize(
          objective,
          column(p.start),
          LBFGSConfig(maxIterations = iterationLimit, tolerance = p.tolerance, historySize = history)
        )
      case Some((lower, upper)) =>
        val box = MatrixBoxBounds.from(column(lower), column(upper)).toOption.get
        LBFGSB.minimize(
          objective,
          box,
          column(p.start),
          LBFGSBConfig(maxIterations = iterationLimit, tolerance = p.tolerance, historySize = history)
        )
    val millis = (System.nanoTime() - started) / 1e6
    solved match
      case Left(error)     => fail(s"${p.name}: gale returned $error")
      case Right(solution) =>
        val x = Array.tabulate(p.n)(solution.primal(_, 0))
        val certified = solution.certificate.primalResidual <=
          solution.certificate.settings.tolerance.threshold(solution.certificate.settings.primalResidualScale)
        Run(
          x,
          solution.objective,
          solution.status == FirstOrderStoppingStatus.Converged && certified,
          solution.status.toString,
          counted.calls,
          millis
        )

  private case object GaleResidualConverged extends ConvergenceReason:
    def reason: String = "gale-equivalent residual converged"

  private def diffFunction(counted: Counted): DiffFunction[BDV[Double]] = new DiffFunction[BDV[Double]]:
    def calculate(x: BDV[Double]): (Double, BDV[Double]) =
      val (value, g) = counted(x.toArray)
      (value, BDV(g))

  /** Gale's rule as a Breeze `ConvergenceCheck`; residual first so a converged last iterate is reported as such. */
  private def galeEquivalentCheck(p: Problem, tau: Double): ConvergenceCheck[BDV[Double]] =
    ConvergenceCheck.fromPartialFunction[BDV[Double]] {
      case s if residual(p, s.x.toArray, s.grad.toArray) <= tau =>
        GaleResidualConverged
    } || FirstOrderMinimizer.maxIterationsReached[BDV[Double]](iterationLimit) ||
      FirstOrderMinimizer.searchFailed[BDV[Double]]

  /** Breeze L-BFGS-B with its built-in `‖r‖∞ < 1e-5` check bypassed: the public iterator is stopped by our rule. */
  private final class GaleCheckedLBFGSB(lower: BDV[Double], upper: BDV[Double])
      extends BreezeLBFGSB(lower, upper, maxIter = iterationLimit, m = history):
    def solveWith(f: DiffFunction[BDV[Double]], init: BDV[Double], stop: State => Option[String]): (State, String) =
      val states = infiniteIterations(f, initialState(adjustFunction(f), init))
      var state = states.next()
      var reason = stop(state)
      while reason.isEmpty do
        state = states.next()
        reason = stop(state)
      (state, reason.get)

  private def runBreeze(p: Problem): Run =
    val tau = threshold(p)
    val counted = Counted(p.oracle)
    val f = diffFunction(counted)
    val init = BDV(p.start.clone())
    val started = System.nanoTime()
    val (x, value, reason) = p.bounds match
      case None =>
        val state = new BreezeLBFGS[BDV[Double]](galeEquivalentCheck(p, tau), history).minimizeAndReturnState(f, init)
        (state.x, state.value, state.convergenceReason.map(_.reason).getOrElse("none"))
      case Some((lower, upper)) =>
        val solver = GaleCheckedLBFGSB(BDV(lower.clone()), BDV(upper.clone()))
        val (state, why) = solver.solveWith(
          f,
          init,
          s =>
            if residual(p, s.x.toArray, s.grad.toArray) <= tau then Some(GaleResidualConverged.reason)
            else if s.searchFailed then Some(FirstOrderMinimizer.SearchFailed.reason)
            else if s.iter >= iterationLimit then Some("max iterations reached")
            else None
        )
        (state.x, state.value, why)
    val millis = (System.nanoTime() - started) / 1e6
    Run(x.toArray, value, reason == GaleResidualConverged.reason, reason, counted.calls, millis)

  // ---------------------------------------------------------------------------
  // Assertions and the informational table
  // ---------------------------------------------------------------------------

  private val report = ArrayBuffer.empty[String]
  private val floorReport = ArrayBuffer.empty[String]

  private def norm2(v: Array[Double]): Double = math.sqrt(v.map(x => x * x).sum)
  private def distance(a: Array[Double], b: Array[Double]): Double = norm2(Array.tabulate(a.length)(i => a(i) - b(i)))

  /** Error bound implied by a residual of at most `tau` per solver: `‖r‖₂ / μ` unconstrained, `(1 + L) / μ · ‖G‖₂`
    * boxed. The box bound is loose by about `κ`; [[assertFreeCoordinates]] tightens it once the active set is known.
    */
  private def singleSolverBound(p: Problem, tau: Double): Double =
    val gNorm = math.sqrt(p.n.toDouble) * tau
    val base = if p.bounds.isEmpty then gNorm / p.mu else (1.0 + p.lipschitz) / p.mu * gNorm
    if p.convex then base else 10.0 * base

  private def valueBound(p: Problem, at: Array[Double], other: Array[Double], value: Double): Double =
    val d = distance(at, other)
    val (_, g) = p.oracle(at)
    norm2(g) * d + 0.5 * p.lipschitz * d * d + 1e-13 * math.max(1.0, math.abs(value))

  /** Coordinates within `delta` of a bound. With `delta` the per-solver error bound, every coordinate active at `x*`
    * qualifies, and a coordinate free at `x*` qualifies only if `x*` lies within `2 delta` of a bound.
    */
  private def active(p: Problem, x: Array[Double], delta: Double): Set[Int] =
    p.bounds match
      case None                 => Set.empty
      case Some((lower, upper)) =>
        x.indices.filter(i => x(i) - lower(i) <= delta || upper(i) - x(i) <= delta).toSet

  /** Tight check on the free coordinates of a box quadratic with a known active set `A` and free set `F`.
    *
    * Free coordinates sit more than `τ` inside the box, so a unit projected-gradient residual `≤ τ` gives
    * `‖g_F‖∞ ≤ τ`. With `g*_F = 0`, `g_F = H_FF d_F + H_FA d_A`, hence `‖d_F‖₂ ≤ (√|F| τ + L ‖d_A‖₂) / λmin(H_FF)`.
    */
  private def assertFreeCoordinates(p: Problem, tau: Double, label: String, x: Array[Double]): Unit =
    for
      h <- p.hessian
      activeSet <- p.knownActive
      xStar <- p.known
      if p.bounds.isDefined
    do
      val free = (0 until p.n).filterNot(activeSet).toArray
      if free.nonEmpty then
        val restricted = BDM.tabulate(free.length, free.length)((i, j) => h(free(i))(free(j)))
        val muFree = eigSym.justEigenvalues(restricted).toArray.min
        val g = p.oracle(x)._2
        val freeGradient = free.map(i => math.abs(g(i))).max
        assert(freeGradient <= tau, s"${p.name}: $label free gradient $freeGradient > $tau")
        val activeOffset = math.sqrt(activeSet.toSeq.map(i => (x(i) - xStar(i)) * (x(i) - xStar(i))).sum)
        val freeError = math.sqrt(free.map(i => (x(i) - xStar(i)) * (x(i) - xStar(i))).sum)
        val bound = (math.sqrt(free.length.toDouble) * tau + p.lipschitz * activeOffset) / muFree
        assert(freeError <= bound, s"${p.name}: $label ‖x_F − x*_F‖ = $freeError > $bound")

  private def compare(p: Problem): Unit =
    val tau = threshold(p)
    val gale = runGale(p)
    val breeze = runBreeze(p)
    // Independent recomputation of the stopping rule at each returned point.
    val galeResidual = residual(p, gale.x, p.oracle(gale.x)._2)
    val breezeResidual = residual(p, breeze.x, p.oracle(breeze.x)._2)
    assert(gale.converged, s"${p.name}: gale did not converge (${gale.status}, residual $galeResidual, f ${gale.value})")
    assert(
      breeze.converged,
      s"${p.name}: breeze did not reach the shared rule (${breeze.status}, residual $breezeResidual, f ${breeze.value})"
    )
    assert(galeResidual <= tau, s"${p.name}: gale exit residual $galeResidual > $tau")
    assert(breezeResidual <= tau, s"${p.name}: breeze exit residual $breezeResidual > $tau")
    p.bounds.foreach: (lower, upper) =>
      for x <- Seq(gale.x, breeze.x); i <- x.indices do
        assert(x(i) >= lower(i) && x(i) <= upper(i), s"${p.name}: infeasible coordinate $i")

    val single = singleSolverBound(p, tau)
    val dx = distance(gale.x, breeze.x)
    assert(dx <= 2.0 * single, s"${p.name}: ‖x_gale − x_breeze‖ = $dx > ${2.0 * single}")
    p.known.foreach: xStar =>
      assert(distance(gale.x, xStar) <= single, s"${p.name}: gale ‖x − x*‖ = ${distance(gale.x, xStar)} > $single")
      assert(
        distance(breeze.x, xStar) <= single,
        s"${p.name}: breeze ‖x − x*‖ = ${distance(breeze.x, xStar)} > $single"
      )
    val df = math.abs(gale.value - breeze.value)
    val fBound = valueBound(p, breeze.x, gale.x, breeze.value)
    assert(df <= fBound, s"${p.name}: |f_gale − f_breeze| = $df > $fBound")

    val galeActive = active(p, gale.x, single)
    if p.bounds.isDefined then
      assertEquals(galeActive, active(p, breeze.x, single), s"${p.name}: active sets differ")
      p.knownActive.foreach(expected => assertEquals(galeActive, expected, s"${p.name}: wrong active set"))
      assertFreeCoordinates(p, tau, "gale", gale.x)
      assertFreeCoordinates(p, tau, "breeze", breeze.x)

    report += f"| ${p.name} | ${if p.bounds.isDefined then s"B (${galeActive.size} active)" else "U"} | $tau%.1e | $dx%.1e | ${2.0 * single}%.1e " +
      f"| $df%.1e | ${gale.callbacks} | ${breeze.callbacks} | ${gale.millis}%.1f | ${breeze.millis}%.1f |"

  // ---------------------------------------------------------------------------
  // Tests
  // ---------------------------------------------------------------------------

  test("unconstrained quadratics, well- to ill-conditioned (κ ≤ 1e6), agree with each other and with x*") {
    Seq(
      unconstrainedQuadratic(10, 10.0, 11L),
      unconstrainedQuadratic(50, 1e3, 12L),
      unconstrainedQuadratic(100, 1e4, 13L),
      unconstrainedQuadratic(200, 1e6, 14L)
    ).foreach(compare)
  }

  test("chained Rosenbrock n = 2, 10, 100 reaches the known optimum 1 in both libraries") {
    Seq(2, 10, 100).map(rosenbrockProblem).foreach(compare)
  }

  test("L2-regularized logistic regression agrees on weights and loss") {
    Seq(logisticProblem(20, 500, 1e-2, 21L), logisticProblem(100, 2000, 1e-2, 22L)).foreach(compare)
  }

  test("box quadratics with known solution: mixed, all-inactive, and all-active bounds") {
    Seq(
      boxQuadratic(30, 1e2, 31L, i => Seq(-1, 0, 1)(i % 3)),
      boxQuadratic(100, 1e4, 32L, i => if i % 4 == 0 then -1 else if i % 4 == 1 then 1 else 0),
      boxQuadratic(40, 1e3, 33L, _ => 0),
      boxQuadratic(20, 1e2, 34L, i => if i % 2 == 0 then -1 else 1)
    ).foreach(compare)
  }

  test("bounded logistic regression agrees on weights, loss, and active set") {
    compare(logisticProblem(50, 1000, 1e-2, 41L, box = Some(0.25)))
  }

  test("edge cases: start at the optimum, start on a bound, and n = 1") {
    val atOptimum = unconstrainedQuadratic(10, 1e2, 51L)
    compare(atOptimum.copy(name = "quadratic, start = x*", start = atOptimum.known.get.clone()))
    val boxed = boxQuadratic(12, 1e2, 52L, i => Seq(-1, 0, 1)(i % 3))
    compare(boxed.copy(name = "box quadratic, start = x*", start = boxed.known.get.clone()))
    compare(boxQuadratic(12, 1e2, 53L, i => Seq(-1, 0, 1)(i % 3), lower => lower.clone()).copy(name = "box, start = l"))
    compare(unconstrainedQuadratic(1, 1.0, 54L).copy(name = "quadratic n=1"))
    compare(boxQuadratic(1, 1.0, 55L, _ => 1).copy(name = "box n=1, upper active"))
  }

  test("relative tolerance: the threshold scales with max(1, ‖r(x₀)‖∞) in both libraries") {
    val relative = FirstOrderTolerance.from(1e-10, 1e-11).toOption.get
    val cases = Seq(
      unconstrainedQuadratic(50, 1e3, 12L).copy(name = "quadratic n=50 κ=1e3, rel 1e-11", tolerance = relative),
      // The unit projected step is capped by the distance to the far bound, so start at a corner where it can exceed 1.
      boxQuadratic(30, 1e2, 31L, i => Seq(-1, 0, 1)(i % 3), lower => lower.clone())
        .copy(name = "box n=30 κ=1e2, start = l, rel 1e-11", tolerance = relative)
    )
    for p <- cases do
      assert(threshold(p) > relative.threshold(1.0), s"${p.name}: ‖r(x₀)‖∞ ≤ 1, scaled branch not exercised")
      compare(p)
  }

  /** The same problem evaluated as `½ xᵀAx − bᵀx`, so `|f*|` is O(100) instead of 0. */
  private def unshifted(p: Problem, kappa: Double, seed: Long): Problem =
    val a = spdWithCondition(p.n, kappa, seed)
    val xStar = p.known.get
    val gStar = p.oracle(xStar)._2
    val ax = matVec(a, xStar)
    p.copy(
      name = s"${p.name}, unshifted",
      oracle = unshiftedQuadraticOracle(a, Array.tabulate(p.n)(i => ax(i) - gStar(i)))
    )

  test("roundoff floor: with |f*| ≈ 1e2 either library may stop on a failed line search; neither claims τ falsely") {
    // Strict sufficient decrease must resolve f differences near eps·|f|. Contract asserted for each library: it either
    // reports convergence with residual ≤ τ, or reports a failed line search. Whether a given case reaches τ is not
    // asserted. The 10·sqrt(eps·|f*|·L) limit below is a heuristic guard against a stall far from the minimizer, not a
    // derived bound.
    val mixed = (i: Int) => Seq(-1, 0, 1)(i % 3)
    val cases = Seq(
      unshifted(unconstrainedQuadratic(10, 10.0, 11L), 10.0, 11L),
      unshifted(boxQuadratic(30, 1e2, 31L, mixed), 1e2, 31L),
      unshifted(boxQuadratic(12, 1e2, 53L, mixed, lower => lower.clone()).copy(name = "box, start = l"), 1e2, 53L)
    )
    for p <- cases do
      val tau = threshold(p)
      val gale = runGale(p)
      val breeze = runBreeze(p)
      val galeResidual = residual(p, gale.x, p.oracle(gale.x)._2)
      val breezeResidual = residual(p, breeze.x, p.oracle(breeze.x)._2)
      if gale.converged then assert(galeResidual <= tau, s"${p.name}: gale claimed convergence at $galeResidual")
      else assertEquals(gale.status, "LineSearchFailed", s"${p.name}: gale stop")
      if breeze.converged then assert(breezeResidual <= tau, s"${p.name}: breeze stopped at $breezeResidual")
      else assertEquals(breeze.status, FirstOrderMinimizer.SearchFailed.reason, s"${p.name}: breeze stop")
      val floor = 10.0 * math.sqrt(2.2e-16 * math.max(1.0, math.abs(gale.value)) * p.lipschitz)
      assert(galeResidual <= floor, s"${p.name}: gale residual $galeResidual above heuristic guard $floor")
      assert(breezeResidual <= floor, s"${p.name}: breeze residual $breezeResidual above heuristic guard $floor")
      floorReport += f"| ${p.name} | ${gale.value}%.2f | $floor%.1e | ${gale.status} | $galeResidual%.1e " +
        f"| ${breeze.status} | $breezeResidual%.1e |"
  }

  test("infeasible starts are rejected by both libraries; the projected start then agrees") {
    val p = boxQuadratic(8, 10.0, 61L, i => Seq(-1, 0, 1)(i % 3))
    val (lower, upper) = p.bounds.get
    val infeasible = Array.tabulate(p.n)(i => if i % 2 == 0 then 3.0 else -3.0)
    val objective = DifferentiableObjective(p.n)(x =>
      val (v, g) = p.oracle(Array.tabulate(p.n)(x(_, 0)))
      Right(ObjectiveEvaluation(v, column(g)))
    ).toOption.get
    val box = MatrixBoxBounds.from(column(lower), column(upper)).toOption.get
    LBFGSB.minimize(objective, box, column(infeasible), LBFGSBConfig(tolerance = sharedTolerance)) match
      case Left(FirstOrderError.InvalidConfiguration(_)) => ()
      case other                                         => fail(s"gale accepted an infeasible start: $other")
    intercept[IllegalArgumentException]:
      BreezeLBFGSB(BDV(lower), BDV(upper)).minimize(diffFunction(Counted(p.oracle)), BDV(infeasible))
    val projected = Array.tabulate(p.n)(i => math.max(lower(i), math.min(upper(i), infeasible(i))))
    compare(p.copy(name = "box, projected infeasible start", start = projected))
  }

  test("informational: Breeze default stopping rules versus the shared rule") {
    // Not assertions about agreement: they record how far default-configured Breeze stops from x*.
    val quad = unconstrainedQuadratic(100, 1e4, 13L)
    val defaultLbfgs = new BreezeLBFGS[BDV[Double]]()
    val lbfgsCounted = Counted(quad.oracle)
    val lbfgsState = defaultLbfgs.minimizeAndReturnState(diffFunction(lbfgsCounted), BDV(quad.start.clone()))
    val lbfgsX = lbfgsState.x.toArray
    val box = boxQuadratic(100, 1e4, 32L, i => if i % 4 == 0 then -1 else if i % 4 == 1 then 1 else 0)
    val (lower, upper) = box.bounds.get
    val boxCounted = Counted(box.oracle)
    val boxState = BreezeLBFGSB(BDV(lower), BDV(upper))
      .minimizeAndReturnState(diffFunction(boxCounted), BDV(box.start.clone()))
    val lbfgsError = distance(lbfgsX, quad.known.get)
    val boxError = distance(boxState.x.toArray, box.known.get)
    val boxResidual = residual(box, boxState.x.toArray, boxState.grad.toArray)
    report += ""
    report += "Breeze with default settings (informational):"
    report += ""
    report += "| Problem | Breeze call | ‖x − x*‖₂ | exit residual | stop reason | callbacks |"
    report += "| --- | --- | ---: | ---: | --- | ---: |"
    report += f"| ${quad.name} | `LBFGS()` (maxIter −1, m 7, tol 1e-9) | $lbfgsError%.1e | " +
      f"${residual(quad, lbfgsX, quad.oracle(lbfgsX)._2)}%.1e | " +
      s"${lbfgsState.convergenceReason.map(_.reason).getOrElse("none")} | ${lbfgsCounted.calls} |"
    report += f"| ${box.name} | `LBFGSB(l, u)` (maxIter 100, m 5, tol 1e-8) | $boxError%.1e | $boxResidual%.1e | " +
      s"${boxState.convergenceReason.map(_.reason).getOrElse("none")} | ${boxCounted.calls} |"
    assert(lbfgsError.isFinite && boxError.isFinite)
  }

  override def afterAll(): Unit =
    val header = Seq(
      "# Optimizer parity: Gale vs Breeze 2.1 (informational timings)",
      "",
      s"Shared rule: ‖r‖∞ ≤ $absoluteTolerance + $relativeTolerance · max(1, ‖r(x₀)‖∞); history m = $history " +
        "(rows marked rel use 1e-10 + 1e-11 · max(1, ‖r(x₀)‖∞)). " +
        "U = LBFGS, B = LBFGSB. Times are single cold runs in one forked JVM; do not read them as benchmarks.",
      "",
      "| Problem | Kind | τ | ‖Δx‖₂ | Δx bound | \\|Δf\\| | Gale callbacks | Breeze callbacks | Gale ms | Breeze ms |",
      "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
    )
    val floor = Seq(
      "",
      "Roundoff floor with unshifted quadratics (|f*| ≈ 1e2; τ = 1e-8):",
      "",
      "| Problem | f* | 10·sqrt(eps·\\|f*\\|·L) (heuristic guard) | Gale status | Gale residual | Breeze status | Breeze residual |",
      "| --- | ---: | ---: | --- | ---: | --- | ---: |"
    )
    val text = (header ++ report ++ floor ++ floorReport).mkString("\n") + "\n"
    println(text)
    val out = Paths.get("target", "optimizer-parity.md")
    Files.createDirectories(out.getParent)
    Files.write(out, text.getBytes(StandardCharsets.UTF_8))
    ()
