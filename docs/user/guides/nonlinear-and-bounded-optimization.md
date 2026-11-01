# Nonlinear fitting and box-constrained minimization

Use `LevenbergMarquardt.minimize` when a model naturally supplies a residual
vector and Jacobian. Use `LBFGSB.minimize` for a smooth scalar objective with
elementwise parameter bounds. Both methods run in portable `gale-core` on JVM
and Scala.js and reuse `SolverControl`, typed failures, and final-point
certificates. These are local solvers; starting points matter for nonconvex
problems.

## Minimize a smooth function

A `DifferentiableObjective` returns the value and the gradient at a point in
one callback. Parameters are a one-column matrix. This is the Rosenbrock
function, whose minimum is at (1, 1):

```scala mdoc
import gale.linalg.{DMat, DVec}
import gale.optim.*

val rosenbrock = DifferentiableObjective(2)(x => {
  val (a, b) = (x(0, 0), x(1, 0))
  val value = (1.0 - a) * (1.0 - a) + 100.0 * (b - a * a) * (b - a * a)
  val gradient = DMat.dense(2, 1, Seq(
    -2.0 * (1.0 - a) - 400.0 * a * (b - a * a),
    200.0 * (b - a * a)
  ))
  Right(ObjectiveEvaluation(value, gradient))
}).toOption.get

val smooth = LBFGS.minimize(rosenbrock, DMat.dense(2, 1, Seq(-1.2, 1.0))).toOption.get
(smooth.status, smooth.primal(0, 0), smooth.primal(1, 0))
assert(smooth.status == FirstOrderStoppingStatus.Converged)
assert(math.abs(smooth.primal(0, 0) - 1.0) < 1e-5)
```

`status` says whether the run converged or stopped for another reason, such as
the iteration limit. Check it before using `primal`.

## Fit a nonlinear residual model

The objective is `0.5 * sum(residual_i²)`. A Jacobian has one row per residual
and one column per parameter. The model below fits exponential decay to
synthetic observations; model construction does not evaluate callbacks.

```scala mdoc
val times = Vector(0.0, 0.5, 1.0, 2.0, 3.0)
val observed = times.map(t => 2.0 * math.exp(-0.4 * t))
val model = LeastSquaresObjective(2, times.size)(
  x => Right(DVec.tabulate(times.size)(i =>
    x(0) * math.exp(-x(1) * times(i)) - observed(i))),
  x => Right(DMat.tabulate(times.size, 2)((i, j) =>
    if j == 0 then math.exp(-x(1) * times(i))
    else -x(0) * times(i) * math.exp(-x(1) * times(i))))
).toOption.get

val fit = LevenbergMarquardt.minimize(model, DVec.fromSeq(Seq(1.0, 0.1))).toOption.get
assert(fit.status == FirstOrderStoppingStatus.Converged)
assert(math.abs(fit.parameters(0) - 2.0) < 1e-5)
assert(math.abs(fit.parameters(1) - 0.4) < 1e-5)
```

`LeastSquaresObjective.finiteDifferences(parameters, residuals)(callback)`
explicitly selects central differences when analytic derivatives are unavailable.
Every perturbation invokes the residual callback through the shared evaluation
budget. The default relative perturbation is `6e-6 * max(1, abs(parameter))`;
choose an appropriate step for differently scaled parameters. Analytic
Jacobians are preferable when available.

`model.weighted(weights)` applies fixed nonnegative observation weights by
multiplying residual and Jacobian rows by `sqrt(weight)`. The returned residuals
and Jacobian include those weights. Zero weights exclude rows; weights never
adapt during a solve.

LM solves the damped linearized problem using pivoted QR of an augmented,
column-scaled Jacobian. It does not form `JᵀJ`. Default scaling tracks column
norms; `parameterScale` supplies positive characteristic parameter magnitudes
instead. The acceptance test compares actual objective reduction with the
prediction for the actual floating-point displacement. Damping grows after a
rejection and falls after a well-predicted improvement.

`fit.solution.certificate.primalResidual` is the raw infinity norm of the final
`Jᵀr`. The solver uses a fixed initial reference for relative tolerance.
`scaledGradientNorm` is a separate diagnostic. Severe cancellation in the
stationarity dot products triggers exact accumulation of their rounded products;
this does not make multiplication exact. Small steps and small objective changes
alone never establish convergence. A nonzero residual at convergence is normal
for noisy data and does not establish a root or a unique fit.

`maxIterations` counts proposed steps, including rejected trials.
`rejectedSteps`, final damping, the final residual/Jacobian, and callback counts
remain inspectable. On evaluation-budget exhaustion after an accepted point
exists, the result retains that fully evaluated point. If the initial point
cannot be evaluated within the budget, the solver returns a typed error.

## Minimize with parameter bounds

L-BFGS-B accepts the same fused value/gradient callback as L-BFGS. The example
has an optimum on one bound and a free optimum in the other coordinate.

```scala mdoc
val objective = DifferentiableObjective(2)(x => {
  val a = x(0, 0) - 2.0
  val b = x(1, 0) - 0.3
  Right(ObjectiveEvaluation(0.5 * (a*a + b*b), DMat.dense(2, 1, Seq(a, b))))
}).toOption.get
val bounds = MatrixBoxBounds.uniform(2, 0.0, 1.0).toOption.get
val bounded = LBFGSB.minimize(objective, bounds, DMat.zeros(2, 1)).toOption.get
assert(bounded.status == FirstOrderStoppingStatus.Converged)
assert(math.abs(bounded.primal(0, 0) - 1.0) < 1e-7)
assert(math.abs(bounded.primal(1, 0) - 0.3) < 1e-7)
```

`MatrixBoxBounds.from(lower, upper)` supplies different bounds for each matrix entry.
These matrix-shaped bounds are distinct from the finite vector `BoxBounds` used
by `BoxQuasiNewton`.
Finite equal bounds fix a coordinate; infinite open sides are supported. The
initial point must already be feasible. Bounds own their input snapshots, and
callbacks are evaluated only at feasible points. A whole matrix remains one
optimization variable, so the objective may couple columns.

The implementation follows L-BFGS-B's generalized Cauchy point, free-variable
quadratic minimization, and feasible line search. The subspace solve uses
conjugate gradients with a model-decrease safeguard. Strong-Wolfe search permits
Armijo-only acceptance at the actual feasible endpoint, where the curvature
condition may be impossible. Unsuitable curvature pairs are skipped.

CG products visit only free coordinates and reuse local scratch space. The
active-coordinate contributions are zero, so this preserves the accumulation
order of the nonzero terms in the full product.

The residual is `||x - clamp(x - gradient, lower, upper)||∞` with a unit step.
It is evaluated using gradient and bound distances to avoid rounded-away
subtraction at large coordinates. Convergence is local box stationarity; global
optimality additionally requires convexity.

History storage is `O(historySize * parameterCount)`. Rebuilding the Hessian's
rank-two factors costs `O(historySize² * parameterCount)`. The generalized Cauchy
path uses sorted breakpoints and incremental updates; severe cancellation can
trigger aggregate recomputation. This is an iterative-subspace implementation,
not an instruction-for-instruction port of SciPy/MINPACK or the Fortran L-BFGS-B
package. The repository reports at `docs/verification/optim-extensions/README.md`
and `docs/verification/optim-performance/README.md` record measured workloads,
the subsequent performance pass, and limitations.

## Use projected BFGS with a vector callback

Choose `BoxQuasiNewton` for a small or moderate finite-box problem when your
callback writes into a supplied gradient array and an n-by-n inverse metric
is affordable. Choose the L-BFGS-B route above when you need limited-memory
history or infinite open bounds. Their bounds, callbacks, and result types
are distinct.

`gale.optim.BoxQuasiNewton` provides portable full-memory projected BFGS for
finite boxes of dimension 1 through 1024. It accepts a value/gradient oracle,
uses a projected Armijo line search, resets the inverse metric when the active
set changes, and falls back to projected steepest descent for a non-descent
direction. It can traverse regions with indefinite objective curvature.

The scalar objective below is minimized at `x = 0.25`, inside `[0, 1]`:

```scala mdoc
val vectorObjective = new BoxDifferentiableObjective:
  val dimension = 1
  def evaluate(x: Array[Double], gradient: Array[Double]): Either[BoxOptimizationError, Double] =
    val delta = x(0) - 0.25
    gradient(0) = delta
    Right(0.5 * delta * delta)

val vectorBounds = BoxBounds.from(Vector(0.0), Vector(1.0)).toOption.get
val vectorConfig = BoxQuasiNewtonConfig.from(
  maxIterations = 20, maxEvaluations = 100, maxLineSearch = 20,
  projectedGradientTolerance = 1e-10
).toOption.get
val vectorFit = BoxQuasiNewton.minimize(
  vectorObjective, vectorBounds, Vector(0.75), vectorConfig
).toOption.get
assert(vectorFit.status == BoxOptimizationStatus.Stationary)
assert(math.abs(vectorFit.point.get.coordinates.head - 0.25) < 1e-10)
```

The checked constructors return `Either`; the example unwraps them because
its inputs are fixed. In reusable code, retain the errors and inspect the
stopping status before using the point. If the initial callback refuses,
`point` is `None`; a later refusal preserves the last finite point.

Construct `BoxBounds` and `BoxQuasiNewtonConfig` through their checked `from`
methods, then call `BoxQuasiNewton.minimize` with an in-box initial vector.
The oracle borrows a read-only coordinate array and must fill every entry of
the supplied gradient array. Invalid input returns an error before an oracle call;
runtime oracle refusal preserves the last finite point and charges the attempt.
Results own their coordinates and gradient and report evaluations, accepted
iterations, rejected proposals, metric resets and the stopping reason.

`Stationary` means the projected gradient meets the declared tolerance at the
returned point. It does not certify positive curvature or a global minimum.
Iteration, evaluation and line-search limits remain distinct from convergence.
Coordinate/energy scaling, independent starts and application-specific acceptance
checks belong to the caller. Solver scratch is one n-by-n inverse metric,
nine length-n Double arrays and two length-n Boolean arrays; it does not grow
with the iteration count. This implementation is not L-BFGS-B.

## Scope and interpretation

LM currently supports unconstrained dense least squares and fixed weights.
Robust adaptive losses, sparse/operator Jacobians, statistical standard errors,
and bounds on the LM parameters are separate extensions. Use the existing
[proximal methods](first-order-optimization.md) for nonsmooth penalties.
None of these solvers supplies automatic differentiation or certifies that callbacks
implement a consistent objective and derivative.

Method references: [L-BFGS-B, Byrd et al.](https://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf),
[nonlinear least squares, Madsen et al.](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf).
