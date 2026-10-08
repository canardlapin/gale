# Smooth optimization with general constraints

`AugmentedLagrangian.minimize` minimizes a smooth scalar objective subject to
nonlinear equalities, inequalities, and optional box bounds. It runs in portable
`gale-core` on JVM and Scala.js. Use it when constraints couple parameters and
you have analytic objective gradients and constraint Jacobians.

The problem is `min f(x)` subject to `h(x) = 0`, `g(x) <= 0`, and optionally
`lower <= x <= upper`. Variables are a single-column `DMat`. Equality rows come
first in the constraint callback; inequality rows follow. The Jacobian has one
row per constraint and one column per parameter.

## Fit subject to an equality

The closest point to `(3, 2)` on `x + y = 1` is `(1, 0)`:

```scala mdoc
import gale.linalg.{DMat, DVec}
import gale.optim.*

val objective = DifferentiableObjective(2)(x => {
  val dx = x(0, 0) - 3.0
  val dy = x(1, 0) - 2.0
  Right(ObjectiveEvaluation(0.5 * (dx * dx + dy * dy), DMat.dense(2, 1, Seq(dx, dy))))
}).toOption.get

val constraints = NonlinearConstraints(2, equalityCount = 1, inequalityCount = 0)(x =>
  Right(ConstraintEvaluation(
    DVec.fromSeq(Seq(x(0, 0) + x(1, 0) - 1.0)),
    DMat.dense(1, 2, Seq(1.0, 1.0))
  ))
).toOption.get

val fit = AugmentedLagrangian.minimize(objective, constraints, DMat.zeros(2, 1)).toOption.get
assert(fit.status == AugmentedLagrangianStatus.Converged)
assert(math.abs(fit.primal(0, 0) - 1.0) < 1e-7)
assert(math.abs(fit.primal(1, 0)) < 1e-7)
assert(math.abs(fit.multipliers(0) - 2.0) < 1e-7)
fit.diagnostics
```

The multiplier sign follows `L = f + lambda*h + mu*g`, with inequality
multipliers nonnegative. Optional `initialMultipliers` use these original
constraint units. Supplying a previous `fit.primal` and `fit.multipliers` warm
starts a related solve; curvature history and the penalty are initialized anew.

## Accuracy and scaling

`ConstrainedDiagnostics` records raw feasibility, scaled feasibility,
stationarity, complementarity, and dual feasibility at the returned point with
the returned multipliers. Convergence requires **all** of these checks:

- Raw and scaled constraint violations satisfy their respective absolute
  tolerances. Box violations contribute to both.
- The infinity norm of the unit box-projected Lagrangian gradient satisfies
  `stationarityTolerance`.
- `max(abs(mu_i * g_i))` satisfies `complementarityTolerance` and all inequality
  multipliers are nonnegative.

Small steps, small objective changes, and convergence of an inner solve are
insufficient. These diagnostics establish approximate first-order conditions;
they do not establish a local minimum, a global minimum, or constraint
qualification. Redundant/degenerate constraints can make multipliers nonunique.
Stationarity uses the supplied parameter coordinates, with no automatic
parameter or objective rescaling.

`constraintScales`, when supplied, divide constraint values and Jacobian rows
by finite positive characteristic magnitudes. Returned constraint values and
multipliers retain original units. `feasibilityTolerance` applies to raw values;
`scaledFeasibilityTolerance` applies after row scaling. Set both deliberately
when constraints have different units. Complementarity is invariant to the
row-scale/multiplier conversion.

The gradient accumulation compensates rounding and uses exact summation of
rounded products when severe cancellation could hide nonstationarity. Products
and caller derivatives still use floating-point arithmetic.

## Execution and failure

The initial point must satisfy any box bounds, and every subsequent callback
point satisfies them. Nonlinear constraints can be violated during a solve;
callbacks must remain finite and defined at those points. Callbacks must be
deterministic for a given point because their last complete evaluation is cached.
An empty constraint set skips the constraint callback entirely.

The solver uses a safeguarded Powell–Hestenes–Rockafellar augmented Lagrangian,
with L-BFGS inside, or L-BFGS-B when bounds are supplied. By default, inner
tolerances adapt to the scaled constraint/complementarity update residual and
never loosen or fall below the requested stationarity tolerance. Zero-step
solves and feasible but nonstationary points force further tightening. Set
`adaptiveInnerTolerance = false` for fixed geometric tightening.

Correction-pair history restarts for each subproblem. A finite positive scalar
inverse-curvature estimate carries across subproblems while the penalty remains
unchanged; it resets on penalty changes and starts anew for every call to
`minimize`. Set `reuseCurvatureScale = false` to disable this reuse. The line
search uses safeguarded interpolation and prefers the ordinary sufficient-decrease
and strong-curvature conditions. This interpolation policy is specific to ALM
inner solves; public L-BFGS/B calls retain their established bisection policy.

Near objective-value resolution, `roundoffAwareLineSearch` permits an approximate
Wolfe step: observed and predicted objective changes must both be within eight
floating-point spacings at the current value, and the approximate-Wolfe derivative
bounds plus the existing strong-curvature test must hold. The allowance is intended
to tolerate rounding; a true increase within that bound can also pass. It never
bypasses curvature at a box boundary or changes the final KKT tolerances. The eight-spacing allowance is
an empirical safeguard, not a certified bound on arbitrary callback error. Set
`roundoffAwareLineSearch = false` for strict sufficient decrease throughout.
The derivative criterion follows [Hager and Zhang, Section 4](https://people.clas.ufl.edu/hager/files/cg_descent.pdf);
Gale retains its own bracketing algorithm.

The penalty grows only when feasibility/complementarity progress fails the configured contraction test. Updated multipliers are used
for diagnostics before clipping the estimates for the next subproblem.

`SolverControl.maxEvaluations` caps **actual attempted user callbacks** across
every inner solve. A fused objective call increments `callbacks`, `values`, and
`gradients`; a fused constraint value/Jacobian call increments `callbacks` and
`jacobians`. Cached recomputation of the augmented objective adds no user call.
Progress and bounded traces occur at checked outer iterates. Cancellation is
also checked inside the inner solver and before each user callback.

Budget or cancellation returns a checked endpoint when one is available. If an
interrupted inner solve's endpoint needs a fresh evaluation beyond the budget,
the previous checked outer endpoint is retained. Before the first complete
evaluation, stopping returns a typed `ExecutionStopped` error. Callback failures
remain typed errors. `lastInnerStatus` and total `innerIterations` aid diagnosis.
With `control.traceCapacity > 0`, `innerTrace` retains that many attempted inner
solves, including interruptions. Each entry records the penalty, inner tolerance,
carried scalar, iterations, augmented requests, and actual objective/constraint
calls including endpoint rechecks. `approximateSteps` counts steps accepted through
the roundoff safeguard. The initial outer evaluation is counted only
in the result totals. `endpointChecked = false` means the previous checked outer
point was retained. Scale fields report the carried scalar (1 when reuse is
disabled), not the inner solver's complete curvature model. Tracing does not add
user callbacks.

`PenaltyLimit`, `LineSearchFailed`, `InnerIterationLimit`, and `IterationLimit`
are unsuccessful terminations, not proofs that constraints are infeasible.
Dense Jacobians require storage proportional to parameters times constraints;
there is no sparse Jacobian or Hessian callback in this version.

For other objective structures, see [first-order optimization](first-order-optimization.md),
[nonlinear least squares and bounds](nonlinear-and-bounded-optimization.md), and
[CMA-ES](cma-es.md). The repository's `docs/verification/optim-alm-robust/README.md`
records independent accuracy checks, policy ablations, and measured comparisons
with SciPy. Nearly dependent constraints remain sensitive to scaling and inner
policy; unsuccessful line searches retain checked diagnostics rather than
implying convergence.

Method reference: [Birgin and Martínez, Practical Augmented Lagrangian Methods](https://www.ime.usp.br/~egbirgin/publications/bmsurveyal.pdf).
