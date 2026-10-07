# Optimization

`gale.optim` provides smooth unconstrained, convex proximal/projected, and
operator-based composite optimization, plus a specialized constrained-Rayleigh
helper. It owns iteration, stopping, typed failures, and final-point diagnostics.
Models and application-specific solver plans belong to the caller.

The same implementation runs on the JVM and Scala.js. Variables are `DMat`
values. The entire matrix is one variable: callbacks may couple rows and
columns. Columns represent independent problems only when the supplied
objective and proximal/projection callbacks are separable across columns.

## Choose the method from the objective structure

| Objective structure | Gale method |
| --- | --- |
| smooth unconstrained objective, no global bound available | `LBFGS.minimize` |
| convex smooth term plus convex proximal term or set | `AcceleratedProximal.minimize` |
| smooth term plus a directly proximable term with a known bound | `proximalGradient` |
| smooth term plus a projection set | `projectedGradient` |
| proximable primal term plus `g(Kx)` | `linearCompositePrimalDual` |
| smooth term plus direct proximal term plus `g(Kx)` | `smoothCompositePrimalDual` |
| exact linear null-space parameterization | `ExactLinearReduction.verify` |

`FirstOrderCapabilities.select` lets a host library separate model lowering
from the methods available in a runtime. Automatic selection is deterministic;
a required but unavailable method is a typed `MissingCapability`.

## Proximal-gradient example

This example minimizes

```text
0.5 ||x - c||² + 0.25 ||x||₁
```

whose independent solution is soft-thresholding `c` by `0.25`.

```scala mdoc
import gale.linalg.DMat
import gale.optim.*

val center = DMat.dense(2, 1, Seq(1.0, -2.0))

val smooth = new SmoothObjective:
  val variableRows = 2
  val lipschitz = 1.0

  def value(at: DMat): Either[FirstOrderError, Double] =
    var squared = 0.0
    var row = 0
    while row < at.rows do
      val difference = at(row, 0) - center(row, 0)
      squared += difference * difference
      row += 1
    Right(0.5 * squared)

  def gradient(at: DMat): Either[FirstOrderError, DMat] =
    Right(DMat.tabulate(at.rows, at.cols): (row, column) =>
      at(row, column) - center(row, column)
    )

val l1 = ProximalTerms.l1(2, 0.25).toOption.get

val fit = FirstOrderSolvers
  .proximalGradient(smooth, l1, DMat.zeros(2, 1))
  .fold(error => throw new IllegalStateException(error.message), identity)

fit.primal(0, 0)
fit.primal(1, 0)
fit.status
```

Objective and proximal callbacks return `Either[FirstOrderError, _]`. A
non-finite callback value, mismatched result shape, invalid configuration, or
failed Gale operator application remains a typed failure.

## Linear operators and norm bounds

Primal-dual methods accept a `BoundedLinearOperator`: a
`DoubleLinearOperator` paired with a finite non-negative upper bound on its
induced norm. Dense `DMat` values already implement `DoubleLinearOperator`;
matrix-free, sparse, block, and composed operators use the same contract.

The bound is evidence supplied by the caller, not estimated by the solver. It
controls the safe primal and dual steps, so an invalid underestimate voids the
method's convergence premise. Construction validates the representation of the
claim, while the model layer remains responsible for how the bound was derived.

Operator application returns `FirstOrderError.OperatorFailure` with the
original `LinAlgError` retained as `cause`. Gale never flattens that failure
into an untyped exception string.

## Read the result and certificate

A successful numerical execution returns `FirstOrderSolution`, including:

- the owned primal value and optional dual value;
- the final objective;
- `Converged`, `IterationLimit`, `EvaluationLimit`, `Cancelled`,
  `NumericalStagnation`, or `LineSearchFailed`;
- fixed-point residuals, objective change, iteration count, and exact settings;
- diagnostic value summaries and separate exact-content binding to the returned matrices;
- numerical-resolution estimates for the fixed-point residuals.

Iteration exhaustion is a successful execution with
`FirstOrderStoppingStatus.IterationLimit`, not a false convergence claim.
`FirstOrderCertificate.binds` compares the bound matrix contents and detects
substitution of either returned value, including permutations with identical
summaries. A certificate constructed from summaries alone does not bind a result.
Residuals describe the implemented fixed-point equations at the reported
settings; application-specific optimality or statistical validity remains the
host library's responsibility.

`ExactLinearReduction.verify` separately checks that a proposed basis lies in
the declared constraint null space. Its `LinearReductionCertificate` records
the basis image, constraint image, measured maximum residual, and tolerance
threshold.

## Numerical and domain assumptions

Proximal-gradient and convex primal-dual guarantees require convex objectives,
exact proximal maps, consistent adjoints, and valid upper bounds. Projected
gradient also accepts nonconvex projections for existing local workflows; its
fixed-point diagnostic does not establish a global minimum. Both primal-dual
methods require `extrapolation = 1`; other values are rejected before callbacks.
For the smooth splitting, the steps satisfy the coupled condition
`1/tau - sigma * normBound² > lipschitz/2` in exact arithmetic. Derived steps
must also be finite and strictly positive in floating-point arithmetic.

Objective components must be finite at every point at which the solver evaluates
them. Direct indicator terms can be used with a feasible initial point and a
feasibility-preserving proximal map. General extended-valued composite indicators
are outside this contract: an infeasible `Kx` yielding positive infinity is a
typed failure, not a finite objective or convergence claim.

Residuals are computed at the returned point. Verification uses
`min(iterationStep, 1)` so an arbitrarily large step cannot dilute a constrained
mapping. `primalResidualStep` and `dualResidualStep` expose those verification
steps, separately from the steps used to advance the solver. This is a
coordinate-dependent fixed-point convention, not a scale-invariant KKT bound. Their relative scale comes from
the first checked residual (at initialization for single-variable methods, after
the first primal-dual update for splitting), with a minimum reference scale of
one, rather than the absolute parameter coordinates. Numerical-resolution
estimates compare coordinate ulps with the chosen step; they are conservative
safeguards against rounded-away steps, not proved bounds on errors inside user
callbacks. A small residual with unresolved numerical precision cannot certify
convergence. Objective change is a diagnostic; it is never a substitute for
stationarity. Optimizer variables must have at least one row and column.

`ExactLinearReduction.verify` establishes only containment in a null space. A
zero or incomplete basis can pass; independence, rank, and completeness require
separate evidence from the caller.

## Smooth minimization without a Lipschitz bound

L-BFGS uses a fused value/analytic-gradient callback, bounded strong-Wolfe line
search, and limited curvature history. It supports a whole matrix variable;
`DifferentiableObjective.vector` adapts a single-column `DVec` callback.
It is unconstrained: use a proximal/projected method for bounds. For nonconvex
objectives its gradient test establishes local stationarity, not a global optimum.

```scala mdoc
val differentiable = DifferentiableObjective.fromSmooth(smooth)
val unconstrained = LBFGS.minimize(differentiable, DMat.zeros(2, 1)).toOption.get
assert(unconstrained.status == FirstOrderStoppingStatus.Converged)
assert(math.abs(unconstrained.primal(0, 0) - 1.0) < 1e-6)

val gradientCheck = GradientCheck.directional(
  differentiable, DMat.zeros(2, 1), DMat.dense(2, 1, Seq(0.3, -0.7))
).toOption.get
assert(gradientCheck.passed)
```

For an expensive callback, implement `DifferentiableObjective(rows)(callback)`
directly and return `ObjectiveEvaluation(value, gradient)` to share intermediate
work. The `fromSmooth` adapter invokes both existing callbacks. Directional
finite differences are an explicit debugging aid; solvers never silently replace
an analytic gradient with finite differences.

`DenseObjectives.leastSquares(design, response)` supplies the objective
`0.5 / m * ||design * x - response||²`. `DenseObjectives.logistic(design,
labels, l2)` supplies mean logistic loss plus `0.5 * l2 * ||x||²`, with labels
exactly `-1` or `+1`. Both use Gale's dense kernels, copy the input data when
constructed, and require a single-column parameter. The least-squares objective
retains an additional contiguous copy of the transposed design for repeated
gradient products; its design storage is about twice the input matrix size. Logistic evaluation stays
finite at saturated finite margins; nonfinite intermediate arithmetic is a
typed failure. Construction belongs outside repeated solve timings.

## Convex acceleration and standard terms

Accelerated proximal gradient backtracks to a local majorization condition and
checks that condition at the point used for convergence. An initial Lipschitz
estimate is a starting guess, not a certified global bound. Optional objective
and gradient restart are enabled by default; `restart = false` disables both.
Accepted objective/gradient evaluations are reused, and the convergence check
does not require another trial objective evaluation. The method requires convex smooth and proximal
components; this interface does not verify convexity or promise an accelerated
rate for arbitrary callbacks/restarts.

```scala mdoc
val accelerated = AcceleratedProximal.minimize(
  differentiable, l1, DMat.zeros(2, 1),
  AcceleratedConfig(initialLipschitz = 0.1)
).toOption.get
assert(accelerated.status == FirstOrderStoppingStatus.Converged)
assert(math.abs(accelerated.primal(0, 0) - 0.75) < 1e-6)

val box = ProjectionSets.box(2, -0.5, 0.5).toOption.get
val bounded = AcceleratedProximal.minimize(
  differentiable, box.indicator, DMat.zeros(2, 1)
).toOption.get
assert(box.violation(bounded.primal).toOption.get == 0.0)
```

`ProximalTerms` supplies zero, scalar/elementwise weighted L1, squared L2, and a
conjugate-prox adapter. Catalog terms use direct stable conjugate formulas.
The general Moreau adapter rejects reconstructions whose cancellation cannot
resolve a useful result; provide a direct `proximalConjugate` in that case.
`ProjectionSets` supplies finite box bounds, nonnegativity, and a simplex per
column. Matrix weights/bounds require exact shape; no implicit broadcasting.
Indicator terms require a feasible start, within a tolerance of `8 * ulp(max(1, maxAbs(x))) * rows`. `violation` measures the actual bound or simplex error.

## Work limits, progress, and warm starts

```scala mdoc
val limited = LBFGS.minimize(
  differentiable, DMat.zeros(2, 1),
  LBFGSConfig(control = SolverControl(maxEvaluations = 20, traceCapacity = 4))
).toOption.get
limited.evaluations
limited.trace.size
```

`EvaluationCounts` counts attempted callback invocations, including rejected
trials and final checks. One fused callback increments callbacks, values, and
gradients once each. An adapter's internal work is part of that callback.
Forward and adjoint counts cover operator calls, not scalar entries or FLOPs.
`maxEvaluations` caps the total callback count. Progress runs once per checked
iterate; returning false or setting `cancelled()` requests cancellation.
Exceptions from hooks become typed `OracleFailure` values. Traces are off by
default and retain only the last `traceCapacity` events.

A stopped run returns its last fully checked point with a nonconverged status.
If cancellation or a budget prevents even the first check, the result is a typed
`ExecutionStopped` error. Callback errors remain errors. `LineSearchFailed`
means the bounded search found no acceptable step; inspect scaling/gradients or
adjust the search budget. `NumericalStagnation` calls for rescaling or a realistic
tolerance, not relabeling the point as converged.

Both primal-dual methods accept `initialDual = Some(previousDual)` alongside
the supplied primal start. Shapes and finiteness are checked before iteration.
Warm starts reuse values only, not a hidden solver state or stale certificate.
Results and certificates retain owned immutable matrix values; exact binding is
content comparison (linear in variable size), not a hash or an optimality proof.
The default disabled trace does not retain the iterate history. L-BFGS storage
is proportional to `historySize * variableSize`.

For all methods, `absolute + relative * reference` defines the residual
threshold. L-BFGS uses the initial gradient infinity norm as its reference,
with a minimum of one. The proximal methods use their first checked mapping
norm. Use a zero relative tolerance when comparing an absolute KKT target.
The certificate records the reference, algorithm options, and evaluation budget;
the solution records the steps actually used at the endpoint.

## Supported guarantees

| Method | Required structure | Meaning of `Converged` |
| --- | --- | --- |
| L-BFGS | differentiable finite objective, consistent analytic gradient | final gradient meets tolerance; local stationarity |
| Accelerated proximal | convex differentiable smooth term and exact convex prox | resolved final mapping after majorization check |
| Proximal/projected gradient | valid positive global smoothness bound, exact prox/projection | resolved final mapping; convex optimality only with convex components |
| Both primal-dual methods | convex components, true adjoint, valid operator bound, extrapolation 1 | resolved primal and dual fixed-point residuals |
| Constrained Rayleigh | symmetric numerator and symmetric positive-definite denominator | positive-cone stationarity; no global extremality claim |
| Exact linear reduction | linear basis and constraint operators | null-space containment only |

The repository protocol at `docs/verification/optim-python/README.md`
records exact inputs, independent endpoint checks, timing boundaries, and measured
wins and losses. These are bounded workload results, not a universal speed claim.
