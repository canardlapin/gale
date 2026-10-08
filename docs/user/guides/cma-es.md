# Derivative-free optimization with CMA-ES

Use `CMAES` when you can evaluate a continuous objective but cannot supply useful
derivatives. It adapts a Gaussian search distribution using ranked candidate
values, a full covariance model, and cumulative step-size adaptation. This is
useful for simulation-based, nonsmooth, or noisy objectives in modest dimensions.
The implementation runs in portable `gale-core` on JVM and Scala.js.

For smooth problems with reliable derivatives, start with
[L-BFGS or accelerated proximal gradient](first-order-optimization.md), or
[LM and L-BFGS-B](nonlinear-and-bounded-optimization.md). CMA-ES typically needs
more objective evaluations. Its dense covariance storage is quadratic in the
number of free coordinates, with periodic cubic eigendecompositions.

## Minimize an objective

Supply an initial mean, a meaningful initial step size in your parameter units,
and an evaluation budget. Scale coordinates before calling the solver when their
natural units differ greatly. The initial covariance is identity.

```scala mdoc
import gale.linalg.DVec
import gale.optim.*

def energy(x: DVec): Double =
  math.abs(x(0) - 0.25) + 3.0 * math.abs(x(1) + 0.5)

val config = CMAESConfig(
  initialStepSize = 0.8,
  seed = 42L,
  maxEvaluations = 10000,
  targetValue = Some(1e-7)
)
val fit = CMAES.minimize(energy, DVec.fromSeq(Seq(2.0, -2.0)), config).toOption.get
assert(fit.stopReason == CMAESStop.TargetReached)
assert(energy(fit.best.get.coordinates) <= 1e-7)
assert(fit.work.evaluations <= config.maxEvaluations)
```

`best` is the best finite observation, with its evaluated coordinates and value.
It can be absent if cancellation, sampling failure, or an objective failure occurs
before a finite observation. The initial mean is not evaluated separately.
`distribution.mean` is an adapted search parameter and need not have been evaluated.

`TargetReached` means an observed value meets the supplied target. With noisy
objectives, assess the returned point using independent repeated measurements.
There is no noise estimator or automatic reevaluation policy in this version.

## Separate sampling from evaluation

The ask/tell interface lets you evaluate a generation using your own scheduler.
Evaluate points concurrently if appropriate, then supply their scores in the
original batch order. The session itself has one owner and is not thread-safe.

```scala mdoc
val search = CMAES.start(
  DVec.fromSeq(Seq(1.0, 1.0)),
  CMAESConfig(seed = 7L, maxEvaluations = 40, stagnationGenerations = 0)
).toOption.get
while search.stopReason.isEmpty do
  search.ask().toOption.get.foreach { batch =>
    val values = batch.points.map(energy)
    search.tell(batch, values).toOption.get
  }
assert(search.work.evaluations == 40)
assert(search.bestPoint.nonEmpty)
```

Only one batch may be outstanding. A foreign, stale, duplicate, wrong-length,
or non-finite receipt is rejected without consuming it, so a caller can correct
an invalid receipt. The coordinates in the batch are owned immutable values.

The last budget-limited batch may be smaller than the normal population. Its
finite observations contribute to `best` and evaluation counts; a partial
population never adapts the mean, covariance, or evolution paths. To cancel an
outstanding batch, call `search.cancel(batch, completedFitness)` with the finite
scores for the completed prefix. If parallel work finishes out of order, account
for that work in your scheduler before submitting a prefix receipt.

Ask/tell counts scores in accepted receipts; it cannot count objective work done
outside that protocol. `minimize` counts every attempted callback, including a
callback that throws or returns a non-finite value. Sampling rejections do not
invoke the objective. Work also records complete generations, sampled proposals,
and rejected bound proposals.

## Bounds and fixed parameters

`MatrixBoxBounds` must have one column and match the vector length. Starts must
be feasible. Equal finite bounds remove fixed coordinates from the covariance
model; an all-fixed box needs one objective evaluation.

```scala mdoc
val box = MatrixBoxBounds.uniform(2, -1.0, 1.0).toOption.get
val bounded = CMAES.minimize(
  energy,
  DVec.zeros(2),
  config.copy(seed = 18L),
  bounds = Some(box)
).toOption.get
assert(bounded.best.get.coordinates.toSeq.forall(v => v >= -1.0 && v <= 1.0))
assert(bounded.best.get.value <= 1e-7)
```

Bounds use rejection sampling over complete vectors of free coordinates.
Candidates are never clipped to a face. `maxSamplingAttempts` caps total proposal
draws per generation; exhaustion yields `SamplingLimit` and discards that
unevaluated batch. Very narrow boxes or high-dimensional corners can exhaust
this cap. Infinite open sides are supported.

The reported covariance is the dimensionless adapted unconstrained model on
`distribution.freeCoordinates`. Sampling uses the most recently decomposed model,
which may lag the reported matrix. Bounds further condition the accepted
samples. `eigenUpdatePeriod = 1` refreshes every generation; the default selects
an amortized interval from the learning rates.

## Work controls, restarts, and stopping

`CMAESControl` supplies cancellation, progress, and bounded traces. Cancellation
is checked before sampling a generation and before each objective evaluation.
Returning false from a progress callback cancels the entire solve, including any
planned restarts. Callback exceptions produce `ControlFailure` or
`ObjectiveFailure` with preceding finite observations preserved.

Restarts are disabled by default. Pass `CMAESRestarts(maxRestarts = 2)` to restart
a run after stagnation, small steps, excessive condition number, unrepresentable
steps, or sampling exhaustion. Each restart uses the original mean and step
size, identity covariance, zero paths, a distinct derived seed, and a larger
population (twice as large by default). Evaluation and generation limits remain
global. `result.runs` retains each run's seed, population, work, and stopping
reason. Restarts increase exploration but do not certify a global optimum.

`FitnessStagnation` counts generations without a strict best-value improvement;
`stagnationGenerations = 0` disables it. `StepTolerance` compares the largest
coordinate standard deviation of the adapted model against `stepTolerance`;
zero disables it. Condition checks occur when the covariance is decomposed.
Invalid or non-finite adaptation states stop explicitly. No eigenvalue floor or
silent covariance reset is used to hide a numerical failure.

These reasons describe why sampling stopped. Small steps, stagnation, and a
small observed value do not certify a zero gradient or optimality. CMA-ES uses
its own result and progress types rather than first-order certificates.

Seeds reproduce a trajectory within a runtime and version. JVM and Scala.js
can differ in transcendental functions and eigendecomposition rounding, so
cross-platform bitwise trajectories are not promised. Tie scores retain batch
order. Explicit restarts are repeatable for the same configuration.

## Numerical qualification

The implementation uses positive logarithmic recombination weights, rank-one
and rank-mu covariance updates, and cumulative step-size adaptation. Active
negative covariance weights, surrogate models, integer parameters, and general
nonlinear constraints are outside this version.

The local qualification checks analytic objectives, update parity with a pinned
independent Python implementation, bounded and partial-generation behavior, and
multiple seeded solves against `pycma`. Compare attained objective quality and
success rates before comparing runtimes; an inexpensive failed solve is not a
speed advantage. See the repository's `docs/verification/optim-cmaes/README.md`
for exact configurations, adapters, raw samples, and limitations.
