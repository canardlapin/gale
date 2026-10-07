# A bounded plan for excellent optimization support

Review date: 2026-10-07. Source revision: `cf1cbe0dcfbc5208e0599de7a2adaef61b890508`.
Status: implemented locally and qualified on 2026-10-07; integrated release acceptance
remains open for downstream status-adapter migration and exact-candidate hosted CI.
The historical review and plan below describe the original baseline.

Implementation/evidence: [optimization qualification](verification/optim-python/README.md).
Mote epic: `bd-01M4BQS2FZZVD18DZTWYK52VBP`; granular tasks OPTIM-01 through OPTIM-19.

Gale should provide a small collection of dependable numerical optimizers:
smooth unconstrained minimization, convex proximal/projected minimization, and
operator-based composite minimization. Keep the existing specialized Rayleigh
helper. The priorities are trustworthy termination, convenient use, and speed
at a verified solution accuracy.

Keep this in `gale.optim` within `gale-core`, portable across JVM and Scala.js.
Add two algorithms: L-BFGS and accelerated proximal gradient with backtracking.
Do not add a modeling language, solver plugin architecture, automatic
differentiation, stochastic/global optimization, general nonlinear constraints,
or an assortment of competing optimizers. L-BFGS-B, trust-region/Newton methods,
nonlinear least-squares solvers, and derivative-free methods can wait for a
concrete consumer requirement.

## What exists and what needs attention

The implementation is compact: four first-order solver entry points, typed
callback/operator failures, caller-supplied smoothness/operator bounds,
capability selection, null-space containment verification, and nonnegative
generalized-Rayleigh optimization. There are 8 first-order tests and 3 Rayleigh
tests. The existing operator abstraction and owned matrix results are useful
foundations; a rewrite would discard working integration unnecessarily.

The [review fixture](verification/optim-review-20261007/README.md) records six
observations reproduced on both JVM and Scala.js:

| Observation | Consequence |
| --- | --- |
| A rounded-away proximal-gradient step reports `Converged` for a quadratic whose actual gradient exceeds the stopping threshold. | Numerical stagnation can be mistaken for stationarity. |
| After one step on `f(x) = x²/2`, the returned point is approximately `0.01` but the certificate residual is `1`, from the preceding point. | The certificate does not consistently describe the returned iterate. |
| Different matrices with permuted entries pass `certificate.binds`. | Sum, squared norm, and maximum magnitude do not provide exact value binding. |
| Finite accepted options can produce a zero step and an uncaught constructor exception. | The numerical failure boundary is not completely represented by `Either`. |
| An asymmetric Rayleigh numerator can produce a zero reported stationarity residual at a point with a positive feasible directional derivative. | Symmetry assumptions need enforcement or explicit symmetric-part semantics. |
| A zero, incomplete basis passes `ExactLinearReduction.verify`. | This checks null-space containment, not rank or completeness of a parameterization. |

The first five are repair targets. The last is primarily a contract boundary:
do not turn a containment check into a costly general null-space construction
unless a consumer needs that separate operation.

Source inspection also shows repeated full-matrix temporaries and repeated
oracle work. Each current primal-dual iteration uses three forward operator
applications and two adjoint applications, including the residual check.
Some are duplicates that can be reused without weakening final diagnostics.
There are no optimization-specific benchmarks in the current benchmark tree.
These are performance opportunities, not measured speedups.

Real consumers constrain the design. Multivar uses the first-order contracts
through its solver compilers and uses the Rayleigh helper for constrained
canonical methods. ScalaFIM's `StructuredPatternOptimizer.TargetProblem` couples
columns and supplies a nonconvex polar projection. A matrix is therefore one
optimization variable unless separability is explicitly supplied; the guide's
claim that columns are independent is too broad. Preserve these low-level uses
and distinguish local fixed-point evidence from convex convergence guarantees.

## 1. Repair termination and mathematical contracts first

- Convert the reproduced defects into permanent regression tests asserting the
  corrected behavior. Compute reported residuals at the returned primal/dual
  point. Keep objective change as supporting evidence, never a replacement for
  stationarity. Distinguish iteration/evaluation exhaustion, numerical
  stagnation, and convergence.
- Detect zero, infinite, or ineffective steps and overflow in solver-owned
  calculations. Return typed errors or a nonconverged result as appropriate.
  For smooth unconstrained problems, check the actual gradient; when a tiny
  proximal displacement cannot resolve stationarity, do not certify convergence.
  Define residual scaling explicitly and test translated/rescaled problems.
- Replace summary-based exact binding with owned-value binding or exact content
  comparison. Keep summaries as diagnostics. Ensure the solution's objective,
  status, settings, and certificate agree.
- Write down the convexity, exact-prox/projection, adjoint, and bound assumptions
  for each guaranteed method. Audit the actual update order and the supported
  extrapolation range against the algorithm reference. For Condat--Vu, include
  the coupled step condition, not just two separate upper bounds. Preserve the
  existing nonconvex projection path with explicitly narrower claims.
- Decide how extended-valued convex terms behave: indicator objectives may be
  `+Infinity` outside their domain, which is different from `NaN` or numerical
  overflow. Start with clear feasible-start rules for built-in projections;
  support composite indicator domains deliberately before claiming general
  extended-valued support.
- For Rayleigh, enforce/document symmetric matrices and an admissible positive
  denominator geometry, check finite arithmetic, and retain the distinction
  between stationary and globally best solutions. Keep the helper specialized.
  Clarify that exact reduction verification establishes containment only.

Acceptance: the new regressions pass on both platforms; independently
recomputed final residuals agree; scaled, boundary, exhausted-budget, and
malformed-oracle cases return truthful outcomes. Review API changes against
`docs/api-stability.md` and compile representative existing consumers.

## 2. Make common problems easy to express

- Add small constructors for callback objectives, including fused value and
  gradient evaluation. A caller using line search should not have to invent a
  global Lipschitz bound. Retain the existing bound-based interface and adapt it
  rather than imposing a new hierarchy on every consumer.
- Supply a deliberately small numerical catalog: zero term, weighted L1,
  squared L2, box projection, nonnegative projection, and simplex projection.
  Provide a checked conjugate-prox adapter using Moreau's identity where valid.
  Do not build an arbitrary penalty-expression language or compose proximal
  maps as though every sum had an exact composed prox.
- Provide a convenient vector-facing entry point backed by the same numerical
  implementation. Define whether simplex/group operations act on the whole
  variable or on individual columns. Keep coupled matrix objectives supported.
- Add evaluation counts, bounded evaluation budgets, final accepted step sizes,
  and a lightweight progress/cancellation hook. Traces remain opt-in and bounded.
  Count extra work spent checking convergence as well as advancing iterates.
- Add a small directional gradient-check utility for debugging analytic
  callbacks. Finite-difference gradients should not silently become the solver
  default. Allow primal/dual warm starts for repeated composite solves without
  introducing persistent solver sessions.

Acceptance: common L1 and box-constrained examples require no handwritten
proximal/projection loops; ownership, dimensions, cancellation, and evaluation
counts are tested. Old callbacks remain usable, including coupled columns.

## 3. Add two complementary algorithms

**Accelerated proximal gradient (FISTA), with backtracking.** Reuse the repaired
proximal kernel and standard terms. Retain a supplied-bound fast path, add
bounded backtracking when a useful bound is unavailable, and use an explicitly
documented restart/monotonicity policy. Acceleration is for the supported convex
problem class; it should not silently apply to arbitrary nonconvex projections.
Always check the accepted primal iterate, not the extrapolated point.

Acceptance: analytic quadratic/L1 solutions, ill-conditioned convex examples,
backtracking/restart behavior, and active constraints. Record evaluation counts
and time to the same independently checked accuracy against ordinary proximal
gradient. The [FISTA paper](https://www.tau.ac.il/~becka/FISTA.pdf) supplies the
algorithmic basis; performance gains in Gale still need measurement.

**L-BFGS for smooth unconstrained objectives.** Use bounded memory, fused
value/gradient callbacks, a safeguarded strong-Wolfe line search, curvature-pair
checks, and descent-direction recovery. Make failed line search, evaluation
limits, non-finite trial points, and stagnation explicit. Allow local nonconvex
use with a stationary-point claim. Do not implement L-BFGS-B by clipping an
unconstrained step.

Acceptance: rotated and ill-conditioned SPD quadratics with direct-solve
oracles, Rosenbrock, stable logistic-loss examples, and adversarial line-search
and curvature cases. Compare objective/gradient accuracy with an independent
reference; identical iteration trajectories are unnecessary. The original
[L-BFGS implementation and references](https://users.iems.northwestern.edu/~nocedal/lbfgs.html)
provide the bounded smooth-optimization model to follow.

## 4. Measure and remove avoidable cost

Capture a benchmark baseline before changing hot loops, alongside stage 1.
After semantics are settled:

- Fuse affine updates and residual reductions instead of allocating matrices
  for every `add`, `scale`, and `subtract`. Reuse internal scratch storage while
  preserving owned results and the existing safe callback API.
- Cache operator images, adjoint results, and gradients at unchanged points.
  Use linearity to update extrapolated operator images where beneficial. Verify
  reduced call counts directly. Keep a full check at the returned point.
- Add destination-writing callback/operator paths only if profiles show they
  are necessary after internal fusion and caching. Reuse Gale's established
  destination/workspace conventions; avoid exposing mutable backing arrays.
- Cover small cheap objectives and larger operator-dominated problems, one and
  several columns, and dense, sparse, and matrix-free operators. Include the
  specialized Rayleigh method at a modest size; cache its fixed matrix scales
  before considering more ambitious algorithmic changes.

Acceptance: JMH measurements with allocation profiling and an optimized Node
benchmark report time, bytes where measurable, evaluations, residual quality,
and termination together. Fix inputs, accuracy targets, warmup, hardware, and
runtime versions. Show lower allocations/callback work and a repeatable benefit
on representative cases, with regressions explained. Never obtain a favorable
timing by accepting a worse solution or skipping final diagnostics. Memory
should remain bounded with iteration count when tracing is disabled.

## 5. Qualify the supported surface and finish the documentation

Maintain a compact, fixed acceptance set: analytic soft-threshold/box/simplex
cases; rotated quadratic systems; Rosenbrock; logistic loss; a small TV/fused
lasso problem; and the existing Rayleigh angle oracle. Add deterministic
scaling, active-set, empty/degenerate-shape, multi-column, failure, and ownership
variants where they expose a distinct risk. Store reference provenance and
tolerances. External reference tools belong in fixture generation, not runtime
dependencies. For convex examples, check KKT conditions or a duality gap when
available, not just agreement with another optimizer's success flag.

Publish one method-selection table and executable examples for smooth,
L1-regularized, box-constrained, and operator-composite problems. Explain
stationarity versus objective-error guarantees and convex versus local use.
Correct the Breeze comparison guide, which still lists optimization wholesale
as outside Gale's scope. Include benchmark recipes and a short failure guide.

Use the smallest relevant checks per change, then the repository gates for the
integrated result:

```sh
sbt 'coreJVM/testOnly gale.optim.*' 'coreJS/testOnly gale.optim.*'
sbt scalafmtCheckAll compileAll testAllFull parityTest interopBreezeTest docsCheck
sbt benchCompile
```

Benchmark compilation alone is not performance evidence: run the new named JMH
cases and optimized JS workloads and retain results. Required CI uses JDK 21
and Node 22; optional native/Wasm lanes are not prerequisites for this portable
scope. Check source/binary compatibility as required by the current release
boundary and test the selected Multivar/ScalaFIM call sites against the exact
candidate. Do not call a local check a hosted-CI or downstream qualification.

## Delivery boundary

Aim for roughly 6–8 reviewable changes: contract repairs, common API/helpers,
FISTA, L-BFGS, measured performance work, and final qualification/docs, splitting
the correctness work when needed. Stage 1 comes first; baseline measurements
start early; each algorithm brings its own tests and examples.

Stop when these methods solve the declared problem classes reliably, common
tasks are concise, and representative JVM/JS timing and allocation evidence is
available. Additional optimizer families need a concrete unmet use case. The
next implementation action is the rounded-step false-convergence regression
and repair, followed by final-point residual and certificate consistency.

For the primal-dual audit, use the assumptions and coupled step condition in
[Condat's original paper](https://lcondat.github.io/publis/Condat-optim-JOTA-2013.pdf).
The present review does not establish that every currently accepted
extrapolation value has that convergence guarantee.
