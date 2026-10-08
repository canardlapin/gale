# ALM roundoff safeguard: robustness follow-up

Local qualification on 2026-10-08. This follow-up addresses the two optional
fixed-tolerance/no-curvature-reuse failures disclosed in the
[previous review](../optim-alm-review/README.md). Historical evidence is preserved.

## Result

The previously failing dense equality and nonlinear ball cases now converge.
All four inner policies pass all 16 fixtures on both JVM and fully optimized
Node: **128 independently checked endpoints**. The refreshed comparison also qualifies
**480 timed Gale endpoints and 480 SciPy endpoints (960 total)**.

| Fixed-tolerance/no-reuse case | Before | After | Objective calls after | Independent stationarity |
|---|---|---|---:|---:|
| equality-dense-32-s0 | LineSearchFailed | Converged | 460 | 9.735e-8 |
| ball-64-s1 | LineSearchFailed | Converged | 63 | 1.249e-10 |

The recovered dense case has raw feasibility 3.445e-9. The ball has zero raw
violation and complementarity 1.626e-9. JVM and Node agree on the statuses and
callback counts. Final KKT thresholds, callback budgets, and iteration limits
were not relaxed.

## Change and its limits

The failing searches rejected useful curvature-qualified steps whose reported
objective differences were within 16 floating-point spacings. A fresh-history
bisection restart did not recover either case and was removed; its failed
experiment remains under `probes/`.

ALM now has a default-on `roundoffAwareLineSearch` safeguard. Ordinary Wolfe
acceptance is tried first. An approximate step requires all of:

- Absolute reported objective change <= eight ULPs of the origin objective.
- Absolute predicted change `step * initialDirectionalDerivative` <= the same
  allowance, with no unit-sized absolute floor.
- The Hager–Zhang approximate-Wolfe derivative interval, with `c1 < 0.5`.
- The existing strong-curvature condition, including at a box boundary.

The derivative criterion comes from [Hager and Zhang, Section 4, equations
4.1–4.3](https://people.clas.ufl.edu/hager/files/cg_descent.pdf). Gale uses its own
bracketing algorithm; this adaptation does not inherit the paper's complete
line-search analysis or a universal convergence guarantee.

The allowance is intended to tolerate rounding, and can admit a true increase
within eight ULPs. It is an empirical safeguard, not a certified callback-error
bound. Overflowing value differences fail closed. Final original-problem
feasibility, stationarity, complementarity, and dual-sign checks remain required.
Set `roundoffAwareLineSearch = false` for strict sufficient decrease throughout;
that opt-out retains the known strict-policy limitations. `innerTrace.approximateSteps`
records how often the safeguard actually accepts a step. Across the 16 fixtures,
JVM recorded 5/1/2/0 approximate steps for fixed/adaptive/reuse/default policies;
Node recorded 4/1/2/0. Counts differ slightly by runtime; independently verified
endpoint accuracy agrees. The default fixtures
never invoked approximate acceptance.

Public L-BFGS/B explicitly disable approximate acceptance and preserve their
existing search policy. The standalone comparison checks endpoint accuracy and
compares statuses, accepted iterations, and callback counts against the previous
qualified source; see `checks/standalone.json`.

## Accuracy and performance evidence

The same independently constructed fixtures, analytic solutions, scalar Python
`math.fsum` diagnostics, common endpoint criteria, and calibrated SciPy protocol
are used as in the previous report. The fixture files are byte-identical. The
four-policy check now uses `--require-all-converged`; failures cannot silently
become qualified results.

The default objective-call total is still **1,187** for one solve of each of the
16 fixtures. This change improves recovery, with no measured reduction in the
default's total work. Timings below are fresh measurements, not a claim of a
causal speed improvement over the previous run.

| Runtime | Reference | Lower median time | Geometric mean reference/Gale |
|---|---|---:|---:|
| JVM | SciPy SLSQP | 14/16 | 3.05x |
| Node | SciPy SLSQP | 11/16 | 1.79x |
| JVM | SciPy trust-constr | 16/16 | 36.40x |
| Node | SciPy trust-constr | 16/16 | 21.32x |

[Per-case results](comparison.md) retain losses and timing ranges. Nonisolated
host timing varies, particularly near ties; the work counts and convergence
checks are the stronger evidence for this follow-up. These synthetic callbacks
are inexpensive; no expensive-callback speed claim is made.

The runtime protocol remains 50 warmups and 15 timed repeats, sequential
runtimes, complete solve costs, JDK 21.0.12.1, Node 22.18.0/full-optimized
Scala.js, Python 3.14.7, NumPy 2.4.3, SciPy 1.17.1 on Apple M3 Max/macOS 14.3.
BLAS thread environment limits were set to 1, but no BLAS pools were enumerated
for independent confirmation. Source SHA-256 hashes, raw rows, calibration data,
and environment records are retained.

## What remains

- Nearly dependent constraints still show policy sensitivity. The default solves
  the tested rows `(1, +/-0.001)` accurately, but some optional policies still
  terminate unsuccessfully. The safeguard does not solve arbitrary conditioning
  or rank-deficiency problems. Tests continue to require honest diagnostics.
- ALM still uses more user callbacks than SLSQP on these fixtures. Reducing that
  gap substantially would require further algorithmic work and measurement with
  expensive callbacks. The present change does not establish that result.
- First-order KKT checks do not prove a local/global minimum or infeasibility.
  This fixture set is not industrial-scale or sparse-constraint qualification.

## Checks and reproduction

931 JVM and 918 Node core tests, 54 laws tests per runtime, 85 parity tests,
29 Breeze interop tests, and 52 affected optimizer tests under fully optimized
Node pass. Formatting, `compileAll`, `docsCheck`, and `benchCompile` pass with
only the existing unrelated documentation warnings. The new tests exercise
positive/negative offsets from 1e-200 to 1e200, analytic minima under one-ULP
value noise, resolved increases, insufficient/negative curvature, predicted
change limits, ordinary-Wolfe preference, budgets/cancellation, box boundaries,
and ALM trace/accounting integration. Eight verifier mutation probes reject
corrupt evidence. Fresh independent numerical review found no blocking issue.

Using the pinned runtime environment:

```sh
python tools/optim-alm/run_comparison.py --output /tmp/gale-alm-robust --warmups 50 --repeats 15
python tools/optim-alm/check_verifier.py /tmp/gale-alm-robust
```

For each of `fixed`, `adaptive`, `scale`, and `combined`, capture an instrumented
benchmark run (replace both occurrences of `fixed`):

```sh
sbt 'benchmarksJVM/runMain gale.bench.OptimizationALMBench 5 1 docs/verification/optim-alm-robust/fixtures.tsv fixed fixed trace'
```

Combine the four raw logs and verify:

```sh
python tools/optim-alm/check_ablations.py --fixtures docs/verification/optim-alm-robust/fixtures.json --log /tmp/ablations.log --output /tmp/ablations.json --require-all-converged
```

The same main and arguments are used through the fully optimized Scala.js main
initializer for the Node evidence. `checks/qualify.py` records the complete local
sequence. No commit or push was made for this follow-up.
