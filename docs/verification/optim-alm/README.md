# Augmented Lagrangian: accuracy and SciPy comparison

Local qualification on 2026-10-07, Apple M3 Max, 36 GiB RAM, macOS 14.3 arm64.
The measured implementation provides smooth mixed nonlinear equality/inequality
constraints and optional box bounds on JVM and Scala.js. It uses safeguarded
Powell–Hestenes–Rockafellar multiplier updates and existing L-BFGS/B inner solves.

## Measured result

All **960 timed endpoints** passed the independent common accuracy gate:
240 each for Gale JVM, Gale Node, SciPy SLSQP, and SciPy trust-constr.
There are eight problem families, two starts per family, and 15 timed repeats
per case. Dimensions range from 2 to 128 parameters.

| Gale runtime | Reference | Cases with lower median time | Geometric mean of reference/Gale median time |
|---|---|---:|---:|
| JVM | calibrated SLSQP | 11/16 | 2.15x |
| Node | calibrated SLSQP | 10/16 | 1.16x |
| JVM | calibrated trust-constr | 16/16 | 24.91x |
| Node | calibrated trust-constr | 16/16 | 13.45x |

These are ratios across these synthetic cases, not a general speed guarantee or
a statistical significance claim. [Per-case timings](comparison.md) show the
losses as well as the wins. JVM remains slower than SLSQP on the tiny mixed
linear/bound problem, one 32-parameter dense equality case, and both 64-parameter
ball cases (one almost tied). Node also loses both dense 32-parameter cases.

**Callback cost matters.** Gale JVM uses approximately 2.4–44.5 times as many
objective callbacks as calibrated SLSQP here. Its lower per-iteration overhead
can outweigh this for these cheap objectives; expensive callbacks can reverse
the timing result. Counts and timing ranges are retained in
[comparison.json](comparison.json), including all individual observations.

## Independent accuracy

The fixtures have independent known solutions:

- SPD quadratics with linear equalities: the construction satisfies a linear
  KKT system with a known primal and multiplier. Both diagonal and dense
  Hessians are covered, with condition numbers 10 and 100.
- Mixed active/inactive constraints and a bound: the solution is `(1, 0, 0)`,
  objective 2.5, equality multiplier 2, and inequality multipliers `(2, 0)`.
- Euclidean ball projection: the solution is the radial projection of the
  center, with an analytic positive multiplier.
- Parabolic equalities: each pair satisfies `y=x^2`; the unique positive root
  of `2*x^3-x-2=0` gives the global solution of the reduced quartic. Bisection
  supplies this oracle independently of every tested optimizer.

The verifier recomputes objectives, constraints, Jacobians and Lagrangian
gradients using scalar Python `math.fsum`, outside timing. It uses the analytic
reference multipliers to check every solver's returned point under the same
conditions. It additionally checks Gale's reported diagnostics with Gale's own
returned multipliers and rejects false convergence.

Common acceptance requires raw feasibility <=1e-8, unit box-projected
Lagrangian-gradient infinity norm <=1e-6, complementarity <=1e-8, parameter
error <=1e-5, objective error <=1e-7 times max(1, absolute reference objective),
exact box feasibility, and nonnegative inequality multipliers. Gale itself
requires stationarity <=1e-7 and both raw/scaled feasibility <=1e-8.

| Runtime / solver | Accurate timed endpoints | Largest raw violation | Largest common stationarity | Largest complementarity |
|---|---:|---:|---:|---:|
| Gale JVM | 240/240 | 9.03e-9 | 9.41e-8 | 9.10e-9 |
| Gale Node | 240/240 | 9.03e-9 | 9.41e-8 | 9.10e-9 |
| SLSQP | 240/240 | 7.68e-10 | 9.35e-7 | 7.59e-10 |
| trust-constr | 240/240 | 5.24e-10 | 9.16e-7 | 1.82e-9 |

These tests establish accuracy for the supplied fixtures. They do not establish
global convergence for general nonlinear constraints or robustness on a broad
industrial optimization collection.

## Timing protocol and SciPy calibration

JVM uses Temurin JDK 21.0.12.1+1; Node uses 22.18.0 and fully optimized
Scala.js; Scala is 3.7.4. Python is 3.14.7, NumPy 2.4.3, SciPy 1.17.1, and
threadpoolctl 3.6.0. `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and
`VECLIB_MAXIMUM_THREADS` are set to 1. `threadpoolctl` enumerated no BLAS pools,
so the actual BLAS thread count was not independently confirmed.

Data, model callbacks, initial points, and bounds are prepared outside timing.
Each timing includes a complete solve, solver state/cache construction, user
callbacks, and solver-produced final diagnostics. Gale performs 50 untimed
solves of **every** family before measuring any family, allowing the shared
JIT paths to warm before early cases. SciPy performs 50 warmup solves at its
selected setting. Runtimes execute sequentially; the host is not OS-isolated.

All solvers receive analytic first derivatives, without supplied Hessians.
SciPy uses nonlinear constraint callback interfaces even for linear fixtures,
matching the general callback information supplied to Gale. Its constraint
values and Jacobian share a one-point cache. The stopping parameters are not
assumed equivalent across algorithms.

Before final timing, each SciPy method receives a per-case calibration sweep:

- SLSQP `ftol`: 1e-6, 1e-8, 1e-10, 1e-12, 1e-14.
- trust-constr `gtol`: 1e-6 through 1e-16 in factors of 100, paired with initial
  barrier parameter/tolerance 0.1 or 1e-9; `xtol=1e-14`, `barrier_tol=1e-12`.

Each candidate gets one calibration warmup and three timed calibration solves.
The candidate with the lowest calibration median **among those passing all
common accuracy checks** is selected before separate warmups and final timings.
Gale uses the same default numerical configuration on all cases, with a 50,000
callback cap. SciPy has a 1,000-iteration cap and the same externally counted
callback budget. The calibration favors the reference by allowing per-case
settings. All candidates, selected settings, solver status flags, and measured
work are retained in [python.json](python.json). A solver's success flag alone
does not establish acceptance.

## Tests, review and provenance

The 15 portable augmented Lagrangian tests cover analytic endpoints and
multipliers, a one-subproblem closed form, multiplier safeguarding, active and
inactive inequalities, fixed bounds, scaling, warm starts, infeasible stationary
subproblems, global work limits, cancellation, callback errors, nonfinite data,
and cancellation among large gradient terms. The last case verifies that
`1 + 1e16 - 1e16` cannot erase a stationarity residual of 1.

An independent numerical review identified the summation hazard and a mismatch
in reference precision; both were corrected and rereviewed without a remaining
concrete blocker. Review was read-only; the primary agent executed the checks.

The archived [gate log](gates.log) and [exit metadata](gates.meta.json) record:

- Formatting and `compileAll` pass.
- `testAllFull`: 914 JVM core tests, 901 JS core tests, 54 laws per runtime pass.
- 85 Breeze parity tests and 29 Breeze interop tests pass.
- `docsCheck` and `benchCompile` pass.
- The 15 new tests also pass under `FullOptStage`.

Existing unrelated Scaladoc link/classpath and mdoc link warnings remain; no
new guide warning was reported. [Verifier mutation probes](verifier-checks.log)
reject a NaN point, incorrect objective, incorrect stationarity diagnostic,
over-budget work, an unknown status, and duplicate evidence.

[source-manifest.json](source-manifest.json) binds the measured optimizer,
harness and test source to SHA-256 hashes plus the base Git commit and runtime
versions. The source was uncommitted at measurement time; the base commit
alone does not identify this implementation. [artifact-manifest.json](artifact-manifest.json)
binds the saved report, fixtures, raw outputs, and logs. This is local numerical
and performance qualification, not hosted CI or a published release.

## Reproduce

Install `tools/optim-alm/requirements.txt` in a Python environment, select JDK 21
and Node 22, then run from the repository root:

```sh
python3 tools/optim-alm/run_comparison.py \
  --python /path/to/pinned-python \
  --output /tmp/gale-alm-comparison --warmups 50 --repeats 15
python3 tools/optim-alm/check_verifier.py /tmp/gale-alm-comparison --repeats 15
```

Use a new output directory to preserve this snapshot. The harness retains
inputs and raw output and distinguishes valid endpoints from accuracy failures.
No reported timing ratio should include a failed accuracy case.
