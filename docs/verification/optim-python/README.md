# Optimization qualification, 2026-10-07

This report and its v1 manifest retain the original local qualification before
the subsequent performance changes. The current implementation, including dense
objective helpers and proximal evaluation/restart changes, is described and
fingerprinted in the [performance follow-up](followup/README.md). Use that report
for current-source checks; v1 timing and source hashes apply to the earlier state.

Gale now provides a bounded optimization portfolio: L-BFGS for smooth
unconstrained objectives, accelerated proximal gradient with backtracking for
convex composite objectives, the four existing fixed-step/operator solvers,
standard proximal terms and convex sets, work controls and explicit warm starts.
It remains part of portable `gale-core`. No modeling framework, AD, stochastic
optimizer, general nonlinear constraints, L-BFGS-B, or global-optimum guarantee
has been added.

The Mote epic is `bd-01M4BQS2FZZVD18DZTWYK52VBP` (OPTIM-01–19). Implementation
and local qualification are complete; OPTIM-19 remains open for downstream
adapter migration and exact-candidate hosted CI before publication. No commit,
release, native/Wasm qualification, or full ScalaFIM build is claimed.

## Evidence and reproducibility

The base revision is `cf1cbe0dcfbc5208e0599de7a2adaef61b890508`. The tested
implementation is an uncommitted worktree change; `source-manifest.json`
identifies it by per-file SHA-256 hashes rather than pretending it has a new
commit SHA. The manifest also fingerprints the materialized inputs and reports.

- [JVM comparison](gale-jvm-comparison-v1.md) and [optimized Node comparison](gale-node-comparison-v1.md)
  independently recompute all nine endpoints in Python.
- `gale-jvm-v1.tsv` and `gale-node-v1.tsv` retain coordinates, status, iterations,
  callback counts and every measured time sample.
- `python-baseline-v1.json` retains the standard comparator endpoints, work,
  times and Python/NumPy/SciPy/scikit-learn/BLAS environment.
- `jmh-v1.json` measures the old two-pass affine expression versus the fused
  expression. This is an isolated kernel comparison, not a historical whole-solver speedup.
- `jmh-solves-v1.json` measures full solves and allocations for diagonal
  quadratic, Rosenbrock and Lasso cases. `jmh-proximal-v1.json` compares
  accelerated and fixed-step proximal methods at independently checked accuracy.
- Final validation logs and exit metadata are retained in `logs/`.

On JDK 21 with Node 22 and sbt available:

```sh
python3 -m venv /tmp/gale-optim-python
/tmp/gale-optim-python/bin/pip install -r tools/optim-python/requirements.txt
python3 tools/optim-python/run_comparison.py --python /tmp/gale-optim-python/bin/python
```

The runner executes JVM and fully optimized Node solves, then the Python
baseline, then the independent verifier. It performs no publication or global
configuration. Use `--output /tmp/optim-rerun` to preserve checked-in results.
The default fixture/reference artifacts are materialized in this directory;
`optim_protocol.py all` regenerates them using the pinned Python packages.

Allocation/full-solve profiling is separately reproducible:

```sh
sbt 'benchmarksJVM/Jmh/run -i 3 -wi 2 -r 1s -w 1s -f 2 -t 1 -prof gc -foe true gale.bench.Optimization.*Bench.*'
sbt 'benchmarksJVM/runMain gale.bench.ProximalComparison' \
  'set benchmarksJS / Compile / mainClass := Some("gale.bench.ProximalComparison")' \
  'set benchmarksJS / scalaJSStage := FullOptStage' 'benchmarksJS/run'
```

JMH uses two forks, two one-second warmup iterations and three one-second
measurement iterations per fork. Its warm steady-state means and allocation
measurements are separate from the standalone medians; do not combine their
ratios. Every full-solve JMH setup checks the fixture endpoint. The proximal
comparison additionally reconstructs analytic KKT conditions.

## Independent numerical protocol

`fixtures-v1.json` contains exact starts, formulae, budgets, bounds and
materialized float64 arrays. Its `.decimal.tsv` companion uses 17 significant
digits and flattened row-major arrays, which both JVM and Scala.js parse.
Hex-float companions are retained as an additional lossless representation.
Gale never regenerates the arrays through a different random-number generator.

The nine workloads are condition-1000 diagonal quadratics (16 and 512 variables),
a box-constrained diagonal quadratic, a rotated 32-variable quadratic,
Rosenbrock in 2 and 32 variables, L2 logistic regression (1024 by 16), and
Lasso (256 by 32 and 1024 by 128). The Lasso objective is
`0.5/m * ||Xx-y||² + lambda * ||x||1`, matching scikit-learn with no intercept.

References use analytic minima for quadratic/Rosenbrock fixtures and strict
external solves for logistic/Lasso. Strict reference failures remain visible:
four SciPy quadratic/Rosenbrock solves did not meet the requested 1e-9 reference
metric; their analytic minima provide the independent expected objective.

`compare_gale.py` uses scalar-loop formulas independent of the Scala callbacks
and the NumPy objective implementation. It checks finite data, dimensions,
feasibility, agreement with the reported objective, objective error at most
`1e-7*(1+abs(reference))`, and gradient/projected-unit-gradient/Lasso KKT infinity
norm at most `1e-6`. It rejects missing/duplicate fixtures and invalid timings.
All 18 final Gale endpoints pass. It never trusts a solver's convergence flag.
The core tests also exercise primal-dual analytic solutions, nonconvex local
projection, Rayleigh angle/geometry checks, malformed callbacks, empty/invalid
inputs, budgets, cancellation, ownership and floating-point adversaries.

## Timing boundaries and interpretation

Hardware: Apple M3 Max, arm64, 36 GiB RAM. Runtimes: Adoptium JDK 21.0.12.1,
Node 22.18.0, Scala 3.7.4; Python package versions and threadpool metadata are
recorded in `python-baseline-v1.json`. These are local single-machine
measurements, without OS isolation or a guarantee that other user processes
were idle. Near-equal timings are not evidence of a meaningful advantage.

The standalone Gale runner uses 100 untimed warmup solves and 11 timed solves;
Python uses 10 warmups and 11 timed solves. JVM warmup is deliberately longer.
Parsing, random generation, array conversion and objective construction are
outside the timing boundary. Solver state setup, callbacks, rejected trials,
final solver diagnostics and Python's endpoint checks remain inside. Python
array conversion was initially inside timing; that exploratory baseline was
replaced with the fair prepared-array baseline retained here.

SciPy uses `minimize(method="L-BFGS-B")`, without bounds for the six smooth
unconstrained cases and with bounds for the box case. Gale uses L-BFGS or its
convex proximal method. Lasso compares Gale's general accelerated proximal
method to scikit-learn's specialized coordinate descent. The algorithms are
appropriate standard implementations, not identical iteration sequences.
Gale uses an absolute 1e-7 stopping target; the timed SciPy runs request
`gtol=1e-7, ftol=0`, and timed scikit-learn runs use `tol=1e-8`. All attained
endpoints are separately checked; option values alone do not establish parity.
No standard Python comparator was replaced with slow Python scalar loops for timing.

The final retained JVM run wins five of nine cases. Logistic regression is
a near tie (0.313 ms Gale versus 0.296 ms Python), with losses for the box and
Lasso cases. Earlier repetitions favored Gale on logistic regression; that
variation is why no meaningful speed advantage is claimed for this case. Specialized bounded/coordinate solvers and
native array kernels remain stronger on those workloads. The Node results also
include near ties and losses. This does not motivate adding an entire solver
catalog merely to win a benchmark.

On the condition-1000, 64-variable quadratic-plus-L1 comparison, JMH reports
about 3.07 ms accelerated versus 20.84 ms fixed-step (6.8x), and 6.67 MB versus
41.97 MB allocated per solve (84% less). Both reach independently recomputed
KKT residual below 1e-7. The fused affine kernel reduces allocations from
328/8264/65608 to 184/4152/32824 bytes per operation at sizes 16/512/4096;
these are approximate JMH measurements. Primal-dual instrumented tests verify
one initial forward plus two per iteration, reduced from three per iteration;
final-point checks are retained.

## Integration limits and remaining release gate

Final local checks pass: formatting, `compileAll`, `testAllFull`, Breeze parity,
interop, `docsCheck`, benchmark compilation and explicit optimized-JS execution
of all 53 optimization tests. The full test invocations report 835 JVM core,
823 JS core, 54 laws on each platform, 85 parity and 29 Breeze tests. Existing
Scaladoc/link warnings outside this optimization work remain.

An isolated build compiled current Multivar sources against this Gale worktree
and passed 23 selected optimization/compiler/canonical tests. This is source
integration evidence, not published-artifact compatibility. Seven exhaustive
status matches in four Multivar files still need explicit handling for
`NumericalStagnation`, `EvaluationLimit`, `Cancelled` and `LineSearchFailed`:

- `family/glrm/FittedLatentEncoder.scala` (two matches)
- `solver/CompositePenaltyCompiler.scala` (one)
- `solver/VariationalFunctionalFrames.scala` (one)
- `solver/VariationalSolverCompiler.scala` (three)

These statuses must retain a nonconverged/unresolved interpretation and must
never acquire an optimality claim. Current happy-path consumer tests do not
exercise every new stop. ScalaFIM's coupled-column and nonconvex-projection
contracts were preserved and tested in Gale; a full ScalaFIM build was not run.
Both downstream migration and hosted exact-candidate CI stay in OPTIM-19.

Method references: [FISTA](https://www.tau.ac.il/~becka/FISTA.pdf),
[Condat splitting](https://lcondat.github.io/publis/Condat-optim-JOTA-2013.pdf),
[L-BFGS](https://users.iems.northwestern.edu/~nocedal/lbfgs.html),
[SciPy L-BFGS-B](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html),
[scikit-learn Lasso](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html).
