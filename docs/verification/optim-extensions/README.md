# Nonlinear least squares and L-BFGS-B qualification

This addition implements dense unconstrained Levenberg–Marquardt and L-BFGS-B
with an iterative subspace solve. The [executable guide](../../user/guides/nonlinear-and-bounded-optimization.md)
defines the public contracts. Mote issue: `bd-01M4BZE6050SETTP6ANWH8ANM6`.

LM uses residual/Jacobian callbacks, explicit central differences when selected,
fixed nonnegative weights, column scaling, augmented pivoted QR, gain-ratio
damping, and a raw final `||Jᵀr||∞` stationarity certificate. L-BFGS-B uses a
generalized Cauchy point, limited-memory Hessian factors, a CG subspace solve,
and feasible Wolfe search with an Armijo endpoint exception. Neither method
uses a small step or small objective change as sufficient convergence evidence.

## Independent review and adversaries

A focused independent numerical review covered algorithm design and source.
It identified and prompted permanent regressions for cancellation in LM's
gradient accumulation, cancellation when a large Cauchy-path coordinate freezes,
and stale proposal counts on evaluation-budget interruption. Exact objective
ties in Wolfe search receive derivative-based treatment; points already
satisfying both Wolfe conditions are accepted before bracket-order checks.
No Armijo slack or relaxed endpoint tolerance was introduced.

Tests include analytic active/free/fixed bounds, linear boundary optima,
translated variables, independently constructed dense BFGS matrices after
history eviction, analytic piecewise Cauchy paths, weighted and nonzero-residual
fits, rank deficiency, parameter scaling, finite-difference work counts,
callback failures, rejected trials, cancellation, and unrepresentable LM steps.

## Materialized comparison protocol

Ten fixtures in `fixtures.json` and losslessly round-trippable `fixtures.tsv`
are shared by Scala and Python. They cover exponential fits with and without
noise/offset, Rosenbrock residuals, scaled linear least squares, rank deficiency,
16/512-variable bounded diagonal quadratics, a rotated bounded quadratic,
bounded Rosenbrock, and bounded logistic regression.

## Final local results

All 10 JVM and all 10 optimized Node endpoints pass both independent quality
and certificate checks. SciPy meets the common gate on 8 of 10 timed cases;
its 512-variable diagonal and 32-variable rotated box cases terminate with
projected residuals about `2.52e-6` and `1.04e-6`, above the `1e-6` gate.
Gale reaches its requested `1e-8` target on both. Their timing ratios remain
unqualified; a faster insufficiently accurate solve is not a comparable win.

| Fixture | JVM ms | Node ms | SciPy ms | Both endpoints qualified |
| --- | ---: | ---: | ---: | :---: |
| Exponential, 64 observations | 0.162 | 0.181 | 0.110 | yes |
| Exponential with offset/noise, 256 | 0.349 | 0.594 | 0.181 | yes |
| Rosenbrock residuals | 0.076 | 0.058 | 0.185 | yes |
| Scaled linear least squares, 8 parameters | 0.135 | 0.063 | 0.059 | yes |
| Rank-deficient least squares, 3 parameters | 0.110 | 0.068 | 0.052 | yes |
| Box diagonal quadratic, 16 | 0.130 | 0.148 | 0.204 | yes |
| Box diagonal quadratic, 512 | 136.256 | 153.688 | 9.737 | no |
| Box rotated quadratic, 32 | 0.896 | 1.669 | 0.736 | no |
| Box Rosenbrock, 2 | 0.202 | 0.238 | 0.600 | yes |
| Box logistic, 1024 by 16 | 0.191 | 0.281 | 0.278 | yes |

The final standalone medians favor Gale on 4/8 qualified JVM cases and 3/8 Node
cases. Node logistic is a near tie. JVM small-rank/logistic timings varied
substantially across the first run and reproduction; these measurements are
not OS-isolated and do not establish stable small-margin speed rankings.
SciPy is faster on most of these dense LM fixtures. The larger box problem is
expensive in this initial implementation: 300 iterations and repeated
`O(history²*n)` Hessian-factor reconstruction. These results establish useful
coverage and truthful convergence, not a general speed advantage.

Detailed [JVM](jvm-comparison.md) and [Node](node-comparison.md) reports retain
accuracy, status and ratios. TSV files retain coordinates, all 11 time samples,
iterations and work counts; `python.json` retains SciPy statuses, references,
work, samples and environment. The source and input identity is in
`source-manifest.json`. Hardware is Apple M3 Max/36 GiB; runtimes are Adoptium
21.0.12.1, Node 22.18.0, Scala 3.7.4, SciPy 1.17.1 and NumPy 2.4.3.

Final gates pass: formatting, `compileAll`, `testAllFull` (864 core JVM,
852 core JS, 54 laws per platform), 85 Breeze parity and 29 interop tests,
`docsCheck`, `benchCompile`, and all 82 optimization tests under FullOptStage.
The new solvers contribute 21 shared adversarial tests. Existing unrelated
Scaladoc/link warnings remain. The full reproduction script ran successfully.
Four verifier mutation probes reject NaN/mismatched residuals, an unknown
status, and a NaN objective. Raw logs and actual exit metadata are retained
in `logs/`; exploratory failures are retained alongside final successful runs.

## Timing and acceptance details

SciPy uses `least_squares(method="lm")` with analytic Jacobians and
`minimize(method="L-BFGS-B", jac=True)` with bounds. Callback/data preparation
is outside all timed solves. Timings include solver setup, callbacks, rejected
trials, final diagnostics, and result creation. Gale warms up 100 solves and
Python 10, then each records 11 timed solves. Optimized Node uses FullOptStage.
Python uses vectorized NumPy callbacks; the independent scalar verifier is
outside timing for both implementations. Thread limits are one.

Gale requests absolute stationarity `1e-8`. SciPy LM uses `gtol=1e-8` and
`ftol=xtol=1e-12`; its `gtol` is a cosine-based test, not Gale's raw gradient
criterion. SciPy L-BFGS-B uses `gtol=1e-8`, `ftol=0`, memory 10, and 40 line-search
trials. A common independent endpoint gate requires stationarity at most `1e-6`,
exact feasibility, and objective error at most `1e-7*(1+abs(reference))`.
Gale certificates must additionally be finite/nonnegative and agree with the
independent metric to `1e-10*max(1,metric)`. A `Converged` result must report at
most `1e-8`, with independent recomputation at most `1.01e-8`. Objective agreement
uses `1e-10*(1+abs(value))`. All raw statuses and callback counts are retained.

References use strict external solves, analytic clipped centers for diagonal
quadratics, and an independently solved/checked active-set linear KKT system
for the rotated quadratic. Failed SciPy endpoints remain visible. Speed ratios
and win counts are qualified only when both endpoints meet the common gate.

Reproduce with JDK 21, Node 22, sbt and a Python environment containing the
versions in `tools/optim-python/requirements.txt`:

```sh
python3 tools/optim-extensions/run_comparison.py \
  --python /path/to/pinned-python --output /tmp/optim-extensions-rerun
```

`compare.py fixtures` regenerates the inputs using that pinned Python
environment; reproduction normally reads the materialized inputs unchanged.
The result is local qualification of uncommitted source. Publication, hosted
CI, downstream adapter migration, native/Wasm qualification, robust losses,
and sparse/operator Jacobians are separate work.

Primary methods: [Byrd et al., L-BFGS-B](https://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf),
[Madsen et al., nonlinear least squares](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf),
[SciPy least squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html).
