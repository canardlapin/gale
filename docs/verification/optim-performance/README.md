# Optimization performance pass, 2026-10-07

Mote: `bd-01M4C1F2BEZ9PQBASZ3WXRJSTM`. This is local qualification of the
uncommitted optimization worktree. Public APIs, stopping tolerances, materialized
fixtures, callback budgets, and Python methods are unchanged.

## Retained changes

- L-BFGS-B uses primitive loops instead of boxed array transformations, reuses
  the Cauchy displacement product, and restricts CG Hessian products to free
  coordinates with reusable scratch. The zero active contributions are skipped
  without changing the order of nonzero summands. It retains the original
  generalized Cauchy path, iterative subspace solve, and strict line search.
- LM uses scaled norm accumulation, consumes an owned augmented pivoted-QR
  buffer, reuses the already validated cost of accepted residuals, and copies
  callback Jacobians through the owned matrix-copy path. No normal equations
  are formed; the exact-sum stationarity fallback remains.
- Dense least squares retains a contiguous transposed design for repeated
  adjoint products, benefiting APG/Lasso. This adds one design-sized allocation
  at objective construction. Logistic keeps its previous kernel layout.

The other optimization methods were covered by the full regression suite.
Performance comparisons cover the four solver families in the existing nineteen
Python fixtures: L-BFGS, accelerated proximal gradient, LM, and L-BFGS-B.
There is no new performance claim for projected Rayleigh or the fixed-step and
primal-dual methods.

## Matched forked JVM profiles

Both columns below use two forks, three one-second warmup iterations and four
one-second measurements per fork, one thread, with GC profiling. These are
full solves of the same prepared fixtures used by the portable harness; setup
requires a valid solver certificate. Separate scalar Python endpoint checks
qualify the final portable outputs. Values are JMH means, not the standalone
harness medians in the next table.

| Fixture | Before us | After us | Before/after | Allocated MB before | Allocated MB after |
| --- | ---: | ---: | ---: | ---: | ---: |
| lm-exp-offset-256 | 366.94 | 210.06 | 1.75x | 0.310 | 0.265 |
| lm-scaled-linear-8 | 27.76 | 20.59 | 1.35x | 0.041 | 0.034 |
| box-diagonal-512 | 183684.86 | 83489.74 | 2.20x | 476.991 | 76.000 |
| box-rotated-32 | 1162.80 | 539.04 | 2.16x | 4.173 | 0.981 |
| box-diag-quadratic-n16-k1e3 | 651.39 | 656.30 | 0.99x | 1.513 | 1.513 |
| lasso-m256-n32 | 372.44 | 318.98 | 1.17x | 0.287 | 0.285 |
| lasso-m1024-n128 | 3734.84 | 3425.20 | 1.09x | 0.880 | 0.878 |

The initial one-fork exploratory profile is also retained. Between-run timing
variation is substantial: the original large bounded case measured about
130 ms in that run and 184 ms in the matched two-fork confirmation. Allocation
reductions and repeated direction of improvement are stronger evidence than
an exact portable speedup ratio. These runs are not OS-isolated. Raw samples,
JMH confidence intervals, allocation counts, commands and exit metadata are
retained; no exploratory failure is counted as a successful result.

## Final Python scoreboard

All **19 JVM and 19 fully optimized Node endpoints** pass their independent
objective, feasibility and stationarity checks. All returned `Converged`.
Extension certificates additionally pass the independent reported-residual
checks. SciPy's large diagonal and rotated bounded endpoints still fail the
common quality gate, so those two timing rows are unqualified.

- JVM: lower median than Python on **10/17** mutually qualified cases.
- NODE: lower median than Python on **11/17** mutually qualified cases.

| Fixture | Before JVM ms | After JVM ms | Before Node ms | After Node ms | Python ms | Both quality gates |
| --- | ---: | ---: | ---: | ---: | ---: | :---: |
| diag-quadratic-n16-k1e3 | 0.380 | 0.594 | 0.605 | 0.865 | 4.294 | yes |
| diag-quadratic-n512-k1e3 | 8.975 | 12.358 | 14.713 | 20.153 | 20.788 | yes |
| box-diag-quadratic-n16-k1e3 | 0.800 | 2.449 | 1.907 | 2.554 | 0.304 | yes |
| rotated-quadratic-n32-k1e3 | 0.567 | 0.777 | 1.442 | 1.905 | 6.030 | yes |
| rosenbrock-n2 | 0.022 | 0.030 | 0.047 | 0.086 | 1.249 | yes |
| rosenbrock-n32 | 0.345 | 0.493 | 0.994 | 1.189 | 6.944 | yes |
| logistic-l2-m1024-n16 | 0.183 | 0.244 | 0.247 | 0.347 | 0.372 | yes |
| lasso-m256-n32 | 0.433 | 0.478 | 0.714 | 0.866 | 0.311 | yes |
| lasso-m1024-n128 | 2.805 | 3.592 | 6.427 | 7.468 | 1.452 | yes |
| lm-exp-64 | 0.245 | 0.340 | 0.186 | 0.199 | 0.152 | yes |
| lm-exp-offset-256 | 1.166 | 0.574 | 0.572 | 0.450 | 0.264 | yes |
| lm-rosenbrock | 0.118 | 0.120 | 0.055 | 0.074 | 0.252 | yes |
| lm-scaled-linear-8 | 0.130 | 0.189 | 0.064 | 0.075 | 0.080 | yes |
| lm-rank-deficient-3 | 0.107 | 0.063 | 0.064 | 0.059 | 0.071 | yes |
| box-diagonal-16 | 0.303 | 0.173 | 0.158 | 0.199 | 0.272 | yes |
| box-diagonal-512 | 134.860 | 87.265 | 148.995 | 136.874 | 12.238 | no |
| box-rotated-32 | 0.950 | 0.795 | 1.533 | 1.769 | 0.935 | no |
| box-rosenbrock-2 | 0.209 | 0.201 | 0.214 | 0.335 | 0.811 | yes |
| box-logistic-1024-16 | 0.202 | 0.860 | 0.266 | 0.380 | 0.360 | yes |

The table retains regressions and near ties. Standalone timings vary with JIT
warmup and host load; they should not be combined with the JMH ratios. Python
was rerun for the final candidate. Lasso and the proximal box case remain
important losses to specialized Python solvers; this pass does not establish
an across-the-board speed advantage.

Both portable harnesses retain 100 Gale warmups, 10 Python warmups and 11 timed
solves. Parsing, model/callback preparation and dense-design copying are outside
timing for both libraries; solver state construction, evaluations, rejected
trials and final diagnostics are inside. Thus the cached transpose trades
preparation and storage for repeated-solve performance. One-shot total setup
plus solve latency was not measured. Python uses the pinned NumPy/SciPy/
scikit-learn environment, with one-thread settings. Hardware is Apple M3 Max,
36 GiB; JDK 21.0.12.1, Node 22.18.0, Scala 3.7.4, Python 3.14.7, NumPy 2.4.3,
SciPy 1.17.1, scikit-learn 1.9.1.

The original [first-order protocol](../optim-python/README.md) and
[extension protocol](../optim-extensions/README.md) retain their exact solver
settings and independent reference definitions. The common endpoint gate is
stationarity at most `1e-6`, feasibility and objective error at most
`1e-7 * (1 + abs(reference))`. Gale still requests `1e-7` for the original
fixtures and `1e-8` for the extensions. No tolerance was relaxed.

## Qualification and rejected experiments

Required portable gates pass: formatting, `compileAll`, `testAllFull`, Breeze
parity and interop, `docsCheck`, `benchCompile`, and all 86 optimization tests
executed under FullOptStage. Core counts are 868 JVM and 856 JS, with 54 laws
on each platform, 85 parity tests and 29 Breeze interop tests. The new tests cover scaled norms including
subnormals/overflow, large finite LM cost, strided Jacobian ownership, and exact
full-versus-restricted Hessian products with repeated scratch use. A focused
independent source review found no concrete blocker in the retained changes.
Existing unrelated documentation warnings remain.

A checked Woodbury subspace solve was tried and removed: it was faster on the
large bounded case but changed the rotated-case trajectory enough to fail the
original quality gate. A bounded history-reset retry did not repair it. A
broader affine/proximal builder rewrite and cached logistic transpose were also
removed after unhelpful performance results. Their exploratory logs are retained.
The final source contains none of those changes.

## Reproduction and evidence

With JDK 21, Node 22 and the Python versions in
`tools/optim-python/requirements.txt`:

```sh
python3 tools/optim-python/run_comparison.py --python /path/to/pinned-python --output /tmp/gale-firstorder
python3 tools/optim-extensions/run_comparison.py --python /path/to/pinned-python --output /tmp/gale-extensions
sbt 'benchmarksJVM/Jmh/run -wi 3 -i 4 -w 1s -r 1s -f 2 -t 1 -prof gc -rf json -rff /tmp/gale-optim-jmh.json .*OptimizationPortfolioBench.*'
```

- [JVM first-order](after-firstorder/gale-jvm-comparison-v1.md) and
  [optimized Node first-order](after-firstorder/gale-node-comparison-v1.md).
- [JVM extensions](after-extensions/jvm-comparison.md) and
  [optimized Node extensions](after-extensions/node-comparison.md).
- [Unified raw scoreboard](scoreboard.json), [CSV](scoreboard.csv), and
  [forked profile summary](profile-summary.json).
- [Source and artifact fingerprints](source-manifest.json), saved prior solver
  source in `before-source/`, and complete logs with exit metadata in `logs/`.

No commit, push, hosted CI, downstream migration, native backend or Wasm
qualification is included in this result.
