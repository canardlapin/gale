# Four-case performance follow-up, 2026-10-07

This preserves the qualification before the subsequent
[LM/L-BFGS-B extension](../../optim-extensions/README.md). Its source manifest
identifies that earlier implementation; use the extension report for the newer
candidate's tests and measurements.

The user asked whether the four workloads still slower than Python could be
improved without expanding the solver portfolio. This follow-up retains the
existing accuracy criteria, fixtures, timing boundaries, warmup policy and pinned
Python comparator. It adds reusable `DenseObjectives.leastSquares` and
`DenseObjectives.logistic`, reuses accepted proximal evaluations, removes trial
objective work from fixed-point checks, and adds a normalized gradient restart
alongside objective restart. No additional solver family was added.

All nine endpoints on JVM and all nine on fully optimized Node pass the
independent Python objective, feasibility and stationarity/KKT checks. Current
JVM results favor Gale on six of nine workloads; the box case and both Lasso
cases remain slower than specialized SciPy/scikit-learn methods.

| Workload | Original JVM ms | Current JVM ms | Current Node ms | Current Python ms |
| --- | ---: | ---: | ---: | ---: |
| Box quadratic, 16 variables | 0.649 | 1.412 | 1.668 | 0.233 |
| Logistic, 1024 by 16 | 0.313 | 0.176 | 0.249 | 0.298 |
| Lasso, 256 by 32 | 0.436 | 0.628 | 0.679 | 0.242 |
| Lasso, 1024 by 128 | 8.337 | 2.842 | 6.412 | 1.032 |

These are medians of 11 solves after 100 Gale/10 Python warmups on the same
Apple M3 Max, JDK 21.0.12.1 and Node 22.18.0 as v1. The original column is the
preserved earlier run, not an interleaved before/after experiment. The larger
Lasso change is substantial in both runtimes (Node was 12.009 ms in v1);
logistic's JVM improvement is also visible in the earlier exploratory follow-up
run (0.172 ms). Exact speedup ratios are not portable guarantees.

The changes do not improve every case. The JVM small Lasso and box regressions
are retained explicitly. Dense kernels and restart decisions have size- and
problem-dependent costs: small Lasso uses the same 22 iterations with fewer
callbacks (126 to 107), while the box case takes 426 rather than 358 iterations
(2133 rather than 2162 callbacks). Reduced callback work alone does not imply
lower wall time. We did not tune a special solver or fixture-specific stopping
rule to erase these losses. The optimized Node small-Lasso measurement improves
from 0.792 to 0.679 ms; its smaller margin should be interpreted cautiously.

The separate forked JMH small-Lasso profile is faster than the earlier profile
(268.6 versus 389.5 microseconds/solve), but allocates more (285.8 versus
135.2 kB/solve). Its warm steady-state timings differ from the standalone harness;
do not combine the ratios or treat the standalone regression as universal.
The current accelerated/fixed-step proximal profile measures 2.565 versus
21.229 ms/solve and 5.510 versus 41.974 MB/solve, with the harness independently
checking KKT conditions. Both profiles use two forks, two one-second warmup
iterations and three one-second measurement iterations per fork, with GC
profiling. JSON reports are `jmh-lasso-v2.json` and `jmh-proximal-v2.json`.

Current-source checks passed on JDK 21/Node 22: formatting, `compileAll`,
`testAllFull` (843 core JVM, 831 core JS, 54 laws on each platform), 85 parity
and 29 Breeze tests, `docsCheck`, and `benchCompile`. All 61 optimization tests
also executed under fully optimized Scala.js. The guide edit was followed by
another successful `docs/mdoc`. Existing unrelated documentation/link warnings
remain. Complete logs are in `logs/`; `resume-gates.log` ends with a JMH result
path error after the portable checks pass. Both profiles subsequently pass
using absolute paths in `resume-profile.log`; no failed profile is counted as
successful evidence.

Detailed current results retain coordinates and all timing samples:

- [JVM accuracy/timing report](gale-jvm-comparison-v1.md)
- [Optimized Node accuracy/timing report](gale-node-comparison-v1.md)
- [Python comparator endpoints and environment](python-baseline-v1.json)
- [Current source and artifact manifest](source-manifest.json)

Reproduce using the original runner with `--output /tmp/optim-followup` and the
pinned Python environment. `previous-source-manifest.json` preserves the v1
identity; the manifest in this directory identifies the current uncommitted
source. This is local evidence, without publication, hosted CI, downstream
adapter migration, native/Wasm qualification, or a general speed claim.
