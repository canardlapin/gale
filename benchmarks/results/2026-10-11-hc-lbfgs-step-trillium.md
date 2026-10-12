# hc-lbfgs-step paired run on trillium (Zen 5, lanes A and B)

`d11d15f` (breeze/hc-lbfgs-step: scaled first step, fixed-budget decrease floor, shared tolerance config, L-BFGS preflight) against base `61c82cb`. Slurm job 2539904, one exclusive compute node (`tri0564`, 2 × AMD EPYC 9655, Zen 5, AVX-512).

## Protocol

- **Bundles:** base `61c82cb` and new `d11d15f`, both compiled locally from one clean bench worktree (branch `bench/lbs-d11d15f`) with `benchmarksJVM/clean` run before each compile; sbt 1.11.7, JDK 25.0.1, Scala 3.7.4. No sbt ran on the cluster.
- **Harness:** `docs/verification/w21-simd-spike/x86-avx512/lb.sh` (unchanged from the hc-lbfgs receipt), which is the reviewed `run.sh` lane paths restricted to `gale.bench.LbfgsBreezeJmh.*`. Lane A uses `-jvmArgs` (no Vector module) and `-p backend=pure`. Lane B uses the inherited `--add-modules` and `-p backend=pure,vector`. Both use `-Xms4g -Xmx4g`, a bundle-bound `gale.bench.sourceSha`, `-foe true`, a fresh netlib sidecar per run, and the annotation defaults: 2 forks, 5 × 500 ms warmup, 5 × 500 ms measurement, average time in µs/op. Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, `numactl --physcpubind=0-7 --membind=0`.
- **Order:** interleaved: laneA-base-r1, laneA-new-r1, laneB-base-r1, laneB-new-r1, laneA-base-r2, laneA-new-r2, laneB-base-r2, laneB-new-r2.
- **Load at start:** A-base-r1: 0.87; A-new-r1: 1.13; B-base-r1: 1.20; B-new-r1: 1.11; A-base-r2: 1.35; A-new-r2: 1.25; B-base-r2: 1.38; B-new-r2: 1.19 (1-minute).
- **Netlib** (blas, lapack, vectorModule): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False), ('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- **Validation:** the new bundle (`d11d15f`) has an `LbfgsBreezeJmh` preflight. Every new run wrote validation records bound to `d11d15f`, and all passed: 16 per lane-A run and 20 per lane-B run. The new runs' scoreboards (`--strict --validation <run>/netlib.jsonl.validation.jsonl`) exited 0, validating 8 (lane A) and 10 (lane B) optimizer result rows with 0 unpaired. The base (`61c82cb`) has no `LbfgsBreezeJmh` preflight, so the base runs were scored `--strict` without `--validation` (legacy receipts) and also exited 0. The fixed-budget harness guard (exactly 20 iterations on both sides) passed in every fork.

## gale / Breeze (speed ratio, >1 means gale is faster; verdict from 99.9% CIs)

| lane | problem | budget | gale backend | base r1 | base r2 | new r1 | new r2 |
|---|---|---|---|---|---|---|---|
| A | Logistic | fixed | pure | 1.14x ahead | 1.14x ahead | 1.20x ahead | 1.20x ahead |
| A | Logistic | tolerance | pure | 0.93x behind | 0.93x behind | 1.32x ahead | 1.30x ahead |
| A | Rosenbrock | fixed | - | 1.74x ahead | 1.73x ahead | 1.90x ahead | 1.90x ahead |
| A | Rosenbrock | tolerance | - | 1.52x ahead | 1.51x ahead | 1.74x ahead | 1.50x ahead |
| B | Logistic | fixed | pure | 0.69x behind | 0.68x behind | 0.71x behind | 0.70x behind |
| B | Logistic | fixed | vector | 0.96x behind | 0.94x behind | 0.99x tie | 0.99x tie |
| B | Logistic | tolerance | pure | 0.55x behind | 0.55x behind | 0.79x behind | 0.79x behind |
| B | Logistic | tolerance | vector | 0.78x behind | 0.77x behind | 1.13x ahead | 1.15x ahead |
| B | Rosenbrock | fixed | - | 1.73x ahead | 1.75x ahead | 1.89x ahead | 1.89x ahead |
| B | Rosenbrock | tolerance | - | 1.56x ahead | 1.54x ahead | 1.50x ahead | 1.61x ahead |

## Times (µs/op, mean ± 99.9% CI half-width)

| lane | problem | budget | series | base r1 | new r1 | base r2 | new r2 | gale new/base speed r1 / r2 |
|---|---|---|---|---:|---:|---:|---:|---|
| A | Logistic | fixed | gale pure | 810.8 ± 9.0 | 761.5 ± 7.9 | 814.6 ± 9.8 | 763.0 ± 8.2 | 1.06x* / 1.07x* improved |
| A | Logistic | fixed | Breeze | 927.1 ± 16.0 | 912.7 ± 19.6 | 927.4 ± 16.6 | 911.9 ± 20.6 | |
| A | Logistic | tolerance | gale pure | 1,895.6 ± 19.1 | 1,302.8 ± 14.2 | 1,888.6 ± 19.1 | 1,305.2 ± 16.0 | 1.46x* / 1.45x* improved |
| A | Logistic | tolerance | Breeze | 1,763.9 ± 23.9 | 1,719.9 ± 68.9 | 1,760.5 ± 37.1 | 1,691.3 ± 43.9 | |
| A | Rosenbrock | fixed | gale - | 44.9 ± 0.7 | 41.4 ± 0.2 | 44.7 ± 0.3 | 41.6 ± 0.4 | 1.08x* / 1.08x* improved |
| A | Rosenbrock | fixed | Breeze | 78.1 ± 0.3 | 78.8 ± 0.3 | 77.5 ± 0.7 | 78.8 ± 1.2 | |
| A | Rosenbrock | tolerance | gale - | 1,226.8 ± 7.0 | 1,245.7 ± 9.6 | 1,227.6 ± 2.5 | 1,245.9 ± 8.3 | 0.98x* / 0.99x* regression |
| A | Rosenbrock | tolerance | Breeze | 1,865.4 ± 11.1 | 2,169.7 ± 10.1 | 1,853.6 ± 5.8 | 1,864.2 ± 10.7 | |
| B | Logistic | fixed | gale pure | 813.3 ± 10.3 | 763.6 ± 8.3 | 813.5 ± 9.7 | 763.4 ± 8.8 | 1.07x* / 1.07x* improved |
| B | Logistic | fixed | Breeze | 558.1 ± 3.3 | 539.1 ± 3.5 | 550.4 ± 11.3 | 536.4 ± 3.9 | |
| B | Logistic | fixed | gale vector | 584.3 ± 7.1 | 543.3 ± 22.5 | 583.9 ± 12.4 | 542.3 ± 15.1 | 1.08x* / 1.08x* improved |
| B | Logistic | fixed | Breeze | 558.1 ± 3.3 | 539.1 ± 3.5 | 550.4 ± 11.3 | 536.4 ± 3.9 | |
| B | Logistic | tolerance | gale pure | 1,903.0 ± 18.0 | 1,313.0 ± 12.0 | 1,892.3 ± 20.1 | 1,307.9 ± 14.9 | 1.45x* / 1.45x* improved |
| B | Logistic | tolerance | Breeze | 1,044.2 ± 16.9 | 1,031.8 ± 9.4 | 1,048.6 ± 5.8 | 1,039.3 ± 4.7 | |
| B | Logistic | tolerance | gale vector | 1,345.2 ± 11.1 | 911.4 ± 13.4 | 1,368.8 ± 15.7 | 901.7 ± 12.1 | 1.48x* / 1.52x* improved |
| B | Logistic | tolerance | Breeze | 1,044.2 ± 16.9 | 1,031.8 ± 9.4 | 1,048.6 ± 5.8 | 1,039.3 ± 4.7 | |
| B | Rosenbrock | fixed | gale - | 44.8 ± 0.3 | 41.7 ± 0.2 | 44.4 ± 0.2 | 41.7 ± 0.4 | 1.08x* / 1.06x* improved |
| B | Rosenbrock | fixed | Breeze | 77.8 ± 0.3 | 78.8 ± 0.9 | 77.8 ± 0.9 | 78.7 ± 0.7 | |
| B | Rosenbrock | tolerance | gale - | 1,203.6 ± 33.1 | 1,247.4 ± 20.0 | 1,208.7 ± 21.8 | 1,245.5 ± 5.7 | 0.96x / 0.97x*  |
| B | Rosenbrock | tolerance | Breeze | 1,877.2 ± 11.8 | 1,865.9 ± 15.9 | 1,866.2 ± 11.6 | 2,006.6 ± 246.0 | |

`*` marks a new-vs-base gale change with disjoint 99.9% CIs. The new/base column is the gale speed ratio (base time / new time).

## Evaluation counts (from `[breeze-lbfgs]` trial-setup lines, all forks)

| budget | problem | gale iterations / evaluations / status | Breeze iterations / evaluations / reason | seen in runs |
|---|---|---|---|---|
| fixed | logistic | 20 / 21 / IterationLimit | 20 / 22 / max | 4/8 (A-new-r1, A-new-r2, B-new-r1, B-new-r2) |
| fixed | logistic | 20 / 23 / IterationLimit | 20 / 22 / max | 4/8 (A-base-r1, A-base-r2, B-base-r1, B-base-r2) |
| fixed | rosenbrock | 20 / 23 / IterationLimit | 20 / 25 / max | 4/8 (A-new-r1, A-new-r2, B-new-r1, B-new-r2) |
| fixed | rosenbrock | 20 / 36 / IterationLimit | 20 / 25 / max | 4/8 (A-base-r1, A-base-r2, B-base-r1, B-base-r2) |
| tolerance | logistic | 34 / 36 / Converged | 39 / 42 / function | 4/8 (A-new-r1, A-new-r2, B-new-r1, B-new-r2) |
| tolerance | logistic | 46 / 53 / Converged | 39 / 42 / function | 4/8 (A-base-r1, A-base-r2, B-base-r1, B-base-r2) |
| tolerance | rosenbrock | 513 / 629 / Converged | 522 / 614 / function | 4/8 (A-base-r1, A-base-r2, B-base-r1, B-base-r2) |
| tolerance | rosenbrock | 522 / 607 / Converged | 522 / 614 / function | 4/8 (A-new-r1, A-new-r2, B-new-r1, B-new-r2) |

## Summary

**Evaluation counts:** they are identical across rounds and lanes within each bundle, and change only between bundles:

| budget | problem | gale base (iterations / evaluations) | gale new | Breeze |
|---|---|---|---|---|
| fixed | Rosenbrock | 20 / 36 | 20 / **23** | 20 / 25 |
| fixed | logistic | 20 / 23 | 20 / **21** | 20 / 22 |
| tolerance | Rosenbrock | 513 / 629 | 522 / 607 | 522 / 614 |
| tolerance | logistic | 46 / 53 | **34 / 36** | 39 / 42 |

The scaled first step brings fixed-budget Rosenbrock from 36 to 23 evaluations, which is now fewer than Breeze's 25.

**gale new vs base (CI-separated in both rounds unless noted):**
- **Logistic tolerance:** 1.45–1.52× faster, on every lane and backend.
- **Logistic fixed:** 1.06–1.08× faster.
- **Rosenbrock fixed:** 1.06–1.08× faster.
- **Rosenbrock tolerance:** **0.98–0.99×** in lane A, a small CI-separated slowdown in both rounds. Lane B measured 0.96× / 0.97×, separated in round 2 only. This is consistent with the new path taking more iterations (522 vs 513) despite fewer evaluations (607 vs 629).

**gale vs Breeze, new bundle:**
- **Lane A:** every pair is ahead. Logistic fixed is 1.20×, Rosenbrock fixed 1.90×. The tolerance pairs are withheld by the scoreboard (each library stops on its own convergence test), but logistic tolerance moves from 0.93× to 1.30–1.32×.
- **Lane B, vector:** logistic fixed moves from behind (0.94–0.96×) to a tie (0.99×). The withheld logistic tolerance pair moves from 0.77–0.78× to 1.13–1.15×. The lane B headline is now 1 ahead / 1 tie / 0 behind (base: 1 / 0 / 1).
- **Lane B, pure:** logistic stays behind VectorBLAS Breeze, at 0.70–0.71× (fixed) and 0.79× (tolerance).

**Breeze noise:** two Breeze Rosenbrock-tolerance measurements in new runs are outliers: lane A r1 (2,170 ± 10 µs, against about 1,860 elsewhere) and lane B r2 (2,007 ± 246). The Rosenbrock-tolerance gale/Breeze ratios in those two runs (1.74× and 1.61×) are inflated by them. Gale's own new/base ratios are not affected.

Raw data: `2026-10-11-hc-lbfgs-step-trillium.json`. It holds the per-run JMH JSON, netlib records, validation records (new runs), the parsed evaluation counts, source SHA and load, keyed by lane, variant and round.
