# hc-lbfgs paired run on trillium (Zen 5, lanes A and B)

`16d8712` (breeze/hc-lbfgs: lower L-BFGS per-iteration overhead without changing a bit) against the Milestone 2 baseline. Slurm job 2537627, one exclusive compute node (`tri0006`, 2 × AMD EPYC 9655, Zen 5, AVX-512).

## Protocol

- **Bundles:** base is the Milestone 2 bundle compiled from `0ff6440`, which is source-identical to `5d3950e`, the branch base. New is `16d8712`, compiled locally in a clean worktree (local branch `bench/hclb-16d8712`) with sbt 1.11.7, JDK 25.0.1 and Scala 3.7.4. No sbt ran on the cluster.
- **Harness:** `docs/verification/w21-simd-spike/x86-avx512/lb.sh`, which is the reviewed `run.sh` lane paths restricted to `gale.bench.LbfgsBreezeJmh.*`. Lane A uses `-jvmArgs` (no Vector module) and `-p backend=pure`. Lane B uses the inherited `--add-modules` and `-p backend=pure,vector`. Both use `-Xms4g -Xmx4g`, a bundle-bound `gale.bench.sourceSha`, `-foe true`, a fresh netlib sidecar per run, and the annotation defaults: 2 forks, 5 × 500 ms warmup, 5 × 500 ms measurement, average time in µs/op. Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, `numactl --physcpubind=0-7 --membind=0`.
- **Order:** interleaved: laneA-base-r1, laneA-new-r1, laneB-base-r1, laneB-new-r1, laneA-base-r2, laneA-new-r2, laneB-base-r2, laneB-new-r2.
- **Load at start:** A-base-r1: 0.00; A-new-r1: 0.98; B-base-r1: 1.31; B-new-r1: 1.24; A-base-r2: 1.08; A-new-r2: 1.17; B-base-r2: 1.39; B-new-r2: 1.14 (1-minute).
- **Netlib** (blas, lapack, vectorModule): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False), ('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- **Validation:** `LbfgsBreezeJmh` is not a preflight-covered class. `BreezeBenchmarkValidation` covers the dense decomposition, factorization, least-squares, multi-RHS, small-dense and symmetric-eigen classes (`PREFLIGHT_CLASSES` in `tools/bench/breeze_scoreboard.py`). Every fork still carried `-Dgale.bench.sourceSha` bound to its bundle (`16d8712` or `0ff6440`), but no run wrote a validation sidecar, because there were no covered checks. Each run's scoreboard was run as `breeze_scoreboard.py --lane {A,B} --strict --validation <empty sidecar>` and exited 0, reporting "0 decomposition/factorization result rows validated" with 0 unpaired rows. The fixed-budget harness's own guard (it throws unless both libraries ran exactly 20 iterations) passed in every fork. The evaluation counts below are the functional check.

## gale / Breeze (speed ratio, >1 means gale is faster; verdict from 99.9% CIs)

| lane | problem | budget | gale backend | base r1 | base r2 | new r1 | new r2 |
|---|---|---|---|---|---|---|---|
| A | Logistic | fixed | pure | 1.11x ahead | 1.10x ahead | 1.14x ahead | 1.14x ahead |
| A | Logistic | tolerance | pure | 0.90x behind | 0.90x behind | 0.93x behind | 0.93x behind |
| A | Rosenbrock | fixed | - | 0.74x behind | 0.75x behind | 1.78x ahead | 1.75x ahead |
| A | Rosenbrock | tolerance | - | 0.73x behind | 0.73x behind | 1.51x ahead | 1.51x ahead |
| B | Logistic | fixed | pure | 0.66x behind | 0.66x behind | 0.68x behind | 0.68x behind |
| B | Logistic | fixed | vector | 0.89x behind | 0.88x behind | 0.93x behind | 0.94x behind |
| B | Logistic | tolerance | pure | 0.54x behind | 0.54x behind | 0.55x behind | 0.55x behind |
| B | Logistic | tolerance | vector | 0.73x behind | 0.73x behind | 0.76x behind | 0.76x behind |
| B | Rosenbrock | fixed | - | 0.74x behind | 0.76x behind | 1.74x ahead | 1.74x ahead |
| B | Rosenbrock | tolerance | - | 0.79x behind | 0.79x behind | 1.52x ahead | 1.54x ahead |

## Times (µs/op, mean ± 99.9% CI half-width)

| lane | problem | budget | series | base r1 | new r1 | base r2 | new r2 | gale new/base speed r1 / r2 |
|---|---|---|---|---:|---:|---:|---:|---|
| A | Logistic | fixed | gale pure | 842.4 ± 7.6 | 813.6 ± 9.8 | 843.5 ± 8.4 | 816.3 ± 9.2 | 1.04x* / 1.03x* improved |
| A | Logistic | fixed | Breeze | 931.7 ± 15.6 | 930.3 ± 12.5 | 929.1 ± 16.0 | 930.4 ± 14.4 | |
| A | Logistic | tolerance | gale pure | 1,964.3 ± 19.4 | 1,893.9 ± 18.2 | 1,956.2 ± 15.6 | 1,898.8 ± 24.4 | 1.04x* / 1.03x* improved |
| A | Logistic | tolerance | Breeze | 1,763.3 ± 33.0 | 1,766.7 ± 27.8 | 1,769.7 ± 37.6 | 1,760.0 ± 26.6 | |
| A | Rosenbrock | fixed | gale - | 105.0 ± 2.5 | 44.7 ± 0.2 | 105.4 ± 2.0 | 44.5 ± 0.2 | 2.35x* / 2.37x* improved |
| A | Rosenbrock | fixed | Breeze | 78.0 ± 0.4 | 79.5 ± 1.2 | 78.9 ± 0.5 | 78.1 ± 0.6 | |
| A | Rosenbrock | tolerance | gale - | 2,581.0 ± 48.1 | 1,233.0 ± 4.7 | 2,581.6 ± 48.8 | 1,232.1 ± 3.2 | 2.09x* / 2.10x* improved |
| A | Rosenbrock | tolerance | Breeze | 1,875.3 ± 36.3 | 1,860.0 ± 18.9 | 1,877.4 ± 18.8 | 1,861.4 ± 12.5 | |
| B | Logistic | fixed | gale pure | 842.1 ± 6.2 | 817.8 ± 9.7 | 843.1 ± 6.7 | 816.4 ± 4.7 | 1.03x* / 1.03x* improved |
| B | Logistic | fixed | Breeze | 555.5 ± 3.0 | 559.6 ± 6.6 | 555.1 ± 3.6 | 556.7 ± 5.9 | |
| B | Logistic | fixed | gale vector | 626.0 ± 13.9 | 600.9 ± 9.8 | 628.9 ± 8.9 | 593.7 ± 9.1 | 1.04x* / 1.06x* improved |
| B | Logistic | fixed | Breeze | 555.5 ± 3.0 | 559.6 ± 6.6 | 555.1 ± 3.6 | 556.7 ± 5.9 | |
| B | Logistic | tolerance | gale pure | 1,963.2 ± 17.4 | 1,903.4 ± 25.9 | 1,964.7 ± 11.8 | 1,897.7 ± 18.2 | 1.03x* / 1.04x* improved |
| B | Logistic | tolerance | Breeze | 1,054.9 ± 13.9 | 1,048.7 ± 12.1 | 1,054.6 ± 6.2 | 1,047.5 ± 10.6 | |
| B | Logistic | tolerance | gale vector | 1,436.8 ± 18.7 | 1,376.5 ± 22.4 | 1,448.5 ± 13.2 | 1,381.3 ± 25.3 | 1.04x* / 1.05x* improved |
| B | Logistic | tolerance | Breeze | 1,054.9 ± 13.9 | 1,048.7 ± 12.1 | 1,054.6 ± 6.2 | 1,047.5 ± 10.6 | |
| B | Rosenbrock | fixed | gale - | 104.1 ± 0.3 | 45.2 ± 0.2 | 104.1 ± 0.6 | 44.7 ± 0.2 | 2.30x* / 2.33x* improved |
| B | Rosenbrock | fixed | Breeze | 77.4 ± 0.4 | 78.6 ± 0.9 | 78.7 ± 1.8 | 77.7 ± 0.9 | |
| B | Rosenbrock | tolerance | gale - | 2,569.7 ± 15.9 | 1,238.6 ± 10.6 | 2,569.2 ± 12.1 | 1,212.6 ± 34.9 | 2.07x* / 2.12x* improved |
| B | Rosenbrock | tolerance | Breeze | 2,039.9 ± 35.3 | 1,878.9 ± 7.1 | 2,030.7 ± 252.2 | 1,873.0 ± 4.7 | |

`*` marks a new-vs-base gale change with disjoint 99.9% CIs. The new/base column is the gale speed ratio (base time / new time).

## Evaluation counts (from `[breeze-lbfgs]` trial-setup lines, all forks)

| budget | problem | gale iterations / evaluations / status | Breeze iterations / evaluations / reason | seen in runs |
|---|---|---|---|---|
| fixed | logistic | 20 / 23 / IterationLimit | 20 / 22 / max | 8/8 (all) |
| fixed | rosenbrock | 20 / 36 / IterationLimit | 20 / 25 / max | 8/8 (all) |
| tolerance | logistic | 46 / 53 / Converged | 39 / 42 / function | 8/8 (all) |
| tolerance | rosenbrock | 513 / 629 / Converged | 522 / 614 / function | 8/8 (all) |

## Summary

**Evaluation counts are unchanged.** They are identical in all 8 runs, base and new, in both lanes:

| budget | problem | gale (iterations / evaluations) | Breeze (iterations / evaluations) |
|---|---|---|---|
| fixed | Rosenbrock | 20 / 36 | 20 / 25 |
| fixed | logistic | 20 / 23 | 20 / 22 |
| tolerance | Rosenbrock | 513 / 629 | 522 / 614 |
| tolerance | logistic | 46 / 53 | 39 / 42 |

The fixed-budget counts are the expected 36/23 against 25/22, so `16d8712` does the same work as base. This is consistent with its no-bit-change claim at the level this benchmark observes.

**gale new vs base, CI-separated in both rounds for every series:**
- **Rosenbrock:** 2.30–2.37× faster at the fixed budget (105.0 ± 2.5 → 44.7 ± 0.2 µs in lane A) and 2.07–2.12× faster at the tolerance budget (2,581 ± 48 → 1,233 ± 5 µs).
- **Logistic:** 1.03–1.06× faster in both lanes, on both backends.
- **Regressions:** none.

**gale vs Breeze, new bundle:**
- **Rosenbrock fixed** goes from behind (0.74–0.76×) to **ahead at 1.74–1.78×** in both lanes.
- **Rosenbrock tolerance** goes from 0.73–0.79× to 1.51–1.54×. The scoreboard withholds this pair as not like-for-like, because each library stops on its own convergence test.
- **Logistic fixed, lane A (pure vs scalar Breeze):** ahead, 1.10× → 1.14×.
- **Logistic fixed, lane B (VectorBLAS Breeze):** still behind, at 0.68× pure and 0.93–0.94× vector (from 0.66× and 0.88–0.89×).
- **Logistic tolerance (withheld pair):** 0.93× in lane A; 0.55× pure and 0.76× vector in lane B.

**Breeze noise:** in lane B, Breeze's Rosenbrock tolerance time in the base runs (2,040 ± 35 and 2,031 ± 252 µs) was higher than in the new runs (about 1,875 µs). Part of that row's ratio change therefore comes from Breeze run-to-run variation, not from gale; gale's own new/base change (2.07–2.12×) is measured directly. All other Breeze series are stable across the eight runs.

Raw data: `2026-10-11-hc-lbfgs-trillium.json`. It holds the per-run JMH JSON, netlib records, any validation records, the parsed evaluation counts, source SHA and load, keyed by lane, variant and round.
