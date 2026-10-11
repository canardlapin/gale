# hc-sigmoid round 2: three-way paired run on trillium (Zen 5)

Baseline vs branch-free (`567a110`) vs hybrid (`70fb262`: `1/(1+exp(-x))` for x ≥ -36, `t/(1+t)` below). Slurm job 2535952, one exclusive compute node (`tri0281`, 2 × AMD EPYC 9655, Zen 5, AVX-512), 2026-10-11 00:13–00:27 EDT.

## Protocol

- **Same as round 1** (`2026-10-11-hc-sigmoid-trillium.md`, job 2535855): `docs/verification/w21-simd-spike/x86-avx512/sig.sh` restricted to `.*ReductionBreezeJmh.*[Ss]igmoid.*`, with lane fork arguments, `-Xms4g -Xmx4g`, a bundle-bound `gale.bench.sourceSha`, a fresh netlib sidecar per run, 2 forks with the annotation warmup and measurement, Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, and `numactl --physcpubind=0-7 --membind=0`. No sbt ran on the cluster.
- **Bundles:** base is the Milestone 2 bundle (`0ff6440`, source-identical to `5d3950e`). `567a110` and `70fb262` were compiled locally from the clean breeze/hc-sigmoid worktree.
- **Order:** interleaved per lane as base, 567a110, 70fb262, repeated twice: laneA base r1; laneA b567 r1; laneA h70f r1; laneA base r2; laneA b567 r2; laneA h70f r2; laneB base r1; laneB b567 r1; laneB h70f r1; laneB base r2; laneB b567 r2; laneB h70f r2.
- **Load:** the node's 1-minute load average at job start was 369, left over from the previous allocation. It decayed at the rate of an idle machine: 115 after about 70 s, 37 after about 140 s, 13 before lane A round 2, and about 1 for all of lane B. Round 1 and round 2 of lane A agree closely (see table), so the decaying load did not show up in the measurements.
- **Netlib** (blas, lapack, vectorModule, lane): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False, 'A'), ('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True, 'B')].
- **Not covered:** sigmoid is not covered by numerical preflight, so there is no validation sidecar. Each of the 12 runs' scoreboards ran with `--strict` (legacy mode) and exited 0.
- **Mixed input not run:** `ReductionBreezeJmh` has no input-distribution parameter, and no benchmark code was added, so the uniform-[-750, 750]-style mixed input was not run.

## gale / Breeze ratio (Breeze measured in the same run)

The ratio is followed by gale ops/s, mean ± 99.9% CI half-width, for each round. Ratios >1 mean gale is faster.

| lane | n | round | base | 567a110 | 70fb262 |
|---|---:|---:|---|---|---|
| A | 1024 | 1 | 0.85x (247,110 ± 4,169) | 0.80x (231,959 ± 3,773) | 0.96x (279,937 ± 2,079) |
| A | 1024 | 2 | 0.85x (247,250 ± 4,115) | 0.80x (232,160 ± 3,613) | 0.96x (280,142 ± 2,071) |
| A | 65536 | 1 | 0.53x (2,319 ± 89.9) | 0.81x (3,523 ± 22.7) | 0.96x (4,171 ± 22.4) |
| A | 65536 | 2 | 0.54x (2,354 ± 36.3) | 0.81x (3,502 ± 48.6) | 0.96x (4,173 ± 21.4) |
| A | 1048576 | 1 | 0.71x (133.2 ± 2.6) | 1.20x (223.8 ± 0.8) | 1.42x (265.1 ± 0.9) |
| A | 1048576 | 2 | 0.71x (133.2 ± 2.6) | 1.19x (223.3 ± 0.7) | 1.42x (265.4 ± 1.0) |
| B | 1024 | 1 | 0.86x (249,513 ± 2,192) | 0.80x (233,348 ± 806.3) | 0.97x (279,757 ± 640.4) |
| B | 1024 | 2 | 0.85x (248,421 ± 1,059) | 0.81x (233,550 ± 845.1) | 0.96x (281,118 ± 759.2) |
| B | 65536 | 1 | 0.53x (2,288 ± 53.2) | 0.79x (3,498 ± 54.8) | 0.96x (4,168 ± 30.9) |
| B | 65536 | 2 | 0.53x (2,307 ± 65.8) | 0.80x (3,499 ± 50.2) | 0.97x (4,177 ± 33.2) |
| B | 1048576 | 1 | 0.71x (133.0 ± 2.6) | 1.20x (223.6 ± 1.0) | 1.45x (265.4 ± 1.0) |
| B | 1048576 | 2 | 0.71x (133.3 ± 2.7) | 1.20x (223.7 ± 0.9) | 1.42x (265.3 ± 0.7) |

## Breeze throughput per run (stability check)

| lane | n | round | in base run | in 567a110 run | in 70fb262 run |
|---|---:|---:|---|---|---|
| A | 1024 | 1 | 290,375 ± 2,131 | 290,987 ± 2,800 | 292,320 ± 1,219 |
| A | 1024 | 2 | 291,040 ± 987.1 | 290,571 ± 999.8 | 291,015 ± 992.3 |
| A | 65536 | 1 | 4,385 ± 38.6 | 4,346 ± 84.0 | 4,337 ± 80.6 |
| A | 65536 | 2 | 4,352 ± 60.8 | 4,350 ± 78.7 | 4,328 ± 60.6 |
| A | 1048576 | 1 | 187.0 ± 0.5 | 187.1 ± 0.2 | 187.0 ± 0.3 |
| A | 1048576 | 2 | 186.9 ± 0.3 | 186.9 ± 0.4 | 187.1 ± 0.3 |
| B | 1024 | 1 | 290,841 ± 1,296 | 290,756 ± 1,178 | 289,680 ± 1,582 |
| B | 1024 | 2 | 290,696 ± 1,447 | 289,823 ± 1,464 | 291,358 ± 982.5 |
| B | 65536 | 1 | 4,326 ± 63.5 | 4,402 ± 29.7 | 4,344 ± 78.6 |
| B | 65536 | 2 | 4,362 ± 73.0 | 4,379 ± 58.5 | 4,323 ± 54.8 |
| B | 1048576 | 1 | 186.9 ± 0.3 | 187.0 ± 0.2 | 183.6 ± 5.1 |
| B | 1048576 | 2 | 187.1 ± 0.2 | 187.0 ± 0.2 | 186.2 ± 1.8 |

## gale throughput relative to base (same lane and round)

| lane | n | 567a110 r1 / r2 | 70fb262 r1 / r2 | CI-separated from base? (567 / 70f) |
|---|---:|---|---|---|
| A | 1024 | 0.94x / 0.94x | 1.13x / 1.13x | both / both |
| A | 65536 | 1.52x / 1.49x | 1.80x / 1.77x | both / both |
| A | 1048576 | 1.68x / 1.68x | 1.99x / 1.99x | both / both |
| B | 1024 | 0.94x / 0.94x | 1.12x / 1.13x | both / both |
| B | 65536 | 1.53x / 1.52x | 1.82x / 1.81x | both / both |
| B | 1048576 | 1.68x / 1.68x | 2.00x / 1.99x | both / both |

## Summary

- **70fb262 (hybrid) is the fastest variant at every size in both lanes.** gale throughput relative to base is 1.12–1.13× at n = 1024, 1.77–1.82× at 64K and 1.99–2.00× at 1M, with CIs separated from base in every round.
- **It removes the 567a110 regression at n = 1024.** 567a110 measured 0.94× of base there, as in round 1; the hybrid is 1.13×.
- **Against Breeze, the hybrid leads only at n = 1M.** It is ahead there, at 1.42–1.45×. At n = 1024 and 64K it is still slightly behind, at 0.96–0.97×, with gale's upper CI bound below Breeze's lower bound in every run. For example, at lane A n = 1024 round 1, gale is 279,937 ± 2,079 and Breeze 292,320 ± 1,219 ops/s.
- **Breeze throughput was stable across all 12 runs** (table above). The only wider CI is one lane-B n = 1M row in the hybrid run (183.6 ± 5.1).

Raw data: `2026-10-11-hc-sigmoid-trillium-r2.json`. It holds the per-run JMH JSON, netlib sidecar records, source SHA and load at start, keyed by lane, variant and round.
