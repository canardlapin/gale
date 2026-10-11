# hc-sigmoid paired confirmation on trillium (Zen 5)

Branch-free sigmoid numerator (`567a110`, breeze/hc-sigmoid) against the Milestone 2 baseline. Slurm job 2535855, one exclusive compute node (`tri0119`, 2 × AMD EPYC 9655, Zen 5, AVX-512), 2026-10-10 23:45–23:55 EDT.

## Protocol

- **Bundles:** base is the Milestone 2 bundle compiled from `0ff6440`. Its source is identical to `5d3950e`, which only adds receipts and `supp.sh`. New is `567a110`, compiled locally from a clean worktree with sbt 1.11.7, JDK 25.0.1 and Scala 3.7.4. No sbt ran on the cluster.
- **Harness:** `sig.sh` (`docs/verification/w21-simd-spike/x86-avx512/sig.sh`) is the reviewed `run.sh` (lane fork arguments, `-Xms4g -Xmx4g`, `-Dgale.bench.sourceSha` from the bundle, `-foe true`, fresh netlib sidecar per run) restricted to `.*ReductionBreezeJmh.*[Ss]igmoid.*`. It used the annotation defaults: 2 forks, 5 × 500 ms warmup, 5 × 500 ms measurement. Temurin 25+36; `JAVA_TOOL_OPTIONS` was unset. Pinning was `numactl --physcpubind=0-7 --membind=0`.
- **Order:** interleaved on the same node to cancel drift. Per lane: base r1, new r1, base r2, new r2 (laneA base r1; laneA new r1; laneA base r2; laneA new r2; laneB base r1; laneB new r1; laneB base r2; laneB new r2).
- **Netlib** (blas, lapack, vectorModule, lane) across all runs: [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False, 'A'), ('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True, 'B')]. Lane A ran the `pure` backend; lane B had the Vector module. Sigmoid is backend-insensitive.
- **Not covered:** elementwise ops are not covered by numerical preflight, so there is no validation sidecar. Each run's scoreboard ran with `--strict` (legacy mode) and exited 0.

## Results

Values are mean ± 99.9% CI half-width, in ops/s. The ratio is gale / Breeze (>1 means gale is faster).

| lane | n | round | base gale | new gale | Breeze (base run / new run) | base ratio | new ratio | gale new/base |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 1024 | 1 | 246,789 ± 4,401 | 231,414 ± 2,827 | 291,756 ± 1,537 / 290,083 ± 2,156 | 0.85x | 0.80x | 0.94x |
| A | 1024 | 2 | 246,626 ± 4,363 | 231,211 ± 5,078 | 291,006 ± 760.6 / 290,953 ± 1,203 | 0.85x | 0.79x | 0.94x |
| A | 65536 | 1 | 2,216 ± 30.4 | 3,482 ± 59.9 | 4,282 ± 13.3 / 4,274 ± 23.5 | 0.52x | 0.81x | 1.57x |
| A | 65536 | 2 | 2,275 ± 84.9 | 3,457 ± 10.8 | 4,281 ± 7.0 / 4,287 ± 9.2 | 0.53x | 0.81x | 1.52x |
| A | 1048576 | 1 | 132.8 ± 2.6 | 222.7 ± 0.6 | 186.3 ± 0.2 / 185.9 ± 0.4 | 0.71x | 1.20x | 1.68x |
| A | 1048576 | 2 | 132.8 ± 2.6 | 222.7 ± 0.6 | 184.3 ± 3.1 / 186.5 ± 0.2 | 0.72x | 1.19x | 1.68x |
| B | 1024 | 1 | 247,729 ± 740.1 | 232,501 ± 395.8 | 290,892 ± 1,390 / 290,053 ± 1,137 | 0.85x | 0.80x | 0.94x |
| B | 1024 | 2 | 246,708 ± 2,775 | 232,538 ± 648.3 | 290,970 ± 790.3 / 291,304 ± 1,108 | 0.85x | 0.80x | 0.94x |
| B | 65536 | 1 | 2,214 ± 25.3 | 3,453 ± 11.7 | 4,310 ± 70.1 / 4,283 ± 18.6 | 0.51x | 0.81x | 1.56x |
| B | 65536 | 2 | 2,271 ± 97.3 | 3,452 ± 11.3 | 4,285 ± 16.6 / 4,281 ± 13.1 | 0.53x | 0.81x | 1.52x |
| B | 1048576 | 1 | 132.8 ± 2.6 | 222.7 ± 0.6 | 186.5 ± 0.3 / 184.6 ± 2.8 | 0.71x | 1.21x | 1.68x |
| B | 1048576 | 2 | 132.7 ± 2.6 | 222.8 ± 0.5 | 174.2 ± 17.4 / 173.4 ± 20.4 | 0.76x | 1.29x | 1.68x |

## Summary

- Lane A, n=1024: ratio 0.85/0.85x → 0.80/0.79x; gale throughput ×0.94/×0.94.
- Lane A, n=65536: ratio 0.52/0.53x → 0.81/0.81x; gale throughput ×1.57/×1.52.
- Lane A, n=1048576: ratio 0.71/0.72x → 1.20/1.19x; gale throughput ×1.68/×1.68.
- Lane B, n=1024: ratio 0.85/0.85x → 0.80/0.80x; gale throughput ×0.94/×0.94.
- Lane B, n=65536: ratio 0.51/0.53x → 0.81/0.81x; gale throughput ×1.56/×1.52.
- Lane B, n=1048576: ratio 0.71/0.76x → 1.21/1.29x; gale throughput ×1.68/×1.68.

At n = 64K and 1M the branch-free numerator lifts gale's sigmoid throughput 1.52–1.68× in both lanes, with non-overlapping CIs. n = 1M is now ahead of Breeze (1.19–1.29×). n = 64K is still behind (0.81×). **At n = 1024 it regresses:** gale drops about 6%, from about 246.7K to about 231.4K ops/s, and the CIs do not overlap in either lane, taking the ratio from 0.85× to 0.79–0.80×. Breeze throughput was stable across the interleaved runs, apart from wider CIs on two lane-B n = 1M Breeze rows.

Raw data: `2026-10-11-hc-sigmoid-trillium.json`. It holds the per-run JMH JSON and netlib sidecar records, keyed by lane, variant, round and source SHA.
