# hc-max2: max hill-climb on trillium (Zen 5, lane A)

breeze/hc-max2 `38f1a41` (x86 max candidates `scan4Max` and `keyMax`, with a bit-exact setup self-check). Slurm job 2539883, one exclusive compute node (`tri0182`, 2 × AMD EPYC 9655, Zen 5, AVX-512). 1-minute load at start: 75.47.

## Protocol

- **Build:** `38f1a41`, compiled locally in a clean bench worktree (branch `bench/mx2-38f1a41`), with `benchmarksJVM/clean` run before compiling. sbt 1.11.7, JDK 25.0.1, Scala 3.7.4. No sbt ran on the cluster.
- **Runner:** `docs/verification/w21-simd-spike/x86-avx512/mx2.sh`, a copy of the reviewed `run.sh` setup. Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, `numactl --physcpubind=0-7 --membind=0`. Lane A: `-jvmArgs "-Xms4g -Xmx4g -Dgale.bench.sourceSha=38f1a41… -Dgale.bench.lane=A …"` replaces the host's `--add-modules`, so no fork has the Vector module. `-foe true`.
- **Job 1:** `gale.bench.MaxHillClimbJmh.(curMax|scanMax|scan4Max|keyMax|reduceMax|breezeShapeMax)$ -p n=1024,65536,1048576 -f 8 -wi 5 -i 15 -w 500ms -r 500ms`, throughput.
- **Units:** G elem/s = ops/s × n / 1e9; higher is better. CIs are JMH's 99.9% intervals over all 8 × 15 measurement iterations. Per-fork values are each fork's mean over its 15 iterations, kept so that bimodal JIT outcomes stay visible. The raw per-iteration data is in the JSON.

**Semantics self-check:** `MaxHillClimbJmh`'s trial setup throws on any bit mismatch between candidates. The run used `-foe true`, and all 6 × 3 × 8 = 144 forks completed (job rc = 0), so **no semantics failure** occurred.

## Aggregate (G elem/s, mean [99.9% CI])

| method | n = 1024 | n = 65536 | n = 1048576 |
|---|---:|---:|---:|
| curMax | 5.29 [5.20, 5.38] | 5.86 [5.84, 5.88] | 5.93 [5.92, 5.93] |
| scanMax | 5.31 [5.28, 5.34] | 5.87 [5.87, 5.87] | 5.93 [5.93, 5.93] |
| scan4Max | 5.54 [5.49, 5.58] | 5.84 [5.83, 5.84] | 5.90 [5.90, 5.90] |
| keyMax | 1.50 [1.50, 1.50] | 1.34 [1.29, 1.39] | 1.17 [1.17, 1.17] |
| reduceMax | 3.27 [3.26, 3.28] | 3.66 [3.66, 3.66] | 3.70 [3.70, 3.70] |
| breezeShapeMax | 4.70 [4.67, 4.73] | 4.72 [4.68, 4.77] | 6.08 [5.90, 6.25] |

## Ratios (speed; `*` marks disjoint CIs)

| candidate | vs curMax n = 1024 / 64K / 1M | vs breezeShapeMax n = 1024 / 64K / 1M |
|---|---|---|
| curMax | 1.00× / 1.00× / 1.00× | 1.13×* / 1.24×* / 0.98× |
| scanMax | 1.00× / 1.00× / 1.00× | 1.13×* / 1.24×* / 0.98× |
| scan4Max | 1.05×* / 1.00×* / 1.00×* | 1.18×* / 1.24×* / 0.97× |
| keyMax | 0.28×* / 0.23×* / 0.20×* | 0.32×* / 0.28×* / 0.19×* |
| reduceMax | 0.62×* / 0.62×* / 0.62×* | 0.69×* / 0.77×* / 0.61×* |
| breezeShapeMax | 0.89×* / 0.81×* / 1.03× | 1.00× / 1.00× / 1.00× |

## Per-fork throughput (G elem/s; each fork's mean over 15 iterations, in fork order)

| method | n | fork 1 … 8 | min / max | max/min |
|---|---:|---|---|---:|
| curMax | 1024 | 5.24, 4.99, 5.29, 5.47, 5.31, 5.35, 5.37, 5.30 | 4.99 / 5.47 | 1.10 |
| curMax | 65536 | 5.83, 5.86, 5.87, 5.87, 5.86, 5.86, 5.84, 5.87 | 5.83 / 5.87 | 1.01 |
| curMax | 1048576 | 5.93, 5.93, 5.92, 5.93, 5.93, 5.92, 5.93, 5.93 | 5.92 / 5.93 | 1.00 |
| scanMax | 1024 | 5.31, 5.29, 5.51, 5.24, 5.38, 5.29, 5.23, 5.28 | 5.23 / 5.51 | 1.05 |
| scanMax | 65536 | 5.88, 5.87, 5.88, 5.87, 5.87, 5.87, 5.87, 5.87 | 5.87 / 5.88 | 1.00 |
| scanMax | 1048576 | 5.93, 5.93, 5.93, 5.93, 5.93, 5.93, 5.92, 5.92 | 5.92 / 5.93 | 1.00 |
| scan4Max | 1024 | 5.38, 5.55, 5.59, 5.34, 5.44, 5.62, 5.67, 5.69 | 5.34 / 5.69 | 1.07 |
| scan4Max | 65536 | 5.84, 5.83, 5.84, 5.84, 5.83, 5.83, 5.83, 5.84 | 5.83 / 5.84 | 1.00 |
| scan4Max | 1048576 | 5.90, 5.90, 5.90, 5.89, 5.91, 5.91, 5.90, 5.90 | 5.89 / 5.91 | 1.00 |
| keyMax | 1024 | 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50 | 1.50 / 1.50 | 1.00 |
| keyMax | 65536 | 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34 | 1.34 / 1.34 | 1.00 |
| keyMax | 1048576 | 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17 | 1.17 / 1.17 | 1.00 |
| reduceMax | 1024 | 3.28, 3.27, 3.26, 3.29, 3.27, 3.19, 3.27, 3.29 | 3.19 / 3.29 | 1.03 |
| reduceMax | 65536 | 3.65, 3.66, 3.66, 3.66, 3.65, 3.66, 3.66, 3.66 | 3.65 / 3.66 | 1.00 |
| reduceMax | 1048576 | 3.70, 3.70, 3.70, 3.70, 3.70, 3.70, 3.70, 3.70 | 3.70 / 3.70 | 1.00 |
| breezeShapeMax | 1024 | 4.68, 4.78, 4.65, 4.67, 4.68, 4.77, 4.67, 4.71 | 4.65 / 4.78 | 1.03 |
| breezeShapeMax | 65536 | 4.78, 4.78, 4.79, 4.73, 4.79, 4.79, 4.33, 4.79 | 4.33 / 4.79 | 1.11 |
| breezeShapeMax | 1048576 | 5.97, 6.36, 6.39, 6.40, 6.03, 6.35, 4.80, 6.31 | 4.80 / 6.40 | 1.33 |

## Job 2: profiling at n = 1048576 (`scanMax`, `scan4Max`, `keyMax`, 2 forks)

- **perfasm (job 2539883, same node):** `perf` is available (`perf_event_paranoid = 2`, so user-space events only), but **hsdis is not installed** on the cluster's JDK. `-prof perfasm` ran (rc = 0) and attributes the cycles: 99.96% in the C2 code of each candidate (for example `MaxHillClimbRef$::keyMax`, compile id 1750). The hot regions are printed as raw machine-code bytes, though, not instructions, so there is no instruction-level reading. The output is kept in the JSON.
- **perfnorm fallback (job 2539982, debug node tri0001, exclusive, same runner with `FORCE_PERFNORM=1`):** per-element hardware counters, user space only:

| candidate | G elem/s | cycles/elem | instructions/elem | IPC | branches/elem | branch-misses/elem | L1d loads/elem | L1d misses/elem |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| keyMax | 1.18 | 3.827 | 14.38 | 3.76 | 0.126 | 0.00018 | 3.38 | 0.125 |
| scan4Max | 5.86 | 0.767 | 3.76 | 4.90 | 1.251 | 0.00027 | 1.00 | 0.125 |
| scanMax | 5.89 | 0.763 | 3.75 | 4.92 | 1.251 | 0.00027 | 1.00 | 0.125 |

`scanMax` and `scan4Max` are identical at the instruction level: 3.75 instructions/elem, IPC ≈ 4.9, 1.25 branches/elem, essentially no mispredicts, and 0.76 cycles/elem, which is about 1.3 elem/cycle. `keyMax` takes 0.13 branches/elem but 14.4 instructions and 3.4 L1d loads per element, at 3.8 cycles/elem. Its bit-key transform costs about 4× the instructions of the scan, and lower branch counts do not make up for that.

## Summary

- **`scanMax` ≈ `curMax`:** CIs overlap at every n (1.00×). They run at 5.3 G elem/s at n = 1024 and 5.9 at 64K and 1M.
- **`scan4Max`:** a small but CI-separated gain at n = 1024 (1.05×, 5.54 against 5.29 G elem/s). It ties at 64K and is marginally slower at 1M (0.995×, separated).
- **`keyMax` is not viable on Zen 5:** 0.20–0.28× of `curMax`, at 1.2–1.5 G elem/s.
- **`reduceMax`:** 0.62× of `curMax`.
- **Against `breezeShapeMax`:** `curMax`, `scanMax` and `scan4Max` are ahead at n = 1024 (1.13–1.18×) and 64K (1.24×). At 1M they are within the Breeze-shape CI (0.97–0.98×, overlapping).
- **Bimodal JIT modes:** visible only in `breezeShapeMax`. At n = 1M, fork 7 ran at 4.80 G elem/s against 5.97–6.40 for the others, and at 64K one fork ran at 4.33 against about 4.78. Every gale candidate is unimodal across its 8 forks (max/min ≤ 1.10 at n = 1024, ≤ 1.01 at 64K/1M).
- **No semantics failure:** the bit-exact setup self-check passed in all 144 forks.

Raw data: `2026-10-11-hc-max2-trillium.json`. It holds the full JMH JSON with per-fork `rawData`, plus both profiling outputs (the perfasm run without hsdis, and the perfnorm fallback).
