# hc-max2 facade A/B on trillium (Zen 5, lane A)

Candidate `75f9191` (breeze/hc-max2: dense `DVec.max`/`min` pass a constant offset and stride to the kernel) against base `96c95e5` (Milestone 3, which showed the slow `galeMax` mode). Slurm job 2540120, one exclusive compute node (`tri0098`, 2 × AMD EPYC 9655, Zen 5, AVX-512).

## Protocol

- **Bundles:** `75f9191`, compiled locally in a clean bench worktree (branch `bench/mx2b-75f9191`) with `benchmarksJVM/clean` first. Base `96c95e5` is the Milestone 3 bundle already on the cluster (`/scratch/brad/gale-bench/m3/bundle`, SHA file `96c95e5…`), the same bundle that produced the M3 receipts. sbt 1.11.7, JDK 25.0.1, Scala 3.7.4. No sbt ran on the cluster.
- **Runner:** `docs/verification/w21-simd-spike/x86-avx512/mf.sh`, a copy of the reviewed `run.sh` setup, with lane-A fork flags (`-jvmArgs "-Xms4g -Xmx4g … lane=A …"`, so no Vector module, and `-p backend=pure`). Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, `numactl --physcpubind=0-7 --membind=0`, `-foe true`.
- **Job 1:** a same-jar A/B at `75f9191`: `MaxHillClimbJmh.(curMax|facadeMax|generalMax)$ -p n=65536,1048576 -f 8 -wi 5 -i 15 -w 500ms -r 500ms`.
- **Job 2:** the scoreboard facade, `ReductionBreezeJmh.(galeMax|breezeMax)$`, same settings, base and new interleaved for 2 rounds.
- **perfnorm:** JMH reports `-prof perfnorm` as one aggregate across the forks of an invocation, and the per-fork raw data comes back empty. So every perfnorm configuration ran as **8 separate `-f 1` invocations**, each with the same warmup and measurement as above. That gives per-fork cycles/elem next to that fork's own throughput. The plain `-f 8` runs are the main per-fork distribution.
- **Units:** G elem/s = ops/s × n / 1e9, per fork = the mean of its 15 iterations. A fork is counted **slow** when it runs at ≤ 4.6 G elem/s.
- **Order:** ab-new, abpn-new, fac-base-r1, fac-new-r1, fac-base-r2, fac-new-r2, facpn-base, facpn-new.

## Job 1: same-jar A/B (75f9191), per-fork G elem/s

| method | n | forks 1…8 | min / max | slow forks |
|---|---:|---|---|---:|
| curMax | 65536 | 5.83, 5.83, 5.82, 5.83, 6.04, 5.83, 5.83, 5.83 | 5.82 / 6.04 | 0/8 |
| curMax | 1048576 | 5.88, 5.89, 5.89, 5.89, 5.88, 5.88, 5.89, 5.88 | 5.88 / 5.89 | 0/8 |
| facadeMax | 65536 | 5.82, 5.83, 5.83, 5.83, 5.83, 5.82, 5.82, 5.83 | 5.82 / 5.83 | 0/8 |
| facadeMax | 1048576 | 5.84, 5.88, 5.84, 5.86, 5.85, 5.85, 5.85, 5.85 | 5.84 / 5.88 | 0/8 |
| generalMax | 65536 | 5.83, 5.83, 5.83, 5.83, 5.83, 5.82, 5.83, 5.83 | 5.82 / 5.83 | 0/8 |
| generalMax | 1048576 | 5.84, 5.85, 5.84, 5.84, 5.85, 5.85, 5.85, 5.84 | 5.84 / 5.85 | 0/8 |

| method | n | mean [99.9% CI] |
|---|---:|---|
| curMax | 65536 | 5.86 [5.82, 5.89] |
| curMax | 1048576 | 5.89 [5.89, 5.89] |
| facadeMax | 65536 | 5.82 [5.82, 5.83] |
| facadeMax | 1048576 | 5.85 [5.85, 5.86] |
| generalMax | 65536 | 5.83 [5.83, 5.83] |
| generalMax | 1048576 | 5.84 [5.84, 5.85] |

## Job 2: ReductionBreezeJmh max, base 96c95e5 vs new 75f9191, per-fork G elem/s

| run | method | n | forks 1…8 | min / max | slow forks |
|---|---|---:|---|---|---:|
| fac-base-r1 | galeMax | 65536 | 4.42, 5.83, 5.81, 5.82, 5.86, 5.78, 5.83, 5.83 | 4.42 / 5.86 | 1/8 |
| fac-base-r1 | galeMax | 1048576 | 5.19, 5.85, 4.43, 5.85, 4.43, 5.21, 5.87, 5.83 | 4.43 / 5.87 | 2/8 |
| fac-base-r1 | breezeMax | 65536 | 6.05, 6.20, 5.06, 5.84, 4.44, 6.00, 5.64, 6.08 | 4.44 / 6.20 | 1/8 |
| fac-base-r1 | breezeMax | 1048576 | 6.13, 5.90, 6.12, 5.53, 5.89, 5.86, 5.75, 5.89 | 5.53 / 6.13 | 0/8 |
| fac-new-r1 | galeMax | 65536 | 4.42, 4.59, 4.43, 5.84, 4.42, 4.42, 4.42, 4.42 | 4.42 / 5.84 | 7/8 |
| fac-new-r1 | galeMax | 1048576 | 5.85, 5.84, 5.85, 5.85, 5.85, 5.87, 5.85, 5.85 | 5.84 / 5.87 | 0/8 |
| fac-new-r1 | breezeMax | 65536 | 6.13, 5.98, 6.08, 5.86, 6.26, 6.07, 5.80, 5.82 | 5.80 / 6.26 | 0/8 |
| fac-new-r1 | breezeMax | 1048576 | 5.60, 4.46, 6.16, 5.84, 5.92, 5.92, 5.88, 6.13 | 4.46 / 6.16 | 1/8 |
| fac-base-r2 | galeMax | 65536 | 5.82, 5.83, 5.82, 5.81, 5.83, 5.77, 5.83, 5.82 | 5.77 / 5.83 | 0/8 |
| fac-base-r2 | galeMax | 1048576 | 4.42, 5.85, 5.85, 5.85, 4.43, 5.84, 5.84, 5.84 | 4.42 / 5.85 | 2/8 |
| fac-base-r2 | breezeMax | 65536 | 6.01, 6.03, 5.81, 5.66, 5.78, 6.20, 6.12, 5.81 | 5.66 / 6.20 | 0/8 |
| fac-base-r2 | breezeMax | 1048576 | 6.16, 4.46, 4.65, 5.86, 5.87, 5.13, 5.88, 6.19 | 4.46 / 6.19 | 1/8 |
| fac-new-r2 | galeMax | 65536 | 4.42, 4.46, 5.19, 4.42, 4.42, 4.42, 5.82, 4.59 | 4.42 / 5.82 | 6/8 |
| fac-new-r2 | galeMax | 1048576 | 5.85, 5.85, 5.85, 5.86, 5.85, 5.85, 5.85, 5.85 | 5.85 / 5.86 | 0/8 |
| fac-new-r2 | breezeMax | 65536 | 6.09, 5.45, 6.08, 6.06, 6.13, 5.05, 6.08, 5.57 | 5.05 / 6.13 | 0/8 |
| fac-new-r2 | breezeMax | 1048576 | 5.87, 5.95, 5.90, 5.14, 6.14, 5.58, 5.90, 5.90 | 5.14 / 6.14 | 0/8 |

| run | n | galeMax mean [CI] | breezeMax mean [CI] | gale/Breeze |
|---|---:|---|---|---|
| fac-base-r1 | 65536 | 5.65 [5.50, 5.79] | 5.66 [5.48, 5.85] | 1.00× |
| fac-base-r1 | 1048576 | 5.33 [5.12, 5.54] | 5.88 [5.82, 5.95] | 0.91×* |
| fac-new-r1 | 65536 | 4.62 [4.47, 4.77] | 6.00 [5.93, 6.07] | 0.77×* |
| fac-new-r1 | 1048576 | 5.85 [5.85, 5.86] | 5.74 [5.58, 5.90] | 1.02× |
| fac-base-r2 | 65536 | 5.82 [5.77, 5.86] | 5.93 [5.86, 6.00] | 0.98× |
| fac-base-r2 | 1048576 | 5.49 [5.30, 5.68] | 5.52 [5.29, 5.76] | 0.99× |
| fac-new-r2 | 65536 | 4.72 [4.54, 4.89] | 5.82 [5.68, 5.95] | 0.81×* |
| fac-new-r2 | 1048576 | 5.85 [5.85, 5.85] | 5.80 [5.69, 5.90] | 1.01× |

## perfnorm, per fork (8 × `-f 1`; G elem/s, cycles/elem, instructions/elem, IPC)

Forks shown as `n/a` returned implausible counter values (thousands of cycles/elem). That is a perf collection or multiplexing glitch in that invocation, so their counters are excluded. Their throughput is real and still shown.

| config | method | n | per fork: G/s @ cycles/elem (insn/elem, IPC) |
|---|---|---:|---|
| abpn-new | curMax | 65536 | 5.82 @ 0.772 (3.75, 4.87); 5.82 @ 0.772 (3.75, 4.87); 5.81 @ n/a; 5.82 @ 0.771 (3.76, 4.87); 5.81 @ 0.772 (3.75, 4.87); 5.82 @ 0.771 (3.75, 4.87); 5.81 @ 0.771 (3.76, 4.87); 5.82 @ 0.770 (3.75, 4.87) |
| abpn-new | curMax | 1048576 | 5.87 @ 0.763 (3.75, 4.92); 5.87 @ 0.763 (3.75, 4.92); 5.87 @ 0.763 (3.75, 4.92); 5.87 @ 0.763 (3.75, 4.92); 5.86 @ 0.764 (3.75, 4.91); 5.87 @ 0.763 (3.76, 4.92); 5.87 @ 0.763 (3.75, 4.92); 5.87 @ n/a |
| abpn-new | facadeMax | 65536 | 5.81 @ n/a; 5.81 @ 0.772 (3.75, 4.87); 5.81 @ 0.771 (3.75, 4.87); 5.81 @ 0.772 (3.76, 4.87); 5.81 @ 0.772 (3.76, 4.87); 5.81 @ 0.771 (3.75, 4.87); 5.81 @ 0.772 (3.76, 4.87); 5.82 @ n/a |
| abpn-new | facadeMax | 1048576 | 5.83 @ 0.768 (3.76, 4.89); 5.83 @ 0.768 (3.76, 4.89); 5.83 @ 0.768 (3.75, 4.89); 4.99 @ 0.898 (3.76, 4.19); 5.84 @ 0.768 (3.76, 4.89); 5.83 @ n/a; 5.83 @ 0.768 (3.75, 4.89); 5.83 @ 0.768 (3.75, 4.88) |
| abpn-new | generalMax | 65536 | 5.81 @ 0.771 (4.00, 5.19); 5.81 @ 0.772 (4.00, 5.18); 5.81 @ 0.772 (4.00, 5.19); 5.82 @ 0.772 (4.01, 5.19); 5.81 @ 0.772 (4.00, 5.19); 5.81 @ 0.772 (4.00, 5.19); 5.82 @ 0.771 (4.00, 5.19); 5.81 @ 0.772 (4.00, 5.19) |
| abpn-new | generalMax | 1048576 | 5.83 @ 0.768 (4.00, 5.21); 5.82 @ 0.769 (4.01, 5.21); 5.81 @ 0.771 (4.01, 5.20); 5.83 @ n/a; 5.82 @ 0.770 (4.01, 5.20); 5.82 @ 0.769 (4.01, 5.21); 5.81 @ 0.770 (4.01, 5.20); 5.83 @ 0.768 (4.00, 5.21) |
| facpn-base | breezeMax | 65536 | 5.83 @ 0.769 (3.50, 4.56); 5.51 @ 0.814 (3.50, 4.29); 6.13 @ 0.731 (3.50, 4.79); 6.12 @ n/a; 6.07 @ 0.738 (3.50, 4.75); 6.03 @ 0.744 (3.51, 4.72); 6.10 @ 0.735 (3.51, 4.77); 5.93 @ 0.756 (3.50, 4.63) |
| facpn-base | breezeMax | 1048576 | 5.62 @ 0.797 (3.50, 4.40); 5.61 @ 0.798 (3.50, 4.39); 5.89 @ 0.760 (3.50, 4.61); 5.89 @ 0.760 (3.50, 4.61); 5.89 @ 0.760 (3.50, 4.61); 5.89 @ 0.760 (3.50, 4.61); 5.85 @ 0.765 (3.50, 4.58); 6.03 @ 0.742 (3.50, 4.72) |
| facpn-base | galeMax | 65536 | 5.80 @ n/a; 5.81 @ 0.771 (4.01, 5.20); 5.81 @ 0.772 (4.01, 5.19); 5.80 @ 0.772 (4.01, 5.19); 5.81 @ 0.772 (4.00, 5.19); 5.81 @ n/a; 5.81 @ 0.773 (4.00, 5.18); 5.79 @ 0.773 (4.00, 5.18) |
| facpn-base | galeMax | 1048576 | 5.83 @ 0.768 (4.00, 5.21); 5.77 @ 0.775 (4.00, 5.16); 4.42 @ 1.014 (4.01, 3.95); 4.42 @ 1.015 (4.00, 3.95); 5.82 @ 0.769 (4.00, 5.20); 4.42 @ 1.014 (4.01, 3.95); 4.41 @ 1.014 (4.00, 3.95); 5.85 @ 0.764 (4.00, 5.24) |
| facpn-new | breezeMax | 65536 | 6.13 @ n/a; 6.11 @ 0.734 (3.50, 4.77); 5.74 @ 0.781 (3.50, 4.49); 5.81 @ 0.771 (3.50, 4.54); 5.04 @ 0.889 (3.50, 3.94); 5.81 @ n/a; 5.04 @ 0.890 (3.51, 3.94); 6.17 @ 0.726 (3.50, 4.83) |
| facpn-new | breezeMax | 1048576 | 6.15 @ 0.728 (3.50, 4.81); 5.89 @ 0.760 (3.50, 4.61); 6.11 @ 0.733 (3.50, 4.77); 5.64 @ 0.794 (3.50, 4.41); 5.54 @ 0.807 (3.50, 4.34); 6.03 @ 0.742 (3.51, 4.73); 4.45 @ 1.007 (3.50, 3.48); 6.23 @ 0.718 (3.50, 4.88) |
| facpn-new | galeMax | 65536 | 5.90 @ 0.760 (3.75, 4.94); 5.81 @ 0.772 (3.76, 4.86); 4.65 @ n/a; 5.83 @ 0.769 (3.75, 4.88); 5.81 @ 0.772 (3.76, 4.87); 5.85 @ 0.767 (3.75, 4.89); 4.41 @ 1.017 (3.75, 3.69); 5.82 @ n/a |
| facpn-new | galeMax | 1048576 | 5.83 @ 0.768 (3.75, 4.89); 5.82 @ 0.769 (3.76, 4.89); 5.83 @ 0.768 (3.75, 4.89); 5.83 @ 0.767 (3.75, 4.89); 5.83 @ 0.768 (3.75, 4.89); 5.83 @ 0.768 (3.75, 4.88); 5.84 @ 0.767 (3.75, 4.89); 5.83 @ 0.768 (3.75, 4.89) |

## Summary

**Job 1, same-jar A/B at `75f9191` (`MaxHillClimbJmh`): no slow forks for any method.** `curMax`, `facadeMax` and `generalMax` all run at 5.82–5.89 G elem/s at both n, in all 8 forks (0/8 slow everywhere). In the separate perfnorm forks, one `facadeMax` n = 1M fork ran at 4.99 (0.898 cycles/elem). In this harness, the facade path and the general path are as fast as the kernel.

**Job 2, scoreboard facade (`ReductionBreezeJmh.galeMax`), slow forks (≤ 4.6 G elem/s) per 8:**

| n | base 96c95e5 r1 / r2 | new 75f9191 r1 / r2 | perfnorm forks base / new |
|---:|---|---|---|
| 65536 | 1 / 0 | **7 / 6** | 0 / 2 |
| 1048576 | 2 / 2 | **0 / 0** | 4 / 0 |

The candidate **removes the slow mode at n = 1M** (0 of 16 slow plain forks, against 4 of 16 for base) but **makes it much more frequent at n = 64K** (13 of 16 slow, against 1 of 16 for base). The slow mode itself is the same in both places: about 4.42 G elem/s against about 5.83 fast. Breeze's `breezeMax` also has occasional slow forks (about 4.4–5.1 G elem/s, 0–1 per 8 per run) in both bundles.

**What the slow mode is (perfnorm, per fork):** it is **not clock state**. If the core clock were lower, cycles/elem would stay the same; instead, slow forks spend more cycles per element (1.014–1.017 against 0.760–0.775). It is also **not extra instructions**. Slow and fast forks execute the same instructions/elem (base galeMax 1M: 4.00–4.01 in all forks; new galeMax 64K: 3.75 in all forks) and the same branches/elem (1.25), but IPC drops from about 5.2 to 3.95 (base) or from about 4.9 to 3.7 (new). The same instruction stream runs about 32% slower. That points to a layout or alignment effect in the JIT-compiled loop (for example, loop or code alignment, or a different scheduling of the same instructions), not a different algorithmic code shape. Distinguishing these further needs perfasm with hsdis, which is not available on the cluster.

The candidate does reduce the facade's instructions per element: 4.00 at base (n = 1M) against 3.75 at new, which matches the kernel's 3.75. This is consistent with passing a constant offset and stride.

Raw data: `2026-10-11-hc-max2-facade-trillium.json` (all JMH JSON, including per-fork rawData and the per-invocation perfnorm results).
