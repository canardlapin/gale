# max loop placement on trillium (Zen 5, lane A)

Does the bimodal `galeMax` slow mode follow C2 loop placement? Source: main `88dde85`, using the dot-rerun bundle already on the cluster. Slurm job 2540487, one exclusive compute node (`tri0119`, 2 × AMD EPYC 9655, Zen 5, AVX-512).

## Protocol

- **Runner:** `docs/verification/w21-simd-spike/x86-avx512/mp.sh`, a copy of the reviewed `run.sh` setup. Lane-A fork flags (`-jvmArgs "-Xms4g -Xmx4g … lane=A …"`, so no Vector module, and `-p backend=pure`), Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, `numactl --physcpubind=0-7 --membind=0`, `-foe true`.
- **Benchmarks:** `ReductionBreezeJmh.(galeMax|breezeMax)$ -p n=65536,1048576 -f 8 -wi 5 -i 15 -w 500ms -r 500ms`.
- **Variants** (passed with `-jvmArgsAppend`):
  - **A:** JDK defaults on this x86 build (OptoLoopAlignment = 16, CodeEntryAlignment = 32).
  - **B:** `-XX:OptoLoopAlignment=64`, as requested. JDK 25 refuses to start with that flag alone ("OptoLoopAlignment (64) must be less or equal to CodeEntryAlignment (32)"), so B also sets `-XX:+UnlockExperimentalVMOptions -XX:CodeEntryAlignment=64`. B therefore moves method entries as well as loop heads.
  - **C:** `-XX:OptoLoopAlignment=32` alone, the largest loop alignment possible without changing entry alignment.
- **Order:** A, B, C interleaved for 2 rounds, then per-fork perfnorm for B as 8 separate `-f 1` invocations (JMH aggregates perfnorm across forks within one invocation).
- **Slow fork:** a per-fork mean of ≤ 4.6 G elem/s (G elem/s = ops/s × n / 1e9, the fork's mean over 15 iterations).

## Slow forks per 8

| n | side | A r1 | A r2 | B r1 | B r2 | C r1 | C r2 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 65536 | galeMax | 0/8 | 0/8 | 0/8 | 0/8 | 1/8 | 0/8 |
| 65536 | breezeMax | 0/8 | 1/8 | 2/8 | 1/8 | 1/8 | 0/8 |
| 1048576 | galeMax | 3/8 | 3/8 | 0/8 | 0/8 | 0/8 | 0/8 |
| 1048576 | breezeMax | 0/8 | 0/8 | 0/8 | 0/8 | 0/8 | 0/8 |

## Per-fork G elem/s

| run | method | n | forks 1…8 |
|---|---|---:|---|
| A-r1 | galeMax | 65536 | 5.83, 5.82, 5.84, 5.81, 5.82, 5.82, 5.83, 5.82 |
| A-r1 | galeMax | 1048576 | 5.85, 4.43, 4.43, 5.79, 4.43, 5.85, 5.85, 5.52 |
| A-r1 | breezeMax | 65536 | 6.06, 6.11, 6.00, 6.00, 6.05, 6.00, 6.05, 6.01 |
| A-r1 | breezeMax | 1048576 | 5.88, 5.79, 6.11, 5.94, 6.13, 5.59, 5.90, 5.92 |
| B-r1 | galeMax | 65536 | 5.84, 5.83, 5.83, 5.84, 5.83, 5.83, 5.82, 5.83 |
| B-r1 | galeMax | 1048576 | 5.85, 5.86, 5.86, 5.86, 5.87, 5.86, 5.85, 6.54 |
| B-r1 | breezeMax | 65536 | 4.95, 5.00, 4.57, 4.90, 4.52, 4.93, 4.93, 4.94 |
| B-r1 | breezeMax | 1048576 | 5.74, 5.85, 5.85, 5.72, 5.25, 5.71, 5.85, 5.27 |
| C-r1 | galeMax | 65536 | 5.82, 5.76, 4.42, 5.82, 5.83, 5.83, 5.83, 5.81 |
| C-r1 | galeMax | 1048576 | 5.84, 5.84, 5.84, 5.85, 5.84, 5.84, 5.84, 5.84 |
| C-r1 | breezeMax | 65536 | 6.01, 5.09, 4.45, 5.93, 5.58, 5.99, 6.16, 6.15 |
| C-r1 | breezeMax | 1048576 | 5.15, 5.84, 5.36, 5.32, 6.14, 5.68, 5.90, 5.33 |
| A-r2 | galeMax | 65536 | 5.83, 5.84, 5.83, 5.82, 5.83, 5.83, 5.83, 5.72 |
| A-r2 | galeMax | 1048576 | 5.83, 4.85, 4.70, 5.85, 4.44, 4.43, 5.85, 4.43 |
| A-r2 | breezeMax | 65536 | 6.06, 6.08, 6.15, 4.54, 5.99, 6.15, 6.02, 6.14 |
| A-r2 | breezeMax | 1048576 | 6.15, 6.20, 6.04, 6.17, 6.08, 5.92, 4.89, 6.15 |
| B-r2 | galeMax | 65536 | 5.81, 5.83, 5.87, 5.83, 5.83, 5.83, 5.83, 5.83 |
| B-r2 | galeMax | 1048576 | 5.86, 5.87, 5.87, 5.90, 5.86, 5.85, 5.85, 6.53 |
| B-r2 | breezeMax | 65536 | 4.99, 4.94, 4.58, 4.99, 4.99, 4.97, 4.91, 4.97 |
| B-r2 | breezeMax | 1048576 | 5.40, 5.85, 5.78, 5.80, 5.87, 5.76, 5.39, 5.73 |
| C-r2 | galeMax | 65536 | 5.83, 5.84, 5.82, 5.84, 5.83, 5.84, 5.83, 5.86 |
| C-r2 | galeMax | 1048576 | 5.84, 5.85, 5.84, 5.84, 5.86, 5.85, 5.84, 5.85 |
| C-r2 | breezeMax | 65536 | 6.02, 6.15, 5.96, 5.06, 5.06, 5.84, 5.89, 6.16 |
| C-r2 | breezeMax | 1048576 | 5.88, 6.11, 4.98, 5.54, 5.91, 5.74, 6.23, 6.09 |

## Per-fork perfnorm for B (8 × `-f 1`)

Per fork: G elem/s, cycles/elem, instructions/elem, IPC. Values above 100 cycles/elem are perf collection glitches (as in the facade job) and are shown as n/a.

| method | n | per fork |
|---|---:|---|
| breezeMax | 65536 | 4.92 @ 0.912 (3.51, IPC 3.84); 5.02 @ 0.893 (3.51, IPC 3.92); 5.03 @ 0.892 (3.50, IPC 3.93); 4.55 @ 0.987 (3.50, IPC 3.55); 5.04 @ 0.891 (3.50, IPC 3.93); 5.01 @ 0.895 (3.51, IPC 3.92); 4.98 @ 0.901 (3.51, IPC 3.89); 4.48 @ 1.003 (3.51, IPC 3.50) |
| breezeMax | 1048576 | 5.29 @ 0.847 (3.50, IPC 4.13); 5.84 @ 0.767 (3.50, IPC 4.57); 5.71 @ 0.784 (3.51, IPC 4.48); 5.84 @ 0.768 (3.50, IPC 4.56); 5.86 @ n/a; 5.80 @ 0.773 (3.50, IPC 4.53); 5.77 @ 0.777 (3.51, IPC 4.51); 5.18 @ 0.864 (3.50, IPC 4.05) |
| galeMax | 65536 | 5.79 @ 0.775 (4.01, IPC 5.17); 5.81 @ 0.771 (4.01, IPC 5.19); 5.81 @ 0.772 (4.00, IPC 5.19); 5.82 @ 0.771 (4.00, IPC 5.19); 5.77 @ 0.778 (4.01, IPC 5.15); 5.82 @ 0.771 (4.01, IPC 5.20); 5.81 @ 0.772 (4.01, IPC 5.19); 5.82 @ 0.771 (4.00, IPC 5.19) |
| galeMax | 1048576 | 5.84 @ 0.767 (4.01, IPC 5.22); 6.51 @ n/a; 5.83 @ 0.768 (4.00, IPC 5.21); 5.94 @ 0.754 (4.01, IPC 5.31); 5.86 @ 0.765 (4.01, IPC 5.24); 6.52 @ 0.687 (4.00, IPC 5.82); 6.51 @ n/a; 5.83 @ 0.768 (4.00, IPC 5.21) |

| method | n | fast forks: IPC (count) | slow forks: IPC (count) |
|---|---:|---|---|
| galeMax | 65536 | 5.15–5.20 (8) | — (0) |
| galeMax | 1048576 | 5.21–5.82 (6) | — (0) |
| breezeMax | 65536 | 3.84–3.93 (6) | 3.50–3.55 (2) |
| breezeMax | 1048576 | 4.05–4.57 (7) | — (0) |

## Summary

**Loop placement controls `galeMax`'s slow mode at n = 1M.**

| `galeMax` slow forks (≤ 4.6 G elem/s) | A (default) | B (loop + entry alignment 64) | C (loop alignment 32) |
|---|---|---|---|
| n = 1048576 | **6 / 16** (3 + 3), plus 3 intermediate forks at 4.70–5.52 | **0 / 16** | **0 / 16** |
| n = 65536 | 0 / 16 | 0 / 16 | 1 / 16 (one fork at 4.42 in C r1) |
| B perfnorm forks (8 × `-f 1`) | — | 0 / 8 at both n | — |

- **With the JDK defaults (A):** the slow mode at n = 1M reproduced in both rounds, at about 4.43 G elem/s against about 5.85.
- **With 64-byte alignment (B):** every one of 16 + 16 + 8 forks was fast, at both sizes. Two B forks at n = 1M ran even faster (6.53–6.54 G elem/s).
- **With 32-byte loop alignment only (C):** the n = 1M slow mode is also gone (0/16), but one n = 64K fork was still slow.
- **Conclusion:** the slow mode follows the code placement of C2's compiled loop, not the input or the clock.

**Breeze (`breezeMax`) is not fixed by alignment.** Slow-fork counts per 8 (r1 / r2) at n = 64K are 0 / 1 in A, 2 / 1 in B and 1 / 0 in C; at n = 1M, 0 in every variant. Breeze's occasional slow forks are a different effect.

**IPC:**
- **B, fast forks (from this job):** `galeMax` 5.15–5.20 (64K) and 5.21–5.82 (1M); `breezeMax` 3.84–3.93 (64K) and 4.05–4.57 (1M).
- **B, Breeze slow forks:** two Breeze forks at 64K ran at IPC 3.50–3.55. No gale fork in B was slow, so B gives no gale slow-mode IPC.
- **Default slow mode, from the facade job (`2026-10-11-hc-max2-facade-trillium.md`, same harness, base `96c95e5`):** slow `galeMax` forks run the same instructions/elem (4.00) at IPC 3.95, against 5.2 for fast forks, at 1.01 against 0.77 cycles/elem.
- **Reading:** the slow mode is an IPC collapse on an identical instruction stream, and moving loop heads to 32- or 64-byte boundaries removes it. This is consistent with a loop that straddles a fetch or uop-cache boundary at the default 16-byte alignment.

Raw data: `2026-10-11-max-placement-trillium.json` (all JMH JSON with per-fork rawData, and the per-invocation perfnorm results).
