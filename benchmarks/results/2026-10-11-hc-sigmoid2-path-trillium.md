# hc-sigmoid2: sigmoid path attribution and x86_64 digests on trillium (Zen 5, lane A)

`4cbec0a` (breeze/hc-sigmoid2). Its main code is identical to `96c95e5`: the commits since then add only `SigmoidPathJmh` and `SigmoidDigestSuite`. Slurm job 2539162, one exclusive compute node (`tri0122`, 2 × AMD EPYC 9655, Zen 5, AVX-512), 2026-10-11. 1-minute load at start: 0.91.

## Protocol

- **Build:** `4cbec0a`, compiled locally in a clean worktree (local branch `bench/sg2-4cbec0a`) with sbt 1.11.7, JDK 25.0.1 and Scala 3.7.4. The exported JMH and `coreJVM` test classpaths ran on the cluster; no sbt ran there.
- **Runner:** `docs/verification/w21-simd-spike/x86-avx512/sg2.sh`, a copy of the reviewed `run.sh` setup. It uses Temurin 25+36 (`module load java/25`), unsets `JAVA_TOOL_OPTIONS`, pins with `numactl --physcpubind=0-7 --membind=0`, and binds `-Dgale.bench.sourceSha` from the bundle.
  - `pathA`: `gale.bench.SigmoidPathJmh.*` with the lane-A fork arguments (`-jvmArgs "-Xms4g -Xmx4g … -Dgale.bench.lane=A …"`, so no Vector module) and `-p backend=pure`, using the annotation defaults: 2 forks, 5 × 500 ms warmup, 5 × 500 ms measurement.
  - `digest`: `SigmoidDigestSuite` through JUnitCore on the `coreJVM` test classpath, unforked and without `--add-modules`, as sbt runs core tests. It ran with `GALE_SIGMOID_DIGEST_OUT=<out>/digests-amd64.txt`, then again in compare mode with the stored resource plus the generated lines overlaid ahead of the test classes.
- **Units:** JMH reports throughput per call; this receipt converts it to **ns per element** = 1e9 / (ops/s × n). CIs are the converted 99.9% JMH confidence bounds (lower is better).

## ns/elem, mean [99.9% CI]

| method | n = 1024 | n = 65536 | n = 1048576 |
|---|---:|---:|---:|
| breeze | 3.352 [3.338, 3.366] | 3.556 [3.541, 3.571] | 5.103 [5.093, 5.113] |
| gale | 3.500 [3.466, 3.536] | 3.694 [3.616, 3.775] | 3.591 [3.576, 3.606] |
| rawHybrid | 3.488 [3.448, 3.529] | 3.679 [3.625, 3.734] | 3.595 [3.588, 3.603] |
| rawDirect | 3.335 [3.308, 3.362] | 3.541 [3.519, 3.565] | 3.528 [3.518, 3.538] |
| rawHybrid4 | 3.719 [3.671, 3.768] | 3.773 [3.750, 3.795] | 3.751 [3.739, 3.764] |
| scanMin | 0.313 [0.302, 0.325] | 0.266 [0.264, 0.267] | 0.265 [0.261, 0.270] |
| scanMin4 | 0.480 [0.479, 0.481] | 0.400 [0.371, 0.435] | 0.373 [0.314, 0.459] |
| scanBits | 0.230 [0.230, 0.230] | 0.229 [0.229, 0.229] | 0.234 [0.234, 0.234] |
| prescanMin4 | 3.746 [3.671, 3.824] | 3.927 [3.889, 3.967] | 3.933 [3.899, 3.968] |
| prescanBits | 3.611 [3.552, 3.673] | 3.788 [3.775, 3.800] | 3.761 [3.751, 3.770] |

## Comparisons (relative ns/elem; negative means faster; `*` marks disjoint CIs)

| comparison | n = 1024 | n = 65536 | n = 1048576 |
|---|---:|---:|---:|
| rawDirect vs rawHybrid | -4.4%* | -3.7%* | -1.9%* |
| rawHybrid4 vs rawHybrid | +6.6%* | +2.6%* | +4.3%* |
| prescanMin4 vs rawHybrid | +7.4%* | +6.8%* | +9.4%* |
| prescanBits vs rawHybrid | +3.5%* | +3.0%* | +4.6%* |
| gale vs rawHybrid | +0.4% | +0.4% | -0.1% |
| gale vs breeze | +4.4%* | +3.9%* | -29.6%* |
| rawDirect vs breeze | -0.5% | -0.4% | -30.9%* |

**Reading:**
- **`rawDirect` vs `rawHybrid`:** the direct form `1/(1+exp(-x))` is 4.4% (n = 1024), 3.7% (64K) and 1.9% (1M) faster than the hybrid, with disjoint CIs. That is the branch cost of the `x < -36` test on this data.
- **`gale` vs `rawHybrid`:** the facade adds no measurable overhead on top of the raw hybrid loop.
- **Scans:** the min scans cost 0.27–0.31 ns/elem (`scanMin`) and 0.37–0.48 ns/elem (`scanMin4`); the bit scan costs 0.23 ns/elem (`scanBits`). On this CPU `scanMin4` is slower than `scanMin`, and `scanBits` is the cheapest.
- **Prescan variants:** `prescanMin4` (+6.8–9.4%) and `prescanBits` (+3.0–4.6%) are both slower than `rawHybrid`, with disjoint CIs at every n. A scan-then-direct strategy does not pay for itself on Zen 5 with this data, because the scan costs more than the branch it removes.
- **`rawHybrid4`** (the 4-way unrolled hybrid) is 2.6–6.6% slower than `rawHybrid`.
- **Against Breeze:** gale is 3.9–4.4% slower than Breeze at n = 1024 and 65536 (the lane-A gap). They are 30% faster at n = 1M, where Breeze slows to 5.1 ns/elem. `rawDirect` is within 0.5% of Breeze at 1024 and 64K.

## x86_64 (`amd64`) sigmoid digests

The first `SigmoidDigestSuite` run, with `GALE_SIGMOID_DIGEST_OUT` set, passed 3 tests and wrote 114 `amd64 …` lines. Those lines are appended, unchanged, to `core/jvm/src/test/resources/gale/linalg/sigmoid-digests.txt` (sorted, after the aarch64 block). The compare-mode rerun on the same node, with the merged resource overlaid on the classpath, passed 3 tests: the frozen-reference check, the JIT-tier check, and the per-family digest check for `os.arch=amd64`.

A local negative control showed the overlay is really read in compare mode: corrupting one stored aarch64 digest makes the suite fail with `1 cases changed bits`.

65 of the 114 amd64 digests differ from their aarch64 counterparts. The 24 `input` digests are identical, as expected: platform `Math.exp` differs while the corpus does not. The cross-platform bit check against the frozen 96c95e5 reference passed on both platforms.

Raw data: `2026-10-11-hc-sigmoid2-path-trillium.json` (JMH JSON, netlib records if any, the digest-step log and the generated digest lines).
