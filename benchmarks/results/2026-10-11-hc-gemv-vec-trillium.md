# HC-4 vector gemv on trillium (Zen 5, lane B)

breeze/hc-gemv-vec: an eight-row gemv tile, masked tail and row-sum tree (`7c460eb`), plus the paired `GemvHillClimbJmh` harness (`7a79d4a`). Slurm job 2539802, one exclusive compute node (`tri0745`, 2 × AMD EPYC 9655, Zen 5, AVX-512, 8 double lanes).

## Protocol

- **Bundles:** base `90c74a3` (main, with the Milestone 3 receipts) and new `7a79d4a`, both compiled locally from one clean bench worktree (branch `bench/gv4-7a79d4a`) with sbt 1.11.7, JDK 25.0.1 and Scala 3.7.4. No sbt ran on the cluster.
- **Harness:** `docs/verification/w21-simd-spike/x86-avx512/gv4.sh`, a copy of the reviewed `run.sh` setup: Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, `numactl --physcpubind=0-7 --membind=0`, host `--add-modules=jdk.incubator.vector` inherited by the forks (lane B), `-jvmArgsAppend "-Xms4g -Xmx4g -Dgale.bench.sourceSha=<bundle SHA> -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=…"`, and `-foe true`.
  - `hill`: `GemvHillClimbJmh`, `-f 3 -wi 5 -i 5`, new bundle only. Its setup throws unless VectorBLAS resolves, and every run did resolve it.
  - `pairB`: `BlasL2BreezeJmh.(galeGemvT|breezeGemvT|galeGemv|breezeGemv)$` and `LbfgsBreezeJmh.*[Ll]ogistic`, with `-p backend=vector -p budget=fixed` and the annotation defaults (2 forks).
- **Order:** hill-new, pair-base-r1, pair-new-r1, pair-base-r2, pair-new-r2.
- **Netlib** (blas, lapack, vectorModule): [('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- **Ratios:** speed ratios; >1 means the first-named variant is faster. Brackets give conservative bounds from the two 99.9% CIs, and `*` marks disjoint CIs.

## Part 1: GemvHillClimbJmh (new bundle)

| shape | curGemv | refGemv | noTreeGemv | vblasGemv | curGemvT | vblasGemvT |
|---|---:|---:|---:|---:|---:|---:|
| 256x256 | 2.534e+05 ± 9.4e+02 ops/s | 2.348e+05 ± 3.9e+03 ops/s | 2.486e+05 ± 6.7e+02 ops/s | 2.343e+05 ± 1.2e+04 ops/s | 2.455e+05 ± 1.1e+03 ops/s | 2.542e+05 ± 3.8e+03 ops/s |
| 1024x1024 | 1.347e+04 ± 1.7e+02 ops/s | 1.26e+04 ± 16 ops/s | 1.324e+04 ± 39 ops/s | 1.252e+04 ± 3.2e+02 ops/s | 1.136e+04 ± 5.6e+02 ops/s | 1.331e+04 ± 1.2e+02 ops/s |
| 2048x2048 | 3,111 ± 11 ops/s | 2,705 ± 23 ops/s | 3,091 ± 13 ops/s | 2,845 ± 27 ops/s | 2,888 ± 48 ops/s | 2,937 ± 44 ops/s |
| 1000x50 | 2.94e+05 ± 1.1e+03 ops/s | 2.287e+05 ± 4.3e+03 ops/s | 2.671e+05 ± 6.8e+03 ops/s | 3.217e+05 ± 7.1e+02 ops/s | 2.626e+05 ± 2.3e+03 ops/s | 3.247e+05 ± 2e+03 ops/s |

| shape | curGemv / refGemv | curGemv / noTreeGemv | curGemv / vblasGemv | curGemvT / vblasGemvT |
|---|---|---|---|---|
| 256x256 | 1.08× [1.06, 1.10]* | 1.02× [1.01, 1.03]* | 1.08× [1.02, 1.15]* | 0.97× [0.95, 0.98]* |
| 1024x1024 | 1.07× [1.06, 1.08]* | 1.02× [1.00, 1.03]* | 1.08× [1.04, 1.12]* | 0.85× [0.80, 0.90]* |
| 2048x2048 | 1.15× [1.14, 1.16]* | 1.01× [1.00, 1.01] | 1.09× [1.08, 1.11]* | 0.98× [0.95, 1.02] |
| 1000x50 | 1.29× [1.26, 1.32]* | 1.10× [1.07, 1.13]* | 0.91× [0.91, 0.92]* | 0.81× [0.80, 0.82]* |

## Part 2: paired base vs new (lane B, `backend=vector`, `budget=fixed`)

| class | op | params | base r1 gale | new r1 gale | gale/Breeze base r1 / r2 | gale/Breeze new r1 / r2 | gale new/base r1 / r2 |
|---|---|---|---:|---:|---|---|---|
| BlasL2BreezeJmh | Gemv | backend=vector, n=1024 | 1.262e+04 ± 23 ops/s | 1.325e+04 ± 2.6e+02 ops/s | 0.98× [0.93, 1.04] / 0.99× [0.95, 1.03] | 1.03× [0.97, 1.10] / 1.04× [0.99, 1.10] | 1.05× [1.03, 1.07]* / 1.06× [1.05, 1.07]* |
| BlasL2BreezeJmh | Gemv | backend=vector, n=2048 | 2,669 ± 17 ops/s | 3,050 ± 30 ops/s | 0.97× [0.96, 0.98]* / 0.95× [0.93, 0.97]* | 1.08× [1.07, 1.09]* / 1.10× [1.08, 1.12]* | 1.14× [1.12, 1.16]* / 1.15× [1.14, 1.16]* |
| BlasL2BreezeJmh | Gemv | backend=vector, n=256 | 2.255e+05 ± 3.4e+03 ops/s | 2.468e+05 ± 4.9e+03 ops/s | 0.96× [0.95, 0.98]* / 0.97× [0.95, 0.99]* | 0.99× [0.89, 1.12] / 1.00× [0.89, 1.12] | 1.09× [1.06, 1.13]* / 1.09× [1.06, 1.13]* |
| BlasL2BreezeJmh | GemvT | backend=vector, n=1024 | 1.225e+04 ± 28 ops/s | 1.252e+04 ± 2.6e+02 ops/s | 0.91× [0.88, 0.93]* / 0.95× [0.90, 1.00] | 0.95× [0.93, 0.97]* / 0.94× [0.91, 0.97]* | 1.02× [1.00, 1.05] / 0.99× [0.93, 1.06] |
| BlasL2BreezeJmh | GemvT | backend=vector, n=2048 | 2,931 ± 6 ops/s | 2,918 ± 13 ops/s | 1.01× [0.98, 1.04] / 0.99× [0.98, 1.01] | 1.02× [1.01, 1.03]* / 1.00× [0.98, 1.03] | 1.00× [0.99, 1.00] / 0.99× [0.99, 1.00] |
| BlasL2BreezeJmh | GemvT | backend=vector, n=256 | 2.411e+05 ± 1.2e+03 ops/s | 2.404e+05 ± 2.3e+03 ops/s | 0.98× [0.96, 1.00]* / 0.99× [0.96, 1.01] | 0.94× [0.86, 1.02] / 0.97× [0.94, 1.00]* | 1.00× [0.98, 1.01] / 0.99× [0.98, 1.01] |
| LbfgsBreezeJmh | Logistic | backend=vector, budget=fixed | 587.6 ± 6.2 us/op | 561.2 ± 7.8 us/op | 0.95× [0.94, 0.96]* / 0.95× [0.93, 0.98]* | 0.99× [0.97, 1.01] / 1.00× [0.99, 1.01] | 1.05× [1.02, 1.07]* / 1.05× [1.02, 1.08]* |

## Scoreboard (per pair run, `breeze_scoreboard.py --lane B --strict --validation <empty>`)

None of these classes is preflight-covered, so no validation records exist; each run was scored with an empty validation sidecar and exited 0 with 0 unpaired rows. Only `vector` gale rows were run, so the headline covers every pair.

- `pair-base-r1`: Headline, gale vector vs Breeze VectorBLAS: 0 ahead, 2 tie, 5 behind, 0 n/a, 0 withheld. All rows: 0 ahead, 2 tie, 5 behind, 0 n/a, 0 withheld of 7 pairs.
- `pair-base-r2`: Headline, gale vector vs Breeze VectorBLAS: 0 ahead, 4 tie, 3 behind, 0 n/a, 0 withheld. All rows: 0 ahead, 4 tie, 3 behind, 0 n/a, 0 withheld of 7 pairs.
- `pair-new-r1`: Headline, gale vector vs Breeze VectorBLAS: 2 ahead, 4 tie, 1 behind, 0 n/a, 0 withheld. All rows: 2 ahead, 4 tie, 1 behind, 0 n/a, 0 withheld of 7 pairs.
- `pair-new-r2`: Headline, gale vector vs Breeze VectorBLAS: 1 ahead, 4 tie, 2 behind, 0 n/a, 0 withheld. All rows: 1 ahead, 4 tie, 2 behind, 0 n/a, 0 withheld of 7 pairs.

## Summary

**Part 1, kernel (new bundle, 3 forks):**
- **The row-sum tree and masked tail are intrinsified on x86.** `curGemv` is not slower than `noTreeGemv` at any shape: 1.02× at 256², 1.02× at 1024², 1.01× (CIs overlap, a tie) at 2048², and 1.10× at 1000×50. The rotate tree and masked `fromArray` are not a liability on AVX-512.
- **The HC-4 tile against the pre-HC-4 reference:** `curGemv / refGemv` is 1.07–1.15× on the square shapes and 1.29× on 1000×50, all with disjoint CIs.
- **Against VectorBLAS, gemv:** ahead on the square shapes (`curGemv / vblasGemv` 1.08–1.09×, disjoint) but behind on the short-wide 1000×50 (0.91×, disjoint).
- **Against VectorBLAS, gemvT:** behind or tied (0.97× at 256², 0.85× at 1024², 0.98× tie at 2048², 0.81× at 1000×50). HC-4 did not target gemvT.

**Part 2, receipt pairs:**
- **Gale gemv (vector) new vs base:** faster at every n with disjoint CIs in both rounds: 1.09× at n = 256, 1.05–1.06× at 1024, and 1.14–1.15× at 2048. Against Breeze VectorBLAS it moves from behind (0.95–0.97× at 256 and 2048) to a tie at 256 (0.99–1.00×) and 1024 (1.03–1.04×, CIs overlap), and to ahead at 2048 (1.08–1.10×).
- **gemvT (vector):** unchanged vs base (0.99–1.02×, CIs overlap). n = 1024 stays behind Breeze (0.94–0.95×), and n = 256 / 2048 are ties or near-ties.
- **L-BFGS logistic, fixed budget (vector):** 1.05× faster than base in both rounds (disjoint). Against Breeze it goes from 0.95× (behind) to 0.99–1.00× (tie).
- **Headline counts:** across the 7 pairs, 0 ahead / 2–4 tie / 3–5 behind in the base runs, against 1–2 ahead / 4 tie / 1–2 behind in the new runs.

Raw data: `2026-10-11-hc-gemv-vec-trillium.json`. It holds the per-run JMH JSON, netlib records, source SHA and load.
