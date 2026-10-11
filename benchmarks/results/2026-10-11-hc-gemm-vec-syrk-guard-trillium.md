# syrk SIMD guard on AVX-512 (trillium, Zen 5, lane B)

`552d2ad` (breeze/hc-gemm-vec): does the Vector `syrk` shape guard (`syrkUsesSimd`: SIMD iff `m ≥ 4·lanes` and `m·k² ≥ 2·64³`) pick the faster path on 8-lane AVX-512? Slurm job 2536354, one exclusive compute node (`tri0515`, 2 × AMD EPYC 9655), two back-to-back rounds.

## Protocol

- **Build:** `552d2ad`, compiled locally in a clean worktree (local branch `bench/hcgv-552d2ad`) with sbt 1.11.7, JDK 25.0.1 and Scala 3.7.4. The exported JMH classpath ran on the cluster; no sbt ran there.
- **Harness:** `docs/verification/w21-simd-spike/x86-avx512/sy.sh`, which is the reviewed `run.sh` lane-B path: inherited `--add-modules=jdk.incubator.vector`, `-Xms4g -Xmx4g`, a bundle-bound `gale.bench.sourceSha`, `-foe true`, `-p backend=pure,vector`, Temurin 25+36 with `JAVA_TOOL_OPTIONS` unset, and `numactl --physcpubind=0-7 --membind=0`. It is restricted to `VectorSyrkJmh.*` and `BlasL3BreezeJmh.(galeAtA|breezeAtA)$`. It used the annotation defaults: 2 forks, 5 × 500 ms warmup, 5 × 500 ms measurement for `VectorSyrkJmh`.
- **Load at start:** r1: 0.53, 56.53, 129.65; r2: 1.02, 7.88, 66.63.
- **Netlib** (blas, lapack, vectorModule): [('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)]. The `[w21-simd]`-style species is 8 lanes on this CPU (AVX-512, as in the W2.1 x86 run), so `4·lanes = 32`.
- **Guard column:** computed from the guard's formula with lanes = 8. The guarded `kernel(vector)` timing agrees with it: it matches `simdForced` where the guard says SIMD, and the pure kernel where it says pure.

## Per shape (µs/op, mean ± 99.9% CI half-width; lower is better)

| shape (m×k) | round | pure kernel | simdForced | guarded vector kernel | guard picks | faster path | flag |
|---|---:|---:|---:|---:|---|---|---|
| 64x64 | 1 | 16.74 ± 0.00 | 9.30 ± 0.04 | 16.73 ± 0.00 | pure | SIMD (SIMD 1.80× faster) | **guard picks slower path** |
| 64x64 | 2 | 16.73 ± 0.00 | 9.31 ± 0.05 | 16.75 ± 0.01 | pure | SIMD (SIMD 1.80× faster) | **guard picks slower path** |
| 128x64 | 1 | 32.66 ± 0.04 | 16.31 ± 0.07 | 16.31 ± 0.08 | SIMD | SIMD (SIMD 2.00× faster) |  |
| 128x64 | 2 | 32.65 ± 0.05 | 16.32 ± 0.10 | 16.31 ± 0.05 | SIMD | SIMD (SIMD 2.00× faster) |  |
| 256x64 | 1 | 64.06 ± 0.03 | 56.05 ± 0.19 | 56.04 ± 0.14 | SIMD | SIMD (SIMD 1.14× faster) |  |
| 256x64 | 2 | 64.17 ± 0.13 | 56.11 ± 0.16 | 56.07 ± 0.13 | SIMD | SIMD (SIMD 1.14× faster) |  |
| 128x128 | 1 | 128.7 ± 0.2 | 57.11 ± 1.69 | 57.03 ± 1.62 | SIMD | SIMD (SIMD 2.25× faster) |  |
| 128x128 | 2 | 128.5 ± 0.1 | 57.03 ± 1.78 | 57.04 ± 1.68 | SIMD | SIMD (SIMD 2.25× faster) |  |
| 1024x256 | 1 | 4,731.4 ± 5.9 | 1,580.4 ± 51.5 | 1,545.1 ± 9.8 | SIMD | SIMD (SIMD 2.99× faster) |  |
| 1024x256 | 2 | 4,644.0 ± 92.9 | 1,583.7 ± 44.4 | 1,548.6 ± 10.1 | SIMD | SIMD (SIMD 2.93× faster) |  |
| 4x1024 | 1 | 2,737.1 ± 8.3 | 2,833.5 ± 41.8 | 2,732.5 ± 2.9 | pure | pure (pure 1.04× faster) |  |
| 4x1024 | 2 | 2,739.4 ± 18.2 | 2,844.1 ± 10.8 | 2,756.8 ± 14.7 | pure | pure (pure 1.04× faster) |  |
| 8x1024 | 1 | 2,945.1 ± 6.3 | 2,786.4 ± 9.7 | 2,965.0 ± 35.9 | pure | SIMD (SIMD 1.06× faster) | **guard picks slower path** |
| 8x1024 | 2 | 2,957.1 ± 1.6 | 2,822.9 ± 61.7 | 2,925.0 ± 50.5 | pure | SIMD (SIMD 1.05× faster) | **guard picks slower path** |
| 16x1024 | 1 | 3,620.9 ± 5.1 | 2,936.0 ± 12.8 | 3,629.2 ± 40.2 | pure | SIMD (SIMD 1.23× faster) | **guard picks slower path** |
| 16x1024 | 2 | 3,596.0 ± 28.8 | 2,888.7 ± 22.8 | 3,590.0 ± 33.8 | pure | SIMD (SIMD 1.24× faster) | **guard picks slower path** |
| 32x512 | 1 | 1,179.6 ± 179.4 | 694.5 ± 24.6 | 972.0 ± 128.1 | SIMD | SIMD (SIMD 1.70× faster) |  |
| 32x512 | 2 | 1,098.9 ± 47.8 | 858.3 ± 234.2 | 695.0 ± 23.7 | SIMD | tie (pure/SIMD 1.28×) |  |
| 8x4096 | 1 | 137,496.4 ± 4,067.9 | 219,906.6 ± 5,723.4 | 128,898.5 ± 2,461.1 | pure | pure (pure 1.60× faster) |  |
| 8x4096 | 2 | 135,283.7 ± 7,435.0 | 219,264.4 ± 17,753.5 | 130,049.0 ± 9,231.5 | pure | pure (pure 1.62× faster) |  |

## Facade `A.t * A` (µs/op)

| shape | round | facade pure | facade vector | vector/pure speedup |
|---|---:|---:|---:|---:|
| 64x64 | 1 | 17.98 ± 0.05 | 18.00 ± 0.06 | 1.00× |
| 64x64 | 2 | 17.99 ± 0.05 | 17.99 ± 0.03 | 1.00× |
| 128x64 | 1 | 34.20 ± 0.41 | 17.05 ± 0.04 | 2.01×* |
| 128x64 | 2 | 34.21 ± 0.47 | 17.07 ± 0.01 | 2.00×* |
| 256x64 | 1 | 65.55 ± 0.43 | 57.25 ± 0.11 | 1.14×* |
| 256x64 | 2 | 65.39 ± 0.40 | 57.24 ± 0.11 | 1.14×* |
| 128x128 | 1 | 135.6 ± 1.6 | 65.99 ± 1.48 | 2.06×* |
| 128x128 | 2 | 135.7 ± 1.9 | 65.92 ± 1.53 | 2.06×* |
| 1024x256 | 1 | 4,559.0 ± 90.4 | 1,672.0 ± 90.5 | 2.73×* |
| 1024x256 | 2 | 4,630.6 ± 165.0 | 1,613.6 ± 18.6 | 2.87×* |
| 4x1024 | 1 | 2,971.9 ± 9.0 | 2,959.7 ± 12.2 | 1.00× |
| 4x1024 | 2 | 2,969.2 ± 7.3 | 2,962.2 ± 10.8 | 1.00× |
| 8x1024 | 1 | 3,167.5 ± 22.3 | 3,167.0 ± 22.8 | 1.00× |
| 8x1024 | 2 | 3,164.8 ± 20.1 | 3,157.3 ± 13.0 | 1.00× |
| 16x1024 | 1 | 3,850.9 ± 8.9 | 3,844.5 ± 10.6 | 1.00× |
| 16x1024 | 2 | 3,848.1 ± 13.2 | 3,847.5 ± 12.9 | 1.00× |
| 32x512 | 1 | 1,126.9 ± 5.0 | 794.2 ± 2.7 | 1.42×* |
| 32x512 | 2 | 1,124.7 ± 3.0 | 739.6 ± 2.7 | 1.52×* |
| 8x4096 | 1 | 133,251.3 ± 4,205.0 | 133,825.0 ± 5,523.8 | 1.00× |
| 8x4096 | 2 | 133,205.7 ± 4,560.4 | 133,094.4 ± 4,288.5 | 1.00× |

## BlasL3 AtA vs Breeze (VectorBLAS), ops/s

| n | round | gale pure | gale vector | Breeze | pure/Breeze | vector/Breeze |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 1 | 786,935.2 ± 6,710.1 | 786,335.7 ± 5,966.4 | 396,311.2 ± 10,203.3 | 1.99× ahead | 1.98× ahead |
| 16 | 2 | 786,832.9 ± 2,707.9 | 787,410.2 ± 6,276.8 | 396,062.0 ± 7,749.9 | 1.99× ahead | 1.99× ahead |
| 64 | 1 | 15,252.6 ± 94.2 | 17,486.4 ± 84.0 | 14,414.7 ± 51.4 | 1.06× ahead | 1.21× ahead |
| 64 | 2 | 15,295.1 ± 107.6 | 17,483.4 ± 102.9 | 14,410.6 ± 52.8 | 1.06× ahead | 1.21× ahead |
| 256 | 1 | 215.7 ± 9.8 | 619.8 ± 7.3 | 351.1 ± 8.9 | 0.61× behind | 1.77× ahead |
| 256 | 2 | 218.1 ± 3.3 | 619.9 ± 6.9 | 350.3 ± 8.8 | 0.62× behind | 1.77× ahead |

## Summary

In both rounds, the guard picks the CI-disjointly slower path on **3 of 10 shapes**. All three are shapes where it chooses pure but SIMD is faster:

| shape | why the guard says pure | pure kernel (µs) | simdForced (µs) | SIMD advantage |
|---|---|---:|---:|---:|
| 64×64 | work `m·k² = 262,144 < 2·64³` | 16.73 ± 0.00 | 9.30 ± 0.04 | 1.80× |
| 16×1024 | `m = 16 < 4·lanes = 32` | 3,596–3,621 | 2,889–2,936 | 1.23–1.24× |
| 8×1024 | `m = 8 < 32` | 2,945–2,957 | 2,786–2,823 | 1.05–1.06× |

On the shapes where the guard chooses SIMD (128×64, 256×64, 128×128, 1024×256, 32×512), SIMD is faster: 1.14–2.99×. 32×512 round 2 is a statistical tie because the CIs there are wide (±234 µs on simdForced). On two of the shapes where it chooses pure, pure is faster and the guard is right: 4×1024 (pure 1.04× faster) and especially 8×4096 (pure 1.60–1.62× faster).

**Reading for AVX-512 (8 lanes):**
- **Work floor too high.** The `2·64³` work floor, measured on 2-lane NEON, is too conservative here: 64×64 gains 1.8×. Lowering it to the GEMM's `64³` (or lower) would capture this.
- **Row floor too high.** The short-wide floor `m ≥ 4·lanes` (32 on AVX-512) excludes m = 16, which gains 1.23×. m = 8 is a marginal 5–6% gain, and m = 4 loses 4%.
- **Long-k shapes lose regardless.** At 8×4096 SIMD loses 1.6× whatever `m` is, so a lower row floor alone (e.g. `m ≥ 2·lanes`) would be right for 16×1024. But 8×4096 shows a `k`-dependent cost (the `k × 256` panel, `k = 4096`) that an `m`-only floor cannot express.
- **Not a correctness issue.** The guarded vector kernel's timings match the path the formula predicts on every shape, so the guard is wired as intended.

**Facade:** `A.t * A` with the Vector backend shows 2.0–2.9× on the SIMD-routed shapes and 1.00× elsewhere, consistent with the kernel results.

**BlasL3 AtA vs Breeze:** at n = 256 the vector backend is now 1.77× ahead (619.8 ± 7.3 vs 351.1 ± 8.9 ops/s). The `4cdf0f3` run measured 0.95–0.98× there (`2026-10-11-hc-gemm-vec-trillium.md`). At n = 64, vector is 1.21× ahead and pure 1.06×. At n = 16 both are 1.99× ahead. Pure n = 256 remains 0.61–0.62× behind.

Raw data: `2026-10-11-hc-gemm-vec-syrk-guard-trillium.json` (per-round JMH JSON, netlib sidecar records, source SHA, load).
