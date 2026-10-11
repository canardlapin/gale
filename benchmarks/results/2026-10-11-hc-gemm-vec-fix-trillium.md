# hc-gemm-vec fix check on trillium: 7a90fb7 (via 4bc265d) vs Milestone 2 (Zen 5, lane B)

`4bc265d` is `7a90fb7` (restore the GEMM's own packed tile loop; syrk gets its own `syrkUpperTiles` loop) plus the scratch `TmpGemmBisectJmh` harness, which is not for merge. It is paired against the Milestone 2 baseline bundle (`0ff6440`, source-identical to `5d3950e`). Slurm job 2536395, one exclusive compute node (`tri0299`, 2 × AMD EPYC 9655, Zen 5, AVX-512, 8 double lanes), 2026-10-11 02:43–04:40 EDT.

## Protocol

- **Build:** `4bc265d`, compiled locally in a clean worktree (local branch `bench/hcgv-4bc265d`) with sbt 1.11.7, JDK 25.0.1 and Scala 3.7.4. No sbt ran on the cluster.
- **Harness:** `docs/verification/w21-simd-spike/x86-avx512/bi.sh`, which is the reviewed `run.sh` lane-B path: inherited `--add-modules`, `-Xms4g -Xmx4g`, a bundle-bound `gale.bench.sourceSha`, `-foe true`, a fresh netlib and validation sidecar per run, `-p backend=pure,vector`, the annotation defaults (2 forks), Temurin 25+36 with `JAVA_TOOL_OPTIONS` unset, and `numactl --physcpubind=0-7 --membind=0`. It has three modes:
  - `gvB`: the gemm/gemv/factorization/decomposition/lstsq Breeze regex, run on both bundles.
  - `bisectB`: `TmpGemmBisectJmh`, run on the new bundle only.
  - `syB`: `VectorSyrkJmh`, run on the new bundle only.
- **Execution order (interleaved ×2):** gv base r1, gv new r1, bisect new r1, syrk new r1, then the same again for r2.
- **Validation:** base r1 (0ff6440): {'passed': 162}; new r1 (4bc265d): {'passed': 162}; base r2 (0ff6440): {'passed': 162}; new r2 (4bc265d): {'passed': 162}. Each run's scoreboard (`breeze_scoreboard.py --lane B --strict --validation`) exited 0. The per-run scoreboards are reproduced from the raw JSON.

## Summary

**Yes: 7a90fb7 restores n = 256 vector gemm and gemmTall to base or better, and keeps (in fact extends) the n = 64 wins.** Values are gale ops/s, mean ± 99.9% CI, for base r1 / new r1 / base r2 / new r2.

| op (vector) | base r1 | new r1 | base r2 | new r2 | new / base | vs Breeze (new) |
|---|---:|---:|---:|---:|---:|---:|
| gemm n = 256 | 1,553.2 ± 8.5 | 1,599.8 ± 121.5 | 1,555.3 ± 7.3 | 1,550.6 ± 7.5 | 1.03× / 1.00× (CIs overlap: parity) | 1.36–1.41× ahead |
| gemmTall n = 256 | 459.5 ± 0.8 | 466.8 ± 1.2 | 459.8 ± 2.6 | 465.2 ± 1.6 | 1.02× / 1.01× (CI-separated) | 1.59–1.61× ahead |
| gemm n = 64 | 31,736 ± 319 | 78,938 ± 750 | 31,735 ± 291 | 81,126 ± 5,422 | 2.49× / 2.56× | 2.78–2.86× ahead |
| gemmTall n = 64 | 7,886 ± 64 | 21,914 ± 2,007 | 7,884 ± 67 | 21,550 ± 1,805 | 2.78× / 2.73× | 1.70–1.73× ahead |
| AtA n = 256 | 215.8 ± 10.4 | 650.9 ± 9.9 | 213.9 ± 2.2 | 651.3 ± 7.2 | 3.02× / 3.04× | 1.83× ahead (Breeze 356.5 ± 10.5) |

In the `4cdf0f3` run, gemm n = 256 vector was 0.85–0.88× of base and gemmTall n = 256 was 0.90×; both are now recovered. Pure is unchanged everywhere (1.00×).

**Bisect harness** (new bundle only, gemm via each variant's kernel, ops/s; `*` marks CIs disjoint from `base`):

| shape | variant | r1 ops/s | r2 ops/s | vs `base` r1 / r2 |
|---|---|---:|---:|---:|
| sq256 | base | 1,527.7 ± 24.3 | 1,519.5 ± 9.7 | 1.00x / 1.00x |
| sq256 | baseInline | 1,572.5 ± 130.1 | 1,532.5 ± 30.9 | 1.03x / 1.01x |
| sq256 | head | 1,388.7 ± 21.1 | 1,381.2 ± 22.2 | 0.91x* / 0.91x* |
| sq256 | fix | 1,522.5 ± 9.1 | 1,538.2 ± 21.6 | 1.00x / 1.01x |
| tall256 | base | 459.7 ± 2.2 | 458.4 ± 2.0 | 1.00x / 1.00x |
| tall256 | baseInline | 458.5 ± 2.5 | 459.8 ± 1.7 | 1.00x / 1.00x |
| tall256 | head | 410.6 ± 1.6 | 410.9 ± 1.7 | 0.89x* / 0.90x* |
| tall256 | fix | 467.0 ± 2.5 | 464.7 ± 3.6 | 1.02x* / 1.01x* |

`head` (the `4cdf0f3`-era shared tile loop) reproduces the regression: 0.91× on sq256 and 0.89–0.90× on tall256. `fix` is at parity with `base` on sq256 and slightly faster on tall256 (1.01–1.02×). `baseInline` is indistinguishable from `base`.

**Other changes vs base, CI-separated in both rounds:** QR n = 256 vector is 1.23× faster, AtA n = 64 vector 1.19–1.20× and gemvT n = 2048 / 256 vector 1.06–1.08×. **There are no regressions in both rounds.** Two slowdowns appeared in a single round only: gemvT n = 1024 vector (0.89× in r1, 0.97× in r2) and QR n = 16 vector (0.99× in r1). All factorization and decomposition rows (LU, Cholesky, QR, solve, lstsq, det, inv, pinv, svd, eig) are otherwise unchanged, and all 81 validated rows passed in each of the four Breeze runs.

**syrk guard at 7a90fb7:** the verdict from `552d2ad` (`2026-10-11-hc-gemm-vec-syrk-guard-trillium.md`) still holds, with the SIMD path now faster.
- **Wrong picks (pure chosen, SIMD faster):** 64×64 (SIMD 2.14× faster, up from 1.80×), 16×1024 (1.25–1.26×) and 8×1024 (1.03–1.04×).
- **Right picks:** 128×64, 256×64, 128×128 and 1024×256 (SIMD 1.20–3.03×); 4×1024 is a tie; 8×4096 correctly stays pure (pure 1.53–1.56× faster).
- **32×512 anomaly:** the guarded vector kernel (674 µs) is faster than `simdForced` (834–1,000 µs, noisy), so the bench hook may not exercise exactly the guarded SIMD path at this shape. Worth a look.

## gale / Breeze ratio by backend

Each cell is the speed ratio (>1 means gale is faster), with the verdict from non-overlapping 99.9% CIs, for round 1 / round 2. Breeze is measured in the same run.

| class | op | params | backend | base r1 | base r2 | new r1 | new r2 |
|---|---|---|---|---|---|---|---|
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 0.49x behind | 0.52x behind | 0.49x behind | 0.50x behind |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | 0.96x behind | 1.01x ahead | 0.98x tie | 0.99x tie |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 0.57x behind | 0.57x behind | 0.57x behind | 0.57x behind |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | 0.96x behind | 0.94x behind | 0.96x behind | 0.96x behind |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 0.42x behind | 0.42x behind | 0.43x behind | 0.43x behind |
| BlasL2BreezeJmh | Gemv | n=256 | vector | 0.96x behind | 0.96x behind | 0.97x behind | 0.98x tie |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 0.44x behind | 0.44x behind | 0.44x behind | 0.44x behind |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | 0.99x tie | 0.96x behind | 0.89x behind | 0.97x tie |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 0.49x behind | 0.49x behind | 0.49x behind | 0.49x behind |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | 0.94x behind | 0.95x behind | 1.00x tie | 1.00x tie |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 0.33x behind | 0.32x behind | 0.33x behind | 0.34x behind |
| BlasL2BreezeJmh | GemvT | n=256 | vector | 0.91x behind | 0.87x behind | 0.96x behind | 0.98x tie |
| BlasL3BreezeJmh | AtA | n=16 | pure | 1.99x ahead | 1.99x ahead | 2.00x ahead | 2.00x ahead |
| BlasL3BreezeJmh | AtA | n=16 | vector | 1.99x ahead | 1.99x ahead | 2.00x ahead | 2.00x ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 0.60x behind | 0.59x behind | 0.61x behind | 0.60x behind |
| BlasL3BreezeJmh | AtA | n=256 | vector | 0.60x behind | 0.60x behind | 1.83x ahead | 1.83x ahead |
| BlasL3BreezeJmh | AtA | n=64 | pure | 1.06x ahead | 1.06x ahead | 1.06x ahead | 1.06x ahead |
| BlasL3BreezeJmh | AtA | n=64 | vector | 1.06x ahead | 1.06x ahead | 1.27x ahead | 1.27x ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.16x ahead | 1.17x ahead | 1.17x ahead | 1.16x ahead |
| BlasL3BreezeJmh | Gemm | n=16 | vector | 1.16x ahead | 1.16x ahead | 1.17x ahead | 1.16x ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 0.46x behind | 0.45x behind | 0.46x behind | 0.45x behind |
| BlasL3BreezeJmh | Gemm | n=256 | vector | 1.36x ahead | 1.36x ahead | 1.41x ahead | 1.36x ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x ahead | 1.12x ahead | 1.12x ahead | 1.12x ahead |
| BlasL3BreezeJmh | Gemm | n=64 | vector | 1.12x ahead | 1.12x ahead | 2.78x ahead | 2.86x ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.10x ahead | 1.10x ahead | 1.10x ahead | 1.10x ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | 1.10x ahead | 1.10x ahead | 1.10x ahead | 1.10x ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 0.44x behind | 0.44x behind | 0.45x behind | 0.44x behind |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | 1.58x ahead | 1.58x ahead | 1.61x ahead | 1.59x ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 0.62x behind | 0.62x behind | 0.62x behind | 0.62x behind |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | 0.62x behind | 0.62x behind | 1.73x ahead | 1.70x ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.44x ahead | 1.51x ahead | 1.46x ahead | 1.46x ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | vector | 1.48x ahead | 1.49x ahead | 1.46x ahead | 1.45x ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.84x ahead | 1.84x ahead | 1.94x ahead | 1.83x ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | vector | 1.84x ahead | 1.84x ahead | 1.95x ahead | 1.95x ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 1.51x ahead | 1.51x ahead | 1.51x ahead | 1.51x ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | vector | 1.51x ahead | 1.51x ahead | 1.51x ahead | 1.51x ahead |
| DenseDecompositionBreezeJmh | Eig | n=16 | - | 2.60x ahead | 2.58x ahead | 2.58x ahead | 2.58x ahead |
| DenseDecompositionBreezeJmh | Eig | n=256 | - | 14.34x ahead | 14.39x ahead | 14.31x ahead | 14.45x ahead |
| DenseDecompositionBreezeJmh | Eig | n=64 | - | 1.51x ahead | 1.51x ahead | 1.51x ahead | 1.47x ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 1.37x ahead | 1.37x ahead | 1.36x ahead | 1.38x ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector | 1.36x ahead | 1.38x ahead | 1.41x ahead | 1.37x ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 2.20x ahead | 2.27x ahead | 2.22x ahead | 2.27x ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector | 2.20x ahead | 2.27x ahead | 2.27x ahead | 2.21x ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 1.78x ahead | 1.79x ahead | 1.78x ahead | 1.78x ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector | 1.78x ahead | 1.78x ahead | 1.78x ahead | 1.78x ahead |
| DenseDecompositionBreezeJmh | Pinv | n=16 | - | 2.93x ahead | 2.95x ahead | 2.95x ahead | 2.94x ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | - | 1.06x ahead | 1.06x ahead | 1.06x ahead | 1.06x ahead |
| DenseDecompositionBreezeJmh | Pinv | n=64 | - | 1.27x ahead | 1.27x ahead | 1.28x ahead | 1.27x ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | - | 2.76x ahead | 2.76x ahead | 2.78x ahead | 2.78x ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | - | 1.17x ahead | 1.17x ahead | 1.17x ahead | 1.17x ahead |
| DenseDecompositionBreezeJmh | Svd | n=64 | - | 1.31x ahead | 1.30x ahead | 1.30x ahead | 1.30x ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.89x ahead | 3.88x ahead | 3.91x ahead | 3.91x ahead |
| FactorizationBreezeJmh | Chol | n=16 | vector | 3.89x ahead | 3.88x ahead | 3.92x ahead | 3.91x ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.63x ahead | 1.61x ahead | 1.64x ahead | 1.64x ahead |
| FactorizationBreezeJmh | Chol | n=256 | vector | 1.63x ahead | 1.62x ahead | 1.64x ahead | 1.64x ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.99x ahead | 2.77x ahead | 2.90x ahead | 3.01x ahead |
| FactorizationBreezeJmh | Chol | n=64 | vector | 3.00x ahead | 2.77x ahead | 2.90x ahead | 3.01x ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.38x ahead | 1.37x ahead | 1.35x ahead | 1.37x ahead |
| FactorizationBreezeJmh | Lu | n=16 | vector | 1.36x ahead | 1.36x ahead | 1.35x ahead | 1.36x ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.95x ahead | 1.94x ahead | 1.95x ahead | 1.95x ahead |
| FactorizationBreezeJmh | Lu | n=256 | vector | 1.94x ahead | 1.95x ahead | 1.86x ahead | 1.96x ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 1.48x ahead | 1.47x ahead | 1.48x ahead | 1.48x ahead |
| FactorizationBreezeJmh | Lu | n=64 | vector | 1.48x ahead | 1.48x ahead | 1.48x ahead | 1.48x ahead |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.41x ahead | 1.42x ahead | 1.43x ahead | 1.39x ahead |
| FactorizationBreezeJmh | Qr | n=16 | vector | 1.42x ahead | 1.42x ahead | 1.40x ahead | 1.40x ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.48x ahead | 1.48x ahead | 1.48x ahead | 1.47x ahead |
| FactorizationBreezeJmh | Qr | n=256 | vector | 1.49x ahead | 1.49x ahead | 1.83x ahead | 1.83x ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 1.09x ahead | 1.09x ahead | 1.09x ahead | 1.09x ahead |
| FactorizationBreezeJmh | Qr | n=64 | vector | 1.09x ahead | 1.09x ahead | 1.09x ahead | 1.09x ahead |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.32x ahead | 1.31x ahead | 1.33x ahead | 1.34x ahead |
| FactorizationBreezeJmh | Solve | n=16 | vector | 1.31x ahead | 1.32x ahead | 1.34x ahead | 1.37x ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.78x ahead | 1.77x ahead | 1.77x ahead | 1.68x ahead |
| FactorizationBreezeJmh | Solve | n=256 | vector | 1.68x ahead | 1.68x ahead | 1.68x ahead | 1.68x ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 1.38x ahead | 1.37x ahead | 1.38x ahead | 1.38x ahead |
| FactorizationBreezeJmh | Solve | n=64 | vector | 1.38x ahead | 1.37x ahead | 1.37x ahead | 1.38x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.88x behind | 0.87x behind | 0.87x behind | 0.88x behind |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector | 0.88x behind | 0.87x behind | 0.87x behind | 0.87x behind |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.22x ahead | 1.22x ahead | 1.20x ahead | 1.22x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector | 1.27x ahead | 1.27x ahead | 1.29x ahead | 1.28x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.23x ahead | 1.22x ahead | 1.23x ahead | 1.22x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector | 1.22x ahead | 1.23x ahead | 1.23x ahead | 1.23x ahead |

## gale new / base (same backend, same round)

Ratio of gale speeds. A `*` marks rounds where the 99.9% CIs of new and base do not overlap. `regression` means both rounds are CI-separated and slower; `improved` means both are CI-separated and faster.

| class | op | params | backend | r1 | r2 | call |
|---|---|---|---|---|---|---|
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 0.99x | 1.00x |  |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | 1.00x | 1.00x |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 1.00x | 1.00x |  |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | 1.00x | 1.01x |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 1.01x | 1.01x |  |
| BlasL2BreezeJmh | Gemv | n=256 | vector | 1.00x | 1.02x |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 0.99x | 0.97x |  |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | 0.89x* | 0.97x |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 1.01x | 1.00x |  |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | 1.08x* | 1.07x* | improved |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 1.00x | 1.00x |  |
| BlasL2BreezeJmh | GemvT | n=256 | vector | 1.06x* | 1.06x* | improved |
| BlasL3BreezeJmh | AtA | n=16 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | AtA | n=16 | vector | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | 1.02x | 1.02x* |  |
| BlasL3BreezeJmh | AtA | n=256 | vector | 3.02x* | 3.04x* | improved |
| BlasL3BreezeJmh | AtA | n=64 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | AtA | n=64 | vector | 1.19x* | 1.20x* | improved |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=16 | vector | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=256 | vector | 1.03x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=64 | vector | 2.49x* | 2.56x* | improved |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | 1.02x* | 1.01x* | improved |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | 2.78x* | 2.73x* | improved |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.00x | 0.96x |  |
| DenseDecompositionBreezeJmh | Det | n=16 | vector | 0.98x | 0.97x |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.05x | 1.00x |  |
| DenseDecompositionBreezeJmh | Det | n=256 | vector | 1.06x | 1.06x |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Det | n=64 | vector | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | - | 1.00x | 0.97x |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 0.97x | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector | 1.02x | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 1.01x | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector | 1.03x | 0.98x |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Pinv | n=16 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | - | 1.01x | 1.00x |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | - | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=16 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=256 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=64 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | 0.99x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=16 | vector | 1.00x | 1.01x |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=256 | vector | 0.96x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=64 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.01x | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=16 | vector | 0.99x* | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.00x | 0.99x |  |
| FactorizationBreezeJmh | Qr | n=256 | vector | 1.23x* | 1.23x* | improved |
| FactorizationBreezeJmh | Qr | n=64 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=64 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Solve | n=16 | vector | 1.01x | 1.02x |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.00x | 0.95x |  |
| FactorizationBreezeJmh | Solve | n=256 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Solve | n=64 | vector | 1.00x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.99x | 1.01x* |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector | 0.99x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.00x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector | 1.02x | 1.01x |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.00x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector | 1.00x | 1.00x |  |

## Threshold sanity: vector / pure at each size (gale only)

Ratio of the vector backend's speed to the pure backend's, in the same run. A `*` marks non-overlapping 99.9% CIs.

| class | op | params | base r1 | base r2 | new r1 | new r2 |
|---|---|---|---|---|---|---|
| BlasL2BreezeJmh | Gemv | n=1024 | 1.96x* | 1.96x* | 1.98x* | 1.97x* |
| BlasL2BreezeJmh | Gemv | n=2048 | 1.68x* | 1.65x* | 1.67x* | 1.66x* |
| BlasL2BreezeJmh | Gemv | n=256 | 2.30x* | 2.28x* | 2.27x* | 2.30x* |
| BlasL2BreezeJmh | GemvT | n=1024 | 2.27x* | 2.19x* | 2.03x* | 2.20x* |
| BlasL2BreezeJmh | GemvT | n=2048 | 1.93x* | 1.92x* | 2.07x* | 2.06x* |
| BlasL2BreezeJmh | GemvT | n=256 | 2.76x* | 2.74x* | 2.91x* | 2.90x* |
| BlasL3BreezeJmh | AtA | n=16 | 1.00x | 1.00x | 1.00x | 1.00x |
| BlasL3BreezeJmh | AtA | n=256 | 1.01x | 1.02x* | 3.00x* | 3.06x* |
| BlasL3BreezeJmh | AtA | n=64 | 1.00x | 1.00x | 1.19x* | 1.19x* |
| BlasL3BreezeJmh | Gemm | n=16 | 1.00x | 1.00x | 1.01x | 1.00x |
| BlasL3BreezeJmh | Gemm | n=256 | 2.99x* | 3.00x* | 3.08x* | 2.99x* |
| BlasL3BreezeJmh | Gemm | n=64 | 1.00x | 1.00x | 2.49x* | 2.55x* |
| BlasL3BreezeJmh | GemmTall | n=16 | 1.00x | 1.00x | 1.00x | 1.00x |
| BlasL3BreezeJmh | GemmTall | n=256 | 3.55x* | 3.56x* | 3.61x* | 3.60x* |
| BlasL3BreezeJmh | GemmTall | n=64 | 1.00x | 1.00x | 2.78x* | 2.73x* |
| DenseDecompositionBreezeJmh | Det | n=16 | 1.03x | 0.99x | 1.00x | 0.99x |
| DenseDecompositionBreezeJmh | Det | n=256 | 1.00x | 1.00x | 1.01x | 1.07x |
| DenseDecompositionBreezeJmh | Det | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| DenseDecompositionBreezeJmh | Inv | n=16 | 1.00x | 1.00x | 1.04x* | 1.00x |
| DenseDecompositionBreezeJmh | Inv | n=256 | 1.00x | 1.00x | 1.02x | 0.97x |
| DenseDecompositionBreezeJmh | Inv | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Chol | n=16 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Chol | n=256 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Chol | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Lu | n=16 | 0.99x | 0.99x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Lu | n=256 | 1.00x | 1.01x | 0.96x | 1.01x |
| FactorizationBreezeJmh | Lu | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Qr | n=16 | 1.01x | 1.00x | 0.98x | 1.00x |
| FactorizationBreezeJmh | Qr | n=256 | 1.01x* | 1.01x* | 1.24x* | 1.25x* |
| FactorizationBreezeJmh | Qr | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Solve | n=16 | 1.00x | 1.01x | 1.01x | 1.02x |
| FactorizationBreezeJmh | Solve | n=256 | 0.94x | 0.95x | 0.95x | 1.00x |
| FactorizationBreezeJmh | Solve | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| LeastSquaresBreezeJmh | Lstsq | n=16 | 1.00x | 1.00x | 1.00x | 1.00x* |
| LeastSquaresBreezeJmh | Lstsq | n=256 | 1.05x* | 1.04x* | 1.07x* | 1.05x* |
| LeastSquaresBreezeJmh | Lstsq | n=64 | 1.00x | 1.00x | 1.00x | 1.01x |

## Per shape (µs/op, mean ± 99.9% CI half-width; lower is better)

| shape (m×k) | round | pure kernel | simdForced | guarded vector kernel | guard picks | faster path | flag |
|---|---:|---:|---:|---:|---|---|---|
| 64x64 | 1 | 16.71 ± 0.01 | 7.80 ± 0.03 | 16.71 ± 0.00 | pure | SIMD (SIMD 2.14× faster) | **guard picks slower path** |
| 64x64 | 2 | 16.71 ± 0.02 | 7.82 ± 0.11 | 16.71 ± 0.01 | pure | SIMD (SIMD 2.14× faster) | **guard picks slower path** |
| 128x64 | 1 | 32.57 ± 0.04 | 14.39 ± 0.03 | 14.37 ± 0.04 | SIMD | SIMD (SIMD 2.26× faster) |  |
| 128x64 | 2 | 32.60 ± 0.08 | 14.41 ± 0.03 | 14.37 ± 0.06 | SIMD | SIMD (SIMD 2.26× faster) |  |
| 256x64 | 1 | 63.95 ± 0.07 | 53.44 ± 0.24 | 53.39 ± 0.23 | SIMD | SIMD (SIMD 1.20× faster) |  |
| 256x64 | 2 | 63.93 ± 0.02 | 53.46 ± 0.24 | 53.34 ± 0.19 | SIMD | SIMD (SIMD 1.20× faster) |  |
| 128x128 | 1 | 128.3 ± 0.1 | 51.36 ± 0.19 | 51.31 ± 0.23 | SIMD | SIMD (SIMD 2.50× faster) |  |
| 128x128 | 2 | 128.3 ± 0.1 | 51.45 ± 0.17 | 51.36 ± 0.31 | SIMD | SIMD (SIMD 2.49× faster) |  |
| 1024x256 | 1 | 4,562.8 ± 45.3 | 1,504.8 ± 9.8 | 1,507.2 ± 10.5 | SIMD | SIMD (SIMD 3.03× faster) |  |
| 1024x256 | 2 | 4,521.9 ± 35.2 | 1,507.6 ± 17.4 | 1,501.6 ± 8.7 | SIMD | SIMD (SIMD 3.00× faster) |  |
| 4x1024 | 1 | 2,716.4 ± 34.1 | 2,767.5 ± 27.6 | 2,705.3 ± 17.0 | pure | tie (pure/SIMD 0.98×) |  |
| 4x1024 | 2 | 2,713.3 ± 21.3 | 2,816.4 ± 94.7 | 2,711.4 ± 35.0 | pure | tie (pure/SIMD 0.96×) |  |
| 8x1024 | 1 | 2,901.8 ± 23.7 | 2,799.2 ± 5.1 | 2,955.1 ± 25.1 | pure | SIMD (SIMD 1.04× faster) | **guard picks slower path** |
| 8x1024 | 2 | 2,887.6 ± 6.4 | 2,802.9 ± 22.1 | 2,902.3 ± 28.9 | pure | SIMD (SIMD 1.03× faster) | **guard picks slower path** |
| 16x1024 | 1 | 3,622.1 ± 11.2 | 2,867.7 ± 53.1 | 3,619.2 ± 33.4 | pure | SIMD (SIMD 1.26× faster) | **guard picks slower path** |
| 16x1024 | 2 | 3,596.9 ± 27.5 | 2,868.1 ± 36.2 | 3,628.8 ± 27.8 | pure | SIMD (SIMD 1.25× faster) | **guard picks slower path** |
| 32x512 | 1 | 1,069.8 ± 1.4 | 1,000.3 ± 2.8 | 674.2 ± 2.8 | SIMD | SIMD (SIMD 1.07× faster) |  |
| 32x512 | 2 | 1,173.0 ± 175.5 | 834.0 ± 251.4 | 674.7 ± 1.2 | SIMD | tie (pure/SIMD 1.41×) |  |
| 8x4096 | 1 | 133,016.1 ± 3,734.3 | 208,003.5 ± 3,530.8 | 131,225.7 ± 2,912.5 | pure | pure (pure 1.56× faster) |  |
| 8x4096 | 2 | 142,374.2 ± 4,232.2 | 217,830.3 ± 1,481.2 | 127,776.1 ± 2,884.1 | pure | pure (pure 1.53× faster) |  |

## Facade `A.t * A` (µs/op)

| shape | round | facade pure | facade vector | vector/pure speedup |
|---|---:|---:|---:|---:|
| 64x64 | 1 | 17.96 ± 0.05 | 17.97 ± 0.07 | 1.00× |
| 64x64 | 2 | 17.96 ± 0.03 | 17.96 ± 0.04 | 1.00× |
| 128x64 | 1 | 34.10 ± 0.38 | 15.31 ± 0.04 | 2.23×* |
| 128x64 | 2 | 34.12 ± 0.36 | 15.30 ± 0.03 | 2.23×* |
| 256x64 | 1 | 65.07 ± 0.45 | 54.67 ± 0.08 | 1.19×* |
| 256x64 | 2 | 65.19 ± 0.34 | 54.67 ± 0.08 | 1.19×* |
| 128x128 | 1 | 135.4 ± 1.6 | 59.54 ± 0.18 | 2.27×* |
| 128x128 | 2 | 135.3 ± 1.4 | 59.61 ± 0.20 | 2.27×* |
| 1024x256 | 1 | 4,585.2 ± 81.2 | 1,535.3 ± 18.2 | 2.99×* |
| 1024x256 | 2 | 4,562.6 ± 114.4 | 1,535.2 ± 17.9 | 2.97×* |
| 4x1024 | 1 | 2,963.2 ± 15.7 | 2,964.1 ± 14.0 | 1.00× |
| 4x1024 | 2 | 2,960.6 ± 11.7 | 2,964.5 ± 19.5 | 1.00× |
| 8x1024 | 1 | 3,158.1 ± 19.6 | 3,164.9 ± 16.2 | 1.00× |
| 8x1024 | 2 | 3,157.9 ± 12.1 | 3,159.1 ± 13.4 | 1.00× |
| 16x1024 | 1 | 3,855.0 ± 44.4 | 3,845.2 ± 8.5 | 1.00× |
| 16x1024 | 2 | 3,841.9 ± 10.1 | 3,843.0 ± 10.8 | 1.00× |
| 32x512 | 1 | 1,121.9 ± 4.8 | 737.6 ± 3.0 | 1.52×* |
| 32x512 | 2 | 1,120.7 ± 6.8 | 736.4 ± 1.6 | 1.52×* |
| 8x4096 | 1 | 132,925.0 ± 4,809.3 | 131,035.8 ± 3,697.0 | 1.01× |
| 8x4096 | 2 | 131,395.2 ± 4,375.8 | 134,130.1 ± 4,632.6 | 0.98× |

Raw data: `2026-10-11-hc-gemm-vec-fix-trillium.json`. It holds the Breeze-paired runs (JMH, netlib and validation records), the bisect runs and the syrk runs, with source SHA and load.
