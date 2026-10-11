# hc-gemm-vec paired run on trillium (Zen 5, lane B)

`4cdf0f3` (breeze/hc-gemm-vec: SIMD syrk, gemm/syrk vector crossover lowered from 128³ to 64³, unboxed gemm edges / vectorDot, 8-column gemvT) against the Milestone 2 baseline. Slurm job 2536008, one exclusive compute node (`tri0949`, 2 × AMD EPYC 9655, Zen 5, AVX-512, 8 double lanes).

## Protocol

- **Bundles:** base is the Milestone 2 bundle compiled from `0ff6440`, which is source-identical to `5d3950e`. New is `4cdf0f3`, compiled locally in a clean worktree. The breeze/hc-gemm-vec worktree had uncommitted work and an active sbt session, so it was not used. The clean worktree was on local branch `bench/hcgv-4cdf0f3` at that commit, because sbt-git cannot read a detached-HEAD linked worktree. No sbt ran on the cluster.
- **Harness:** `docs/verification/w21-simd-spike/x86-avx512/gv.sh`, which is the reviewed `run.sh` lane-B path: inherited `--add-modules`, `-Xms4g -Xmx4g`, a bundle-bound `gale.bench.sourceSha`, `-foe true`, a fresh netlib and validation sidecar per run, and `-p backend=pure,vector`. It is restricted to `BlasL3BreezeJmh.(galeAtA|breezeAtA|galeGemmTall|breezeGemmTall|galeGemm|breezeGemm)$`, `BlasL2BreezeJmh.(galeGemvT|breezeGemvT|galeGemv|breezeGemv)$`, `FactorizationBreezeJmh`, `DenseDecompositionBreezeJmh` and `LeastSquaresBreezeJmh`. It used the annotation defaults (2 forks), Temurin 25+36 with `JAVA_TOOL_OPTIONS` unset, and `numactl --physcpubind=0-7 --membind=0`.
- **Order:** interleaved: laneB base r1; laneB new r1; laneB base r2; laneB new r2. Load at start of each run: base r1: 363.16, 613.36, 407.60; new r1: 1.16, 6.00, 86.82; base r2: 1.17, 1.23, 19.20; new r2: 1.27, 1.16, 4.93.
- **Netlib** (blas, lapack, vectorModule): [('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- **Validation:** base r1 (0ff6440): {'passed': 162}; new r1 (4cdf0f3): {'passed': 162}; base r2 (0ff6440): {'passed': 162}; new r2 (4cdf0f3): {'passed': 162}. Each run's scoreboard (`breeze_scoreboard.py --lane B --strict --validation`) exited 0. The per-run scoreboards are reproduced from the raw JSON.

## gale / Breeze ratio by backend

Each cell is the speed ratio (>1 means gale is faster), with the verdict from non-overlapping 99.9% CIs, for round 1 / round 2. Breeze is measured in the same run.

| class | op | params | backend | base r1 | base r2 | new r1 | new r2 |
|---|---|---|---|---|---|---|---|
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 0.49x behind | 0.50x behind | 0.52x behind | 0.52x behind |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | 0.96x behind | 0.98x tie | 1.01x ahead | 1.00x tie |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 0.56x behind | 0.57x behind | 0.57x behind | 0.56x behind |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | 0.94x behind | 0.95x behind | 0.98x tie | 0.94x behind |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 0.42x behind | 0.40x behind | 0.42x behind | 0.42x behind |
| BlasL2BreezeJmh | Gemv | n=256 | vector | 1.14x ahead | 0.91x tie | 0.96x behind | 0.96x behind |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 0.44x behind | 0.44x behind | 0.45x behind | 0.43x behind |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | 0.99x tie | 0.98x tie | 0.95x behind | 0.93x behind |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 0.48x behind | 0.49x behind | 0.48x behind | 0.49x behind |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | 0.92x behind | 0.95x behind | 1.01x tie | 1.00x tie |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 0.32x behind | 0.33x behind | 0.34x behind | 0.32x behind |
| BlasL2BreezeJmh | GemvT | n=256 | vector | 0.88x behind | 0.99x tie | 0.98x tie | 0.94x tie |
| BlasL3BreezeJmh | AtA | n=16 | pure | 2.00x ahead | 1.99x ahead | 2.00x ahead | 1.99x ahead |
| BlasL3BreezeJmh | AtA | n=16 | vector | 2.00x ahead | 1.99x ahead | 1.99x ahead | 1.99x ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 0.60x behind | 0.60x behind | 0.59x behind | 0.62x behind |
| BlasL3BreezeJmh | AtA | n=256 | vector | 0.61x behind | 0.61x behind | 0.95x tie | 0.98x tie |
| BlasL3BreezeJmh | AtA | n=64 | pure | 1.06x ahead | 1.06x ahead | 1.06x ahead | 1.06x ahead |
| BlasL3BreezeJmh | AtA | n=64 | vector | 1.06x ahead | 1.06x ahead | 1.26x ahead | 1.26x ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.17x ahead | 1.17x ahead | 1.17x ahead | 1.17x ahead |
| BlasL3BreezeJmh | Gemm | n=16 | vector | 1.17x ahead | 1.17x ahead | 1.17x ahead | 1.16x ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 0.46x behind | 0.45x behind | 0.45x behind | 0.46x behind |
| BlasL3BreezeJmh | Gemm | n=256 | vector | 1.40x ahead | 1.45x ahead | 1.22x ahead | 1.23x ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x ahead | 1.12x ahead | 1.12x ahead | 1.12x ahead |
| BlasL3BreezeJmh | Gemm | n=64 | vector | 1.12x ahead | 1.12x ahead | 2.34x ahead | 2.34x ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.10x ahead | 1.09x ahead | 1.10x ahead | 1.10x ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | 1.10x ahead | 1.09x ahead | 1.10x ahead | 1.10x ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 0.44x behind | 0.44x behind | 0.45x behind | 0.44x behind |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | 1.57x ahead | 1.57x ahead | 1.42x ahead | 1.42x ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 0.62x behind | 0.62x behind | 0.62x behind | 0.62x behind |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | 0.62x behind | 0.62x behind | 1.40x ahead | 1.42x ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.46x ahead | 1.45x ahead | 1.48x ahead | 1.47x ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | vector | 1.46x ahead | 1.46x ahead | 1.49x ahead | 1.47x ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.94x ahead | 1.94x ahead | 1.83x ahead | 1.83x ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | vector | 1.94x ahead | 1.95x ahead | 1.94x ahead | 1.95x ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 1.51x ahead | 1.51x ahead | 1.51x ahead | 1.51x ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | vector | 1.50x ahead | 1.52x ahead | 1.51x ahead | 1.51x ahead |
| DenseDecompositionBreezeJmh | Eig | n=16 | - | 2.58x ahead | 2.66x ahead | 2.64x ahead | 2.59x ahead |
| DenseDecompositionBreezeJmh | Eig | n=256 | - | 14.31x ahead | 14.34x ahead | 14.36x ahead | 14.36x ahead |
| DenseDecompositionBreezeJmh | Eig | n=64 | - | 1.50x ahead | 1.50x ahead | 1.51x ahead | 1.51x ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 1.39x ahead | 1.36x ahead | 1.37x ahead | 1.40x ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector | 1.37x ahead | 1.40x ahead | 1.35x ahead | 1.39x ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 2.26x ahead | 2.20x ahead | 2.27x ahead | 2.20x ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector | 2.26x ahead | 2.27x ahead | 2.26x ahead | 2.26x ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 1.78x ahead | 1.78x ahead | 1.78x ahead | 1.78x ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector | 1.78x ahead | 1.78x ahead | 1.78x ahead | 1.79x ahead |
| DenseDecompositionBreezeJmh | Pinv | n=16 | - | 2.94x ahead | 2.94x ahead | 2.95x ahead | 2.94x ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | - | 1.06x ahead | 1.06x ahead | 1.06x ahead | 1.06x ahead |
| DenseDecompositionBreezeJmh | Pinv | n=64 | - | 1.27x ahead | 1.28x ahead | 1.28x ahead | 1.27x ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | - | 2.76x ahead | 2.77x ahead | 2.77x ahead | 2.76x ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | - | 1.17x ahead | 1.17x ahead | 1.17x ahead | 1.17x ahead |
| DenseDecompositionBreezeJmh | Svd | n=64 | - | 1.30x ahead | 1.30x ahead | 1.30x ahead | 1.30x ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.72x ahead | 3.75x ahead | 3.89x ahead | 3.90x ahead |
| FactorizationBreezeJmh | Chol | n=16 | vector | 3.72x ahead | 3.75x ahead | 3.89x ahead | 3.87x ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.63x ahead | 1.63x ahead | 1.64x ahead | 1.62x ahead |
| FactorizationBreezeJmh | Chol | n=256 | vector | 1.63x ahead | 1.63x ahead | 1.64x ahead | 1.62x ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.90x ahead | 2.89x ahead | 2.89x ahead | 2.77x ahead |
| FactorizationBreezeJmh | Chol | n=64 | vector | 2.90x ahead | 2.89x ahead | 2.89x ahead | 2.77x ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.35x ahead | 1.57x ahead | 1.35x ahead | 1.33x ahead |
| FactorizationBreezeJmh | Lu | n=16 | vector | 1.35x ahead | 1.55x ahead | 1.35x ahead | 1.36x ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.95x ahead | 1.95x ahead | 1.94x ahead | 1.95x ahead |
| FactorizationBreezeJmh | Lu | n=256 | vector | 1.96x ahead | 1.95x ahead | 1.96x ahead | 1.96x ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 1.48x ahead | 1.48x ahead | 1.48x ahead | 1.48x ahead |
| FactorizationBreezeJmh | Lu | n=64 | vector | 1.48x ahead | 1.48x ahead | 1.53x ahead | 1.48x ahead |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.39x ahead | 1.45x ahead | 1.41x ahead | 1.41x ahead |
| FactorizationBreezeJmh | Qr | n=16 | vector | 1.41x ahead | 1.42x ahead | 1.40x ahead | 1.40x ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.48x ahead | 1.48x ahead | 1.48x ahead | 1.48x ahead |
| FactorizationBreezeJmh | Qr | n=256 | vector | 1.49x ahead | 1.49x ahead | 1.87x ahead | 1.87x ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 1.08x ahead | 1.09x ahead | 1.09x ahead | 1.09x ahead |
| FactorizationBreezeJmh | Qr | n=64 | vector | 1.09x ahead | 1.09x ahead | 1.09x ahead | 1.09x ahead |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.32x ahead | 1.33x ahead | 1.30x ahead | 1.32x ahead |
| FactorizationBreezeJmh | Solve | n=16 | vector | 1.33x ahead | 1.35x ahead | 1.36x ahead | 1.29x ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.68x ahead | 1.68x ahead | 1.67x ahead | 1.68x ahead |
| FactorizationBreezeJmh | Solve | n=256 | vector | 1.67x ahead | 1.68x ahead | 1.67x ahead | 1.67x ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 1.38x ahead | 1.32x ahead | 1.37x ahead | 1.37x ahead |
| FactorizationBreezeJmh | Solve | n=64 | vector | 1.38x ahead | 1.37x ahead | 1.38x ahead | 1.38x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.87x behind | 0.87x behind | 0.87x behind | 0.88x behind |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector | 0.88x behind | 0.88x behind | 0.87x behind | 0.88x behind |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.22x ahead | 1.22x ahead | 1.20x ahead | 1.22x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector | 1.27x ahead | 1.27x ahead | 1.29x ahead | 1.31x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.22x ahead | 1.22x ahead | 1.23x ahead | 1.23x ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector | 1.23x ahead | 1.23x ahead | 1.23x ahead | 1.22x ahead |

## gale new / base (same backend, same round)

Ratio of gale speeds. A `*` marks rounds where the 99.9% CIs of new and base do not overlap. `regression` means both rounds are CI-separated and slower; `improved` means both are CI-separated and faster.

| class | op | params | backend | r1 | r2 | call |
|---|---|---|---|---|---|---|
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 1.00x | 1.01x |  |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | 1.00x | 1.00x |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 1.00x | 0.99x |  |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | 1.02x | 1.00x |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 1.01x | 1.01x* |  |
| BlasL2BreezeJmh | Gemv | n=256 | vector | 0.85x* | 1.00x |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 1.00x | 1.00x |  |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | 0.94x* | 0.97x |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 0.99x | 1.01x |  |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | 1.09x* | 1.07x* | improved |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 0.99x | 1.00x |  |
| BlasL2BreezeJmh | GemvT | n=256 | vector | 1.06x* | 0.98x |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | AtA | n=16 | vector | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | 0.98x | 1.04x |  |
| BlasL3BreezeJmh | AtA | n=256 | vector | 1.56x* | 1.61x* | improved |
| BlasL3BreezeJmh | AtA | n=64 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | AtA | n=64 | vector | 1.19x* | 1.19x* | improved |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=16 | vector | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=256 | vector | 0.88x* | 0.85x* | regression |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | Gemm | n=64 | vector | 2.09x* | 2.09x* | improved |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | 0.90x* | 0.90x* | regression |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 1.00x | 1.00x |  |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | 2.24x* | 2.29x* | improved |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.01x | 1.00x |  |
| DenseDecompositionBreezeJmh | Det | n=16 | vector | 1.01x | 1.00x |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 0.94x | 0.94x |  |
| DenseDecompositionBreezeJmh | Det | n=256 | vector | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Det | n=64 | vector | 1.00x | 0.99x* |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | - | 1.03x | 0.98x |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | - | 1.00x* | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 1.00x | 1.03x* |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector | 1.00x | 0.99x |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 1.01x | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Pinv | n=16 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | - | 1.01x | 1.00x |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | - | 1.00x | 1.00x |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | - | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=16 | vector | 1.00x | 0.99x |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=256 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Chol | n=64 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.00x | 0.98x |  |
| FactorizationBreezeJmh | Lu | n=16 | vector | 1.00x | 1.01x |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | 0.99x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=256 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Lu | n=64 | vector | 1.03x | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | 0.99x | 0.98x |  |
| FactorizationBreezeJmh | Qr | n=16 | vector | 0.98x | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=256 | vector | 1.26x* | 1.26x* | improved |
| FactorizationBreezeJmh | Qr | n=64 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Qr | n=64 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | 0.99x | 1.01x |  |
| FactorizationBreezeJmh | Solve | n=16 | vector | 1.03x | 0.97x |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Solve | n=256 | vector | 1.00x | 1.00x |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | 1.00x | 1.04x |  |
| FactorizationBreezeJmh | Solve | n=64 | vector | 1.00x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 1.00x | 1.01x* |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector | 1.00x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.00x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector | 1.04x | 1.03x |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.00x | 1.00x |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector | 1.00x | 1.00x |  |

## Threshold sanity: vector / pure at each size (gale only)

Ratio of the vector backend's speed to the pure backend's, in the same run. A `*` marks non-overlapping 99.9% CIs.

| class | op | params | base r1 | base r2 | new r1 | new r2 |
|---|---|---|---|---|---|---|
| BlasL2BreezeJmh | Gemv | n=1024 | 1.95x* | 1.97x* | 1.96x* | 1.94x* |
| BlasL2BreezeJmh | Gemv | n=2048 | 1.68x* | 1.66x* | 1.71x* | 1.68x* |
| BlasL2BreezeJmh | Gemv | n=256 | 2.72x* | 2.29x* | 2.28x* | 2.27x* |
| BlasL2BreezeJmh | GemvT | n=1024 | 2.26x* | 2.19x* | 2.12x* | 2.13x* |
| BlasL2BreezeJmh | GemvT | n=2048 | 1.91x* | 1.92x* | 2.10x* | 2.04x* |
| BlasL2BreezeJmh | GemvT | n=256 | 2.74x* | 2.98x* | 2.93x* | 2.90x* |
| BlasL3BreezeJmh | AtA | n=16 | 1.00x | 1.00x | 1.00x | 1.00x |
| BlasL3BreezeJmh | AtA | n=256 | 1.01x | 1.01x | 1.60x* | 1.56x* |
| BlasL3BreezeJmh | AtA | n=64 | 1.00x | 1.00x | 1.19x* | 1.19x* |
| BlasL3BreezeJmh | Gemm | n=16 | 1.00x | 1.00x | 1.00x | 1.00x |
| BlasL3BreezeJmh | Gemm | n=256 | 3.08x* | 3.19x* | 2.69x* | 2.69x* |
| BlasL3BreezeJmh | Gemm | n=64 | 1.00x | 1.00x | 2.10x* | 2.09x* |
| BlasL3BreezeJmh | GemmTall | n=16 | 1.00x | 1.00x | 1.00x | 1.00x |
| BlasL3BreezeJmh | GemmTall | n=256 | 3.56x* | 3.55x* | 3.19x* | 3.20x* |
| BlasL3BreezeJmh | GemmTall | n=64 | 1.00x | 1.00x | 2.25x* | 2.29x* |
| DenseDecompositionBreezeJmh | Det | n=16 | 1.00x | 1.00x | 1.00x | 1.00x |
| DenseDecompositionBreezeJmh | Det | n=256 | 1.00x | 1.01x | 1.06x | 1.06x |
| DenseDecompositionBreezeJmh | Det | n=64 | 1.00x | 1.01x* | 1.00x | 1.00x |
| DenseDecompositionBreezeJmh | Inv | n=16 | 0.99x | 1.03x | 0.99x | 0.99x |
| DenseDecompositionBreezeJmh | Inv | n=256 | 1.00x | 1.03x | 0.99x | 1.03x |
| DenseDecompositionBreezeJmh | Inv | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Chol | n=16 | 1.00x | 1.00x | 1.00x | 0.99x |
| FactorizationBreezeJmh | Chol | n=256 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Chol | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Lu | n=16 | 1.00x | 0.99x | 1.00x | 1.02x |
| FactorizationBreezeJmh | Lu | n=256 | 1.01x | 1.00x | 1.01x | 1.00x |
| FactorizationBreezeJmh | Lu | n=64 | 1.00x | 1.00x | 1.03x | 1.00x |
| FactorizationBreezeJmh | Qr | n=16 | 1.01x | 0.98x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Qr | n=256 | 1.01x* | 1.01x | 1.27x* | 1.27x* |
| FactorizationBreezeJmh | Qr | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Solve | n=16 | 1.01x | 1.02x | 1.05x* | 0.97x |
| FactorizationBreezeJmh | Solve | n=256 | 1.00x | 1.00x | 1.00x | 1.00x |
| FactorizationBreezeJmh | Solve | n=64 | 1.00x | 1.04x | 1.00x | 1.00x |
| LeastSquaresBreezeJmh | Lstsq | n=16 | 1.00x | 1.00x | 1.00x | 0.99x |
| LeastSquaresBreezeJmh | Lstsq | n=256 | 1.04x* | 1.04x* | 1.08x* | 1.07x* |
| LeastSquaresBreezeJmh | Lstsq | n=64 | 1.00x | 1.00x | 1.00x | 1.00x |

## Summary

**Threshold sanity at n = 64 on AVX-512 (8 lanes): yes, vector beats pure clearly now that the crossover is 64³.** Values are gale ops/s, mean ± 99.9% CI, for new r1 / new r2. Base routed n = 64 to pure on both backends, so base vector/pure was 1.00×.

| op, n = 64 | pure | vector | vector / pure | vector / Breeze |
|---|---:|---:|---:|---:|
| gemm | 31,674 ± 289 / 31,661 ± 292 | 66,361 ± 556 / 66,183 ± 1,307 | 2.09–2.10× | 2.34× (Breeze 28,335 ± 163) |
| gemmTall (256×64) | 7,871 ± 64 / 7,867 ± 64 | 17,673 ± 1,011 / 17,995 ± 760 | 2.25–2.29× | 1.40–1.42× (Breeze 12,641 ± 26); pure is 0.62× |
| AtA (syrk) | 15,307 ± 119 / 15,249 ± 76 | 18,184 ± 39 / 18,176 ± 35 | 1.19× | 1.26× (Breeze 14,405 ± 28) |

The factorizations at n = 64 (LU, Cholesky, QR, solve, lstsq) are unchanged by the threshold: vector/pure is 1.00×. Their blocked updates at n = 64 are evidently below the size that reaches the gemm path.

**Improvements over base, CI-separated in both rounds:**

| op | backend | new / base | gale ops/s, base r1 → new r1 |
|---|---|---|---|
| gemm n = 64 | vector | 2.09× | 31,679 ± 307 → 66,361 ± 556 |
| gemmTall n = 64 | vector | 2.24–2.29× | 7,875 ± 55 → 17,673 ± 1,011 |
| AtA n = 64 | vector | 1.19× | 15,246 ± 81 → 18,184 ± 39 |
| AtA n = 256 | vector | 1.56–1.61× | 216.2 ± 6.5 → 336.8 ± 11.1 |
| QR n = 256 | vector | 1.26× | 344.5 ± 1.4 → 433.4 ± 8.4 |
| gemvT n = 2048 | vector | 1.07–1.09× | 2,713 ± 17 → 2,967 ± 32 |

These move AtA n = 256 vector from 0.61× of Breeze (behind) to 0.95–0.98× (tie), and QR n = 256 vector from 1.49× to 1.87×.

**Regressions, CI-separated in both rounds:**

| op | backend | new / base | gale ops/s, base r1 / r2 → new r1 / r2 | vs Breeze, base → new |
|---|---|---|---|---|
| gemm n = 256 | vector | 0.88× / 0.85× | 1,592.7 ± 103.7 / 1,649.7 ± 135.7 → 1,395.1 ± 16.3 / 1,395.9 ± 18.7 | 1.40–1.45× → 1.22–1.23× (still ahead) |
| gemmTall n = 256 | vector | 0.90× / 0.90× | 458.2 ± 1.4 / 457.9 ± 2.2 → 411.6 ± 1.0 / 412.3 ± 1.1 | 1.57× → 1.42× (still ahead) |

Both are large-gemm vector-path slowdowns of 10–15%. The pure path is unchanged. These are candidates for the unboxed-edge/tile change in `4cdf0f3`.

Single-round, CI-separated slowdowns that did not repeat in the other round: gemv n = 256 vector (0.85× in r1, 1.00× in r2; the base r1 Breeze-relative ratio was itself an outlier, 1.14× against 0.91× in r2), gemvT n = 1024 vector (0.94× / 0.97×), and det n = 64 vector (1.00× / 0.99×).

**Factorization regressions: none.** LU, Cholesky, QR, solve, lstsq, det, inv, pinv, svd and eig show no CI-separated slowdown in both rounds, on either backend. All 81 validated decomposition/factorization rows passed in each of the four runs.

Still behind Breeze in the new bundle: pure gemv/gemvT (0.32–0.57×), pure gemm/gemmTall/AtA at n = 256 (0.44–0.62×), pure gemmTall n = 64 (0.62×), lstsq n = 16 (0.87–0.88×), and vector gemvT n = 1024 (0.93–0.95×).

Raw data: `2026-10-11-hc-gemm-vec-trillium.json`. It holds the per-run JMH JSON, netlib and validation sidecar records, source SHA and load, keyed by variant and round.
