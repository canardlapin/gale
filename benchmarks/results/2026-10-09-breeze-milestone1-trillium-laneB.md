# gale vs Breeze scoreboard — lane B (SIMD)

- Lane: B
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: --add-modules=jdk.incubator.vector -Xms4g -Xmx4g -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/m1/out/laneB-2526675/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.VectorBLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `ff03eed494bfa2d246b1c363c59f3425d97fcd67`
- Machine: SciNet trillium tri0017: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512, 8 double lanes), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2526675 exclusive

**132 ahead, 8 tie, 58 behind, 0 n/a, 11 withheld** of 209 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 171467 ± 7725 | 174731 ± 2002 | ops/s | 0.98x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 29236.4 ± 2834 | 29880.2 ± 2818 | ops/s | 0.98x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7326.72 ± 877.4 | 7980.74 ± 415.5 | ops/s | 0.92x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 64115.4 ± 2753 | 67387 ± 11.63 | ops/s | 0.95x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16783.9 ± 194.6 | 17044.8 ± 8.542 | ops/s | 0.98x | behind |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4193.9 ± 36.45 | 4259.31 ± 11.67 | ops/s | 0.98x | behind |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 65704 ± 1763 | 4624.86 ± 15.56 | ops/s | 14.21x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 16864.7 ± 197.5 | 1238.24 ± 272.9 | ops/s | 13.62x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 4215.98 ± 30.3 | 357.052 ± 0.1254 | ops/s | 11.81x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 98171.9 ± 659.4 | 233329 ± 564.7 | ops/s | 0.42x | behind |  |
| BlasL2BreezeJmh | Gemv | n=256 | vector | thrpt | 224157 ± 3930 | 233329 ± 564.7 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6515.93 ± 38.37 | 12765.1 ± 714.9 | ops/s | 0.51x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | thrpt | 12644.8 ± 52.49 | 12765.1 ± 714.9 | ops/s | 0.99x | tie |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1610.94 ± 11.31 | 2775.81 ± 1.69 | ops/s | 0.58x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | thrpt | 2698.72 ± 28.23 | 2775.81 ± 1.69 | ops/s | 0.97x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 82851.5 ± 288.2 | 245601 ± 4313 | ops/s | 0.34x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | vector | thrpt | 226735 ± 932.8 | 245601 ± 4313 | ops/s | 0.92x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5769.75 ± 1.142 | 13513.5 ± 98 | ops/s | 0.43x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | thrpt | 12911.3 ± 54.42 | 13513.5 ± 98 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1420.42 ± 5.824 | 2943.26 ± 12.59 | ops/s | 0.48x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | thrpt | 2741.41 ± 38.5 | 2943.26 ± 12.59 | ops/s | 0.93x | behind |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 787944 ± 1271 | 394331 ± 7131 | ops/s | 2.00x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | vector | thrpt | 785390 ± 2607 | 394331 ± 7131 | ops/s | 1.99x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15232.2 ± 85.03 | 14375.6 ± 27.7 | ops/s | 1.06x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | vector | thrpt | 15289.8 ± 127.9 | 14375.6 ± 27.7 | ops/s | 1.06x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 213.705 ± 7.159 | 355.284 ± 10.16 | ops/s | 0.60x | behind |  |
| BlasL3BreezeJmh | AtA | n=256 | vector | thrpt | 218.664 ± 0.2855 | 355.284 ± 10.16 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.78526e+06 ± 6792 | 1.52531e+06 ± 1.391e+04 | ops/s | 1.17x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | vector | thrpt | 1.78497e+06 ± 8021 | 1.52531e+06 ± 1.391e+04 | ops/s | 1.17x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31616.1 ± 278 | 28276.1 ± 144.2 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | vector | thrpt | 31623.8 ± 342.9 | 28276.1 ± 144.2 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 516.55 ± 3.407 | 1136.12 ± 18.37 | ops/s | 0.45x | behind |  |
| BlasL3BreezeJmh | Gemm | n=256 | vector | thrpt | 1557.69 ± 36.4 | 1136.12 ± 18.37 | ops/s | 1.37x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 435814 ± 1001 | 396669 ± 2486 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | thrpt | 436295 ± 1383 | 396669 ± 2486 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7851.1 ± 70.99 | 12631.4 ± 24.79 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | thrpt | 7859.95 ± 62.62 | 12631.4 ± 24.79 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 128.831 ± 1.679 | 291.51 ± 7.148 | ops/s | 0.44x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | thrpt | 456.798 ± 3.013 | 291.51 ± 7.148 | ops/s | 1.57x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.2037e+06 ± 6090 | 934145 ± 5.425e+04 | ops/s | 1.29x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | vector (gemm-routed only) | thrpt | 1.20716e+06 ± 2713 | 934145 ± 5.425e+04 | ops/s | 1.29x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 30567.8 ± 201.1 | 31530.4 ± 335.5 | ops/s | 0.97x | behind |  |
| DenseDecompositionBreezeJmh | Det | n=64 | vector (gemm-routed only) | thrpt | 30354.9 ± 187.3 | 31530.4 ± 335.5 | ops/s | 0.96x | behind |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 526.754 ± 3.757 | 475.68 ± 2.28 | ops/s | 1.11x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | vector (gemm-routed only) | thrpt | 527.03 ± 3.94 | 475.68 ± 2.28 | ops/s | 1.11x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 288.91 ± 2.797 | 13121.6 ± 149.5 | ops/s | 0.02x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 7.1729 ± 0.1659 | 720.595 ± 5.664 | ops/s | 0.01x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 0.125836 ± 0.002445 | 1.3086 ± 0.0278 | ops/s | 0.10x | behind |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 201560 ± 3021 | 334431 ± 1.078e+04 | ops/s | 0.60x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector (gemm-routed only) | thrpt | 201010 ± 3775 | 334431 ± 1.078e+04 | ops/s | 0.60x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 4147.97 ± 9.062 | 9761.25 ± 46.9 | ops/s | 0.42x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector (gemm-routed only) | thrpt | 4149.44 ± 9.035 | 9761.25 ± 46.9 | ops/s | 0.43x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 52.0448 ± 0.1682 | 164.226 ± 0.07704 | ops/s | 0.32x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector (gemm-routed only) | thrpt | 52.1252 ± 0.07981 | 164.226 ± 0.07704 | ops/s | 0.32x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 40332.8 ± 648 | 17695.5 ± 74.86 | ops/s | 2.28x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 1234.92 ± 11.65 | 1153.39 ± 17.24 | ops/s | 1.07x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 24.6775 ± 0.03105 | 26.5237 ± 0.05716 | ops/s | 0.93x | behind |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 43465.2 ± 433 | 18589.9 ± 135.2 | ops/s | 2.34x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 1391.35 ± 12.63 | 1293.31 ± 10.58 | ops/s | 1.08x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 26.7789 ± 0.02011 | 28.1255 ± 0.04342 | ops/s | 0.95x | behind |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 43263.8 ± 500.2 | 37932.2 ± 190.8 | ops/s | 1.14x | ahead |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 2081.63 ± 27.07 | 1622.34 ± 7.798 | ops/s | 1.28x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 44672 ± 311.1 | 38195.6 ± 387.3 | ops/s | 1.17x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 2071.39 ± 11.81 | 1609.01 ± 7.992 | ops/s | 1.29x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 43953.2 ± 233.6 | 38003 ± 377.5 | ops/s | 1.16x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 2060.71 ± 45.29 | 1614.99 ± 15.02 | ops/s | 1.28x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.45485e+06 ± 1.753e+04 | 464238 ± 2977 | ops/s | 3.13x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | vector (gemm-routed only) | thrpt | 1.44989e+06 ± 2.265e+04 | 464238 ± 2977 | ops/s | 3.12x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 55919.4 ± 125.8 | 26199 ± 92.69 | ops/s | 2.13x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | vector (gemm-routed only) | thrpt | 55983.1 ± 141.9 | 26199 ± 92.69 | ops/s | 2.14x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1069.81 ± 11.54 | 766.661 ± 1.81 | ops/s | 1.40x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | vector (gemm-routed only) | thrpt | 1068.95 ± 13.86 | 766.661 ± 1.81 | ops/s | 1.39x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.24502e+06 ± 1.481e+04 | 1.04966e+06 ± 1.37e+04 | ops/s | 1.19x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | vector (gemm-routed only) | thrpt | 1.23867e+06 ± 8795 | 1.04966e+06 ± 1.37e+04 | ops/s | 1.18x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 30153 ± 465.2 | 32240.1 ± 335.4 | ops/s | 0.94x | behind |  |
| FactorizationBreezeJmh | Lu | n=64 | vector (gemm-routed only) | thrpt | 30545.9 ± 257.3 | 32240.1 ± 335.4 | ops/s | 0.95x | behind |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 525.681 ± 4.666 | 476.062 ± 3.576 | ops/s | 1.10x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | vector (gemm-routed only) | thrpt | 525.427 ± 4.182 | 476.062 ± 3.576 | ops/s | 1.10x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 473678 ± 2127 | 365073 ± 3954 | ops/s | 1.30x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | vector (gemm-routed only) | thrpt | 472432 ± 2315 | 365073 ± 3954 | ops/s | 1.29x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 14230 ± 147.9 | 15508.6 ± 39.87 | ops/s | 0.92x | behind |  |
| FactorizationBreezeJmh | Qr | n=64 | vector (gemm-routed only) | thrpt | 14246.5 ± 152.8 | 15508.6 ± 39.87 | ops/s | 0.92x | behind |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 330.884 ± 2.622 | 230.871 ± 1.521 | ops/s | 1.43x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | vector (gemm-routed only) | thrpt | 335.205 ± 1.202 | 230.871 ± 1.521 | ops/s | 1.45x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 971818 ± 2.6e+04 | 800920 ± 2.627e+04 | ops/s | 1.21x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | vector (gemm-routed only) | thrpt | 956870 ± 1.612e+04 | 800920 ± 2.627e+04 | ops/s | 1.19x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 27628.4 ± 186.3 | 29929 ± 160.8 | ops/s | 0.92x | behind |  |
| FactorizationBreezeJmh | Solve | n=64 | vector (gemm-routed only) | thrpt | 27492.8 ± 256.3 | 29929 ± 160.8 | ops/s | 0.92x | behind |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 513.681 ± 3.71 | 468.968 ± 2.332 | ops/s | 1.10x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | vector (gemm-routed only) | thrpt | 513.318 ± 3.655 | 468.968 ± 2.332 | ops/s | 1.09x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.21269 ± 0.01962 | 9.88425 ± 0.08955 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | vector (gemm-routed only) | avgt | 8.19834 ± 0.02576 | 9.88425 ± 0.08955 | ms/op | 1.21x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.4408 ± 0.3385 | 72.9512 ± 0.1903 | ms/op | 1.02x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | vector (gemm-routed only) | avgt | 71.604 ± 0.2069 | 72.9512 ± 0.1903 | ms/op | 1.02x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 197.623 ± 1.734 | 143.616 ± 0.5454 | ms/op | 0.73x | behind |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 1552.39 ± 4.122 | 1075.32 ± 1.269 | ms/op | 0.69x | behind |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 77.7366 ± 0.3136 | 94.8106 ± 0.02082 | ms/op | 1.22x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | vector (gemm-routed only) | avgt | 68.9715 ± 0.4457 | 94.8106 ± 0.02082 | ms/op | 1.37x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 545.18 ± 5.134 | 751.309 ± 1.054 | ms/op | 1.38x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | vector (gemm-routed only) | avgt | 382.103 ± 20.02 | 751.309 ± 1.054 | ms/op | 1.97x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 15.1338 ± 0.01561 | 16.414 ± 0.02333 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | vector (gemm-routed only) | avgt | 15.0929 ± 0.02022 | 16.414 ± 0.02333 | ms/op | 1.09x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 120.796 ± 0.6409 | 129.077 ± 0.1383 | ms/op | 1.07x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | vector (gemm-routed only) | avgt | 120.513 ± 0.1498 | 129.077 ± 0.1383 | ms/op | 1.07x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 26.2058 ± 0.1655 | 36.8563 ± 0.0867 | ms/op | 1.41x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | vector (gemm-routed only) | avgt | 25.5588 ± 0.07674 | 36.8563 ± 0.0867 | ms/op | 1.44x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 196.134 ± 5.393 | 294.222 ± 0.5981 | ms/op | 1.50x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | vector (gemm-routed only) | avgt | 135.619 ± 3.988 | 294.222 ± 0.5981 | ms/op | 2.17x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 15.3423 ± 0.02027 | 16.5309 ± 0.02399 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | vector (gemm-routed only) | avgt | 15.3436 ± 0.02291 | 16.5309 ± 0.02399 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 121.654 ± 0.1187 | 129.454 ± 0.1443 | ms/op | 1.06x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | vector (gemm-routed only) | avgt | 121.649 ± 0.1824 | 129.454 ± 0.1443 | ms/op | 1.06x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 841.091 ± 4.635 | 558.827 ± 5.258 | us/op | 0.66x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | avgt | 614.878 ± 7.607 | 558.827 ± 5.258 | us/op | 0.91x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1966.83 ± 18.54 | 1056.79 ± 4.398 | us/op | 0.54x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | avgt | 1442.66 ± 16.73 | 1056.79 ± 4.398 | us/op | 0.73x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 103.699 ± 0.7714 | 78.1541 ± 0.6746 | us/op | 0.75x | behind |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 2556.43 ± 6.575 | 2172.66 ± 15.87 | us/op | 0.85x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 82732.9 ± 350.4 | 108043 ± 542.3 | ops/s | 0.77x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector (gemm-routed only) | thrpt | 82607 ± 208.9 | 108043 ± 542.3 | ops/s | 0.76x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2388.36 ± 17.94 | 2287.82 ± 8.034 | ops/s | 1.04x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector (gemm-routed only) | thrpt | 2387.97 ± 19.91 | 2287.82 ± 8.034 | ops/s | 1.04x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 44.5603 ± 0.1908 | 37.4905 ± 0.03244 | ops/s | 1.19x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector (gemm-routed only) | thrpt | 46.5033 ± 1.298 | 37.4905 ± 0.03244 | ops/s | 1.24x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 284.142 ± 1.399 | 103.85 ± 0.4011 | ops/s | 2.74x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 2366.89 ± 2.434 | 2816.96 ± 24.67 | ops/s | 0.84x | behind | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 2783.4 ± 69.23 | 190.219 ± 0.6881 | ops/s | 14.63x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4213.39 ± 9.3 | 357.107 ± 0.06952 | ops/s | 11.80x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 244.58 ± 0.731 | 66.3812 ± 0.2812 | ops/s | 3.68x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7422.14 ± 1.476 | 2144.77 ± 0.7573 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 3581.78 ± 30.18 | 2258.72 ± 16.99 | ops/s | 1.59x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7171.2 ± 11.08 | 387.809 ± 10.77 | ops/s | 18.49x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 216.815 ± 2.327 | 319.674 ± 0.7832 | ops/s | 0.68x | behind |  |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 215.538 ± 0.517 | 319.674 ± 0.7832 | ops/s | 0.67x | behind |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 183.403 ± 0.6392 | 247.231 ± 2.624 | ops/s | 0.74x | behind |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 182.869 ± 1.493 | 247.231 ± 2.624 | ops/s | 0.74x | behind |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.33992e+06 ± 1879 | 1.27523e+06 ± 3961 | ops/s | 2.62x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 53770.9 ± 34.36 | 32551.8 ± 4504 | ops/s | 1.65x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3382.5 ± 12.73 | 2135.87 ± 4.346 | ops/s | 1.58x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 342184 ± 1438 | 340418 ± 1985 | ops/s | 1.01x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5040.67 ± 125.4 | 5021.61 ± 128.1 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 323.673 ± 0.7608 | 216.567 ± 0.4701 | ops/s | 1.49x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 293710 ± 428.4 | 300641 ± 1238 | ops/s | 0.98x | behind |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4730.81 ± 6.515 | 4585.72 ± 32.32 | ops/s | 1.03x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 295.668 ± 1.075 | 282.961 ± 0.9844 | ops/s | 1.04x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 3.49729e+06 ± 1.523e+05 | 5.02372e+06 ± 3.928e+05 | ops/s | 0.70x | behind |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 65751.9 ± 1865 | 92124.9 ± 5509 | ops/s | 0.71x | behind |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4252.94 ± 154.1 | 5732.72 ± 195.3 | ops/s | 0.74x | behind |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.49579e+06 ± 5032 | 245516 ± 145.4 | ops/s | 30.53x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118231 ± 101.6 | 3786.84 ± 2.494 | ops/s | 31.22x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7339.48 ± 24.49 | 237.075 ± 0.1331 | ops/s | 30.96x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.35836e+06 ± 1430 | 107542 ± 27.25 | ops/s | 40.53x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68511.2 ± 55.43 | 1673.22 ± 3.122 | ops/s | 40.95x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4268.07 ± 0.6357 | 104.297 ± 0.04968 | ops/s | 40.92x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.86217e+06 ± 2.162e+04 | 2.3204e+06 ± 3500 | ops/s | 1.23x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55528.9 ± 39.33 | 34160.6 ± 35.05 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3498.33 ± 1.709 | 2130.79 ± 1.487 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 247727 ± 1085 | 290176 ± 2209 | ops/s | 0.85x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 2277.83 ± 65.04 | 4322.63 ± 93.18 | ops/s | 0.53x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 132.566 ± 2.664 | 185.654 ± 1.186 | ops/s | 0.71x | behind |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 263421 ± 1035 | 148193 ± 306.9 | ops/s | 1.78x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3961.7 ± 50.03 | 2173.31 ± 259 | ops/s | 1.82x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 251.474 ± 0.8846 | 117.02 ± 0.4938 | ops/s | 2.15x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.66445e+06 ± 5478 | 2.34442e+06 ± 639.4 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 123043 ± 696.9 | 34376.8 ± 32.79 | ops/s | 3.58x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7241.27 ± 21.49 | 2146.1 ± 0.497 | ops/s | 3.37x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.71407e+07 ± 3.274e+05 | 9.18114e+06 ± 9568 | ops/s | 2.96x | ahead |  |
| SmallDenseBreezeJmh | Det | n=3 | vector (gemm-routed only) | thrpt | 2.74974e+07 ± 1.041e+05 | 9.18114e+06 ± 9568 | ops/s | 2.99x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.86123e+07 ± 4.727e+05 | 7.54757e+06 ± 1.099e+05 | ops/s | 2.47x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | vector (gemm-routed only) | thrpt | 1.88215e+07 ± 4.786e+05 | 7.54757e+06 ± 1.099e+05 | ops/s | 2.49x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.13306e+07 ± 7.259e+04 | 2.95426e+07 ± 6.075e+05 | ops/s | 1.06x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | vector | thrpt | 3.09415e+07 ± 2.422e+05 | 2.95426e+07 ± 6.075e+05 | ops/s | 1.05x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.75141e+07 ± 1.854e+06 | 2.0309e+07 ± 5.051e+04 | ops/s | 1.85x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | vector | thrpt | 3.86284e+07 ± 1.773e+05 | 2.0309e+07 ± 5.051e+04 | ops/s | 1.90x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 1.03477e+07 ± 1.328e+05 | 3.51465e+06 ± 1.566e+04 | ops/s | 2.94x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=3 | vector (gemm-routed only) | thrpt | 1.04586e+07 ± 2.238e+04 | 3.51465e+06 ± 1.566e+04 | ops/s | 2.98x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 6.07367e+06 ± 7.476e+04 | 2.80229e+06 ± 9.174e+04 | ops/s | 2.17x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=4 | vector (gemm-routed only) | thrpt | 6.07146e+06 ± 1.742e+05 | 2.80229e+06 ± 9.174e+04 | ops/s | 2.17x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.54391e+07 ± 5.538e+04 | 5.92901e+06 ± 1.338e+05 | ops/s | 2.60x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | vector (gemm-routed only) | thrpt | 1.55864e+07 ± 2.195e+04 | 5.92901e+06 ± 1.338e+05 | ops/s | 2.63x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.10649e+07 ± 9.035e+04 | 5.12372e+06 ± 1.362e+05 | ops/s | 2.16x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | vector (gemm-routed only) | thrpt | 1.11912e+07 ± 1.325e+05 | 5.12372e+06 ± 1.362e+05 | ops/s | 2.18x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.230039 ± 0.0003179 | 1.28141 ± 0.007354 | ms/op | 5.57x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.5708 ± 0.09491 | 122.453 ± 0.1651 | ms/op | 4.98x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.8528 ± 0.005986 | 12.3081 ± 0.01701 | ms/op | 6.64x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 221.526 ± 0.7678 | 1231.41 ± 55.75 | ms/op | 5.56x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00578383 ± 1.211e-05 | 0.0404691 ± 0.0002067 | ms/op | 7.00x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.594982 ± 0.00108 | 3.95564 ± 0.2425 | ms/op | 6.65x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0462841 ± 0.0002374 | 0.374016 ± 0.009222 | ms/op | 8.08x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.1767 ± 0.005013 | 37.8557 ± 0.06841 | ms/op | 6.13x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140824 ± 0.0004131 | 1.29845 ± 0.06676 | ms/op | 9.22x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.4937 ± 0.2682 | 122.384 ± 0.1118 | ms/op | 8.44x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19848 ± 0.007478 | 12.2794 ± 0.02028 | ms/op | 10.25x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 148.769 ± 2.842 | 1231.58 ± 51.48 | ms/op | 8.28x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00498109 ± 2.459e-05 | 0.0403379 ± 7.067e-05 | ms/op | 8.10x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.546417 ± 0.001226 | 3.82004 ± 0.002585 | ms/op | 6.99x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0530993 ± 0.0003278 | 0.372412 ± 0.0007137 | ms/op | 7.01x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.63859 ± 0.01553 | 37.8959 ± 0.02394 | ms/op | 6.72x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 210434 ± 678 | 198427 ± 1114 | ops/s | 1.06x | ahead |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 17733.5 ± 45.73 | 16200.4 ± 297.4 | ops/s | 1.09x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.08047e+06 ± 1.077e+04 | 314534 ± 1441 | ops/s | 6.61x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 208739 ± 984.8 | 32542.8 ± 1238 | ops/s | 6.41x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.13104e+06 ± 6.587e+04 | 324329 ± 433.8 | ops/s | 3.49x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 40658.8 ± 2.364e+04 | 35106.2 ± 105.9 | ops/s | 1.16x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.26326e+06 ± 829.6 | 2.25101e+06 ± 3165 | ops/s | 1.01x | ahead |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 223786 ± 43.82 | 223571 ± 243 | ops/s | 1.00x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 76430.7 ± 2458 | 33957.4 ± 598 | ops/s | 2.25x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 2334.65 ± 12.68 | 1766.45 ± 31.19 | ops/s | 1.32x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 335.811 ± 1.339 | 312.994 ± 2.921 | ops/s | 1.07x | ahead |  |

## Delta against the d102321 baseline

The baseline is `2026-10-09-breeze-two-lane-baseline-trillium-laneB.md` (source 17813d0, committed in d102321). It used the same lane, cluster, CPU model, JDK, pinning and method, on a different node of the same type.

Ratios are gale speed / Breeze speed (>1 means gale is faster). `change` is new ratio / old ratio. Old is baseline 17813d0 (d102321); new is milestone 1 ff03eed.

| class | op | params | gale backend | old ratio | new ratio | change | old verdict | new verdict |
|---|---|---|---|---:|---:|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | 0.43x | 0.92x | 2.14x | behind | tie |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | 0.44x | 0.98x | 2.23x | behind | tie |
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | 0.32x | 0.98x | 3.06x | behind | tie |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | 0.96x | 0.98x | 1.02x | behind | behind |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | 0.99x | 0.98x | 0.99x | tie | behind |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | 0.90x | 0.95x | 1.06x | behind | behind |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | 2.94x | 11.81x | 4.02x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | 3.39x | 13.62x | 4.02x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | 3.58x | 14.21x | 3.97x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 0.32x | 0.51x | 1.59x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | 0.97x | 0.99x | 1.02x | behind | tie |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 0.36x | 0.58x | 1.61x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | 1.01x | 0.97x | 0.96x | tie | behind |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 0.27x | 0.42x | 1.56x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=256 | vector | 1.04x | 0.96x | 0.92x | ahead | behind |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 0.18x | 0.43x | 2.39x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | 0.01x | 0.96x | 96.00x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 0.20x | 0.48x | 2.40x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | 0.01x | 0.93x | 93.00x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 0.18x | 0.34x | 1.89x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=256 | vector | 0.06x | 0.92x | 15.33x | behind | behind |
| BlasL3BreezeJmh | AtA | n=16 | pure | 2.00x | 2.00x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=16 | vector | 2.00x | 1.99x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 0.62x | 0.60x | 0.97x | behind | behind |
| BlasL3BreezeJmh | AtA | n=256 | vector | 0.60x | 0.62x | 1.03x | behind | behind |
| BlasL3BreezeJmh | AtA | n=64 | pure | 1.06x | 1.06x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=64 | vector | 1.07x | 1.06x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.20x | 1.17x | 0.97x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | vector | 1.19x | 1.17x | 0.98x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 0.46x | 0.45x | 0.98x | behind | behind |
| BlasL3BreezeJmh | Gemm | n=256 | vector | 1.37x | 1.37x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | vector | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.10x | 1.10x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | 1.10x | 1.10x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 0.44x | 0.44x | 1.00x | behind | behind |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | 1.58x | 1.57x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 0.62x | 0.62x | 1.00x | behind | behind |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | 0.62x | 0.62x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.27x | 1.29x | 1.02x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | vector (gemm-routed only) | 1.28x | 1.29x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.11x | 1.11x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | vector (gemm-routed only) | 1.11x | 1.11x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 0.97x | 0.97x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Det | n=64 | vector (gemm-routed only) | 0.96x | 0.96x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | 0.02x | 0.02x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | 0.10x | 0.10x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | 0.01x | 0.01x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 0.60x | 0.60x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector (gemm-routed only) | 0.60x | 0.60x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 0.32x | 0.32x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector (gemm-routed only) | 0.32x | 0.32x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 0.42x | 0.42x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector (gemm-routed only) | 0.42x | 0.43x | 1.02x | behind | behind |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | 1.76x | 2.28x | 1.30x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | 0.11x | 0.93x | 8.45x | behind | behind |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | 0.66x | 1.07x | 1.62x | behind | ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | 1.94x | 2.34x | 1.21x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | 0.11x | 0.95x | 8.64x | behind | behind |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | 0.74x | 1.08x | 1.46x | behind | ahead |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | 0.45x | 1.28x | 2.84x | behind | ahead |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | 0.86x | 1.14x | 1.33x | behind | ahead |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | 0.76x | 1.29x | 1.70x | behind | ahead |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | 0.41x | 1.17x | 2.85x | behind | ahead |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | 0.43x | 1.28x | 2.98x | behind | ahead |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | 0.82x | 1.16x | 1.41x | behind | ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.14x | 3.13x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=16 | vector (gemm-routed only) | 3.13x | 3.12x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.39x | 1.40x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | vector (gemm-routed only) | 1.39x | 1.39x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.12x | 2.13x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | vector (gemm-routed only) | 2.12x | 2.14x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.16x | 1.19x | 1.03x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | vector (gemm-routed only) | 1.15x | 1.18x | 1.03x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.10x | 1.10x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | vector (gemm-routed only) | 1.11x | 1.10x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 0.94x | 0.94x | 1.00x | behind | behind |
| FactorizationBreezeJmh | Lu | n=64 | vector (gemm-routed only) | 0.94x | 0.95x | 1.01x | behind | behind |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.25x | 1.30x | 1.04x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=16 | vector (gemm-routed only) | 1.26x | 1.29x | 1.02x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.44x | 1.43x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | vector (gemm-routed only) | 1.45x | 1.45x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 0.91x | 0.92x | 1.01x | behind | behind |
| FactorizationBreezeJmh | Qr | n=64 | vector (gemm-routed only) | 0.91x | 0.92x | 1.01x | behind | behind |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.20x | 1.21x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=16 | vector (gemm-routed only) | 1.18x | 1.19x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.10x | 1.10x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | vector (gemm-routed only) | 1.09x | 1.09x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 0.92x | 0.92x | 1.00x | behind | behind |
| FactorizationBreezeJmh | Solve | n=64 | vector (gemm-routed only) | 0.93x | 0.92x | 0.99x | behind | behind |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | 1.01x | 1.02x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=1024 | vector (gemm-routed only) | 1.01x | 1.02x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | 1.20x | 1.20x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | vector (gemm-routed only) | 1.20x | 1.21x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | 0.04x | 0.69x | 17.25x | behind | behind |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | 0.04x | 0.73x | 18.25x | behind | behind |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | 1.43x | 1.38x | 0.97x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | vector (gemm-routed only) | 2.08x | 1.97x | 0.95x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | 1.22x | 1.22x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | vector (gemm-routed only) | 1.36x | 1.37x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | 1.07x | 1.07x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | vector (gemm-routed only) | 1.07x | 1.07x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | 1.09x | 1.08x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | vector (gemm-routed only) | 1.09x | 1.09x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.55x | 1.50x | 0.97x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | vector (gemm-routed only) | 2.15x | 2.17x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.40x | 1.41x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | vector (gemm-routed only) | 1.46x | 1.44x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | 1.06x | 1.06x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | vector (gemm-routed only) | 1.06x | 1.06x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | 1.08x | 1.08x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | vector (gemm-routed only) | 1.08x | 1.08x | 1.00x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | 0.48x | 0.66x | 1.38x | behind | behind |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | 0.37x | 0.91x | 2.46x | behind | behind |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | 0.40x | 0.54x | 1.35x | withheld | withheld |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | 0.30x | 0.73x | 2.43x | withheld | withheld |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | 0.75x | 0.75x | 1.00x | behind | behind |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | 0.80x | 0.85x | 1.06x | withheld | withheld |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.76x | 0.77x | 1.01x | behind | behind |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector (gemm-routed only) | 0.76x | 0.76x | 1.00x | behind | behind |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.17x | 1.19x | 1.02x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector (gemm-routed only) | 1.24x | 1.24x | 1.00x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.02x | 1.04x | 1.02x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector (gemm-routed only) | 1.02x | 1.04x | 1.02x | ahead | ahead |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | 2.65x | 2.74x | 1.03x | ahead | ahead |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | 0.81x | 0.84x | 1.04x | behind | behind |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | 11.04x | 14.63x | 1.33x | ahead | ahead |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | 11.61x | 11.80x | 1.02x | ahead | ahead |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | 3.04x | 3.68x | 1.21x | ahead | ahead |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.46x | 3.46x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | 1.84x | 1.59x | 0.86x | ahead | ahead |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | 17.94x | 18.49x | 1.03x | ahead | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | 0.66x | 0.68x | 1.03x | behind | behind |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | vector (gemm-routed only) | 0.66x | 0.67x | 1.02x | behind | behind |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | 0.74x | 0.74x | 1.00x | behind | behind |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | vector (gemm-routed only) | 0.74x | 0.74x | 1.00x | behind | behind |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | 2.64x | 2.62x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | 1.59x | 1.58x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | 1.66x | 1.65x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | 0.99x | 1.01x | 1.02x | behind | tie |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | 1.47x | 1.49x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | 0.98x | 1.00x | 1.02x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | 0.96x | 0.98x | 1.02x | behind | behind |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | 1.03x | 1.04x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | 1.01x | 1.03x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | 0.70x | 0.70x | 1.00x | behind | behind |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | 0.60x | 0.74x | 1.23x | behind | behind |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | 0.59x | 0.71x | 1.20x | behind | behind |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | 30.56x | 30.53x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | 30.93x | 30.96x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | 31.21x | 31.22x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | 40.50x | 40.53x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | 40.95x | 40.92x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | 40.96x | 40.95x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | 1.23x | 1.23x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | 1.64x | 1.64x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | 1.63x | 1.63x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | 0.89x | 0.85x | 0.96x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | 0.75x | 0.71x | 0.95x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | 0.54x | 0.53x | 0.98x | behind | behind |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | 1.47x | 1.78x | 1.21x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | 1.76x | 2.15x | 1.22x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | 1.51x | 1.82x | 1.21x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.25x | 3.27x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | 3.38x | 3.37x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | 3.58x | 3.58x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | pure | 3.02x | 2.96x | 0.98x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | vector (gemm-routed only) | 3.04x | 2.99x | 0.98x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | pure | 2.51x | 2.47x | 0.98x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | vector (gemm-routed only) | 2.47x | 2.49x | 1.01x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | 1.05x | 1.06x | 1.01x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | vector | 1.04x | 1.05x | 1.01x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | 1.91x | 1.85x | 0.97x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | vector | 1.90x | 1.90x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | pure | 3.06x | 2.94x | 0.96x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | vector (gemm-routed only) | 3.07x | 2.98x | 0.97x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | pure | 2.27x | 2.17x | 0.96x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | vector (gemm-routed only) | 2.27x | 2.17x | 0.96x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | pure | 2.57x | 2.60x | 1.01x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | vector (gemm-routed only) | 2.58x | 2.63x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | pure | 2.16x | 2.16x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | vector (gemm-routed only) | 2.17x | 2.18x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | 5.56x | 5.57x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | 5.02x | 4.98x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | 6.65x | 6.64x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | 5.51x | 5.56x | 1.01x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | 6.99x | 7.00x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | 6.42x | 6.65x | 1.04x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | 8.10x | 8.08x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | 6.17x | 6.13x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | 9.10x | 9.22x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | 8.33x | 8.44x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | 10.25x | 10.25x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | 8.42x | 8.28x | 0.98x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | 8.18x | 8.10x | 0.99x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | 7.06x | 6.99x | 0.99x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | 7.12x | 7.01x | 0.98x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | 6.72x | 6.72x | 1.00x | withheld | withheld |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | 0.88x | 1.06x | 1.20x | behind | ahead |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | 1.02x | 1.09x | 1.07x | tie | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | 6.59x | 6.61x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | 6.37x | 6.41x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | 3.42x | 3.49x | 1.02x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | 1.17x | 1.16x | 0.99x | tie | tie |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | 1.01x | 1.01x | 1.00x | tie | ahead |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | tie |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | 0.54x | 1.07x | 1.98x | behind | ahead |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | 2.07x | 2.25x | 1.09x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | 0.95x | 1.32x | 1.39x | behind | ahead |

## Provenance

- Date: 2026-10-09, Slurm job 2526675, start/end/elapsed: 2026-10-09T16:11:55 2026-10-09T18:25:08 02:13:13 
- Benchmarked source: `ff03eed494bfa2d246b1c363c59f3425d97fcd67` (breeze/integration). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster.
- Node: SciNet trillium `tri0017`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5), AVX-512, SMT off, 192 CPUs, NPS4 (8 NUMA nodes of 24 cores), boost enabled. Load average at job start: 176.33, 149.11, 75.37.
- JDK: Temurin 25+36 (`module load java/25`). The site default `JAVA_TOOL_OPTIONS=-Xmx2g` was unset before launch.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the CCD of cpu0 and its NUMA node). It applies to the host JVM and is inherited by every fork.
- Netlib (blas, lapack, vectorModule) over all sidecar records (`2026-10-09-breeze-milestone1-trillium-laneB.netlib.jsonl`): [('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)]. No native BLAS/LAPACK is visible to the loader on these nodes.
- Command, equivalent to `breezeLaneB`: `numactl … java --add-modules=jdk.incubator.vector -cp <classpath> org.openjdk.jmh.Main -jvmArgsAppend "-Xms4g -Xmx4g -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure,vector -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings come from the annotations. The runner is `docs/verification/w21-simd-spike/x86-avx512/run.sh` (with `BUNDLE` pointing at the milestone bundle).
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane B --strict`. It exited 0, with 0 unpaired rows.
