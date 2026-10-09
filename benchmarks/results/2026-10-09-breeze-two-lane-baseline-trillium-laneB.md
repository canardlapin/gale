# gale vs Breeze scoreboard — lane B (SIMD)

- Lane: B
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: --add-modules=jdk.incubator.vector -Xms4g -Xmx4g -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/out/laneB-2524301/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.VectorBLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `17813d0d761467187ba0764f1081fd30747ad18e`
- Machine: SciNet trillium tri0732: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512, 8 double lanes), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2524301 exclusive

**120 ahead, 7 tie, 71 behind, 0 n/a, 11 withheld** of 209 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 54290.6 ± 1740 | 170063 ± 5389 | ops/s | 0.32x | behind |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 13097.6 ± 1551 | 30000.5 ± 2736 | ops/s | 0.44x | behind |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 3223.1 ± 473.8 | 7535.46 ± 767.8 | ops/s | 0.43x | behind |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 61311.2 ± 205.7 | 68070.4 ± 610.6 | ops/s | 0.90x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16482.5 ± 285 | 16673.1 ± 723.4 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4104.22 ± 130.5 | 4270.93 ± 6.507 | ops/s | 0.96x | behind |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 16725.6 ± 50.67 | 4667.32 ± 9.028 | ops/s | 3.58x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 4192.92 ± 9.907 | 1238.41 ± 282.5 | ops/s | 3.39x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 1054.94 ± 3.424 | 358.63 ± 0.1208 | ops/s | 2.94x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 64663.1 ± 183.6 | 235212 ± 1204 | ops/s | 0.27x | behind |  |
| BlasL2BreezeJmh | Gemv | n=256 | vector | thrpt | 244544 ± 7432 | 235212 ± 1204 | ops/s | 1.04x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 4149.52 ± 1.732 | 13091.9 ± 13.81 | ops/s | 0.32x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | thrpt | 12654.6 ± 172 | 13091.9 ± 13.81 | ops/s | 0.97x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1026.14 ± 22.5 | 2816.01 ± 42.6 | ops/s | 0.36x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | thrpt | 2845.16 ± 5.306 | 2816.01 ± 42.6 | ops/s | 1.01x | tie |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 46537.6 ± 123 | 257148 ± 1.852e+04 | ops/s | 0.18x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | vector | thrpt | 14736.3 ± 160.5 | 257148 ± 1.852e+04 | ops/s | 0.06x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 2416.7 ± 2.063 | 13127.6 ± 142.6 | ops/s | 0.18x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | thrpt | 188.724 ± 3.094 | 13127.6 ± 142.6 | ops/s | 0.01x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 601.926 ± 0.3351 | 2970.79 ± 16.47 | ops/s | 0.20x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | thrpt | 32.4651 ± 0.2676 | 2970.79 ± 16.47 | ops/s | 0.01x | behind |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 788760 ± 4768 | 395365 ± 6831 | ops/s | 2.00x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | vector | thrpt | 791351 ± 1839 | 395365 ± 6831 | ops/s | 2.00x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15282.5 ± 96.27 | 14411.6 ± 25.31 | ops/s | 1.06x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | vector | thrpt | 15359.7 ± 136.8 | 14411.6 ± 25.31 | ops/s | 1.07x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 216.799 ± 7.856 | 350.873 ± 11.66 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | AtA | n=256 | vector | thrpt | 210.738 ± 2.419 | 350.873 ± 11.66 | ops/s | 0.60x | behind |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.79656e+06 ± 6000 | 1.49783e+06 ± 1.793e+04 | ops/s | 1.20x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | vector | thrpt | 1.78962e+06 ± 9894 | 1.49783e+06 ± 1.793e+04 | ops/s | 1.19x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31772.2 ± 304.3 | 28416.5 ± 158.1 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | vector | thrpt | 31803.1 ± 353.1 | 28416.5 ± 158.1 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 518.945 ± 4.612 | 1140.51 ± 19.08 | ops/s | 0.46x | behind |  |
| BlasL3BreezeJmh | Gemm | n=256 | vector | thrpt | 1563.53 ± 35.67 | 1140.51 ± 19.08 | ops/s | 1.37x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 438040 ± 972.5 | 397180 ± 1834 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | thrpt | 438627 ± 1278 | 397180 ± 1834 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7901.41 ± 57.98 | 12676.9 ± 23.15 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | thrpt | 7901.97 ± 56.03 | 12676.9 ± 23.15 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 129.417 ± 1.592 | 290.969 ± 5.033 | ops/s | 0.44x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | thrpt | 458.91 ± 2.327 | 290.969 ± 5.033 | ops/s | 1.58x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.20743e+06 ± 1.527e+04 | 950045 ± 1.012e+04 | ops/s | 1.27x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | vector (gemm-routed only) | thrpt | 1.21277e+06 ± 5889 | 950045 ± 1.012e+04 | ops/s | 1.28x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 30854.8 ± 168.1 | 31683.7 ± 310.1 | ops/s | 0.97x | behind |  |
| DenseDecompositionBreezeJmh | Det | n=64 | vector (gemm-routed only) | thrpt | 30259.1 ± 396.2 | 31683.7 ± 310.1 | ops/s | 0.96x | behind |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 529.699 ± 3.786 | 477.271 ± 2.674 | ops/s | 1.11x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | vector (gemm-routed only) | thrpt | 529.305 ± 3.69 | 477.271 ± 2.674 | ops/s | 1.11x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 293.713 ± 2.273 | 13105.2 ± 64.05 | ops/s | 0.02x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 7.13966 ± 0.181 | 724.147 ± 6.621 | ops/s | 0.01x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 0.127147 ± 0.001859 | 1.30877 ± 0.005049 | ops/s | 0.10x | behind |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 202357 ± 2941 | 336060 ± 1.155e+04 | ops/s | 0.60x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector (gemm-routed only) | thrpt | 201789 ± 2480 | 336060 ± 1.155e+04 | ops/s | 0.60x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 4159.29 ± 11.25 | 9825.59 ± 28.38 | ops/s | 0.42x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector (gemm-routed only) | thrpt | 4168.27 ± 9.359 | 9825.59 ± 28.38 | ops/s | 0.42x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 52.3171 ± 0.06258 | 165.058 ± 0.04533 | ops/s | 0.32x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector (gemm-routed only) | thrpt | 52.2861 ± 0.1961 | 165.058 ± 0.04533 | ops/s | 0.32x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 31139.6 ± 1156 | 17692.8 ± 107.2 | ops/s | 1.76x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 760.88 ± 5.367 | 1155.95 ± 14.65 | ops/s | 0.66x | behind |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 3.06425 ± 0.02611 | 26.6617 ± 0.07948 | ops/s | 0.11x | behind |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 36205.7 ± 1384 | 18655.9 ± 96.14 | ops/s | 1.94x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 965.309 ± 11.24 | 1296.32 ± 9.864 | ops/s | 0.74x | behind |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 3.24444 ± 0.01992 | 28.2712 ± 0.0482 | ops/s | 0.11x | behind |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 32654.1 ± 904.9 | 38139.8 ± 392 | ops/s | 0.86x | behind |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 718.287 ± 1.35 | 1610.67 ± 12.8 | ops/s | 0.45x | behind |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 15709.6 ± 874.8 | 38003.3 ± 345.3 | ops/s | 0.41x | behind |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 1232.36 ± 34.84 | 1627.28 ± 9.919 | ops/s | 0.76x | behind |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 30962.9 ± 1212 | 37946.9 ± 452.1 | ops/s | 0.82x | behind |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 701.397 ± 22.29 | 1620.44 ± 23.1 | ops/s | 0.43x | behind |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.46525e+06 ± 1.306e+04 | 467348 ± 2999 | ops/s | 3.14x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | vector (gemm-routed only) | thrpt | 1.4619e+06 ± 9835 | 467348 ± 2999 | ops/s | 3.13x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 56174.7 ± 153.1 | 26473.7 ± 618.7 | ops/s | 2.12x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | vector (gemm-routed only) | thrpt | 56197.6 ± 132 | 26473.7 ± 618.7 | ops/s | 2.12x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1075.08 ± 12.58 | 770.922 ± 3.043 | ops/s | 1.39x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | vector (gemm-routed only) | thrpt | 1074.57 ± 12.63 | 770.922 ± 3.043 | ops/s | 1.39x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.23163e+06 ± 2.602e+04 | 1.06431e+06 ± 1.303e+04 | ops/s | 1.16x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | vector (gemm-routed only) | thrpt | 1.22762e+06 ± 2.175e+04 | 1.06431e+06 ± 1.303e+04 | ops/s | 1.15x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 30431 ± 572 | 32361.3 ± 384.4 | ops/s | 0.94x | behind |  |
| FactorizationBreezeJmh | Lu | n=64 | vector (gemm-routed only) | thrpt | 30275.1 ± 449 | 32361.3 ± 384.4 | ops/s | 0.94x | behind |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 526.957 ± 3.983 | 477.418 ± 2.344 | ops/s | 1.10x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | vector (gemm-routed only) | thrpt | 527.703 ± 3.247 | 477.418 ± 2.344 | ops/s | 1.11x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 462631 ± 2785 | 370258 ± 9049 | ops/s | 1.25x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | vector (gemm-routed only) | thrpt | 466127 ± 4091 | 370258 ± 9049 | ops/s | 1.26x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 14269.6 ± 137.1 | 15605 ± 47.89 | ops/s | 0.91x | behind |  |
| FactorizationBreezeJmh | Qr | n=64 | vector (gemm-routed only) | thrpt | 14273.8 ± 151.6 | 15605 ± 47.89 | ops/s | 0.91x | behind |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 332.434 ± 2.136 | 231.594 ± 1.514 | ops/s | 1.44x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | vector (gemm-routed only) | thrpt | 334.782 ± 2.356 | 231.594 ± 1.514 | ops/s | 1.45x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 982669 ± 9889 | 819016 ± 4112 | ops/s | 1.20x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | vector (gemm-routed only) | thrpt | 964074 ± 5451 | 819016 ± 4112 | ops/s | 1.18x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 27736.3 ± 203.5 | 30078 ± 230.7 | ops/s | 0.92x | behind |  |
| FactorizationBreezeJmh | Solve | n=64 | vector (gemm-routed only) | thrpt | 27826.8 ± 243.8 | 30078 ± 230.7 | ops/s | 0.93x | behind |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 516.402 ± 3.654 | 470.939 ± 2.794 | ops/s | 1.10x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | vector (gemm-routed only) | thrpt | 515.266 ± 4.288 | 470.939 ± 2.794 | ops/s | 1.09x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.17854 ± 0.004431 | 9.83435 ± 0.08992 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | vector (gemm-routed only) | avgt | 8.16352 ± 0.01779 | 9.83435 ± 0.08992 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.3188 ± 0.183 | 72.1026 ± 0.3753 | ms/op | 1.01x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | vector (gemm-routed only) | avgt | 71.3668 ± 0.1918 | 72.1026 ± 0.3753 | ms/op | 1.01x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 3231.32 ± 50.01 | 143.055 ± 0.391 | ms/op | 0.04x | behind |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 28179.2 ± 851.5 | 1068.53 ± 0.8438 | ms/op | 0.04x | behind |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 77.2762 ± 0.373 | 94.3705 ± 0.07346 | ms/op | 1.22x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | vector (gemm-routed only) | avgt | 69.4804 ± 1.763 | 94.3705 ± 0.07346 | ms/op | 1.36x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 521.941 ± 9.679 | 747.491 ± 0.4967 | ms/op | 1.43x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | vector (gemm-routed only) | avgt | 360.036 ± 8.213 | 747.491 ± 0.4967 | ms/op | 2.08x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 15.0551 ± 0.01691 | 16.3372 ± 0.02213 | ms/op | 1.09x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | vector (gemm-routed only) | avgt | 15.0213 ± 0.01855 | 16.3372 ± 0.02213 | ms/op | 1.09x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 120.091 ± 0.8363 | 128.343 ± 0.1538 | ms/op | 1.07x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | vector (gemm-routed only) | avgt | 119.851 ± 0.07905 | 128.343 ± 0.1538 | ms/op | 1.07x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 26.1975 ± 0.1659 | 36.6783 ± 0.1231 | ms/op | 1.40x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | vector (gemm-routed only) | avgt | 25.1388 ± 0.07442 | 36.6783 ± 0.1231 | ms/op | 1.46x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 188.423 ± 2.411 | 292.889 ± 0.4931 | ms/op | 1.55x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | vector (gemm-routed only) | avgt | 135.991 ± 16.7 | 292.889 ± 0.4931 | ms/op | 2.15x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 15.2676 ± 0.01583 | 16.4562 ± 0.02972 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | vector (gemm-routed only) | avgt | 15.2726 ± 0.009201 | 16.4562 ± 0.02972 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 120.996 ± 0.1139 | 128.641 ± 0.1572 | ms/op | 1.06x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | vector (gemm-routed only) | avgt | 120.944 ± 0.1336 | 128.641 ± 0.1572 | ms/op | 1.06x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 1123.07 ± 19.76 | 542.692 ± 3.322 | us/op | 0.48x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | avgt | 1478 ± 6.224 | 542.692 ± 3.322 | us/op | 0.37x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 2602.98 ± 26.59 | 1045.24 ± 2.835 | us/op | 0.40x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | avgt | 3454.06 ± 9.409 | 1045.24 ± 2.835 | us/op | 0.30x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 103.982 ± 1.3 | 78.1037 ± 0.7904 | us/op | 0.75x | behind |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 2550.5 ± 14.38 | 2028.27 ± 13.39 | us/op | 0.80x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 82379.3 ± 193.9 | 108204 ± 895.5 | ops/s | 0.76x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector (gemm-routed only) | thrpt | 82094.9 ± 254 | 108204 ± 895.5 | ops/s | 0.76x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2347.74 ± 18.67 | 2294.86 ± 7.199 | ops/s | 1.02x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector (gemm-routed only) | thrpt | 2337.44 ± 29.05 | 2294.86 ± 7.199 | ops/s | 1.02x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 44.6873 ± 0.4079 | 38.1828 ± 0.8562 | ops/s | 1.17x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector (gemm-routed only) | thrpt | 47.2147 ± 1.352 | 38.1828 ± 0.8562 | ops/s | 1.24x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 275.483 ± 0.2063 | 104.111 ± 0.2047 | ops/s | 2.65x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 2317.96 ± 12.12 | 2866.05 ± 4.505 | ops/s | 0.81x | behind | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 2096.96 ± 156.1 | 189.965 ± 3.472 | ops/s | 11.04x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4164.08 ± 119.6 | 358.619 ± 0.07076 | ops/s | 11.61x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 201.459 ± 0.4449 | 66.3399 ± 0.676 | ops/s | 3.04x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7453.04 ± 4.117 | 2153.41 ± 0.6173 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 4156.25 ± 12.6 | 2258.45 ± 0.7565 | ops/s | 1.84x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7210.05 ± 2.412 | 401.957 ± 13.21 | ops/s | 17.94x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 216.571 ± 0.6721 | 326.37 ± 7.838 | ops/s | 0.66x | behind |  |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 216.137 ± 0.4021 | 326.37 ± 7.838 | ops/s | 0.66x | behind |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 184.215 ± 1.006 | 248.301 ± 2.464 | ops/s | 0.74x | behind |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 183.323 ± 1.161 | 248.301 ± 2.464 | ops/s | 0.74x | behind |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.39123e+06 ± 6.111e+04 | 1.28464e+06 ± 743.4 | ops/s | 2.64x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 54031.6 ± 61.36 | 32615.2 ± 4842 | ops/s | 1.66x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3403.69 ± 3.799 | 2146.87 ± 3.111 | ops/s | 1.59x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 337590 ± 2366 | 341703 ± 913 | ops/s | 0.99x | behind |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 4994.74 ± 86.96 | 5112.46 ± 59.28 | ops/s | 0.98x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 318.682 ± 0.6607 | 217.023 ± 1.883 | ops/s | 1.47x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 291832 ± 342.2 | 302799 ± 839.8 | ops/s | 0.96x | behind |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4678.3 ± 3.135 | 4623.26 ± 10.2 | ops/s | 1.01x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 292.472 ± 0.1424 | 284.663 ± 0.8715 | ops/s | 1.03x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 3.42713e+06 ± 8945 | 4.92359e+06 ± 3.375e+05 | ops/s | 0.70x | behind |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 53938 ± 25.92 | 91947 ± 4828 | ops/s | 0.59x | behind |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 3403.37 ± 1.412 | 5666.39 ± 19.92 | ops/s | 0.60x | behind |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.52907e+06 ± 1632 | 246380 ± 40.18 | ops/s | 30.56x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118593 ± 27.32 | 3799.62 ± 0.8097 | ops/s | 31.21x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7359.02 ± 9.741 | 237.936 ± 0.07566 | ops/s | 30.93x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.3731e+06 ± 9265 | 107989 ± 19.82 | ops/s | 40.50x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68794.6 ± 14.87 | 1679.57 ± 3.386 | ops/s | 40.96x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4287.89 ± 2.203 | 104.715 ± 0.01894 | ops/s | 40.95x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.87298e+06 ± 1.867e+04 | 2.3293e+06 ± 1101 | ops/s | 1.23x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55834.2 ± 62.66 | 34257.4 ± 29.89 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3516.08 ± 1.165 | 2138.31 ± 0.319 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 259841 ± 714.2 | 292081 ± 1169 | ops/s | 0.89x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 2346.65 ± 26.57 | 4354.37 ± 79.76 | ops/s | 0.54x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 131.627 ± 2.587 | 174.986 ± 17.31 | ops/s | 0.75x | behind |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 218534 ± 855.1 | 148760 ± 323.5 | ops/s | 1.47x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3294.45 ± 51.89 | 2188.2 ± 264 | ops/s | 1.51x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 207.517 ± 0.5822 | 117.73 ± 0.4218 | ops/s | 1.76x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.64862e+06 ± 7.295e+04 | 2.35121e+06 ± 699.5 | ops/s | 3.25x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 123300 ± 517 | 34481.7 ± 32.2 | ops/s | 3.58x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7270.98 ± 5.947 | 2153.02 ± 0.6126 | ops/s | 3.38x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.74167e+07 ± 1.5e+05 | 9.09022e+06 ± 3.07e+05 | ops/s | 3.02x | ahead |  |
| SmallDenseBreezeJmh | Det | n=3 | vector (gemm-routed only) | thrpt | 2.75995e+07 ± 1.364e+05 | 9.09022e+06 ± 3.07e+05 | ops/s | 3.04x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.90654e+07 ± 1.542e+05 | 7.60575e+06 ± 4.4e+04 | ops/s | 2.51x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | vector (gemm-routed only) | thrpt | 1.8801e+07 ± 3.646e+05 | 7.60575e+06 ± 4.4e+04 | ops/s | 2.47x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.14781e+07 ± 8.619e+04 | 3.00426e+07 ± 1.679e+05 | ops/s | 1.05x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | vector | thrpt | 3.12082e+07 ± 6.432e+04 | 3.00426e+07 ± 1.679e+05 | ops/s | 1.04x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.89354e+07 ± 2.192e+05 | 2.04133e+07 ± 7.739e+04 | ops/s | 1.91x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | vector | thrpt | 3.87913e+07 ± 5.514e+04 | 2.04133e+07 ± 7.739e+04 | ops/s | 1.90x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 1.04844e+07 ± 1.443e+04 | 3.42875e+06 ± 8.946e+04 | ops/s | 3.06x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=3 | vector (gemm-routed only) | thrpt | 1.05184e+07 ± 1.344e+04 | 3.42875e+06 ± 8.946e+04 | ops/s | 3.07x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 6.17428e+06 ± 5.792e+04 | 2.72385e+06 ± 1.558e+05 | ops/s | 2.27x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=4 | vector (gemm-routed only) | thrpt | 6.18273e+06 ± 5.531e+04 | 2.72385e+06 ± 1.558e+05 | ops/s | 2.27x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.56639e+07 ± 6.637e+04 | 6.08403e+06 ± 2.208e+04 | ops/s | 2.57x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | vector (gemm-routed only) | thrpt | 1.57264e+07 ± 8.309e+04 | 6.08403e+06 ± 2.208e+04 | ops/s | 2.58x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.12448e+07 ± 1.799e+04 | 5.20805e+06 ± 4.481e+04 | ops/s | 2.16x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | vector (gemm-routed only) | thrpt | 1.12793e+07 ± 4041 | 5.20805e+06 ± 4.481e+04 | ops/s | 2.17x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.229061 ± 0.0005384 | 1.27307 ± 0.002637 | ms/op | 5.56x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.3335 ± 0.1573 | 122.212 ± 0.398 | ms/op | 5.02x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.84186 ± 0.007743 | 12.2441 ± 0.05885 | ms/op | 6.65x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 222.075 ± 0.9531 | 1222.79 ± 58.81 | ms/op | 5.51x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00576826 ± 2.515e-05 | 0.0403479 ± 0.0002156 | ms/op | 6.99x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.593055 ± 0.001398 | 3.80535 ± 0.004259 | ms/op | 6.42x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0459478 ± 0.0001125 | 0.372178 ± 0.001926 | ms/op | 8.10x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.15903 ± 0.003871 | 38.0043 ± 0.3938 | ms/op | 6.17x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140101 ± 0.0001985 | 1.27479 ± 0.004094 | ms/op | 9.10x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.6762 ± 0.2172 | 122.283 ± 0.6931 | ms/op | 8.33x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.1914 ± 0.006137 | 12.2085 ± 0.02433 | ms/op | 10.25x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 145.755 ± 0.3036 | 1226.74 ± 66.1 | ms/op | 8.42x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00494639 ± 9.873e-06 | 0.0404518 ± 0.0001852 | ms/op | 8.18x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.541461 ± 0.003074 | 3.82061 ± 0.03363 | ms/op | 7.06x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0528843 ± 0.0004903 | 0.376319 ± 0.009984 | ms/op | 7.12x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.62056 ± 0.0118 | 37.7612 ± 0.2688 | ms/op | 6.72x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 175139 ± 3146 | 199352 ± 1077 | ops/s | 0.88x | behind |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 16397.7 ± 452.1 | 16078.5 ± 37.68 | ops/s | 1.02x | tie |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.07504e+06 ± 3.718e+04 | 314989 ± 2021 | ops/s | 6.59x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 209447 ± 2126 | 32897 ± 1372 | ops/s | 6.37x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.11301e+06 ± 2.445e+04 | 325576 ± 297.4 | ops/s | 3.42x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 41310.6 ± 2.537e+04 | 35253.9 ± 35.53 | ops/s | 1.17x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.27016e+06 ± 966 | 2.25513e+06 ± 1.435e+04 | ops/s | 1.01x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 224972 ± 70.91 | 224868 ± 45.25 | ops/s | 1.00x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 70619.9 ± 2445 | 34158.8 ± 601.2 | ops/s | 2.07x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 1670.01 ± 11.91 | 1756.63 ± 16.03 | ops/s | 0.95x | behind |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 170.065 ± 1.925 | 314.257 ± 2.839 | ops/s | 0.54x | behind |  |

## Provenance

- Date: 2026-10-09, Slurm job 2524301, start 11:52:43, end 14:14:06 EDT (elapsed 2h21m23s)
- Benchmarked source: `17813d0d761467187ba0764f1081fd30747ad18e` (breeze/integration). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster.
- Node: SciNet trillium `tri0732`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5, family 26 model 2), AVX-512, SMT off, 192 CPUs, NPS4 (8 NUMA nodes of 24 cores), boost enabled. Load average at job start: 0.66, 0.22, 11.55.
- JDK: Temurin 25+36 (`module load java/25`, java/25.36). The site default `JAVA_TOOL_OPTIONS=-Xmx2g` was unset before launch, so no fork inherited it.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the L3/CCD domain of cpu0 plus its NUMA node). It applies to the JMH host JVM and, by inheritance, to every fork.
- Netlib (sidecar `2026-10-09-breeze-two-lane-baseline-trillium-laneB.netlib.jsonl`, (blas, lapack, vectorModule) over all records): [('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)]. No native BLAS/LAPACK library is visible to the loader on these nodes (`ldconfig -p` lists none; `LD_LIBRARY_PATH=/opt/slurm/lib64`), so nothing had to be forced.
- Command, equivalent to the `breezeLaneB` alias: `numactl … java --add-modules=jdk.incubator.vector -cp <classpath> org.openjdk.jmh.Main -jvmArgsAppend "-Xms4g -Xmx4g -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure,vector -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings come from the annotations. As under sbt, the host JVM carries `--add-modules`: lane A's `-jvmArgs` replaces it in the forks, and lane B's forks inherit it. The runner script is `docs/verification/w21-simd-spike/x86-avx512/run.sh`.
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane B --strict --commit 17813d0d… --netlib 2026-10-09-breeze-two-lane-baseline-trillium-laneB.netlib.jsonl 2026-10-09-breeze-two-lane-baseline-trillium-laneB.json`. It exited 0, with 0 unpaired rows.
