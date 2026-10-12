# gale vs Breeze scoreboard — lane B (SIMD)

- Lane: B
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: --add-modules=jdk.incubator.vector -Dgale.bench.sourceSha=96c95e577f6b437c3de995a081b8ba075a6054e4 -Xms4g -Xmx4g -Dgale.bench.sourceSha=96c95e577f6b437c3de995a081b8ba075a6054e4 -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/m3/out/laneB-2538727/laneB/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.VectorBLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `96c95e577f6b437c3de995a081b8ba075a6054e4`
- Numerical preflight: 151 decomposition/factorization result rows validated on this SHA.
- Machine: SciNet trillium tri0501: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512, 8 double lanes), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2538727 exclusive

**Headline, gale vector vs Breeze VectorBLAS: 118 ahead, 14 tie, 9 behind, 0 n/a, 10 withheld** of 151 pairs (`vector` rows plus `backend-insensitive` rows, gale's best configuration in this lane).

Secondary, gale pure vs Breeze VectorBLAS: 46 ahead, 0 tie, 11 behind, 0 n/a, 1 withheld of 58 pairs.

All rows: 164 ahead, 14 tie, 20 behind, 0 n/a, 11 withheld of 209 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 186230 ± 3.68e+04 | 188835 ± 3.623e+04 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 29574.1 ± 2365 | 29616.6 ± 3327 | ops/s | 1.00x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7688.42 ± 974.2 | 7824.41 ± 691.1 | ops/s | 0.98x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 65853.8 ± 647.3 | 68035.1 ± 560.3 | ops/s | 0.97x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16949.7 ± 111.8 | 16835.2 ± 455.5 | ops/s | 1.01x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4199.25 ± 45.13 | 4238.06 ± 42.41 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 66986.1 ± 1641 | 4660.13 ± 11.74 | ops/s | 14.37x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 16962.4 ± 130.8 | 1242.47 ± 274 | ops/s | 13.65x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 4236.66 ± 72.46 | 358.496 ± 0.4085 | ops/s | 11.82x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 98494.1 ± 476.7 | 234667 ± 946.5 | ops/s | 0.42x | behind |  |
| BlasL2BreezeJmh | Gemv | n=256 | vector | thrpt | 225633 ± 4394 | 234667 ± 946.5 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6518.39 ± 12.64 | 13287.1 ± 43.12 | ops/s | 0.49x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | thrpt | 12607.6 ± 41 | 13287.1 ± 43.12 | ops/s | 0.95x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1609.73 ± 15.19 | 2843.01 ± 22.01 | ops/s | 0.57x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | thrpt | 2724.76 ± 14.22 | 2843.01 ± 22.01 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 83452.6 ± 502.9 | 258440 ± 1.802e+04 | ops/s | 0.32x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | vector | thrpt | 253835 ± 1.673e+04 | 258440 ± 1.802e+04 | ops/s | 0.98x | tie |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5790.73 ± 32.04 | 13224.9 ± 126.6 | ops/s | 0.44x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | thrpt | 12645.8 ± 35.44 | 13224.9 ± 126.6 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1425.41 ± 10.74 | 2973.64 ± 9.758 | ops/s | 0.48x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | thrpt | 2956.65 ± 14.89 | 2973.64 ± 9.758 | ops/s | 0.99x | tie |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 790587 ± 1328 | 396547 ± 1.024e+04 | ops/s | 1.99x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | vector | thrpt | 791483 ± 1324 | 396547 ± 1.024e+04 | ops/s | 2.00x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15316.1 ± 86.42 | 14433 ± 36.82 | ops/s | 1.06x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | vector | thrpt | 18290.8 ± 42.33 | 14433 ± 36.82 | ops/s | 1.27x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 214.672 ± 5.491 | 356.925 ± 9.868 | ops/s | 0.60x | behind |  |
| BlasL3BreezeJmh | AtA | n=256 | vector | thrpt | 651.932 ± 8.057 | 356.925 ± 9.868 | ops/s | 1.83x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.79407e+06 ± 3895 | 1.53348e+06 ± 1.719e+04 | ops/s | 1.17x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | vector | thrpt | 1.793e+06 ± 3747 | 1.53348e+06 ± 1.719e+04 | ops/s | 1.17x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31740.8 ± 307 | 28411.9 ± 162.5 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | vector | thrpt | 79272.5 ± 4097 | 28411.9 ± 162.5 | ops/s | 2.79x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 518.193 ± 3.501 | 1143.61 ± 20.08 | ops/s | 0.45x | behind |  |
| BlasL3BreezeJmh | Gemm | n=256 | vector | thrpt | 1583.42 ± 111.5 | 1143.61 ± 20.08 | ops/s | 1.38x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 438103 ± 1055 | 401227 ± 2194 | ops/s | 1.09x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | thrpt | 437617 ± 1394 | 401227 ± 2194 | ops/s | 1.09x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7893.88 ± 54.32 | 12469.2 ± 337.9 | ops/s | 0.63x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | thrpt | 21864.6 ± 1193 | 12469.2 ± 337.9 | ops/s | 1.75x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 129.443 ± 1.533 | 294.196 ± 8.022 | ops/s | 0.44x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | thrpt | 467.454 ± 2.112 | 294.196 ± 8.022 | ops/s | 1.59x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.3847e+06 ± 1.632e+04 | 937285 ± 1.298e+04 | ops/s | 1.48x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | vector (gemm-routed only) | thrpt | 1.37347e+06 ± 6967 | 937285 ± 1.298e+04 | ops/s | 1.47x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 47610.7 ± 87.1 | 31642.8 ± 274.9 | ops/s | 1.50x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | vector (gemm-routed only) | thrpt | 47734.6 ± 145.7 | 31642.8 ± 274.9 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 922.977 ± 9.72 | 476.973 ± 2.664 | ops/s | 1.94x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | vector (gemm-routed only) | thrpt | 926.268 ± 9.074 | 476.973 ± 2.664 | ops/s | 1.94x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 34139.8 ± 1446 | 12914 ± 103 | ops/s | 2.64x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1086.64 ± 1.056 | 721.682 ± 5.371 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 19.0087 ± 0.02898 | 1.32442 ± 0.004474 | ops/s | 14.35x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 467374 ± 1.949e+04 | 334949 ± 1.106e+04 | ops/s | 1.40x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector (gemm-routed only) | thrpt | 459008 ± 9216 | 334949 ± 1.106e+04 | ops/s | 1.37x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 17482.5 ± 52.11 | 9822.52 ± 36.24 | ops/s | 1.78x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector (gemm-routed only) | thrpt | 17495.5 ± 42.62 | 9822.52 ± 36.24 | ops/s | 1.78x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 371.883 ± 13.87 | 165.001 ± 0.04207 | ops/s | 2.25x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector (gemm-routed only) | thrpt | 371.747 ± 14.87 | 165.001 ± 0.04207 | ops/s | 2.25x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 45448.4 ± 1059 | 15412.1 ± 60.26 | ops/s | 2.95x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 1373.1 ± 9.27 | 1081.57 ± 13.67 | ops/s | 1.27x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 27.7125 ± 0.08407 | 26.195 ± 0.05551 | ops/s | 1.06x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 51440.6 ± 855.7 | 18539 ± 135.9 | ops/s | 2.77x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 1688.04 ± 5.17 | 1297.3 ± 10.78 | ops/s | 1.30x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 33.0832 ± 0.04683 | 28.9025 ± 1.02 | ops/s | 1.14x | ahead |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 44311.6 ± 953.5 | 37696.1 ± 430.1 | ops/s | 1.18x | ahead |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 2062.65 ± 25.51 | 1611.21 ± 16.83 | ops/s | 1.28x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 43998.2 ± 734.6 | 37100.7 ± 1624 | ops/s | 1.19x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 2090.66 ± 12.6 | 1610.62 ± 13.65 | ops/s | 1.30x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 43706.5 ± 540.4 | 38282.2 ± 249.5 | ops/s | 1.14x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 2069.85 ± 35.41 | 1623.89 ± 18.86 | ops/s | 1.27x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.4654e+06 ± 4161 | 377360 ± 4556 | ops/s | 3.88x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | vector (gemm-routed only) | thrpt | 1.4665e+06 ± 4510 | 377360 ± 4556 | ops/s | 3.89x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 56138.9 ± 146.3 | 20323.3 ± 241.6 | ops/s | 2.76x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | vector (gemm-routed only) | thrpt | 56169 ± 220.3 | 20323.3 ± 241.6 | ops/s | 2.76x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1075.28 ± 14.22 | 672.901 ± 2.166 | ops/s | 1.60x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | vector (gemm-routed only) | thrpt | 1073.92 ± 12.67 | 672.901 ± 2.166 | ops/s | 1.60x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.46937e+06 ± 9.914e+04 | 1.04525e+06 ± 1.117e+04 | ops/s | 1.41x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | vector (gemm-routed only) | thrpt | 1.41228e+06 ± 9164 | 1.04525e+06 ± 1.117e+04 | ops/s | 1.35x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 47717.1 ± 122 | 32439.2 ± 310.8 | ops/s | 1.47x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | vector (gemm-routed only) | thrpt | 47894 ± 98.79 | 32439.2 ± 310.8 | ops/s | 1.48x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 928.263 ± 10.46 | 477.042 ± 2.182 | ops/s | 1.95x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | vector (gemm-routed only) | thrpt | 931.402 ± 11.3 | 477.042 ± 2.182 | ops/s | 1.95x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 774484 ± 4710 | 363120 ± 2684 | ops/s | 2.13x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | vector (gemm-routed only) | thrpt | 764953 ± 1.092e+04 | 363120 ± 2684 | ops/s | 2.11x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 16839 ± 181.7 | 15582.5 ± 41.66 | ops/s | 1.08x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | vector (gemm-routed only) | thrpt | 16827 ± 199.8 | 15582.5 ± 41.66 | ops/s | 1.08x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 342.128 ± 1.027 | 231.757 ± 1.594 | ops/s | 1.48x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | vector (gemm-routed only) | thrpt | 423.416 ± 4.936 | 231.757 ± 1.594 | ops/s | 1.83x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 1.06103e+06 ± 1.629e+04 | 779875 ± 7611 | ops/s | 1.36x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | vector (gemm-routed only) | thrpt | 1.06268e+06 ± 1.814e+04 | 779875 ± 7611 | ops/s | 1.36x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 41312.7 ± 276.2 | 30035.9 ± 170.2 | ops/s | 1.38x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | vector (gemm-routed only) | thrpt | 41352.3 ± 234.6 | 30035.9 ± 170.2 | ops/s | 1.38x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 790.708 ± 3.661 | 471.019 ± 2.461 | ops/s | 1.68x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | vector (gemm-routed only) | thrpt | 790.448 ± 4.645 | 471.019 ± 2.461 | ops/s | 1.68x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.18052 ± 0.01093 | 10.449 ± 0.1307 | ms/op | 1.28x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | vector (gemm-routed only) | avgt | 8.18183 ± 0.02374 | 10.449 ± 0.1307 | ms/op | 1.28x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.2984 ± 0.2514 | 73.8615 ± 0.2383 | ms/op | 1.04x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | vector (gemm-routed only) | avgt | 71.3004 ± 0.1881 | 73.8615 ± 0.2383 | ms/op | 1.04x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 102.884 ± 1.423 | 143.612 ± 0.1457 | ms/op | 1.40x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 773.433 ± 6.133 | 1070.48 ± 1.345 | ms/op | 1.38x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 76.8518 ± 0.3561 | 94.3768 ± 0.02929 | ms/op | 1.23x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | vector (gemm-routed only) | avgt | 68.9561 ± 0.4432 | 94.3768 ± 0.02929 | ms/op | 1.37x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 538.616 ± 8.857 | 748.972 ± 1.112 | ms/op | 1.39x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | vector (gemm-routed only) | avgt | 372.883 ± 6.578 | 748.972 ± 1.112 | ms/op | 2.01x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 8.79223 ± 0.01517 | 16.3525 ± 0.025 | ms/op | 1.86x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | vector (gemm-routed only) | avgt | 8.81022 ± 0.01321 | 16.3525 ± 0.025 | ms/op | 1.86x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 69.2717 ± 0.1419 | 128.435 ± 0.1573 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | vector (gemm-routed only) | avgt | 69.3488 ± 0.175 | 128.435 ± 0.1573 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 25.9326 ± 0.03904 | 36.6736 ± 0.1423 | ms/op | 1.41x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | vector (gemm-routed only) | avgt | 24.5713 ± 0.2196 | 36.6736 ± 0.1423 | ms/op | 1.49x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 193.146 ± 4.484 | 292.606 ± 0.2591 | ms/op | 1.51x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | vector (gemm-routed only) | avgt | 125.911 ± 1.576 | 292.606 ± 0.2591 | ms/op | 2.32x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 9.06662 ± 0.01408 | 16.4562 ± 0.02272 | ms/op | 1.82x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | vector (gemm-routed only) | avgt | 9.11453 ± 0.04584 | 16.4562 ± 0.02272 | ms/op | 1.81x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 70.3804 ± 0.1757 | 128.927 ± 0.2338 | ms/op | 1.83x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | vector (gemm-routed only) | avgt | 70.2763 ± 0.2271 | 128.927 ± 0.2338 | ms/op | 1.83x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 814.5 ± 10.11 | 556.114 ± 2.054 | us/op | 0.68x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | avgt | 579.51 ± 7.848 | 556.114 ± 2.054 | us/op | 0.96x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1900.06 ± 16.79 | 1047.9 ± 5.737 | us/op | 0.55x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | avgt | 1358.84 ± 30.91 | 1047.9 ± 5.737 | us/op | 0.77x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 44.6954 ± 0.4005 | 77.7522 ± 0.4836 | us/op | 1.74x | ahead |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 1237.37 ± 8.69 | 1855.32 ± 7.69 | us/op | 1.50x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 134762 ± 1864 | 106795 ± 245.6 | ops/s | 1.26x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector (gemm-routed only) | thrpt | 135491 ± 1204 | 106795 ± 245.6 | ops/s | 1.27x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2954.33 ± 11.5 | 2295.78 ± 9.597 | ops/s | 1.29x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector (gemm-routed only) | thrpt | 2952.07 ± 15.55 | 2295.78 ± 9.597 | ops/s | 1.29x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 45.7056 ± 0.3232 | 37.6405 ± 0.04867 | ops/s | 1.21x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector (gemm-routed only) | thrpt | 47.5737 ± 1.436 | 37.6405 ± 0.04867 | ops/s | 1.26x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 292.875 ± 0.2703 | 103.94 ± 0.1347 | ops/s | 2.82x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 3697.52 ± 76.97 | 2868.43 ± 10.44 | ops/s | 1.29x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 4643.34 ± 16.28 | 191.37 ± 0.9409 | ops/s | 24.26x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4234.73 ± 12.88 | 358.889 ± 0.2253 | ops/s | 11.80x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 254.763 ± 0.5683 | 66.4301 ± 0.8187 | ops/s | 3.84x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7460.4 ± 3.49 | 2150.31 ± 2.15 | ops/s | 3.47x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 3588.1 ± 41.8 | 2256.97 ± 0.9079 | ops/s | 1.59x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7215.54 ± 2.224 | 389.015 ± 10.91 | ops/s | 18.55x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 660.158 ± 8.254 | 321.166 ± 0.8021 | ops/s | 2.06x | ahead |  |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 660.859 ± 10.09 | 321.166 ± 0.8021 | ops/s | 2.06x | ahead |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 621.074 ± 4.555 | 242.183 ± 0.04753 | ops/s | 2.56x | ahead |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 622.741 ± 4.759 | 242.183 ± 0.04753 | ops/s | 2.57x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.37352e+06 ± 6.889e+04 | 1.28317e+06 ± 4885 | ops/s | 2.63x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 53970.5 ± 90.88 | 32623.5 ± 4855 | ops/s | 1.65x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3399.45 ± 1.455 | 2091.69 ± 158.2 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 342333 ± 1594 | 341468 ± 1264 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5061.57 ± 118.2 | 5041.31 ± 102.4 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 324.985 ± 0.3623 | 217.672 ± 0.5476 | ops/s | 1.49x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 298826 ± 3593 | 302869 ± 656.9 | ops/s | 0.99x | tie |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4838.2 ± 2.068 | 4627.8 ± 19.41 | ops/s | 1.05x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 302.635 ± 0.07788 | 284.754 ± 0.3953 | ops/s | 1.06x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 5.11082e+06 ± 1.766e+04 | 5.18955e+06 ± 1.34e+05 | ops/s | 0.98x | tie |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 89352.6 ± 56.75 | 87719.5 ± 9906 | ops/s | 1.02x | tie |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4921.85 ± 1082 | 5457.77 ± 531 | ops/s | 0.90x | tie |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.52403e+06 ± 2549 | 246433 ± 40.11 | ops/s | 30.53x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118741 ± 56.53 | 3798.66 ± 1.095 | ops/s | 31.26x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7388.52 ± 11.38 | 237.823 ± 0.2153 | ops/s | 31.07x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.37603e+06 ± 2718 | 108065 ± 30.06 | ops/s | 40.49x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68803.5 ± 11.37 | 1680.3 ± 3.179 | ops/s | 40.95x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4287.62 ± 4.761 | 104.764 ± 0.04665 | ops/s | 40.93x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.88791e+06 ± 2.073e+04 | 2.32349e+06 ± 4320 | ops/s | 1.24x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55775.9 ± 31.93 | 34245 ± 32.56 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3507.38 ± 10.75 | 2136.88 ± 0.5148 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 280270 ± 666.6 | 291829 ± 2122 | ops/s | 0.96x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 4097.21 ± 58.95 | 4303.39 ± 18.45 | ops/s | 0.95x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 265.665 ± 0.9974 | 186.74 ± 0.4201 | ops/s | 1.42x | ahead |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 272569 ± 696.8 | 148501 ± 347.2 | ops/s | 1.84x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3983.53 ± 18.79 | 2173.26 ± 265.9 | ops/s | 1.83x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 256.453 ± 1.302 | 116.076 ± 0.4586 | ops/s | 2.21x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.69132e+06 ± 1.181e+04 | 2.3491e+06 ± 2962 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 122686 ± 669 | 34467.7 ± 37.13 | ops/s | 3.56x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7272.86 ± 9.022 | 2150.4 ± 2.586 | ops/s | 3.38x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.8197e+07 ± 3.082e+05 | 7.89684e+06 ± 1.267e+05 | ops/s | 3.57x | ahead |  |
| SmallDenseBreezeJmh | Det | n=3 | vector (gemm-routed only) | thrpt | 2.83447e+07 ± 3.475e+05 | 7.89684e+06 ± 1.267e+05 | ops/s | 3.59x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.91049e+07 ± 1.397e+06 | 6.81691e+06 ± 1.494e+05 | ops/s | 2.80x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | vector (gemm-routed only) | thrpt | 1.83921e+07 ± 1.746e+05 | 6.81691e+06 ± 1.494e+05 | ops/s | 2.70x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.14452e+07 ± 1.46e+05 | 2.5441e+07 ± 5.922e+05 | ops/s | 1.24x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | vector | thrpt | 3.13221e+07 ± 1.157e+05 | 2.5441e+07 ± 5.922e+05 | ops/s | 1.23x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.87113e+07 ± 7.096e+05 | 1.82792e+07 ± 1.32e+05 | ops/s | 2.12x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | vector | thrpt | 3.86386e+07 ± 2.892e+05 | 1.82792e+07 ± 1.32e+05 | ops/s | 2.11x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 9.82481e+06 ± 1.361e+05 | 3.38114e+06 ± 6.616e+04 | ops/s | 2.91x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | vector (gemm-routed only) | thrpt | 9.90264e+06 ± 1.324e+04 | 3.38114e+06 ± 6.616e+04 | ops/s | 2.93x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 8.50452e+06 ± 4.483e+04 | 2.78891e+06 ± 4.537e+04 | ops/s | 3.05x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=4 | vector (gemm-routed only) | thrpt | 8.51073e+06 ± 6.11e+04 | 2.78891e+06 ± 4.537e+04 | ops/s | 3.05x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.60498e+07 ± 2.752e+04 | 5.40575e+06 ± 1.316e+05 | ops/s | 2.97x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | vector (gemm-routed only) | thrpt | 1.61321e+07 ± 5.566e+04 | 5.40575e+06 ± 1.316e+05 | ops/s | 2.98x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.12336e+07 ± 2.874e+05 | 4.55008e+06 ± 4.172e+05 | ops/s | 2.47x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | vector (gemm-routed only) | thrpt | 1.09307e+07 ± 1.405e+05 | 4.55008e+06 ± 4.172e+05 | ops/s | 2.40x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.229715 ± 0.00124 | 1.27366 ± 0.00526 | ms/op | 5.54x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.3753 ± 0.1795 | 122.024 ± 0.3486 | ms/op | 5.01x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.84234 ± 0.00386 | 12.2145 ± 0.01924 | ms/op | 6.63x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 221.542 ± 1.75 | 1225.75 ± 55.73 | ms/op | 5.53x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00577443 ± 1.249e-05 | 0.0403759 ± 0.000149 | ms/op | 6.99x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.593474 ± 0.001609 | 3.82006 ± 0.0554 | ms/op | 6.44x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0461864 ± 0.0004242 | 0.370837 ± 0.0009733 | ms/op | 8.03x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.17314 ± 0.01141 | 40.6557 ± 5.636 | ms/op | 6.59x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140046 ± 0.0003984 | 1.27823 ± 0.008975 | ms/op | 9.13x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.5458 ± 0.09704 | 121.97 ± 0.4268 | ms/op | 8.39x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19333 ± 0.007258 | 12.2149 ± 0.007124 | ms/op | 10.24x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 145.81 ± 2.87 | 1228.01 ± 50.63 | ms/op | 8.42x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00494833 ± 3.096e-06 | 0.0403577 ± 0.000279 | ms/op | 8.16x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.543913 ± 0.001658 | 3.80498 ± 0.02886 | ms/op | 7.00x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0529125 ± 0.000249 | 0.370746 ± 0.001161 | ms/op | 7.01x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.64577 ± 0.002769 | 37.7333 ± 0.06931 | ms/op | 6.68x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 211045 ± 1808 | 191310 ± 1.063e+04 | ops/s | 1.10x | ahead |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 17743.7 ± 50.58 | 16084.4 ± 42.25 | ops/s | 1.10x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.0878e+06 ± 8709 | 315968 ± 815.9 | ops/s | 6.61x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 208514 ± 5149 | 32628 ± 1274 | ops/s | 6.39x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.13879e+06 ± 5.052e+04 | 325643 ± 390.4 | ops/s | 3.50x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 41171.1 ± 2.501e+04 | 35306.5 ± 25.54 | ops/s | 1.17x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.26814e+06 ± 2776 | 2.25023e+06 ± 939.2 | ops/s | 1.01x | ahead |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 214957 ± 1693 | 224836 ± 38.42 | ops/s | 0.96x | behind |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97388.5 ± 880.6 | 34219.5 ± 547.8 | ops/s | 2.85x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 3713.49 ± 9.262 | 1771.29 ± 10.04 | ops/s | 2.10x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 596.561 ± 18.91 | 315.506 ± 3.565 | ops/s | 1.89x | ahead |  |

**Milestone 2 under the same scoreboard:** re-scoring the committed Milestone 2 lane-B receipt with this tool version gives a headline of 114 ahead, 16 tie, 11 behind, 10 withheld; a secondary of 45/0/12/1; and all rows 159/16/23/11.

## Delta against Milestone 2

Milestone 2 is `2026-10-10-breeze-milestone2-trillium-laneB.md` (source 0ff6440, receipts 5d3950e). It used the same lane, CPU model, JDK, pinning, runner and preflight, on a different node of the same type. This is an unpaired, cross-run comparison: Breeze is re-measured in each run, so read the gale/Breeze ratio changes rather than gale's absolute times.

Ratios are gale speed / Breeze speed (>1 means gale is faster). `change` is new ratio / old ratio. Old is milestone 2 0ff6440; new is milestone 3 96c95e5.

| class | op | params | gale backend | old ratio | new ratio | change | old verdict | new verdict |
|---|---|---|---|---:|---:|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | 0.96x | 0.98x | 1.02x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | 0.94x | 1.00x | 1.06x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | 1.04x | 0.99x | 0.95x | ahead | tie |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | 0.99x | 0.99x | 1.00x | tie | tie |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | 0.92x | 1.01x | 1.10x | tie | tie |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | 0.96x | 0.97x | 1.01x | behind | behind |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | 11.76x | 11.82x | 1.01x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | 13.69x | 13.65x | 1.00x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | 14.44x | 14.37x | 1.00x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 0.51x | 0.49x | 0.96x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | 0.99x | 0.95x | 0.96x | tie | behind |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 0.57x | 0.57x | 1.00x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | 0.97x | 0.96x | 0.99x | tie | behind |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 0.42x | 0.42x | 1.00x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=256 | vector | 1.05x | 0.96x | 0.91x | tie | behind |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 0.44x | 0.44x | 1.00x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | 1.00x | 0.96x | 0.96x | tie | behind |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 0.49x | 0.48x | 0.98x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | 0.94x | 0.99x | 1.05x | behind | tie |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 0.34x | 0.32x | 0.94x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=256 | vector | 1.01x | 0.98x | 0.97x | tie | tie |
| BlasL3BreezeJmh | AtA | n=16 | pure | 1.99x | 1.99x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=16 | vector | 1.98x | 2.00x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 0.61x | 0.60x | 0.98x | behind | behind |
| BlasL3BreezeJmh | AtA | n=256 | vector | 0.61x | 1.83x | 3.00x | behind | ahead |
| BlasL3BreezeJmh | AtA | n=64 | pure | 1.07x | 1.06x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=64 | vector | 1.06x | 1.27x | 1.20x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.17x | 1.17x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | vector | 1.17x | 1.17x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 0.45x | 0.45x | 1.00x | behind | behind |
| BlasL3BreezeJmh | Gemm | n=256 | vector | 1.36x | 1.38x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | vector | 1.12x | 2.79x | 2.49x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.10x | 1.09x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | 1.10x | 1.09x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 0.44x | 0.44x | 1.00x | behind | behind |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | 1.59x | 1.59x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 0.62x | 0.63x | 1.02x | behind | behind |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | 0.62x | 1.75x | 2.82x | behind | ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.46x | 1.48x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | vector (gemm-routed only) | 1.46x | 1.47x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.93x | 1.94x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | vector (gemm-routed only) | 1.94x | 1.94x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 1.51x | 1.50x | 0.99x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | vector (gemm-routed only) | 1.51x | 1.51x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | 2.61x | 2.64x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | 14.34x | 14.35x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | 1.50x | 1.51x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 1.41x | 1.40x | 0.99x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector (gemm-routed only) | 1.35x | 1.37x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 2.21x | 2.25x | 1.02x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector (gemm-routed only) | 2.28x | 2.25x | 0.99x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 1.78x | 1.78x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector (gemm-routed only) | 1.78x | 1.78x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | 2.95x | 2.95x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | 1.06x | 1.06x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | 1.27x | 1.27x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | 3.70x | 2.77x | 0.75x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | 1.17x | 1.14x | 0.97x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | 1.30x | 1.30x | 1.00x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | 1.30x | 1.28x | 0.98x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | 1.17x | 1.18x | 1.01x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | 1.30x | 1.30x | 1.00x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | 1.16x | 1.19x | 1.03x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | 1.28x | 1.27x | 0.99x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | 1.19x | 1.14x | 0.96x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.91x | 3.88x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=16 | vector (gemm-routed only) | 3.91x | 3.89x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.63x | 1.60x | 0.98x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | vector (gemm-routed only) | 1.63x | 1.60x | 0.98x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.78x | 2.76x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | vector (gemm-routed only) | 2.77x | 2.76x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.35x | 1.41x | 1.04x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | vector (gemm-routed only) | 1.36x | 1.35x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.95x | 1.95x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | vector (gemm-routed only) | 1.95x | 1.95x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 1.47x | 1.47x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=64 | vector (gemm-routed only) | 1.47x | 1.48x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.54x | 2.13x | 1.38x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=16 | vector (gemm-routed only) | 1.51x | 2.11x | 1.40x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.48x | 1.48x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | vector (gemm-routed only) | 1.49x | 1.83x | 1.23x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 1.09x | 1.08x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=64 | vector (gemm-routed only) | 1.09x | 1.08x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.37x | 1.36x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=16 | vector (gemm-routed only) | 1.37x | 1.36x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.68x | 1.68x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | vector (gemm-routed only) | 1.77x | 1.68x | 0.95x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 1.37x | 1.38x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=64 | vector (gemm-routed only) | 1.38x | 1.38x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | 1.04x | 1.04x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=1024 | vector (gemm-routed only) | 1.04x | 1.04x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | 1.28x | 1.28x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | vector (gemm-routed only) | 1.28x | 1.28x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | 1.38x | 1.38x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | 1.39x | 1.40x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | 1.38x | 1.39x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | vector (gemm-routed only) | 1.89x | 2.01x | 1.06x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | 1.23x | 1.23x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | vector (gemm-routed only) | 1.36x | 1.37x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | 1.85x | 1.85x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | vector (gemm-routed only) | 1.85x | 1.85x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | 1.85x | 1.86x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | vector (gemm-routed only) | 1.85x | 1.86x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.45x | 1.51x | 1.04x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | vector (gemm-routed only) | 2.10x | 2.32x | 1.10x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.39x | 1.41x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | vector (gemm-routed only) | 1.42x | 1.49x | 1.05x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | 1.83x | 1.83x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | vector (gemm-routed only) | 1.83x | 1.83x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | 1.81x | 1.82x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | vector (gemm-routed only) | 1.81x | 1.81x | 1.00x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | 0.66x | 0.68x | 1.03x | behind | behind |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | 0.90x | 0.96x | 1.07x | behind | behind |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | 0.53x | 0.55x | 1.04x | withheld | withheld |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | 0.73x | 0.77x | 1.05x | withheld | withheld |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | 0.75x | 1.74x | 2.32x | behind | ahead |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | 0.82x | 1.50x | 1.83x | withheld | withheld |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.87x | 1.26x | 1.45x | behind | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector (gemm-routed only) | 0.88x | 1.27x | 1.44x | behind | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.22x | 1.21x | 0.99x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector (gemm-routed only) | 1.28x | 1.26x | 0.98x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.23x | 1.29x | 1.05x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector (gemm-routed only) | 1.23x | 1.29x | 1.05x | ahead | ahead |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | 2.83x | 2.82x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | 1.29x | 1.29x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | 24.53x | 24.26x | 0.99x | ahead | ahead |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | 11.81x | 11.80x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | 3.81x | 3.84x | 1.01x | ahead | ahead |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.46x | 3.47x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | 1.59x | 1.59x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | 18.46x | 18.55x | 1.00x | ahead | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | 2.06x | 2.06x | 1.00x | ahead | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | vector (gemm-routed only) | 2.06x | 2.06x | 1.00x | ahead | ahead |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | 2.57x | 2.56x | 1.00x | ahead | ahead |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | vector (gemm-routed only) | 2.58x | 2.57x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | 2.68x | 2.63x | 0.98x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | 1.62x | 1.63x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | 1.66x | 1.65x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | tie |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | 1.51x | 1.49x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | 0.99x | 1.00x | 1.01x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | 0.99x | 0.99x | 1.00x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | 1.07x | 1.06x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | 1.05x | 1.05x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | 1.06x | 0.98x | 0.92x | tie | tie |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | 0.76x | 0.90x | 1.18x | behind | tie |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | 0.84x | 1.02x | 1.21x | tie | tie |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | 30.48x | 30.53x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | 30.98x | 31.07x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | 31.20x | 31.26x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | 40.51x | 40.49x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | 40.93x | 40.93x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | 40.96x | 40.95x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | 1.24x | 1.24x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | 1.64x | 1.64x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | 1.63x | 1.63x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | 0.85x | 0.96x | 1.13x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | 0.72x | 1.42x | 1.97x | behind | ahead |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | 0.54x | 0.95x | 1.76x | behind | behind |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | 1.83x | 1.84x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | 2.20x | 2.21x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | 1.83x | 1.83x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.27x | 3.27x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | 3.38x | 3.38x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | 3.57x | 3.56x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | pure | 3.52x | 3.57x | 1.01x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | vector (gemm-routed only) | 3.53x | 3.59x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | pure | 2.74x | 2.80x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | vector (gemm-routed only) | 2.72x | 2.70x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | 1.22x | 1.24x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | vector | 1.21x | 1.23x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | 2.14x | 2.12x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | vector | 2.13x | 2.11x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | pure | 3.11x | 2.91x | 0.94x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | vector (gemm-routed only) | 3.12x | 2.93x | 0.94x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | pure | 3.06x | 3.05x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | vector (gemm-routed only) | 3.07x | 3.05x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | pure | 2.90x | 2.97x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | vector (gemm-routed only) | 2.88x | 2.98x | 1.03x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | pure | 2.27x | 2.47x | 1.09x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | vector (gemm-routed only) | 2.24x | 2.40x | 1.07x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | 5.55x | 5.54x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | 5.02x | 5.01x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | 6.63x | 6.63x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | 5.57x | 5.53x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | 7.01x | 6.99x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | 6.44x | 6.44x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | 8.00x | 8.03x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | 6.11x | 6.59x | 1.08x | ahead | ahead |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | 9.08x | 9.13x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | 8.28x | 8.39x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | 10.24x | 10.24x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | 8.46x | 8.42x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | 8.12x | 8.16x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | 7.01x | 7.00x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | 7.02x | 7.01x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | 6.86x | 6.68x | 0.97x | withheld | withheld |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | 1.06x | 1.10x | 1.04x | ahead | ahead |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | 1.11x | 1.10x | 0.99x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | 6.50x | 6.61x | 1.02x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | 6.41x | 6.39x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | 3.49x | 3.50x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | 1.17x | 1.17x | 1.00x | tie | tie |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | 1.01x | 1.01x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | 1.00x | 0.96x | 0.96x | tie | behind |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | 1.82x | 1.89x | 1.04x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | 2.83x | 2.85x | 1.01x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | 2.10x | 2.10x | 1.00x | ahead | ahead |

## Supplementary: controlled n=16 eigen recheck

Run in the same Slurm job, immediately after the sweep, with `supp.sh eig16B` (`docs/verification/w21-simd-spike/x86-avx512/supp.sh`, unchanged from Milestone 2). Settings: 4 forks × 10 warmup × 10 measurement iterations of 1 s, `-p n=16`, `DenseDecompositionBreezeJmh.(gale|breeze)Eig` and `SymEigenBreezeJmh.(gale|breeze)EigSym`. The scoreboard ran with `--validation --strict`. Files: `2026-10-11-breeze-milestone3-trillium-laneB-eig16.{json,netlib.jsonl,validation.jsonl}`.

- Numerical preflight: 4 decomposition/factorization result rows validated on this SHA.
**Headline, gale vector vs Breeze VectorBLAS: 2 ahead, 0 tie, 0 behind, 0 n/a, 0 withheld** of 2 pairs (`vector` rows plus `backend-insensitive` rows, gale's best configuration in this lane).
| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33388.1 ± 98.88 | 12974.5 ± 32.78 | ops/s | 2.57x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 98066.4 ± 64.86 | 34707.5 ± 32.65 | ops/s | 2.83x | ahead |  |

## Provenance

- Date: Slurm job 2538727, start/end/elapsed: 2026-10-11T16:30:59 2026-10-11T18:52:26 02:21:27 
- Benchmarked source: `96c95e577f6b437c3de995a081b8ba075a6054e4` (breeze/integration = main, clean worktree `bench/m3-96c95e5`). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster. Every fork received `-Dgale.bench.sourceSha=96c95e577f6b437c3de995a081b8ba075a6054e4` from the bundle's `SHA` file.
- Node: SciNet trillium `tri0501`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5), AVX-512, SMT off, 192 CPUs, NPS4, boost enabled. Load average at job start: 1.40, 44.23, 87.98.
- JDK: Temurin 25+36 (`module load java/25`). `JAVA_TOOL_OPTIONS`, `JDK_JAVA_OPTIONS` and `_JAVA_OPTIONS` were unset.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the CCD of cpu0 and its NUMA node). It applies to the host JVM and is inherited by every fork.
- Netlib (blas, lapack, vectorModule) over all sidecar records (`2026-10-11-breeze-milestone3-trillium-laneB.netlib.jsonl`): [('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- Numerical preflight (`2026-10-11-breeze-milestone3-trillium-laneB.validation.jsonl`): {'passed': 302} sourceSha: ['96c95e577f6b437c3de995a081b8ba075a6054e4'].
- Runner: the reviewed `docs/verification/w21-simd-spike/x86-avx512/run.sh laneB` (unmodified, `GALE_BENCH_ROOT=/scratch/brad/gale-bench/m3`), i.e. `numactl … java --add-modules=jdk.incubator.vector -Dgale.bench.sourceSha=<sha> -cp <classpath> org.openjdk.jmh.Main -foe true -jvmArgsAppend "-Xms4g -Xmx4g -Dgale.bench.sourceSha=<sha> -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure,vector -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings (2 forks etc.) come from the annotations.
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane B --strict --commit 96c95e577f6b437c3de995a081b8ba075a6054e4 --netlib 2026-10-11-breeze-milestone3-trillium-laneB.netlib.jsonl --validation 2026-10-11-breeze-milestone3-trillium-laneB.validation.jsonl 2026-10-11-breeze-milestone3-trillium-laneB.json`. It exited 0, with 0 unpaired rows.
