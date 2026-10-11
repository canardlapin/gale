# gale vs Breeze scoreboard — lane A (out-of-box scalar)

- Lane: A
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: -Xms4g -Xmx4g -Dgale.bench.sourceSha=0ff64405ea981f5a6dbf026a12e2bcc81617f412 -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/m2/out/laneA-2535351/laneA/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.Java11BLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `0ff64405ea981f5a6dbf026a12e2bcc81617f412`
- Numerical preflight: 110 decomposition/factorization result rows validated on this SHA.
- Machine: SciNet trillium tri0464: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512, 8 double lanes), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2535351 exclusive

**122 ahead, 12 tie, 7 behind, 0 n/a, 10 withheld** of 151 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 196797 ± 4.504e+04 | 174963 ± 6.343e+04 | ops/s | 1.12x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 28278.1 ± 1436 | 29829.6 ± 2117 | ops/s | 0.95x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7858.98 ± 503.5 | 7573.07 ± 337 | ops/s | 1.04x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 65000.5 ± 1192 | 68149.2 ± 461 | ops/s | 0.95x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16802.3 ± 173.2 | 17119.3 ± 24.17 | ops/s | 0.98x | behind |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4170.29 ± 80.54 | 4085.77 ± 241.6 | ops/s | 1.02x | tie |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 67060.3 ± 1351 | 4655.65 ± 16.55 | ops/s | 14.40x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 16887.7 ± 60.13 | 1242.51 ± 262.2 | ops/s | 13.59x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 4216.85 ± 48.01 | 358.463 ± 0.3301 | ops/s | 11.76x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 99169.8 ± 2666 | 79321.1 ± 2051 | ops/s | 1.25x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6449.49 ± 124.3 | 5071.74 ± 23.46 | ops/s | 1.27x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1611.72 ± 15.09 | 1232.22 ± 5.777 | ops/s | 1.31x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 82237.5 ± 835 | 69179.2 ± 32.25 | ops/s | 1.19x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5809.37 ± 0.9238 | 4245.49 ± 21.95 | ops/s | 1.37x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1427.06 ± 7.83 | 1046.29 ± 3.551 | ops/s | 1.36x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 780132 ± 1.891e+04 | 399332 ± 1.013e+04 | ops/s | 1.95x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15351.8 ± 112.4 | 6242.5 ± 57.31 | ops/s | 2.46x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 215.897 ± 1.234 | 113.825 ± 0.06888 | ops/s | 1.90x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.79637e+06 ± 5757 | 1.51291e+06 ± 1.584e+04 | ops/s | 1.19x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31825.7 ± 226 | 28411.2 ± 201.1 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 519.109 ± 3.058 | 421.836 ± 4.055 | ops/s | 1.23x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 439389 ± 2007 | 396138 ± 1.055e+04 | ops/s | 1.11x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7885.69 ± 53.68 | 6121.34 ± 32.7 | ops/s | 1.29x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 129.161 ± 1.571 | 105.817 ± 0.5105 | ops/s | 1.22x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.34326e+06 ± 2.124e+05 | 945791 ± 2467 | ops/s | 1.42x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 47697.6 ± 325.4 | 31690.9 ± 201.6 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 926.6 ± 9.891 | 476.883 ± 2.167 | ops/s | 1.94x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33344.2 ± 321.6 | 12975.4 ± 47.45 | ops/s | 2.57x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1087.59 ± 2.131 | 723.846 ± 6.738 | ops/s | 1.50x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 19.0021 ± 0.04544 | 1.32243 ± 0.001098 | ops/s | 14.37x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 463639 ± 9934 | 340370 ± 8001 | ops/s | 1.36x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 17478.2 ± 103.3 | 9825.94 ± 30.35 | ops/s | 1.78x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 373.149 ± 14.03 | 164.699 ± 0.4988 | ops/s | 2.27x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 45606.8 ± 1250 | 15320.6 ± 239.4 | ops/s | 2.98x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 1371.59 ± 6.136 | 1082.1 ± 13.9 | ops/s | 1.27x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 27.753 ± 0.05201 | 24.3211 ± 0.0481 | ops/s | 1.14x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 51710.9 ± 969.2 | 18283.6 ± 420.6 | ops/s | 2.83x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 1692.27 ± 10.04 | 1297.12 ± 9.248 | ops/s | 1.30x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 33.0792 ± 0.04289 | 28.2595 ± 0.05579 | ops/s | 1.17x | ahead |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 43730.5 ± 389.2 | 37925.8 ± 1054 | ops/s | 1.15x | ahead |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 2088.71 ± 12.47 | 1614.89 ± 10.44 | ops/s | 1.29x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 44062.4 ± 815.5 | 37770.3 ± 454.4 | ops/s | 1.17x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 2038.6 ± 31.32 | 1621.36 ± 19.29 | ops/s | 1.26x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 44051.4 ± 404.5 | 37819.2 ± 791.1 | ops/s | 1.16x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 2072.78 ± 14.43 | 1609.14 ± 15.45 | ops/s | 1.29x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.46384e+06 ± 3706 | 409670 ± 2578 | ops/s | 3.57x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 56018 ± 749.6 | 19606.7 ± 1029 | ops/s | 2.86x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1077.7 ± 12.53 | 658.328 ± 5.232 | ops/s | 1.64x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.41705e+06 ± 1.014e+04 | 1.04475e+06 ± 2948 | ops/s | 1.36x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 47706.9 ± 346.5 | 32360.3 ± 223.2 | ops/s | 1.47x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 927.103 ± 11.53 | 477.509 ± 2.312 | ops/s | 1.94x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 515877 ± 6462 | 354608 ± 2.108e+04 | ops/s | 1.45x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 16912.9 ± 122.5 | 15495.5 ± 198.6 | ops/s | 1.09x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 342.372 ± 1.241 | 231.506 ± 1.464 | ops/s | 1.48x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 1.04619e+06 ± 2.947e+04 | 797421 ± 2.339e+04 | ops/s | 1.31x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 41331.2 ± 335.1 | 30097.2 ± 249.1 | ops/s | 1.37x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 788.922 ± 7.303 | 471.289 ± 2.484 | ops/s | 1.67x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.18232 ± 0.00627 | 10.5362 ± 0.0104 | ms/op | 1.29x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.3666 ± 0.1427 | 74.2708 ± 0.8549 | ms/op | 1.04x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 100.026 ± 4.332 | 143.271 ± 0.1191 | ms/op | 1.43x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 774.302 ± 2.801 | 1071.04 ± 0.7084 | ms/op | 1.38x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 76.8628 ± 0.4734 | 94.4813 ± 0.0538 | ms/op | 1.23x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 539.98 ± 10.95 | 736.045 ± 20.41 | ms/op | 1.36x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 8.80329 ± 0.03701 | 16.342 ± 0.01517 | ms/op | 1.86x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 69.3742 ± 0.1433 | 128.455 ± 0.1147 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 25.9264 ± 0.05224 | 36.6881 ± 0.1315 | ms/op | 1.42x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 194.286 ± 2.932 | 293.052 ± 0.4685 | ms/op | 1.51x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 9.05961 ± 0.01337 | 16.4452 ± 0.02371 | ms/op | 1.82x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 70.4215 ± 0.231 | 128.774 ± 0.2046 | ms/op | 1.83x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 838.775 ± 6.808 | 929.391 ± 18.05 | us/op | 1.11x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1949.96 ± 16.12 | 1753.8 ± 32.3 | us/op | 0.90x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 104.514 ± 1.928 | 78.0533 ± 1.124 | us/op | 0.75x | behind |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 2568.96 ± 55.18 | 1869.66 ± 28.18 | us/op | 0.73x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 92169.9 ± 726.7 | 106727 ± 825.1 | ops/s | 0.86x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2835.74 ± 28.4 | 2297.38 ± 7.123 | ops/s | 1.23x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 45.613 ± 0.4607 | 37.5949 ± 0.03103 | ops/s | 1.21x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 292.481 ± 0.499 | 103.783 ± 1.188 | ops/s | 2.82x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 3638.34 ± 99.47 | 2852.75 ± 24.26 | ops/s | 1.28x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 4629.51 ± 29.25 | 192.364 ± 0.4993 | ops/s | 24.07x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4236.86 ± 31.21 | 358.438 ± 0.1147 | ops/s | 11.82x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 254.644 ± 0.5534 | 65.9613 ± 2.552 | ops/s | 3.86x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7442.59 ± 18.88 | 2152.82 ± 1.085 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 3580.78 ± 60.83 | 2276.37 ± 7.18 | ops/s | 1.57x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7207.67 ± 4.614 | 405.067 ± 9.903 | ops/s | 17.79x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 661.616 ± 8.332 | 320.476 ± 0.2523 | ops/s | 2.06x | ahead |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 622.65 ± 5.286 | 242.204 ± 0.08328 | ops/s | 2.57x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.43642e+06 ± 3032 | 1.28475e+06 ± 1708 | ops/s | 2.67x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 54028.7 ± 37.7 | 32598.6 ± 4872 | ops/s | 1.66x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3398.55 ± 8.833 | 2076.89 ± 127.8 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 342227 ± 1489 | 341622 ± 1584 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5016.41 ± 80.21 | 5106.79 ± 54.76 | ops/s | 0.98x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 325.019 ± 0.5893 | 212.142 ± 8.539 | ops/s | 1.53x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 301623 ± 475.1 | 302468 ± 1066 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4841.36 ± 6.634 | 4631 ± 5.085 | ops/s | 1.05x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 302.676 ± 0.06706 | 284.781 ± 0.3795 | ops/s | 1.06x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 5.42166e+06 ± 9.781e+05 | 5.12084e+06 ± 8.568e+04 | ops/s | 1.06x | tie |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 78494.5 ± 1.71e+04 | 83669.9 ± 1.483e+04 | ops/s | 0.94x | tie |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 5385.8 ± 719.6 | 5641.39 ± 46.68 | ops/s | 0.95x | tie |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.51687e+06 ± 1980 | 246439 ± 139 | ops/s | 30.50x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118618 ± 17.96 | 3799.3 ± 0.8621 | ops/s | 31.22x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7375.42 ± 12.93 | 238.058 ± 0.06779 | ops/s | 30.98x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.37982e+06 ± 3269 | 108064 ± 57.31 | ops/s | 40.53x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68795.1 ± 64.23 | 1680.53 ± 3.257 | ops/s | 40.94x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4288.45 ± 3.214 | 104.792 ± 0.01659 | ops/s | 40.92x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.88665e+06 ± 2.198e+04 | 2.32252e+06 ± 4740 | ops/s | 1.24x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55795.3 ± 48.46 | 34252 ± 31.48 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3511.62 ± 5.605 | 2137.45 ± 1.084 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 247736 ± 3553 | 292445 ± 752.9 | ops/s | 0.85x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 2286.96 ± 99.37 | 4362.63 ± 69.52 | ops/s | 0.52x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 132.753 ± 2.615 | 187.094 ± 0.1716 | ops/s | 0.71x | behind |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 271149 ± 3607 | 148510 ± 1983 | ops/s | 1.83x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3994.01 ± 21.68 | 2177.43 ± 253.7 | ops/s | 1.83x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 256.846 ± 0.6488 | 117.655 ± 0.394 | ops/s | 2.18x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.6864e+06 ± 7852 | 2.35241e+06 ± 341.8 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 122725 ± 488 | 34481.2 ± 32.08 | ops/s | 3.56x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7272.17 ± 10.91 | 2152.66 ± 0.6129 | ops/s | 3.38x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.80875e+07 ± 4.784e+05 | 7.82538e+06 ± 1.591e+04 | ops/s | 3.59x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.90118e+07 ± 1.158e+06 | 6.62193e+06 ± 4.846e+04 | ops/s | 2.87x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.14974e+07 ± 5.683e+04 | 2.48778e+07 ± 6.201e+05 | ops/s | 1.27x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.78768e+07 ± 1.01e+06 | 1.78826e+07 ± 2.537e+05 | ops/s | 2.12x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 9.84749e+06 ± 9.1e+04 | 3.36771e+06 ± 5.531e+04 | ops/s | 2.92x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 8.00143e+06 ± 5.702e+05 | 2.78635e+06 ± 1.574e+04 | ops/s | 2.87x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.60642e+07 ± 1.914e+04 | 5.24759e+06 ± 9.664e+04 | ops/s | 3.06x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.10108e+07 ± 3.482e+05 | 4.73843e+06 ± 3.054e+04 | ops/s | 2.32x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.229293 ± 0.0004866 | 1.3157 ± 0.07045 | ms/op | 5.74x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.3988 ± 0.1452 | 133.234 ± 0.1339 | ms/op | 5.46x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.8428 ± 0.005202 | 11.9914 ± 0.01804 | ms/op | 6.51x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 221.524 ± 1.877 | 1269.41 ± 79.56 | ms/op | 5.73x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00576559 ± 5.474e-06 | 0.040701 ± 0.0002936 | ms/op | 7.06x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.594764 ± 0.001198 | 3.77395 ± 0.01399 | ms/op | 6.35x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0463834 ± 0.0002344 | 0.369448 ± 0.002175 | ms/op | 7.97x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.17693 ± 0.007608 | 37.4113 ± 0.01761 | ms/op | 6.06x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.139935 ± 0.0003058 | 1.30266 ± 0.0222 | ms/op | 9.31x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.4413 ± 0.2119 | 133.499 ± 0.2668 | ms/op | 9.24x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19138 ± 0.006529 | 12.2207 ± 0.2305 | ms/op | 10.26x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 147.862 ± 2.808 | 1270.58 ± 74.93 | ms/op | 8.59x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00494939 ± 7.917e-06 | 0.0404615 ± 0.0002915 | ms/op | 8.18x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.542698 ± 0.0009153 | 3.78147 ± 0.004706 | ms/op | 6.97x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0529026 ± 0.0001326 | 0.368493 ± 0.001699 | ms/op | 6.97x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.60811 ± 0.02089 | 37.5977 ± 0.7771 | ms/op | 6.70x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 209670 ± 2452 | 198336 ± 1170 | ops/s | 1.06x | ahead |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 17732.5 ± 37.25 | 16068.3 ± 27.12 | ops/s | 1.10x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.05644e+06 ± 4.188e+04 | 315997 ± 1063 | ops/s | 6.51x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 207817 ± 5186 | 32683 ± 1287 | ops/s | 6.36x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.12526e+06 ± 3.508e+04 | 325707 ± 552.8 | ops/s | 3.45x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 41102.5 ± 2.474e+04 | 35260.6 ± 48.95 | ops/s | 1.17x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.27332e+06 ± 8749 | 2.25564e+06 ± 8501 | ops/s | 1.01x | ahead |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 219357 ± 8404 | 225022 ± 45.79 | ops/s | 0.97x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97654.9 ± 1109 | 34112.4 ± 563.7 | ops/s | 2.86x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 3711.02 ± 18.88 | 1770.13 ± 15.8 | ops/s | 2.10x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 575.504 ± 14.73 | 315.919 ± 3.013 | ops/s | 1.82x | ahead |  |

## Delta against Milestone 1

Milestone 1 is `2026-10-09-breeze-milestone1-trillium-laneA.md` (source ff03eed). It used the same lane, CPU model, JDK, pinning and runner protocol, on a different node of the same type. Milestone 1 was a legacy receipt with no numerical preflight. The `Pinv` fixture changed between the runs (seed 2200 to seed 2100; see `benchmarks/README.md`), so its rows carry no change ratio.

Ratios are gale speed / Breeze speed (>1 means gale is faster). `change` is new ratio / old ratio. Old is milestone 1 ff03eed; new is milestone 2 0ff6440.

| class | op | params | gale backend | old ratio | new ratio | change | old verdict | new verdict |
|---|---|---|---|---:|---:|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | 0.97x | 1.04x | 1.07x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | 0.99x | 0.95x | 0.96x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | 1.05x | 1.12x | 1.07x | tie | tie |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | 0.97x | 1.02x | 1.05x | tie | tie |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | 0.98x | 0.98x | 1.00x | behind | behind |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | 0.97x | 0.95x | 0.98x | tie | behind |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | 11.76x | 11.76x | 1.00x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | 13.64x | 13.59x | 1.00x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | 14.28x | 14.40x | 1.01x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 1.26x | 1.27x | 1.01x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 1.29x | 1.31x | 1.02x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 1.26x | 1.25x | 0.99x | ahead | ahead |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 1.36x | 1.37x | 1.01x | ahead | ahead |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 1.36x | 1.36x | 1.00x | ahead | ahead |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 1.21x | 1.19x | 0.98x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=16 | pure | 1.94x | 1.95x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 1.87x | 1.90x | 1.02x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=64 | pure | 2.45x | 2.46x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.20x | 1.19x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 1.22x | 1.23x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.12x | 1.11x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 1.22x | 1.22x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 1.29x | 1.29x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.29x | 1.42x | 1.10x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.11x | 1.94x | 1.75x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 0.96x | 1.51x | 1.57x | behind | ahead |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | 0.02x | 2.57x | 128.50x | behind | ahead |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | 0.10x | 14.37x | 143.70x | behind | ahead |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | 0.01x | 1.50x | 150.00x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 0.61x | 1.36x | 2.23x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 0.32x | 2.27x | 7.09x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 0.42x | 1.78x | 4.24x | behind | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | 2.27x | 2.98x | fixture changed | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | 1.01x | 1.14x | fixture changed | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | 1.08x | 1.27x | fixture changed | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | 2.34x | 2.83x | 1.21x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | 0.95x | 1.17x | 1.23x | behind | ahead |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | 1.07x | 1.30x | 1.21x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | 1.28x | 1.29x | 1.01x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | 1.15x | 1.15x | 1.00x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | 1.28x | 1.26x | 0.98x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | 1.15x | 1.17x | 1.02x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | 1.28x | 1.29x | 1.01x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | 1.16x | 1.16x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.13x | 3.57x | 1.14x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.40x | 1.64x | 1.17x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.12x | 2.86x | 1.35x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.19x | 1.36x | 1.14x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.11x | 1.94x | 1.75x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 0.93x | 1.47x | 1.58x | behind | ahead |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.23x | 1.45x | 1.18x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.44x | 1.48x | 1.03x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 0.93x | 1.09x | 1.17x | behind | ahead |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.21x | 1.31x | 1.08x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.09x | 1.67x | 1.53x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 0.93x | 1.37x | 1.47x | behind | ahead |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | 1.01x | 1.04x | 1.03x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | 1.20x | 1.29x | 1.08x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | 0.69x | 1.38x | 2.00x | behind | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | 0.73x | 1.43x | 1.96x | behind | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | 1.41x | 1.36x | 0.96x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | 1.20x | 1.23x | 1.03x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | 1.07x | 1.85x | 1.73x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | 1.08x | 1.86x | 1.72x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.52x | 1.51x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.41x | 1.42x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | 1.06x | 1.83x | 1.73x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | 1.08x | 1.82x | 1.69x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | 1.10x | 1.11x | 1.01x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | 0.90x | 0.90x | 1.00x | withheld | withheld |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | 0.75x | 0.75x | 1.00x | behind | behind |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | 0.73x | 0.73x | 1.00x | withheld | withheld |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.77x | 0.86x | 1.12x | behind | behind |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.19x | 1.21x | 1.02x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.04x | 1.23x | 1.18x | ahead | ahead |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | 2.74x | 2.82x | 1.03x | ahead | ahead |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | 0.76x | 1.28x | 1.68x | behind | ahead |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | 14.68x | 24.07x | 1.64x | ahead | ahead |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | 11.79x | 11.82x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | 3.73x | 3.86x | 1.03x | ahead | ahead |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.46x | 3.46x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | 1.58x | 1.57x | 0.99x | ahead | ahead |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | 18.27x | 17.79x | 0.97x | ahead | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | 0.67x | 2.06x | 3.07x | behind | ahead |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | 0.74x | 2.57x | 3.47x | behind | ahead |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | 2.64x | 2.67x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | 1.58x | 1.64x | 1.04x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | 1.66x | 1.66x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | tie |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | 1.50x | 1.53x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | 1.00x | 0.98x | 0.98x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | 0.97x | 1.00x | 1.03x | behind | tie |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | 1.04x | 1.06x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | 1.03x | 1.05x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | 0.70x | 1.06x | 1.51x | behind | tie |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | 0.72x | 0.95x | 1.32x | behind | tie |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | 0.63x | 0.94x | 1.49x | behind | tie |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | 30.51x | 30.50x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | 30.97x | 30.98x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | 31.19x | 31.22x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | 40.53x | 40.53x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | 40.96x | 40.92x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | 40.96x | 40.94x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | 1.23x | 1.24x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | 1.64x | 1.64x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | 1.63x | 1.63x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | 0.85x | 0.85x | 1.00x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | 0.71x | 0.71x | 1.00x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | 0.53x | 0.52x | 0.98x | behind | behind |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | 1.79x | 1.83x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | 2.15x | 2.18x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | 1.80x | 1.83x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.27x | 3.27x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | 3.38x | 3.38x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | 3.57x | 3.56x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | pure | 3.62x | 3.59x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | pure | 2.79x | 2.87x | 1.03x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | 1.23x | 1.27x | 1.03x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | 2.15x | 2.12x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | pure | 3.03x | 2.92x | 0.96x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | pure | 2.22x | 2.87x | 1.29x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | pure | 2.80x | 3.06x | 1.09x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | pure | 2.37x | 2.32x | 0.98x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | 5.73x | 5.74x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | 5.46x | 5.46x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | 6.52x | 6.51x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | 5.73x | 5.73x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | 7.01x | 7.06x | 1.01x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | 6.38x | 6.35x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | 8.07x | 7.97x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | 6.07x | 6.06x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | 9.36x | 9.31x | 0.99x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | 9.21x | 9.24x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | 10.08x | 10.26x | 1.02x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | 8.53x | 8.59x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | 8.17x | 8.18x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | 6.97x | 6.97x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | 6.99x | 6.97x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | 6.66x | 6.70x | 1.01x | withheld | withheld |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | 1.05x | 1.06x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | 1.11x | 1.10x | 0.99x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | 6.65x | 6.51x | 0.98x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | 6.37x | 6.36x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | 3.54x | 3.45x | 0.97x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | 1.16x | 1.17x | 1.01x | tie | tie |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | 1.00x | 1.01x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | 1.00x | 0.97x | 0.97x | tie | tie |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | 1.07x | 1.82x | 1.70x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | 2.24x | 2.86x | 1.28x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | 1.34x | 2.10x | 1.57x | ahead | ahead |

## Supplementary: controlled n=16 eigen recheck

Run separately on debug node tri0035 (Slurm job 2535336, exclusive), same bundle, SHA binding, pinning and lane-A fork arguments, through `supp.sh eig16A` (`docs/verification/w21-simd-spike/x86-avx512/supp.sh`). It used more samples than the sweep: 4 forks × 10 warmup × 10 measurement iterations of 1 s, `-p n=16`, `DenseDecompositionBreezeJmh.(gale|breeze)Eig` and `SymEigenBreezeJmh.(gale|breeze)EigSym`. The scoreboard ran with `--validation --strict` and exited 0 with 4 rows validated. Files: `2026-10-10-breeze-milestone2-trillium-laneA-eig16.{json,netlib.jsonl,validation.jsonl}`.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 34219.3 ± 545.5 | 12902 ± 31.34 | ops/s | 2.65x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97710.3 ± 274.8 | 34174.1 ± 143.9 | ops/s | 2.86x | ahead |  |

## Supplementary: JIT-fair nonsymmetric eig (lane A, `-XX:-DontCompileHugeMethods`)

Breeze's F2J `dgeev` / `dhseqr` path contains methods above HotSpot's huge-method limit (8000 bytecodes). By default these are never JIT-compiled and stay partly interpreted, which inflates the lane-A nonsymmetric-eig margins in the sweep above. Gale's eig now sits below the limit (0ff6440). This supplementary run passes `-XX:-DontCompileHugeMethods` to every fork, on both sides, so Breeze's huge methods can compile too. It used the default annotations (2 forks), `supp.sh hugeA`, in the same Slurm job as the lane-A sweep (2535351, tri0464), immediately after it. The scoreboard ran with `--validation --strict` and exited 0 with 6 rows validated. Files: `2026-10-10-breeze-milestone2-trillium-laneA-hugemethods.{json,netlib.jsonl,validation.jsonl}`. **Report these ratios as the JIT-fair nonsymmetric-eig comparison.**

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 34163.7 ± 1313 | 21814.3 ± 64.58 | ops/s | 1.57x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1086.79 ± 3.516 | 1026.2 ± 9.563 | ops/s | 1.06x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 19.0234 ± 0.0285 | 14.3037 ± 0.8635 | ops/s | 1.33x | ahead |  |

## Provenance

- Date: Slurm job 2535351, start/end/elapsed: 2026-10-10T21:13:01 2026-10-10T23:06:48 01:53:47 
- Benchmarked source: `0ff64405ea981f5a6dbf026a12e2bcc81617f412` (breeze/integration, clean tree). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster. Every fork received `-Dgale.bench.sourceSha=0ff64405ea981f5a6dbf026a12e2bcc81617f412` from the bundle's `SHA` file.
- Node: SciNet trillium `tri0464`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5), AVX-512, SMT off, 192 CPUs, NPS4, boost enabled. Load average at job start: 2258.93, 2296.22, 1728.79.
- JDK: Temurin 25+36 (`module load java/25`). `JAVA_TOOL_OPTIONS`, `JDK_JAVA_OPTIONS` and `_JAVA_OPTIONS` were unset.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the CCD of cpu0 and its NUMA node). It applies to the host JVM and is inherited by every fork.
- Netlib (blas, lapack, vectorModule) over all sidecar records (`2026-10-10-breeze-milestone2-trillium-laneA.netlib.jsonl`): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False)].
- Numerical preflight (`2026-10-10-breeze-milestone2-trillium-laneA.validation.jsonl`): {'passed': 220} sourceSha: ['0ff64405ea981f5a6dbf026a12e2bcc81617f412'].
- Runner: the reviewed `docs/verification/w21-simd-spike/x86-avx512/run.sh laneA` (unmodified, `GALE_BENCH_ROOT=/scratch/brad/gale-bench/m2`), i.e. `numactl … java --add-modules=jdk.incubator.vector -Dgale.bench.sourceSha=<sha> -cp <classpath> org.openjdk.jmh.Main -foe true -jvmArgs "-Xms4g -Xmx4g -Dgale.bench.sourceSha=<sha> -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings (2 forks etc.) come from the annotations.
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane A --strict --commit 0ff64405ea981f5a6dbf026a12e2bcc81617f412 --netlib 2026-10-10-breeze-milestone2-trillium-laneA.netlib.jsonl --validation 2026-10-10-breeze-milestone2-trillium-laneA.validation.jsonl 2026-10-10-breeze-milestone2-trillium-laneA.json`. It exited 0, with 0 unpaired rows.
