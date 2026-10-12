# gale vs Breeze scoreboard — lane A (out-of-box scalar)

- Lane: A
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: -Xms4g -Xmx4g -Dgale.bench.sourceSha=96c95e577f6b437c3de995a081b8ba075a6054e4 -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/m3/out/laneA-2538726/laneA/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.Java11BLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `96c95e577f6b437c3de995a081b8ba075a6054e4`
- Numerical preflight: 110 decomposition/factorization result rows validated on this SHA.
- Machine: SciNet trillium tri0463: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512, 8 double lanes), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2538726 exclusive

**124 ahead, 13 tie, 4 behind, 0 n/a, 10 withheld** of 151 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 178469 ± 6.959e+04 | 187164 ± 3.602e+04 | ops/s | 0.95x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 29410.6 ± 3057 | 29185.5 ± 2766 | ops/s | 1.01x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7172.44 ± 784.6 | 7609.52 ± 763.4 | ops/s | 0.94x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 65460.4 ± 784.2 | 67497.1 ± 76.53 | ops/s | 0.97x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16209.9 ± 979.7 | 16950.7 ± 149.4 | ops/s | 0.96x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4202.93 ± 14.31 | 4243.21 ± 41.85 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 66484.5 ± 1730 | 4647.43 ± 13.88 | ops/s | 14.31x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 16828.8 ± 265.8 | 1239.4 ± 258 | ops/s | 13.58x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 4199.86 ± 24.93 | 357.421 ± 0.1816 | ops/s | 11.75x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 98684.9 ± 477.3 | 79078.4 ± 1944 | ops/s | 1.25x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6512.89 ± 12.99 | 5097.73 ± 13.94 | ops/s | 1.28x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1612.59 ± 19.38 | 1220.66 ± 10.01 | ops/s | 1.32x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 82874.9 ± 195.7 | 68880.8 ± 45.89 | ops/s | 1.20x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5727.07 ± 8.637 | 4219.04 ± 10.26 | ops/s | 1.36x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1420.81 ± 5.403 | 1046.77 ± 1.367 | ops/s | 1.36x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 774071 ± 1.855e+04 | 397292 ± 9261 | ops/s | 1.95x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15278.1 ± 78.97 | 6258.28 ± 43.47 | ops/s | 2.44x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 215.715 ± 5.578 | 113.371 ± 0.0688 | ops/s | 1.90x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.78606e+06 ± 8011 | 1.50336e+06 ± 1.622e+04 | ops/s | 1.19x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31726.5 ± 321.8 | 28332.1 ± 184.8 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 517.783 ± 5.022 | 417.692 ± 3.203 | ops/s | 1.24x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 436170 ± 1043 | 397773 ± 1460 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7857.56 ± 61.99 | 6099.58 ± 36.11 | ops/s | 1.29x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 128.777 ± 1.685 | 105.108 ± 0.9533 | ops/s | 1.23x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.3677e+06 ± 4070 | 944059 ± 1.844e+04 | ops/s | 1.45x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 47651.8 ± 264.4 | 31534.7 ± 211.1 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 870.107 ± 79.62 | 475.344 ± 2.519 | ops/s | 1.83x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33197.7 ± 299.5 | 12878.2 ± 61.74 | ops/s | 2.58x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1084.49 ± 0.941 | 720.489 ± 6.567 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 18.9002 ± 0.05056 | 1.31477 ± 0.02058 | ops/s | 14.38x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 465548 ± 1.24e+04 | 340246 ± 8001 | ops/s | 1.37x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 17386.4 ± 108.2 | 9789.1 ± 24.46 | ops/s | 1.78x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 362.227 ± 4.136 | 164.439 ± 0.04758 | ops/s | 2.20x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 45538.7 ± 1173 | 15310.1 ± 213.6 | ops/s | 2.97x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 1374.92 ± 11.86 | 1076.96 ± 12.82 | ops/s | 1.28x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 27.6434 ± 0.04689 | 24.2384 ± 0.03626 | ops/s | 1.14x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 51389.5 ± 836.1 | 18562.9 ± 179.1 | ops/s | 2.77x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 1686.35 ± 11.75 | 1292.03 ± 11.96 | ops/s | 1.31x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 32.9521 ± 0.0447 | 28.1554 ± 0.07188 | ops/s | 1.17x | ahead |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 44260.3 ± 438.7 | 37666.6 ± 199.2 | ops/s | 1.18x | ahead |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 2063.31 ± 17.33 | 1608.02 ± 22.59 | ops/s | 1.28x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 44195.7 ± 379.1 | 38052.9 ± 416.5 | ops/s | 1.16x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 2044.28 ± 8.556 | 1615.05 ± 24.92 | ops/s | 1.27x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 43610.4 ± 419.4 | 38206.4 ± 95.82 | ops/s | 1.14x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 2034.87 ± 9.786 | 1621.36 ± 10 | ops/s | 1.26x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.46053e+06 ± 3355 | 405552 ± 6558 | ops/s | 3.60x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 55876.9 ± 353.9 | 20205 ± 77.21 | ops/s | 2.77x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1072.65 ± 12.28 | 653.943 ± 6.526 | ops/s | 1.64x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.40703e+06 ± 1.495e+04 | 1.04176e+06 ± 5787 | ops/s | 1.35x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 47588.5 ± 276.3 | 32276.7 ± 333 | ops/s | 1.47x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 924.305 ± 10.29 | 475.784 ± 2.433 | ops/s | 1.94x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 773808 ± 1.205e+04 | 368128 ± 3770 | ops/s | 2.10x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 16777.1 ± 171.3 | 15381.9 ± 248.7 | ops/s | 1.09x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 340.76 ± 0.812 | 230.765 ± 1.982 | ops/s | 1.48x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 1.07547e+06 ± 2.136e+04 | 787011 ± 2.143e+04 | ops/s | 1.37x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 41265 ± 250.3 | 29943.2 ± 202.8 | ops/s | 1.38x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 787.147 ± 5.251 | 469.665 ± 2.42 | ops/s | 1.68x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.2104 ± 0.01123 | 10.4887 ± 0.00906 | ms/op | 1.28x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.6164 ± 0.1905 | 74.2817 ± 0.2921 | ms/op | 1.04x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 104.319 ± 1.854 | 144.215 ± 0.2098 | ms/op | 1.38x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 772.359 ± 5.799 | 1078.17 ± 3.141 | ms/op | 1.40x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 77.5877 ± 0.7946 | 93.1493 ± 2.513 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 535.852 ± 5.149 | 736.376 ± 23.07 | ms/op | 1.37x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 8.84362 ± 0.03005 | 16.4023 ± 0.01208 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 69.6041 ± 0.1938 | 128.798 ± 0.1417 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 25.9314 ± 0.0754 | 36.7898 ± 0.1385 | ms/op | 1.42x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 194.495 ± 3.625 | 289.078 ± 7.962 | ms/op | 1.49x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 9.13491 ± 0.0408 | 16.5078 ± 0.02352 | ms/op | 1.81x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 70.5868 ± 0.3249 | 129.271 ± 0.1598 | ms/op | 1.83x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 815.697 ± 7.908 | 932.429 ± 15.33 | us/op | 1.14x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1894.91 ± 18.61 | 1766.24 ± 32.49 | us/op | 0.93x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 45.1146 ± 0.2034 | 77.7914 ± 0.4456 | us/op | 1.72x | ahead |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 1238.45 ± 4.886 | 1871.51 ± 18.25 | us/op | 1.51x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 134782 ± 2233 | 105840 ± 767.8 | ops/s | 1.27x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2940.3 ± 15.1 | 2292.61 ± 4.808 | ops/s | 1.28x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 45.6575 ± 0.2991 | 38.0479 ± 0.8664 | ops/s | 1.20x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 291.365 ± 0.07458 | 103.76 ± 0.9512 | ops/s | 2.81x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 3574.05 ± 60.5 | 3083.48 ± 435.2 | ops/s | 1.16x | tie | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 4621.02 ± 18.26 | 187.912 ± 1.132 | ops/s | 24.59x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4244.48 ± 19.06 | 357.162 ± 0.09516 | ops/s | 11.88x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 253.436 ± 0.68 | 64.956 ± 2.263 | ops/s | 3.90x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7428.67 ± 4.382 | 2145.95 ± 0.7686 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 3540.26 ± 70.11 | 2269.14 ± 8.674 | ops/s | 1.56x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7187.75 ± 3.312 | 390.895 ± 20.04 | ops/s | 18.39x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 658.509 ± 8.711 | 319.728 ± 0.8647 | ops/s | 2.06x | ahead |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 619.664 ± 3.346 | 245.285 ± 6.426 | ops/s | 2.53x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.40713e+06 ± 5.067e+04 | 1.28014e+06 ± 1630 | ops/s | 2.66x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 53793.2 ± 44.57 | 32555.8 ± 4522 | ops/s | 1.65x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3387.62 ± 3.073 | 2140.71 ± 3.67 | ops/s | 1.58x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 340926 ± 808.1 | 341208 ± 1259 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5021.92 ± 55.77 | 5012.11 ± 56.02 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 323.608 ± 0.4378 | 216.999 ± 0.3323 | ops/s | 1.49x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 300307 ± 930 | 301260 ± 1141 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4828.42 ± 13.91 | 4608.12 ± 11.94 | ops/s | 1.05x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 299.096 ± 4.256 | 283.547 ± 0.8746 | ops/s | 1.05x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 5.62406e+06 ± 5.047e+05 | 5.22969e+06 ± 1.505e+05 | ops/s | 1.08x | tie |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 78077.1 ± 1.686e+04 | 72557.5 ± 7569 | ops/s | 1.08x | tie |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4225.91 ± 3.417 | 5543.38 ± 332 | ops/s | 0.76x | behind |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.49719e+06 ± 1069 | 245562 ± 246.6 | ops/s | 30.53x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118161 ± 21.64 | 3789 ± 0.7185 | ops/s | 31.19x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7352.53 ± 31.08 | 237.264 ± 0.07126 | ops/s | 30.99x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.36118e+06 ± 1297 | 107642 ± 49.98 | ops/s | 40.52x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68566 ± 10.98 | 1674.33 ± 3.147 | ops/s | 40.95x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4273.67 ± 1.5 | 104.361 ± 0.04779 | ops/s | 40.95x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.86964e+06 ± 9796 | 2.31916e+06 ± 4993 | ops/s | 1.24x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55597 ± 17.58 | 34140.7 ± 30.39 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3498.56 ± 1.982 | 2129.72 ± 0.7619 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 278961 ± 1639 | 290778 ± 1031 | ops/s | 0.96x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 4164.28 ± 35.82 | 4289.44 ± 15.44 | ops/s | 0.97x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 264.022 ± 0.6472 | 186.397 ± 0.2735 | ops/s | 1.42x | ahead |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 269738 ± 4456 | 147852 ± 2273 | ops/s | 1.82x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3980.49 ± 24.22 | 2170.46 ± 255.4 | ops/s | 1.83x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 256.029 ± 0.7024 | 117.21 ± 0.4665 | ops/s | 2.18x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.66707e+06 ± 1.741e+04 | 2.3442e+06 ± 4433 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 122890 ± 816.5 | 34389.5 ± 30.77 | ops/s | 3.57x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7241.14 ± 25.36 | 2146.41 ± 0.9092 | ops/s | 3.37x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.8239e+07 ± 7.097e+04 | 7.65217e+06 ± 1.305e+05 | ops/s | 3.69x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.88645e+07 ± 1.958e+06 | 6.63636e+06 ± 1.722e+04 | ops/s | 2.84x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.12498e+07 ± 2.051e+05 | 2.48481e+07 ± 7.515e+05 | ops/s | 1.26x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.81162e+07 ± 1.601e+06 | 1.77847e+07 ± 3.469e+05 | ops/s | 2.14x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 9.87179e+06 ± 1.177e+04 | 3.32575e+06 ± 7.207e+04 | ops/s | 2.97x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 8.28318e+06 ± 2.625e+05 | 2.78387e+06 ± 2.909e+04 | ops/s | 2.98x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.583e+07 ± 2.535e+05 | 5.46582e+06 ± 5.176e+04 | ops/s | 2.90x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.12167e+07 ± 8.052e+04 | 4.75671e+06 ± 1.586e+05 | ops/s | 2.36x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.230103 ± 0.0006735 | 1.32992 ± 0.008756 | ms/op | 5.78x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.3302 ± 0.09494 | 133.904 ± 0.412 | ms/op | 5.50x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.85066 ± 0.005471 | 12.0577 ± 0.03633 | ms/op | 6.52x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 222.827 ± 0.9248 | 1276.58 ± 74.46 | ms/op | 5.73x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00578915 ± 2.347e-05 | 0.0406048 ± 0.0004838 | ms/op | 7.01x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.595366 ± 0.002092 | 3.86725 ± 0.1502 | ms/op | 6.50x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0465821 ± 0.0002329 | 0.370606 ± 0.002 | ms/op | 7.96x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.1913 ± 0.001807 | 37.5547 ± 0.02451 | ms/op | 6.07x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140691 ± 0.0003203 | 1.31426 ± 0.01737 | ms/op | 9.34x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.6044 ± 0.09684 | 133.883 ± 0.1513 | ms/op | 9.17x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19677 ± 0.007083 | 12.0672 ± 0.03205 | ms/op | 10.08x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 145.603 ± 5.146 | 1274.15 ± 82.58 | ms/op | 8.75x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00497075 ± 8.995e-06 | 0.0407935 ± 5.192e-05 | ms/op | 8.21x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.543703 ± 0.002865 | 3.79305 ± 0.004615 | ms/op | 6.98x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0532194 ± 0.0001792 | 0.369588 ± 0.001226 | ms/op | 6.94x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.62393 ± 0.01572 | 37.6054 ± 0.1011 | ms/op | 6.69x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 211949 ± 1532 | 191978 ± 9497 | ops/s | 1.10x | ahead |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 17717.5 ± 51.76 | 15977.3 ± 275.6 | ops/s | 1.11x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.0777e+06 ± 2.306e+04 | 314889 ± 823.1 | ops/s | 6.60x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 207477 ± 6220 | 32537.3 ± 1371 | ops/s | 6.38x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.13297e+06 ± 6.088e+04 | 324106 ± 557 | ops/s | 3.50x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 40373.4 ± 2.278e+04 | 35138.2 ± 53.5 | ops/s | 1.15x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.26188e+06 ± 3923 | 2.25174e+06 ± 804 | ops/s | 1.00x | ahead |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 222646 ± 3242 | 223760 ± 140.9 | ops/s | 1.00x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97212.8 ± 835.8 | 34065.2 ± 517.5 | ops/s | 2.85x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 3706.14 ± 16.28 | 1766.47 ± 8.402 | ops/s | 2.10x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 584.76 ± 16.31 | 314.132 ± 2.936 | ops/s | 1.86x | ahead |  |

## Delta against Milestone 2

Milestone 2 is `2026-10-10-breeze-milestone2-trillium-laneA.md` (source 0ff6440, receipts 5d3950e). It used the same lane, CPU model, JDK, pinning, runner and preflight, on a different node of the same type. This is an unpaired, cross-run comparison: Breeze is re-measured in each run, so read the gale/Breeze ratio changes rather than gale's absolute times.

Ratios are gale speed / Breeze speed (>1 means gale is faster). `change` is new ratio / old ratio. Old is milestone 2 0ff6440; new is milestone 3 96c95e5.

| class | op | params | gale backend | old ratio | new ratio | change | old verdict | new verdict |
|---|---|---|---|---:|---:|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | 1.04x | 0.94x | 0.90x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | 0.95x | 1.01x | 1.06x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | 1.12x | 0.95x | 0.85x | tie | tie |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | 1.02x | 0.99x | 0.97x | tie | tie |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | 0.98x | 0.96x | 0.98x | behind | tie |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | 0.95x | 0.97x | 1.02x | behind | behind |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | 11.76x | 11.75x | 1.00x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | 13.59x | 13.58x | 1.00x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | 14.40x | 14.31x | 0.99x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 1.27x | 1.28x | 1.01x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 1.31x | 1.32x | 1.01x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 1.25x | 1.25x | 1.00x | ahead | ahead |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 1.37x | 1.36x | 0.99x | ahead | ahead |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 1.36x | 1.36x | 1.00x | ahead | ahead |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 1.19x | 1.20x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=16 | pure | 1.95x | 1.95x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 1.90x | 1.90x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=64 | pure | 2.46x | 2.44x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.19x | 1.19x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 1.23x | 1.24x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.11x | 1.10x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 1.22x | 1.23x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 1.29x | 1.29x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.42x | 1.45x | 1.02x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.94x | 1.83x | 0.94x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 1.51x | 1.51x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | 1.57x | 2.58x | 1.64x | ahead | ahead |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | 1.33x | 14.38x | 10.81x | ahead | ahead |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | 1.06x | 1.51x | 1.42x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 1.36x | 1.37x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 2.27x | 2.20x | 0.97x | ahead | ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 1.78x | 1.78x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | 2.98x | 2.97x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | 1.14x | 1.14x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | 1.27x | 1.28x | 1.01x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | 2.83x | 2.77x | 0.98x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | 1.17x | 1.17x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | 1.30x | 1.31x | 1.01x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | 1.29x | 1.28x | 0.99x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | 1.15x | 1.18x | 1.03x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | 1.26x | 1.27x | 1.01x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | 1.17x | 1.16x | 0.99x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | 1.29x | 1.26x | 0.98x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | 1.16x | 1.14x | 0.98x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.57x | 3.60x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.64x | 1.64x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.86x | 2.77x | 0.97x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.36x | 1.35x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.94x | 1.94x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 1.47x | 1.47x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.45x | 2.10x | 1.45x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.48x | 1.48x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 1.09x | 1.09x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.31x | 1.37x | 1.05x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.67x | 1.68x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 1.37x | 1.38x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | 1.04x | 1.04x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | 1.29x | 1.28x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | 1.38x | 1.40x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | 1.43x | 1.38x | 0.97x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | 1.36x | 1.37x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | 1.23x | 1.20x | 0.98x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | 1.85x | 1.85x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | 1.86x | 1.85x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.51x | 1.49x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.42x | 1.42x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | 1.83x | 1.83x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | 1.82x | 1.81x | 0.99x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | 1.11x | 1.14x | 1.03x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | 0.90x | 0.93x | 1.03x | withheld | withheld |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | 0.75x | 1.72x | 2.29x | behind | ahead |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | 0.73x | 1.51x | 2.07x | withheld | withheld |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.86x | 1.27x | 1.48x | behind | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.21x | 1.20x | 0.99x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.23x | 1.28x | 1.04x | ahead | ahead |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | 2.82x | 2.81x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | 1.28x | 1.16x | 0.91x | ahead | tie |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | 24.07x | 24.59x | 1.02x | ahead | ahead |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | 11.82x | 11.88x | 1.01x | ahead | ahead |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | 3.86x | 3.90x | 1.01x | ahead | ahead |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.46x | 3.46x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | 1.57x | 1.56x | 0.99x | ahead | ahead |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | 17.79x | 18.39x | 1.03x | ahead | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | 2.06x | 2.06x | 1.00x | ahead | ahead |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | 2.57x | 2.53x | 0.98x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | 2.67x | 2.66x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | 1.64x | 1.58x | 0.96x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | 1.66x | 1.65x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | tie |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | 1.53x | 1.49x | 0.97x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | 0.98x | 1.00x | 1.02x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | 1.06x | 1.05x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | 1.05x | 1.05x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | 1.06x | 1.08x | 1.02x | tie | tie |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | 0.95x | 0.76x | 0.80x | tie | behind |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | 0.94x | 1.08x | 1.15x | tie | tie |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | 30.50x | 30.53x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | 30.98x | 30.99x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | 31.22x | 31.19x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | 40.53x | 40.52x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | 40.92x | 40.95x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | 40.94x | 40.95x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | 1.24x | 1.24x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | 1.64x | 1.64x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | 1.63x | 1.63x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | 0.85x | 0.96x | 1.13x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | 0.71x | 1.42x | 2.00x | behind | ahead |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | 0.52x | 0.97x | 1.87x | behind | behind |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | 1.83x | 1.82x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | 2.18x | 2.18x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | 1.83x | 1.83x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.27x | 3.27x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | 3.38x | 3.37x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | 3.56x | 3.57x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | pure | 3.59x | 3.69x | 1.03x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | pure | 2.87x | 2.84x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | 1.27x | 1.26x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | 2.12x | 2.14x | 1.01x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | pure | 2.92x | 2.97x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | pure | 2.87x | 2.98x | 1.04x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | pure | 3.06x | 2.90x | 0.95x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | pure | 2.32x | 2.36x | 1.02x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | 5.74x | 5.78x | 1.01x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | 5.46x | 5.50x | 1.01x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | 6.51x | 6.52x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | 5.73x | 5.73x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | 7.06x | 7.01x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | 6.35x | 6.50x | 1.02x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | 7.97x | 7.96x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | 6.06x | 6.07x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | 9.31x | 9.34x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | 9.24x | 9.17x | 0.99x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | 10.26x | 10.08x | 0.98x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | 8.59x | 8.75x | 1.02x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | 8.18x | 8.21x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | 6.97x | 6.98x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | 6.97x | 6.94x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | 6.70x | 6.69x | 1.00x | withheld | withheld |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | 1.06x | 1.10x | 1.04x | ahead | ahead |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | 1.10x | 1.11x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | 6.51x | 6.60x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | 6.36x | 6.38x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | 3.45x | 3.50x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | 1.17x | 1.15x | 0.98x | tie | tie |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | 1.01x | 1.00x | 0.99x | ahead | ahead |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | 0.97x | 1.00x | 1.03x | tie | tie |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | 1.82x | 1.86x | 1.02x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | 2.86x | 2.85x | 1.00x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | 2.10x | 2.10x | 1.00x | ahead | ahead |

## Supplementary: controlled n=16 eigen recheck

Run in the same Slurm job, immediately after the sweep, with `supp.sh eig16A` (`docs/verification/w21-simd-spike/x86-avx512/supp.sh`, unchanged from Milestone 2). Settings: 4 forks × 10 warmup × 10 measurement iterations of 1 s, `-p n=16`, `DenseDecompositionBreezeJmh.(gale|breeze)Eig` and `SymEigenBreezeJmh.(gale|breeze)EigSym`. The scoreboard ran with `--validation --strict`. Files: `2026-10-11-breeze-milestone3-trillium-laneA-eig16.{json,netlib.jsonl,validation.jsonl}`.

- Numerical preflight: 4 decomposition/factorization result rows validated on this SHA.
**2 ahead, 0 tie, 0 behind, 0 n/a, 0 withheld** of 2 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.
| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 34172.8 ± 519.3 | 12913.1 ± 32.47 | ops/s | 2.65x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97869.5 ± 243.1 | 34417.5 ± 79.68 | ops/s | 2.84x | ahead |  |

## Supplementary: JIT-fair nonsymmetric eig (lane A, `-XX:-DontCompileHugeMethods`)

Breeze's F2J `dgeev`/`dhseqr` contains methods above HotSpot's huge-method limit (8000 bytecodes), which stay partly interpreted by default and inflate the sweep's lane-A nonsymmetric-eig margins. This run passes `-XX:-DontCompileHugeMethods` to every fork on both sides, using `supp.sh hugeA` with 2 forks, in the same Slurm job right after the sweep. The scoreboard ran with `--validation --strict`. Files: `2026-10-11-breeze-milestone3-trillium-laneA-hugemethods.{json,netlib.jsonl,validation.jsonl}`. **Report these ratios as the JIT-fair nonsymmetric-eig comparison.**

- Numerical preflight: 6 decomposition/factorization result rows validated on this SHA.
**3 ahead, 0 tie, 0 behind, 0 n/a, 0 withheld** of 3 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.
| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33186.5 ± 297.1 | 21857.1 ± 141.9 | ops/s | 1.52x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1081.71 ± 0.8313 | 1023.36 ± 8.694 | ops/s | 1.06x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 18.8832 ± 0.04402 | 14.2852 ± 0.7325 | ops/s | 1.32x | ahead |  |

## Provenance

- Date: Slurm job 2538726, start/end/elapsed: 2026-10-11T16:30:57 2026-10-11T18:30:14 01:59:17 
- Benchmarked source: `96c95e577f6b437c3de995a081b8ba075a6054e4` (breeze/integration = main, clean worktree `bench/m3-96c95e5`). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster. Every fork received `-Dgale.bench.sourceSha=96c95e577f6b437c3de995a081b8ba075a6054e4` from the bundle's `SHA` file.
- Node: SciNet trillium `tri0463`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5), AVX-512, SMT off, 192 CPUs, NPS4, boost enabled. Load average at job start: 0.26, 2.69, 37.36.
- JDK: Temurin 25+36 (`module load java/25`). `JAVA_TOOL_OPTIONS`, `JDK_JAVA_OPTIONS` and `_JAVA_OPTIONS` were unset.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the CCD of cpu0 and its NUMA node). It applies to the host JVM and is inherited by every fork.
- Netlib (blas, lapack, vectorModule) over all sidecar records (`2026-10-11-breeze-milestone3-trillium-laneA.netlib.jsonl`): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False)].
- Numerical preflight (`2026-10-11-breeze-milestone3-trillium-laneA.validation.jsonl`): {'passed': 220} sourceSha: ['96c95e577f6b437c3de995a081b8ba075a6054e4'].
- Runner: the reviewed `docs/verification/w21-simd-spike/x86-avx512/run.sh laneA` (unmodified, `GALE_BENCH_ROOT=/scratch/brad/gale-bench/m3`), i.e. `numactl … java --add-modules=jdk.incubator.vector -Dgale.bench.sourceSha=<sha> -cp <classpath> org.openjdk.jmh.Main -foe true -jvmArgs "-Xms4g -Xmx4g -Dgale.bench.sourceSha=<sha> -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings (2 forks etc.) come from the annotations.
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane A --strict --commit 96c95e577f6b437c3de995a081b8ba075a6054e4 --netlib 2026-10-11-breeze-milestone3-trillium-laneA.netlib.jsonl --validation 2026-10-11-breeze-milestone3-trillium-laneA.validation.jsonl 2026-10-11-breeze-milestone3-trillium-laneA.json`. It exited 0, with 0 unpaired rows.
