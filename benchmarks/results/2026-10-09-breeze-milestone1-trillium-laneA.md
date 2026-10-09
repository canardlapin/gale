# gale vs Breeze scoreboard — lane A (out-of-box scalar)

- Lane: A
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: -Xms4g -Xmx4g -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/m1/out/laneA-2526674/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.Java11BLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `ff03eed494bfa2d246b1c363c59f3425d97fcd67`
- Machine: SciNet trillium tri0101: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512, 8 double lanes), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2526674 exclusive

**106 ahead, 9 tie, 26 behind, 0 n/a, 10 withheld** of 151 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 182785 ± 6.383e+04 | 173654 ± 7962 | ops/s | 1.05x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 28460.7 ± 1225 | 28846.1 ± 1248 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7376.8 ± 687.2 | 7591.78 ± 750.2 | ops/s | 0.97x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 65975 ± 2389 | 67696.2 ± 42.1 | ops/s | 0.97x | tie |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16746.9 ± 171.2 | 17113.1 ± 15.34 | ops/s | 0.98x | behind |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4121.51 ± 57 | 4232.23 ± 82.83 | ops/s | 0.97x | tie |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 66523.1 ± 2388 | 4656.9 ± 9.655 | ops/s | 14.28x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 16916.8 ± 126.8 | 1240.44 ± 266.3 | ops/s | 13.64x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 4218.12 ± 29.39 | 358.536 ± 0.2293 | ops/s | 11.76x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 99514 ± 1465 | 79109.3 ± 1828 | ops/s | 1.26x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6429.36 ± 154 | 5103.12 ± 2.945 | ops/s | 1.26x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1590.49 ± 22.19 | 1229.06 ± 7.437 | ops/s | 1.29x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 83513.8 ± 556.8 | 69136.8 ± 28.66 | ops/s | 1.21x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5786.57 ± 26.47 | 4241.62 ± 10.26 | ops/s | 1.36x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1424.8 ± 4.46 | 1051.27 ± 2.749 | ops/s | 1.36x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 775229 ± 1.642e+04 | 398687 ± 1.058e+04 | ops/s | 1.94x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15312.5 ± 84.07 | 6252.23 ± 60.86 | ops/s | 2.45x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 212.829 ± 3.934 | 113.625 ± 0.09738 | ops/s | 1.87x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.79552e+06 ± 4453 | 1.50085e+06 ± 1.171e+04 | ops/s | 1.20x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31817.6 ± 300.8 | 28420.1 ± 222 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 518.723 ± 4.414 | 423.811 ± 2.438 | ops/s | 1.22x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 439096 ± 2731 | 393239 ± 9501 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7881.33 ± 57.32 | 6113.77 ± 42.57 | ops/s | 1.29x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 129.181 ± 1.555 | 105.994 ± 0.2335 | ops/s | 1.22x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.21164e+06 ± 4792 | 938630 ± 1.093e+04 | ops/s | 1.29x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 30530.9 ± 161.1 | 31710.3 ± 203 | ops/s | 0.96x | behind |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 529.275 ± 3.444 | 475.791 ± 2.085 | ops/s | 1.11x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 288.861 ± 4.792 | 13130.9 ± 156.4 | ops/s | 0.02x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 7.2989 ± 0.004978 | 719.966 ± 7.673 | ops/s | 0.01x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 0.12878 ± 0.001995 | 1.31535 ± 0.02187 | ops/s | 0.10x | behind |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 203723 ± 3726 | 332717 ± 1.008e+04 | ops/s | 0.61x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 4167.91 ± 9.119 | 9813.31 ± 42.68 | ops/s | 0.42x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 52.344 ± 0.09364 | 164.522 ± 0.2721 | ops/s | 0.32x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 40140.6 ± 930.4 | 17655.7 ± 254.9 | ops/s | 2.27x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 1245.93 ± 14.48 | 1156.49 ± 16.09 | ops/s | 1.08x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 24.788 ± 0.03578 | 24.6594 ± 0.03998 | ops/s | 1.01x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 43391.1 ± 1387 | 18570.7 ± 237.4 | ops/s | 2.34x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 1394.22 ± 15.59 | 1297.87 ± 12.04 | ops/s | 1.07x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 26.8972 ± 0.03919 | 28.2424 ± 0.05302 | ops/s | 0.95x | behind |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 43635.1 ± 735.4 | 37884.2 ± 89.72 | ops/s | 1.15x | ahead |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 2066.91 ± 11.38 | 1621.01 ± 19.76 | ops/s | 1.28x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 43975.9 ± 1373 | 38171.1 ± 137.6 | ops/s | 1.15x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 2071.83 ± 24.78 | 1619.89 ± 17 | ops/s | 1.28x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 44440.3 ± 282.8 | 38278.7 ± 42.04 | ops/s | 1.16x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 2076.82 ± 21.21 | 1625.89 ± 14.22 | ops/s | 1.28x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.46231e+06 ± 7575 | 467494 ± 2161 | ops/s | 3.13x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 56029 ± 372.8 | 26474.9 ± 405.5 | ops/s | 2.12x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1075.06 ± 14.09 | 766.497 ± 8.759 | ops/s | 1.40x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.22778e+06 ± 2.859e+04 | 1.03559e+06 ± 7909 | ops/s | 1.19x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 30044.2 ± 163.2 | 32354.8 ± 274.2 | ops/s | 0.93x | behind |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 528.049 ± 4.511 | 476.551 ± 2.273 | ops/s | 1.11x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 458860 ± 5330 | 371884 ± 1.015e+04 | ops/s | 1.23x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 14349 ± 153.5 | 15502.8 ± 254.3 | ops/s | 0.93x | behind |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 332.795 ± 2.354 | 231.512 ± 1.451 | ops/s | 1.44x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 967209 ± 2163 | 797924 ± 1.715e+04 | ops/s | 1.21x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 27854.5 ± 254.6 | 30049 ± 252 | ops/s | 0.93x | behind |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 515.253 ± 3.423 | 470.66 ± 2.055 | ops/s | 1.09x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.17622 ± 0.01262 | 9.83793 ± 0.1007 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.396 ± 0.2051 | 72.0612 ± 0.2863 | ms/op | 1.01x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 196.186 ± 1.248 | 142.843 ± 0.3629 | ms/op | 0.73x | behind |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 1541.9 ± 1.917 | 1070.95 ± 1.285 | ms/op | 0.69x | behind |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 77.4166 ± 0.3272 | 92.7798 ± 2.56 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 531.19 ± 13.29 | 746.924 ± 0.5791 | ms/op | 1.41x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 15.0828 ± 0.02101 | 16.3533 ± 0.01843 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 120.12 ± 0.6552 | 128.461 ± 0.2197 | ms/op | 1.07x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 26.0555 ± 0.1579 | 36.6569 ± 0.09361 | ms/op | 1.41x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 192.855 ± 3.462 | 292.786 ± 0.7742 | ms/op | 1.52x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 15.2658 ± 0.01581 | 16.4629 ± 0.02456 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 121.166 ± 0.1642 | 128.884 ± 0.139 | ms/op | 1.06x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 841.062 ± 8.318 | 928.286 ± 17.79 | us/op | 1.10x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1953.57 ± 16.02 | 1758.99 ± 28.49 | us/op | 0.90x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 104.036 ± 1.897 | 77.7493 ± 0.4591 | us/op | 0.75x | behind |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 2559.73 ± 46.92 | 1868.05 ± 13.23 | us/op | 0.73x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 81887.1 ± 826.4 | 106154 ± 1531 | ops/s | 0.77x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2401.07 ± 22.34 | 2298.8 ± 6.655 | ops/s | 1.04x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 44.6019 ± 0.4328 | 37.6242 ± 0.0336 | ops/s | 1.19x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 284.88 ± 0.3222 | 104.152 ± 0.2423 | ops/s | 2.74x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 2379.02 ± 1.707 | 3122.15 ± 473.7 | ops/s | 0.76x | behind | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 2797.49 ± 57.94 | 190.613 ± 2.326 | ops/s | 14.68x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4230.31 ± 51.7 | 358.914 ± 0.04091 | ops/s | 11.79x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 245.173 ± 0.9178 | 65.741 ± 1.138 | ops/s | 3.73x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7459.12 ± 2.5 | 2155.45 ± 0.4974 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 3591.55 ± 30.18 | 2277 ± 8.654 | ops/s | 1.58x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7216.82 ± 5.251 | 395.104 ± 8.783 | ops/s | 18.27x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 216.452 ± 0.3207 | 321.524 ± 0.2041 | ops/s | 0.67x | behind |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 183.243 ± 1.408 | 246.519 ± 0.1258 | ops/s | 0.74x | behind |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.39324e+06 ± 6.333e+04 | 1.28343e+06 ± 3146 | ops/s | 2.64x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 53943.4 ± 123.9 | 32585.4 ± 4731 | ops/s | 1.66x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3396.81 ± 6.821 | 2145.39 ± 4.522 | ops/s | 1.58x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 342292 ± 1059 | 340631 ± 1516 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5061.64 ± 53.06 | 5083.48 ± 98.87 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 325.677 ± 0.7616 | 217.121 ± 0.4836 | ops/s | 1.50x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 294652 ± 190.6 | 302453 ± 857.9 | ops/s | 0.97x | behind |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4746.41 ± 21.72 | 4618.27 ± 7.974 | ops/s | 1.03x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 296.761 ± 0.4195 | 284.462 ± 0.405 | ops/s | 1.04x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 3.67805e+06 ± 1.209e+05 | 5.22808e+06 ± 1.481e+05 | ops/s | 0.70x | behind |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 59109.8 ± 1.092e+04 | 93586.2 ± 4546 | ops/s | 0.63x | behind |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4146.37 ± 41.72 | 5798.33 ± 163.7 | ops/s | 0.72x | behind |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.52925e+06 ± 6232 | 246791 ± 187.9 | ops/s | 30.51x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118724 ± 43.49 | 3805.89 ± 0.4567 | ops/s | 31.19x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7385.94 ± 23.9 | 238.477 ± 0.03149 | ops/s | 30.97x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.3787e+06 ± 1193 | 108036 ± 64.27 | ops/s | 40.53x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68808.1 ± 31.23 | 1679.69 ± 3.887 | ops/s | 40.96x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4289.44 ± 5.534 | 104.735 ± 0.03523 | ops/s | 40.96x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.86905e+06 ± 7173 | 2.33006e+06 ± 2566 | ops/s | 1.23x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55769.6 ± 54.74 | 34302.6 ± 33.42 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3512.88 ± 0.632 | 2140.62 ± 0.4757 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 247078 ± 5197 | 291640 ± 2469 | ops/s | 0.85x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 2303.22 ± 123.8 | 4318.33 ± 32.8 | ops/s | 0.53x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 133.315 ± 2.646 | 186.785 ± 0.3002 | ops/s | 0.71x | behind |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 262389 ± 5006 | 146889 ± 311.9 | ops/s | 1.79x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3916.22 ± 26.18 | 2176.83 ± 256.2 | ops/s | 1.80x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 253.181 ± 0.8694 | 117.522 ± 0.285 | ops/s | 2.15x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.69633e+06 ± 8549 | 2.35366e+06 ± 894.2 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 123077 ± 1058 | 34513.2 ± 33.66 | ops/s | 3.57x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7283.27 ± 9.468 | 2155.12 ± 0.9602 | ops/s | 3.38x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.74342e+07 ± 8.841e+04 | 7.58207e+06 ± 6.155e+05 | ops/s | 3.62x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.87652e+07 ± 4.625e+05 | 6.73074e+06 ± 1.158e+05 | ops/s | 2.79x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.14151e+07 ± 1.125e+05 | 2.56168e+07 ± 2.578e+05 | ops/s | 1.23x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.87649e+07 ± 4.957e+05 | 1.79969e+07 ± 5.249e+05 | ops/s | 2.15x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 1.04058e+07 ± 1.052e+05 | 3.42958e+06 ± 4.264e+04 | ops/s | 3.03x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 6.17407e+06 ± 4.257e+04 | 2.78578e+06 ± 6.89e+04 | ops/s | 2.22x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.5575e+07 ± 1.517e+05 | 5.55587e+06 ± 1.197e+05 | ops/s | 2.80x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.12449e+07 ± 5530 | 4.7499e+06 ± 9.681e+04 | ops/s | 2.37x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.229603 ± 0.0007464 | 1.3167 ± 0.01826 | ms/op | 5.73x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.4341 ± 0.1491 | 133.481 ± 0.3833 | ms/op | 5.46x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.84554 ± 0.006337 | 12.0326 ± 0.02392 | ms/op | 6.52x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 222.92 ± 0.3784 | 1278.38 ± 84.94 | ms/op | 5.73x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00576641 ± 1.458e-05 | 0.0404281 ± 0.0001925 | ms/op | 7.01x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.594511 ± 0.001745 | 3.79279 ± 0.01969 | ms/op | 6.38x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0457668 ± 0.0004443 | 0.369162 ± 0.002957 | ms/op | 8.07x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.16554 ± 0.01056 | 37.4458 ± 0.05937 | ms/op | 6.07x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140209 ± 0.0004667 | 1.31289 ± 0.01567 | ms/op | 9.36x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.5332 ± 0.3873 | 133.788 ± 0.5716 | ms/op | 9.21x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19292 ± 0.006251 | 12.026 ± 0.0376 | ms/op | 10.08x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 149.68 ± 0.0436 | 1276.09 ± 62.87 | ms/op | 8.53x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00494863 ± 6.257e-06 | 0.0404281 ± 0.0003362 | ms/op | 8.17x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.54372 ± 0.0006822 | 3.79063 ± 0.01566 | ms/op | 6.97x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0527268 ± 0.0004775 | 0.36865 ± 0.002854 | ms/op | 6.99x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.61942 ± 0.002122 | 37.4052 ± 0.03522 | ms/op | 6.66x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 211045 ± 572.1 | 200192 ± 694.8 | ops/s | 1.05x | ahead |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 17747.4 ± 53.33 | 16032.6 ± 63.14 | ops/s | 1.11x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.08386e+06 ± 9802 | 313187 ± 4802 | ops/s | 6.65x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 208987 ± 5006 | 32824.9 ± 1358 | ops/s | 6.37x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.15301e+06 ± 2.362e+04 | 325286 ± 697.5 | ops/s | 3.54x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 41053.2 ± 2.449e+04 | 35279 ± 55.33 | ops/s | 1.16x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.26792e+06 ± 552.3 | 2.2613e+06 ± 1249 | ops/s | 1.00x | ahead |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 223693 ± 1295 | 222884 ± 3031 | ops/s | 1.00x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 76458.9 ± 2388 | 34097.3 ± 562.1 | ops/s | 2.24x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 2346.08 ± 14.08 | 1754.36 ± 17.19 | ops/s | 1.34x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 336.644 ± 1.229 | 313.479 ± 3.005 | ops/s | 1.07x | ahead |  |

## Delta against the d102321 baseline

The baseline is `2026-10-09-breeze-two-lane-baseline-trillium-laneA.md` (source 17813d0, committed in d102321). It used the same lane, cluster, CPU model, JDK, pinning and method, on a different node of the same type.

Ratios are gale speed / Breeze speed (>1 means gale is faster). `change` is new ratio / old ratio. Old is baseline 17813d0 (d102321); new is milestone 1 ff03eed.

| class | op | params | gale backend | old ratio | new ratio | change | old verdict | new verdict |
|---|---|---|---|---:|---:|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | 0.41x | 0.97x | 2.37x | behind | tie |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | 0.43x | 0.99x | 2.30x | behind | tie |
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | 0.29x | 1.05x | 3.62x | behind | tie |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | 0.98x | 0.97x | 0.99x | tie | tie |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | 1.00x | 0.98x | 0.98x | tie | behind |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | 0.90x | 0.97x | 1.08x | behind | tie |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | 2.93x | 11.76x | 4.01x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | 3.39x | 13.64x | 4.02x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | 3.59x | 14.28x | 3.98x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 0.81x | 1.26x | 1.56x | behind | ahead |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 0.84x | 1.29x | 1.54x | behind | ahead |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 0.82x | 1.26x | 1.54x | behind | ahead |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 0.57x | 1.36x | 2.39x | behind | ahead |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 0.57x | 1.36x | 2.39x | behind | ahead |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 0.67x | 1.21x | 1.81x | behind | ahead |
| BlasL3BreezeJmh | AtA | n=16 | pure | 1.95x | 1.94x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 1.88x | 1.87x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=64 | pure | 2.47x | 2.45x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.18x | 1.20x | 1.02x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 1.23x | 1.22x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.11x | 1.12x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 1.23x | 1.22x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 1.29x | 1.29x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.29x | 1.29x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.11x | 1.11x | 1.00x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 0.97x | 0.96x | 0.99x | behind | behind |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | 0.02x | 0.02x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | 0.10x | 0.10x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | 0.01x | 0.01x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 0.60x | 0.61x | 1.02x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 0.32x | 0.32x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 0.43x | 0.42x | 0.98x | behind | behind |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | 1.77x | 2.27x | 1.28x | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | 0.12x | 1.01x | 8.42x | behind | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | 0.66x | 1.08x | 1.64x | behind | ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | 1.97x | 2.34x | 1.19x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | 0.11x | 0.95x | 8.64x | behind | behind |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | 0.75x | 1.07x | 1.43x | behind | ahead |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | 0.44x | 1.28x | 2.91x | behind | ahead |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | 0.86x | 1.15x | 1.34x | behind | ahead |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | 0.75x | 1.28x | 1.71x | behind | ahead |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | 0.43x | 1.15x | 2.67x | behind | ahead |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | 0.45x | 1.28x | 2.84x | behind | ahead |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | 0.86x | 1.16x | 1.35x | behind | ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.12x | 3.13x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.40x | 1.40x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.12x | 2.12x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.19x | 1.19x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.10x | 1.11x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 0.95x | 0.93x | 0.98x | behind | behind |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.22x | 1.23x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.43x | 1.44x | 1.01x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 0.92x | 0.93x | 1.01x | behind | behind |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.21x | 1.21x | 1.00x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.10x | 1.09x | 0.99x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 0.93x | 0.93x | 1.00x | behind | behind |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | 1.02x | 1.01x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | 1.20x | 1.20x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | 0.04x | 0.69x | 17.25x | behind | behind |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | 0.04x | 0.73x | 18.25x | behind | behind |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | 1.39x | 1.41x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | 1.22x | 1.20x | 0.98x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | 1.07x | 1.07x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | 1.08x | 1.08x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.49x | 1.52x | 1.02x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.41x | 1.41x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | 1.06x | 1.06x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | 1.08x | 1.08x | 1.00x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | 0.83x | 1.10x | 1.33x | behind | ahead |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | 0.67x | 0.90x | 1.34x | withheld | withheld |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | 0.75x | 0.75x | 1.00x | behind | behind |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | 0.72x | 0.73x | 1.01x | withheld | withheld |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.77x | 0.77x | 1.00x | behind | behind |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.19x | 1.19x | 1.00x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.02x | 1.04x | 1.02x | ahead | ahead |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | 2.68x | 2.74x | 1.02x | ahead | ahead |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | 0.82x | 0.76x | 0.93x | behind | behind |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | 11.16x | 14.68x | 1.32x | ahead | ahead |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | 11.64x | 11.79x | 1.01x | ahead | ahead |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | 3.03x | 3.73x | 1.23x | ahead | ahead |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.46x | 3.46x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | 1.83x | 1.58x | 0.86x | ahead | ahead |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | 18.34x | 18.27x | 1.00x | ahead | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | 0.67x | 0.67x | 1.00x | behind | behind |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | 0.74x | 0.74x | 1.00x | behind | behind |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | 2.68x | 2.64x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | 1.61x | 1.58x | 0.98x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | 1.66x | 1.66x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | 0.98x | 1.00x | 1.02x | behind | tie |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | 1.47x | 1.50x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | 0.99x | 1.00x | 1.01x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | 0.96x | 0.97x | 1.01x | behind | behind |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | 1.03x | 1.04x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | 1.01x | 1.03x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | 0.66x | 0.70x | 1.06x | behind | behind |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | 0.61x | 0.72x | 1.18x | behind | behind |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | 0.60x | 0.63x | 1.05x | behind | behind |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | 30.49x | 30.51x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | 30.97x | 30.97x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | 31.22x | 31.19x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | 40.56x | 40.53x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | 40.96x | 40.96x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | 41.01x | 40.96x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | 1.24x | 1.23x | 0.99x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | 1.64x | 1.64x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | 1.62x | 1.63x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | 0.89x | 0.85x | 0.96x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | 0.71x | 0.71x | 1.00x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | 0.53x | 0.53x | 1.00x | behind | behind |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | 1.47x | 1.79x | 1.22x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | 1.76x | 2.15x | 1.22x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | 1.51x | 1.80x | 1.19x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.27x | 3.27x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | 3.39x | 3.38x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | 3.56x | 3.57x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | pure | 3.50x | 3.62x | 1.03x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | pure | 2.89x | 2.79x | 0.97x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | 1.24x | 1.23x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | 2.13x | 2.15x | 1.01x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | pure | 3.08x | 3.03x | 0.98x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | pure | 2.18x | 2.22x | 1.02x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | pure | 2.83x | 2.80x | 0.99x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | pure | 2.39x | 2.37x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | 5.97x | 5.73x | 0.96x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | 5.46x | 5.46x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | 6.51x | 6.52x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | 5.77x | 5.73x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | 6.99x | 7.01x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | 6.35x | 6.38x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | 8.00x | 8.07x | 1.01x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | 6.05x | 6.07x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | 9.24x | 9.36x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | 9.09x | 9.21x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | 10.08x | 10.08x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | 8.80x | 8.53x | 0.97x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | 8.30x | 8.17x | 0.98x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | 6.93x | 6.97x | 1.01x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | 6.97x | 6.99x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | 6.64x | 6.66x | 1.00x | withheld | withheld |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | 0.95x | 1.05x | 1.11x | behind | ahead |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | 1.02x | 1.11x | 1.09x | tie | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | 6.60x | 6.65x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | 6.41x | 6.37x | 0.99x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | 3.49x | 3.54x | 1.01x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | 1.16x | 1.16x | 1.00x | tie | tie |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | ahead |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | tie |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | 0.54x | 1.07x | 1.98x | behind | ahead |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | 2.06x | 2.24x | 1.09x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | 0.97x | 1.34x | 1.38x | tie | ahead |

## Provenance

- Date: 2026-10-09, Slurm job 2526674, start/end/elapsed: 2026-10-09T16:11:45 2026-10-09T18:02:16 01:50:31 
- Benchmarked source: `ff03eed494bfa2d246b1c363c59f3425d97fcd67` (breeze/integration). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster.
- Node: SciNet trillium `tri0101`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5), AVX-512, SMT off, 192 CPUs, NPS4 (8 NUMA nodes of 24 cores), boost enabled. Load average at job start: 3.87, 87.47, 149.08.
- JDK: Temurin 25+36 (`module load java/25`). The site default `JAVA_TOOL_OPTIONS=-Xmx2g` was unset before launch.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the CCD of cpu0 and its NUMA node). It applies to the host JVM and is inherited by every fork.
- Netlib (blas, lapack, vectorModule) over all sidecar records (`2026-10-09-breeze-milestone1-trillium-laneA.netlib.jsonl`): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False)]. No native BLAS/LAPACK is visible to the loader on these nodes.
- Command, equivalent to `breezeLaneA`: `numactl … java --add-modules=jdk.incubator.vector -cp <classpath> org.openjdk.jmh.Main -jvmArgs "-Xms4g -Xmx4g -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings come from the annotations. The runner is `docs/verification/w21-simd-spike/x86-avx512/run.sh` (with `BUNDLE` pointing at the milestone bundle).
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane A --strict`. It exited 0, with 0 unpaired rows.
