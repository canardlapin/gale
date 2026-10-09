# gale vs Breeze scoreboard — lane A (out-of-box scalar)

- Lane: A
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: -Xms4g -Xmx4g -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/out/laneA-2524300/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.Java11BLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `17813d0d761467187ba0764f1081fd30747ad18e`
- Machine: SciNet trillium tri0819: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2524300 exclusive

**85 ahead, 8 tie, 48 behind, 0 n/a, 10 withheld** of 151 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 54436 ± 1824 | 190731 ± 4.109e+04 | ops/s | 0.29x | behind |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 12719.8 ± 2125 | 29731 ± 2184 | ops/s | 0.43x | behind |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 3277.34 ± 396.3 | 7917.67 ± 706.9 | ops/s | 0.41x | behind |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 61172.4 ± 67.56 | 68065.4 ± 149.9 | ops/s | 0.90x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16851.1 ± 189.4 | 16899.9 ± 480.2 | ops/s | 1.00x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4146.61 ± 47.02 | 4233.62 ± 108.9 | ops/s | 0.98x | tie |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 16815.3 ± 43.48 | 4680.85 ± 9.944 | ops/s | 3.59x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 4215.01 ± 6.434 | 1244.86 ± 264.2 | ops/s | 3.39x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 1055.85 ± 2.836 | 360.387 ± 0.073 | ops/s | 2.93x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 64843.6 ± 39.42 | 79555.2 ± 2057 | ops/s | 0.82x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 4160.57 ± 10.65 | 5136.66 ± 3.815 | ops/s | 0.81x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1034.89 ± 20.95 | 1224.99 ± 20.83 | ops/s | 0.84x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 46624.3 ± 38.35 | 69469.5 ± 40.76 | ops/s | 0.67x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 2426.16 ± 1.978 | 4255.37 ± 35.65 | ops/s | 0.57x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 603.531 ± 0.6579 | 1057.84 ± 1.998 | ops/s | 0.57x | behind |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 780219 ± 3.138e+04 | 400022 ± 1.013e+04 | ops/s | 1.95x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15431.2 ± 94.22 | 6235.04 ± 26.15 | ops/s | 2.47x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 214.323 ± 3.714 | 114.154 ± 0.1501 | ops/s | 1.88x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.79062e+06 ± 4522 | 1.51137e+06 ± 1.937e+04 | ops/s | 1.18x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31919 ± 335.9 | 28566.7 ± 186.9 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 521.061 ± 5.445 | 423.025 ± 4.834 | ops/s | 1.23x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 442011 ± 3647 | 399723 ± 1855 | ops/s | 1.11x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7919.49 ± 53.38 | 6142.54 ± 32.47 | ops/s | 1.29x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 129.952 ± 1.585 | 105.801 ± 1.246 | ops/s | 1.23x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.21385e+06 ± 6075 | 942956 ± 9731 | ops/s | 1.29x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 30797.4 ± 259.9 | 31829.1 ± 203.3 | ops/s | 0.97x | behind |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 531.409 ± 3.553 | 478.071 ± 1.824 | ops/s | 1.11x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 294.623 ± 5.224 | 13141.5 ± 188.9 | ops/s | 0.02x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 7.18205 ± 0.03678 | 729.938 ± 6.489 | ops/s | 0.01x | behind |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 0.129336 ± 0.001272 | 1.31776 ± 0.01482 | ops/s | 0.10x | behind |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 204497 ± 3655 | 338219 ± 1.123e+04 | ops/s | 0.60x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 4190.63 ± 9.196 | 9845.78 ± 33.55 | ops/s | 0.43x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 52.6116 ± 0.03919 | 165.543 ± 0.07549 | ops/s | 0.32x | behind | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 31331 ± 1110 | 17700.7 ± 270.2 | ops/s | 1.77x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 763.975 ± 4.515 | 1160.31 ± 16.64 | ops/s | 0.66x | behind |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 3.07684 ± 0.01466 | 24.7979 ± 0.05152 | ops/s | 0.12x | behind |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 36825.6 ± 1207 | 18731.3 ± 226.6 | ops/s | 1.97x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 969.908 ± 11.45 | 1300.53 ± 10.15 | ops/s | 0.75x | behind |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 3.25883 ± 0.0174 | 28.3517 ± 0.05647 | ops/s | 0.11x | behind |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 32795.7 ± 261.4 | 38214.6 ± 228.3 | ops/s | 0.86x | behind |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 713.038 ± 15.12 | 1613.24 ± 31.09 | ops/s | 0.44x | behind |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 16404.8 ± 871.7 | 38122.2 ± 85.26 | ops/s | 0.43x | behind |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 1203.54 ± 29.86 | 1605.89 ± 13.79 | ops/s | 0.75x | behind |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 32945 ± 1448 | 38111.6 ± 48.37 | ops/s | 0.86x | behind |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 720.906 ± 1.246 | 1617.53 ± 26.48 | ops/s | 0.45x | behind |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.46138e+06 ± 2.086e+04 | 468008 ± 3180 | ops/s | 3.12x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 56271.8 ± 431.1 | 26534.5 ± 255 | ops/s | 2.12x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1081.85 ± 16.15 | 773.269 ± 7.361 | ops/s | 1.40x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.24242e+06 ± 5091 | 1.04515e+06 ± 4195 | ops/s | 1.19x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 30744.1 ± 213.5 | 32506.6 ± 254.9 | ops/s | 0.95x | behind |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 530.008 ± 3.624 | 479.829 ± 2.298 | ops/s | 1.10x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 456047 ± 8515 | 374935 ± 1523 | ops/s | 1.22x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 14368.5 ± 139.6 | 15573.7 ± 242.3 | ops/s | 0.92x | behind |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 333.5 ± 1.819 | 232.882 ± 1.532 | ops/s | 1.43x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 971198 ± 1.286e+04 | 802666 ± 1.713e+04 | ops/s | 1.21x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 27940.3 ± 213.2 | 30184.2 ± 150.1 | ops/s | 0.93x | behind |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 518.116 ± 3.705 | 472.851 ± 2.792 | ops/s | 1.10x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.14598 ± 0.009524 | 9.76776 ± 0.04243 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.0324 ± 0.2267 | 72.1979 ± 0.2101 | ms/op | 1.02x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 3206.78 ± 16.83 | 142.381 ± 0.3564 | ms/op | 0.04x | behind |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 28168 ± 532.6 | 1065.88 ± 1.292 | ms/op | 0.04x | behind |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 76.7572 ± 0.2862 | 93.9933 ± 0.08116 | ms/op | 1.22x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 534.404 ± 8.157 | 744.252 ± 0.7137 | ms/op | 1.39x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 15.0033 ± 0.02841 | 16.2744 ± 0.02307 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 119.718 ± 0.6893 | 127.867 ± 0.2155 | ms/op | 1.07x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 25.8661 ± 0.1816 | 36.5238 ± 0.0989 | ms/op | 1.41x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 192.729 ± 4.37 | 286.703 ± 7.74 | ms/op | 1.49x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 15.2001 ± 0.01653 | 16.3858 ± 0.02387 | ms/op | 1.08x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 120.672 ± 0.2147 | 128.217 ± 0.2038 | ms/op | 1.06x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 1110.69 ± 21.58 | 922.27 ± 16.91 | us/op | 0.83x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 2601.82 ± 43 | 1751.17 ± 38.78 | us/op | 0.67x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 103.497 ± 2.492 | 77.2897 ± 0.7576 | us/op | 0.75x | behind |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 2554.35 ± 39.19 | 1842.14 ± 8.54 | us/op | 0.72x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 82191.5 ± 752.8 | 107052 ± 410.4 | ops/s | 0.77x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2357.07 ± 20.15 | 2306.78 ± 8.923 | ops/s | 1.02x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 45.0552 ± 0.2597 | 37.8009 ± 0.05247 | ops/s | 1.19x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 277.191 ± 0.06219 | 103.276 ± 0.3328 | ops/s | 2.68x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 2337.52 ± 3.051 | 2847.98 ± 40.01 | ops/s | 0.82x | behind | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 2121.86 ± 174.4 | 190.193 ± 1.106 | ops/s | 11.16x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4193.32 ± 159.1 | 360.312 ± 0.166 | ops/s | 11.64x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 202.542 ± 0.3812 | 66.8804 ± 1.179 | ops/s | 3.03x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7490.6 ± 9.956 | 2164.86 ± 0.4588 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 4172.28 ± 9.856 | 2278.64 ± 15.65 | ops/s | 1.83x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7249.37 ± 2.064 | 395.309 ± 19.35 | ops/s | 18.34x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 219.55 ± 3.204 | 327.418 ± 8.424 | ops/s | 0.67x | behind |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 184.162 ± 1.127 | 247.522 ± 0.3609 | ops/s | 0.74x | behind |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.45097e+06 ± 4550 | 1.28942e+06 ± 1382 | ops/s | 2.68x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 54238.9 ± 15.96 | 32674.4 ± 5170 | ops/s | 1.66x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3416.05 ± 2.578 | 2121.8 ± 52.8 | ops/s | 1.61x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 337669 ± 2501 | 342892 ± 1178 | ops/s | 0.98x | behind |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5037.53 ± 36.19 | 5073.91 ± 103.2 | ops/s | 0.99x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 319.981 ± 0.6088 | 217.904 ± 1.388 | ops/s | 1.47x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 292552 ± 539.5 | 303251 ± 2236 | ops/s | 0.96x | behind |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4694.12 ± 3.635 | 4641.66 ± 15.22 | ops/s | 1.01x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 293.715 ± 0.07853 | 285.863 ± 0.288 | ops/s | 1.03x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 3.43712e+06 ± 1.64e+04 | 5.21681e+06 ± 1.642e+05 | ops/s | 0.66x | behind |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 54232.8 ± 22.17 | 89719.7 ± 1.308e+04 | ops/s | 0.60x | behind |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 3415.93 ± 3.071 | 5639.94 ± 54.12 | ops/s | 0.61x | behind |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.55478e+06 ± 6719 | 247752 ± 54.58 | ops/s | 30.49x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 119251 ± 34.48 | 3820.01 ± 0.3871 | ops/s | 31.22x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7408.21 ± 17.18 | 239.223 ± 0.07122 | ops/s | 30.97x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.39965e+06 ± 1336 | 108464 ± 83.47 | ops/s | 40.56x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 69157 ± 8.922 | 1686.28 ± 3.583 | ops/s | 41.01x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4308.67 ± 5.656 | 105.202 ± 0.01813 | ops/s | 40.96x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.89678e+06 ± 7544 | 2.34043e+06 ± 435.9 | ops/s | 1.24x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55931.8 ± 199.9 | 34434.7 ± 36.58 | ops/s | 1.62x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3529.42 ± 3.113 | 2149.73 ± 0.7002 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 259896 ± 4168 | 292440 ± 725.5 | ops/s | 0.89x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 2338.63 ± 31.42 | 4380 ± 73.38 | ops/s | 0.53x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 132.326 ± 2.651 | 186.901 ± 1.554 | ops/s | 0.71x | behind |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 217944 ± 4454 | 147860 ± 374.1 | ops/s | 1.47x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3310.12 ± 30.91 | 2189.32 ± 257 | ops/s | 1.51x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 208.173 ± 0.7369 | 118.065 ± 0.3177 | ops/s | 1.76x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.73181e+06 ± 7009 | 2.3639e+06 ± 2029 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 123509 ± 944 | 34675 ± 28.67 | ops/s | 3.56x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7321.01 ± 24.24 | 2162.76 ± 0.5376 | ops/s | 3.39x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.75475e+07 ± 1.384e+05 | 7.86548e+06 ± 2.487e+05 | ops/s | 3.50x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.91987e+07 ± 1.008e+05 | 6.63804e+06 ± 7.031e+04 | ops/s | 2.89x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.16187e+07 ± 6.535e+04 | 2.54046e+07 ± 6.369e+04 | ops/s | 1.24x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.90641e+07 ± 2.258e+05 | 1.83476e+07 ± 2.757e+05 | ops/s | 2.13x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 1.0458e+07 ± 1.05e+05 | 3.40029e+06 ± 3.682e+04 | ops/s | 3.08x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 6.13357e+06 ± 1.719e+05 | 2.81023e+06 ± 2.436e+04 | ops/s | 2.18x | ahead | gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale) |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.57029e+07 ± 6788 | 5.55665e+06 ± 1.623e+04 | ops/s | 2.83x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.12561e+07 ± 7.603e+04 | 4.70698e+06 ± 1.919e+05 | ops/s | 2.39x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.228796 ± 0.000749 | 1.36699 ± 0.05159 | ms/op | 5.97x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.3304 ± 0.03978 | 132.901 ± 0.1712 | ms/op | 5.46x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.84022 ± 0.008242 | 11.9832 ± 0.02227 | ms/op | 6.51x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 220.56 ± 1.762 | 1271.95 ± 82.25 | ms/op | 5.77x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00574123 ± 9.235e-06 | 0.0401443 ± 7.338e-05 | ms/op | 6.99x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.590829 ± 0.002171 | 3.75186 ± 0.008045 | ms/op | 6.35x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0459107 ± 0.0001763 | 0.367259 ± 0.001218 | ms/op | 8.00x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.16235 ± 0.01825 | 37.2654 ± 0.03826 | ms/op | 6.05x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.139552 ± 0.0002959 | 1.28984 ± 0.008326 | ms/op | 9.24x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.6362 ± 0.4427 | 133.017 ± 0.2878 | ms/op | 9.09x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.18894 ± 0.00595 | 11.9852 ± 0.01828 | ms/op | 10.08x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 144.137 ± 4 | 1268.64 ± 70.72 | ms/op | 8.80x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00493563 ± 1.966e-05 | 0.0409616 ± 0.0006312 | ms/op | 8.30x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.541213 ± 0.0006943 | 3.7499 ± 0.004948 | ms/op | 6.93x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.052688 ± 0.00021 | 0.367407 ± 0.001839 | ms/op | 6.97x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.61125 ± 0.007658 | 37.2351 ± 0.06725 | ms/op | 6.64x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 178378 ± 759.2 | 188070 ± 1206 | ops/s | 0.95x | behind |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 16639.7 ± 278.3 | 16305.8 ± 264.9 | ops/s | 1.02x | tie |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.095e+06 ± 8872 | 317467 ± 770.8 | ops/s | 6.60x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 210444 ± 1475 | 32850.8 ± 1329 | ops/s | 6.41x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.14389e+06 ± 6.385e+04 | 327314 ± 568.8 | ops/s | 3.49x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 41157.9 ± 2.44e+04 | 35452.7 ± 53.27 | ops/s | 1.16x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.2841e+06 ± 5911 | 2.27634e+06 ± 3874 | ops/s | 1.00x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 225771 ± 300.2 | 225737 ± 408.8 | ops/s | 1.00x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 70920.9 ± 1802 | 34355 ± 559.6 | ops/s | 2.06x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 1713.55 ± 61.52 | 1764.96 ± 22.33 | ops/s | 0.97x | tie |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 171.425 ± 0.6718 | 315.342 ± 3.177 | ops/s | 0.54x | behind |  |

## Provenance

- Date: 2026-10-09, Slurm job 2524300, start 11:52:00, end 13:50:44 EDT (elapsed 1h58m44s)
- Benchmarked source: `17813d0d761467187ba0764f1081fd30747ad18e` (breeze/integration). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster.
- Node: SciNet trillium `tri0819`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5, family 26 model 2), AVX-512, SMT off, 192 CPUs, NPS4 (8 NUMA nodes of 24 cores), boost enabled. Load average at job start: 1.39, 42.62, 87.08.
- JDK: Temurin 25+36 (`module load java/25`, java/25.36). The site default `JAVA_TOOL_OPTIONS=-Xmx2g` was unset before launch, so no fork inherited it.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the L3/CCD domain of cpu0 plus its NUMA node). It applies to the JMH host JVM and, by inheritance, to every fork.
- Netlib (sidecar `2026-10-09-breeze-two-lane-baseline-trillium-laneA.netlib.jsonl`, (blas, lapack, vectorModule) over all records): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False)]. No native BLAS/LAPACK library is visible to the loader on these nodes (`ldconfig -p` lists none; `LD_LIBRARY_PATH=/opt/slurm/lib64`), so nothing had to be forced.
- Command, equivalent to the `breezeLaneA` alias: `numactl … java --add-modules=jdk.incubator.vector -cp <classpath> org.openjdk.jmh.Main -jvmArgs "-Xms4g -Xmx4g -Dgale.bench.lane=A -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings come from the annotations. As under sbt, the host JVM carries `--add-modules`: lane A's `-jvmArgs` replaces it in the forks, and lane B's forks inherit it. The runner script is `docs/verification/w21-simd-spike/x86-avx512/run.sh`.
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane A --strict --commit 17813d0d… --netlib 2026-10-09-breeze-two-lane-baseline-trillium-laneA.netlib.jsonl 2026-10-09-breeze-two-lane-baseline-trillium-laneA.json`. It exited 0, with 0 unpaired rows.
