# gale vs Breeze scoreboard — lane B (SIMD)

- Lane: B
- JDK: 25 (OpenJDK 64-Bit Server VM)
- Fork JVM args: --add-modules=jdk.incubator.vector -Dgale.bench.sourceSha=0ff64405ea981f5a6dbf026a12e2bcc81617f412 -Xms4g -Xmx4g -Dgale.bench.sourceSha=0ff64405ea981f5a6dbf026a12e2bcc81617f412 -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=/scratch/brad/gale-bench/m2/out/laneB-2535352/laneB/netlib.jsonl
- Breeze netlib BLAS: `dev.ludovic.netlib.blas.VectorBLAS`
- Breeze netlib LAPACK: `dev.ludovic.netlib.lapack.F2jLAPACK`
- Commit: `0ff64405ea981f5a6dbf026a12e2bcc81617f412`
- Numerical preflight: 151 decomposition/factorization result rows validated on this SHA.
- Machine: SciNet trillium tri0468: 2x AMD EPYC 9655 (Zen 5, 96c/socket, AVX-512, 8 double lanes), NPS4; numactl --physcpubind=0-7 --membind=0 (one CCD); Temurin 25+36; Slurm job 2535352 exclusive

**159 ahead, 16 tie, 23 behind, 0 n/a, 11 withheld** of 209 pairs; 0 unpaired. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name a known asymmetry behind a verdict.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 179935 ± 3935 | 172602 ± 3156 | ops/s | 1.04x | ahead |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 28107.5 ± 2353 | 30019.2 ± 2368 | ops/s | 0.94x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7307.55 ± 867.5 | 7599.51 ± 459.3 | ops/s | 0.96x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 64987.6 ± 1133 | 67846.2 ± 51.63 | ops/s | 0.96x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 15673.2 ± 1804 | 17108 ± 18.51 | ops/s | 0.92x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4170.03 ± 80.54 | 4226.48 ± 56.92 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 67156.8 ± 1273 | 4649.75 ± 16.56 | ops/s | 14.44x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 16940.8 ± 156.9 | 1237.49 ± 280.1 | ops/s | 13.69x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 4218.65 ± 29.17 | 358.669 ± 0.1829 | ops/s | 11.76x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 99056.5 ± 666.9 | 234004 ± 517.7 | ops/s | 0.42x | behind |  |
| BlasL2BreezeJmh | Gemv | n=256 | vector | thrpt | 246472 ± 3.375e+04 | 234004 ± 517.7 | ops/s | 1.05x | tie |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6514.56 ± 15.18 | 12793 ± 419.4 | ops/s | 0.51x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | thrpt | 12682.5 ± 6.516 | 12793 ± 419.4 | ops/s | 0.99x | tie |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1600.21 ± 24.34 | 2803.18 ± 58.28 | ops/s | 0.57x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | thrpt | 2707.43 ± 38.32 | 2803.18 ± 58.28 | ops/s | 0.97x | tie |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 83123.8 ± 442.7 | 246191 ± 4334 | ops/s | 0.34x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | vector | thrpt | 247548 ± 3.172e+04 | 246191 ± 4334 | ops/s | 1.01x | tie |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5782.33 ± 14.4 | 13078.5 ± 96.57 | ops/s | 0.44x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | thrpt | 13054.5 ± 265.4 | 13078.5 ± 96.57 | ops/s | 1.00x | tie |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1427.08 ± 4.476 | 2910.15 ± 65.73 | ops/s | 0.49x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | thrpt | 2735.94 ± 40.79 | 2910.15 ± 65.73 | ops/s | 0.94x | behind |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 789366 ± 1590 | 397136 ± 8124 | ops/s | 1.99x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | vector | thrpt | 786927 ± 4881 | 397136 ± 8124 | ops/s | 1.98x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15363.8 ± 93.16 | 14418.7 ± 41.11 | ops/s | 1.07x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | vector | thrpt | 15284.9 ± 87.12 | 14418.7 ± 41.11 | ops/s | 1.06x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 213.961 ± 3.801 | 352.759 ± 13.78 | ops/s | 0.61x | behind |  |
| BlasL3BreezeJmh | AtA | n=256 | vector | thrpt | 214.394 ± 3.91 | 352.759 ± 13.78 | ops/s | 0.61x | behind |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.79621e+06 ± 4338 | 1.53189e+06 ± 1.619e+04 | ops/s | 1.17x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | vector | thrpt | 1.79344e+06 ± 3609 | 1.53189e+06 ± 1.619e+04 | ops/s | 1.17x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31748.9 ± 320.6 | 28395.9 ± 154.8 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | vector | thrpt | 31698.1 ± 311.3 | 28395.9 ± 154.8 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 518.577 ± 4.445 | 1140.94 ± 19.46 | ops/s | 0.45x | behind |  |
| BlasL3BreezeJmh | Gemm | n=256 | vector | thrpt | 1551.52 ± 12.26 | 1140.94 ± 19.46 | ops/s | 1.36x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 436381 ± 7488 | 397731 ± 1416 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | thrpt | 437977 ± 1865 | 397731 ± 1416 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7890.3 ± 69.82 | 12665.2 ± 21.78 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | thrpt | 7887.29 ± 64.61 | 12665.2 ± 21.78 | ops/s | 0.62x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 129.154 ± 1.523 | 290.6 ± 8.356 | ops/s | 0.44x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | thrpt | 461.454 ± 3.93 | 290.6 ± 8.356 | ops/s | 1.59x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.38046e+06 ± 1.042e+04 | 946220 ± 2634 | ops/s | 1.46x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | vector (gemm-routed only) | thrpt | 1.38336e+06 ± 2.228e+04 | 946220 ± 2634 | ops/s | 1.46x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 47821.1 ± 128 | 31607.1 ± 205.1 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | vector (gemm-routed only) | thrpt | 47681.6 ± 204.9 | 31607.1 ± 205.1 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 918.728 ± 12.76 | 477.214 ± 2.298 | ops/s | 1.93x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | vector (gemm-routed only) | thrpt | 927.638 ± 11.15 | 477.214 ± 2.298 | ops/s | 1.94x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33310.2 ± 331.8 | 12934.3 ± 92.7 | ops/s | 2.58x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1086.23 ± 0.7227 | 722.168 ± 5.883 | ops/s | 1.50x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 18.9689 ± 0.02302 | 1.32243 ± 0.0007876 | ops/s | 14.34x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 474167 ± 5287 | 336569 ± 1.084e+04 | ops/s | 1.41x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector (gemm-routed only) | thrpt | 455768 ± 1.207e+04 | 336569 ± 1.084e+04 | ops/s | 1.35x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 17496.1 ± 35.02 | 9825.12 ± 35.34 | ops/s | 1.78x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector (gemm-routed only) | thrpt | 17497.5 ± 35.57 | 9825.12 ± 35.34 | ops/s | 1.78x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 363.079 ± 4.116 | 164.549 ± 0.4792 | ops/s | 2.21x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector (gemm-routed only) | thrpt | 374.895 ± 14.93 | 164.549 ± 0.4792 | ops/s | 2.28x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 45499.9 ± 962.3 | 15449.5 ± 38.91 | ops/s | 2.95x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 1373.24 ± 7.581 | 1077.42 ± 12.88 | ops/s | 1.27x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 27.7404 ± 0.04685 | 26.1644 ± 0.07939 | ops/s | 1.06x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 51508.4 ± 830 | 13911.9 ± 6934 | ops/s | 3.70x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 1688.23 ± 3.552 | 1296.65 ± 11.07 | ops/s | 1.30x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 33.0667 ± 0.04673 | 28.2601 ± 0.06506 | ops/s | 1.17x | ahead |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 44143.6 ± 405.9 | 37572.9 ± 785 | ops/s | 1.17x | ahead |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 2080.22 ± 32.44 | 1602.96 ± 22.48 | ops/s | 1.30x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 44257.5 ± 1039 | 38256.4 ± 280.2 | ops/s | 1.16x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 2086.21 ± 16.67 | 1607.72 ± 9.339 | ops/s | 1.30x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 44734.4 ± 788.5 | 37555.9 ± 93.44 | ops/s | 1.19x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 2068.79 ± 27.99 | 1618.23 ± 29.67 | ops/s | 1.28x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.46436e+06 ± 4493 | 374467 ± 3770 | ops/s | 3.91x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | vector (gemm-routed only) | thrpt | 1.46356e+06 ± 4667 | 374467 ± 3770 | ops/s | 3.91x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 56218.9 ± 137.4 | 20227.3 ± 81.47 | ops/s | 2.78x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | vector (gemm-routed only) | thrpt | 56093.2 ± 172 | 20227.3 ± 81.47 | ops/s | 2.77x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1074.13 ± 13.76 | 660.286 ± 21.35 | ops/s | 1.63x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | vector (gemm-routed only) | thrpt | 1076.18 ± 13.93 | 660.286 ± 21.35 | ops/s | 1.63x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.40501e+06 ± 1.455e+04 | 1.03703e+06 ± 1.064e+04 | ops/s | 1.35x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | vector (gemm-routed only) | thrpt | 1.41131e+06 ± 4419 | 1.03703e+06 ± 1.064e+04 | ops/s | 1.36x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 47586.3 ± 176.6 | 32376.2 ± 288 | ops/s | 1.47x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | vector (gemm-routed only) | thrpt | 47646.1 ± 138.7 | 32376.2 ± 288 | ops/s | 1.47x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 932.032 ± 9.298 | 477.814 ± 2.121 | ops/s | 1.95x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | vector (gemm-routed only) | thrpt | 931.356 ± 10.25 | 477.814 ± 2.121 | ops/s | 1.95x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 524076 ± 3039 | 340971 ± 3.703e+04 | ops/s | 1.54x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | vector (gemm-routed only) | thrpt | 516196 ± 5210 | 340971 ± 3.703e+04 | ops/s | 1.51x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 16930.3 ± 179.3 | 15510.6 ± 49.44 | ops/s | 1.09x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | vector (gemm-routed only) | thrpt | 16907.5 ± 209.7 | 15510.6 ± 49.44 | ops/s | 1.09x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 342.129 ± 1.733 | 231.655 ± 1.526 | ops/s | 1.48x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | vector (gemm-routed only) | thrpt | 345.401 ± 1.746 | 231.655 ± 1.526 | ops/s | 1.49x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 1.07047e+06 ± 2.246e+04 | 780741 ± 6029 | ops/s | 1.37x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | vector (gemm-routed only) | thrpt | 1.07206e+06 ± 1.033e+04 | 780741 ± 6029 | ops/s | 1.37x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 41208.5 ± 257.4 | 30013.7 ± 205.4 | ops/s | 1.37x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | vector (gemm-routed only) | thrpt | 41439 ± 203.7 | 30013.7 ± 205.4 | ops/s | 1.38x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 791.31 ± 5.103 | 471.18 ± 2.433 | ops/s | 1.68x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | vector (gemm-routed only) | thrpt | 834.386 ± 74.42 | 471.18 ± 2.433 | ops/s | 1.77x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.17497 ± 0.009719 | 10.454 ± 0.1303 | ms/op | 1.28x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | vector (gemm-routed only) | avgt | 8.18475 ± 0.006638 | 10.454 ± 0.1303 | ms/op | 1.28x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.3753 ± 0.1638 | 73.9419 ± 0.1005 | ms/op | 1.04x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | vector (gemm-routed only) | avgt | 71.357 ± 0.1834 | 73.9419 ± 0.1005 | ms/op | 1.04x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 103.158 ± 2.061 | 143.36 ± 0.1654 | ms/op | 1.39x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 775.168 ± 3.168 | 1072.08 ± 1.515 | ms/op | 1.38x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 76.8717 ± 0.4553 | 94.4375 ± 0.04944 | ms/op | 1.23x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | vector (gemm-routed only) | avgt | 69.4408 ± 2.672 | 94.4375 ± 0.04944 | ms/op | 1.36x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 542.573 ± 6.922 | 748.911 ± 0.6908 | ms/op | 1.38x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | vector (gemm-routed only) | avgt | 395.626 ± 17.63 | 748.911 ± 0.6908 | ms/op | 1.89x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 8.81747 ± 0.01152 | 16.35 ± 0.01698 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | vector (gemm-routed only) | avgt | 8.81735 ± 0.0218 | 16.35 ± 0.01698 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 69.3747 ± 0.1725 | 128.531 ± 0.1624 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | vector (gemm-routed only) | avgt | 69.3587 ± 0.211 | 128.531 ± 0.1624 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 25.9444 ± 0.09873 | 36.0651 ± 1.011 | ms/op | 1.39x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | vector (gemm-routed only) | avgt | 25.4454 ± 0.1778 | 36.0651 ± 1.011 | ms/op | 1.42x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 194.999 ± 2.744 | 283.013 ± 0.2674 | ms/op | 1.45x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | vector (gemm-routed only) | avgt | 135.088 ± 5.915 | 283.013 ± 0.2674 | ms/op | 2.10x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 9.0894 ± 0.03485 | 16.4597 ± 0.0212 | ms/op | 1.81x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | vector (gemm-routed only) | avgt | 9.09449 ± 0.03835 | 16.4597 ± 0.0212 | ms/op | 1.81x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 70.3673 ± 0.2933 | 128.863 ± 0.273 | ms/op | 1.83x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | vector (gemm-routed only) | avgt | 70.4542 ± 0.2723 | 128.863 ± 0.273 | ms/op | 1.83x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 839.47 ± 5.389 | 556.44 ± 5.123 | us/op | 0.66x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | avgt | 616.763 ± 9.145 | 556.44 ± 5.123 | us/op | 0.90x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1957.32 ± 13.44 | 1044.49 ± 14.4 | us/op | 0.53x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | avgt | 1431.01 ± 14.54 | 1044.49 ± 14.4 | us/op | 0.73x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 103.322 ± 0.379 | 77.9367 ± 0.7797 | us/op | 0.75x | behind |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 2496.08 ± 105.7 | 2045.69 ± 12.69 | us/op | 0.82x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 92762.3 ± 278.4 | 106438 ± 530.7 | ops/s | 0.87x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector (gemm-routed only) | thrpt | 93359.6 ± 374.1 | 106438 ± 530.7 | ops/s | 0.88x | behind |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2819.69 ± 14.11 | 2298.29 ± 8.654 | ops/s | 1.23x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector (gemm-routed only) | thrpt | 2824.13 ± 26.33 | 2298.29 ± 8.654 | ops/s | 1.23x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 45.9113 ± 0.2849 | 37.643 ± 0.03084 | ops/s | 1.22x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector (gemm-routed only) | thrpt | 48.2765 ± 1.189 | 37.643 ± 0.03084 | ops/s | 1.28x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 292.247 ± 0.797 | 103.277 ± 0.729 | ops/s | 2.83x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 3665.92 ± 63.07 | 2842.36 ± 42.89 | ops/s | 1.29x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 4641.57 ± 22.19 | 189.193 ± 0.718 | ops/s | 24.53x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4237.89 ± 71.46 | 358.72 ± 0.1689 | ops/s | 11.81x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 254.594 ± 0.6031 | 66.7532 ± 0.6095 | ops/s | 3.81x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7442.32 ± 17.55 | 2153.96 ± 1.726 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 3603.22 ± 9.61 | 2266.57 ± 14.72 | ops/s | 1.59x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7215.1 ± 3.548 | 390.817 ± 11.35 | ops/s | 18.46x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 659.741 ± 8.136 | 320.941 ± 0.6443 | ops/s | 2.06x | ahead |  |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 660.141 ± 9.178 | 320.941 ± 0.6443 | ops/s | 2.06x | ahead |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 620.689 ± 4.285 | 241.943 ± 0.1827 | ops/s | 2.57x | ahead |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | vector (gemm-routed only) | thrpt | 623.18 ± 3.639 | 241.943 ± 0.1827 | ops/s | 2.58x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.42949e+06 ± 5510 | 1.27984e+06 ± 4309 | ops/s | 2.68x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 53966.6 ± 21.61 | 32580.4 ± 4877 | ops/s | 1.66x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3395.41 ± 21.22 | 2098.86 ± 134.3 | ops/s | 1.62x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 342158 ± 1278 | 342324 ± 1153 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 4987.43 ± 17.07 | 5050.51 ± 89.8 | ops/s | 0.99x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 324.744 ± 0.7442 | 214.919 ± 3.782 | ops/s | 1.51x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 299268 ± 3955 | 302604 ± 921.7 | ops/s | 0.99x | tie |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4839.3 ± 4.149 | 4598.91 ± 45.39 | ops/s | 1.05x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 303.432 ± 1.271 | 283.24 ± 0.5983 | ops/s | 1.07x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 5.37722e+06 ± 6.587e+05 | 5.05622e+06 ± 3.177e+05 | ops/s | 1.06x | tie |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 78451 ± 1.68e+04 | 93640.6 ± 4662 | ops/s | 0.84x | tie |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4240.51 ± 1.766 | 5607.96 ± 308.4 | ops/s | 0.76x | behind |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.52062e+06 ± 3032 | 246725 ± 54.79 | ops/s | 30.48x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118675 ± 81.31 | 3804.15 ± 2.62 | ops/s | 31.20x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7373.16 ± 23.6 | 237.983 ± 0.06531 | ops/s | 30.98x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.37509e+06 ± 4796 | 107990 ± 84.18 | ops/s | 40.51x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68770.9 ± 70.41 | 1679.14 ± 3.591 | ops/s | 40.96x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4285.99 ± 0.8028 | 104.715 ± 0.08087 | ops/s | 40.93x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.88006e+06 ± 8174 | 2.329e+06 ± 6493 | ops/s | 1.24x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55776.8 ± 9.831 | 34264.6 ± 31.26 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3511.33 ± 1.598 | 2139.44 ± 1.091 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 248078 ± 1383 | 291168 ± 863.3 | ops/s | 0.85x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 2324.85 ± 20.24 | 4294.35 ± 16.27 | ops/s | 0.54x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 133.029 ± 2.657 | 185.144 ± 2.892 | ops/s | 0.72x | behind |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 271439 ± 1705 | 148324 ± 320.7 | ops/s | 1.83x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3991.67 ± 25.18 | 2186.01 ± 261.9 | ops/s | 1.83x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 256.87 ± 0.8395 | 116.83 ± 1.044 | ops/s | 2.20x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.69038e+06 ± 1.901e+04 | 2.35199e+06 ± 543.3 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 123182 ± 262.9 | 34485.5 ± 34.4 | ops/s | 3.57x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7279.98 ± 15.11 | 2152.1 ± 5.63 | ops/s | 3.38x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.81129e+07 ± 2.959e+05 | 7.98636e+06 ± 3.701e+04 | ops/s | 3.52x | ahead |  |
| SmallDenseBreezeJmh | Det | n=3 | vector (gemm-routed only) | thrpt | 2.81792e+07 ± 4.717e+05 | 7.98636e+06 ± 3.701e+04 | ops/s | 3.53x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.84721e+07 ± 4.776e+04 | 6.74156e+06 ± 1.51e+05 | ops/s | 2.74x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | vector (gemm-routed only) | thrpt | 1.83463e+07 ± 2.079e+05 | 6.74156e+06 ± 1.51e+05 | ops/s | 2.72x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.14487e+07 ± 7.546e+04 | 2.57817e+07 ± 4.444e+04 | ops/s | 1.22x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | vector | thrpt | 3.12014e+07 ± 5.551e+04 | 2.57817e+07 ± 4.444e+04 | ops/s | 1.21x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.89321e+07 ± 6.076e+04 | 1.81825e+07 ± 3.98e+05 | ops/s | 2.14x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | vector | thrpt | 3.87536e+07 ± 1.145e+05 | 1.81825e+07 ± 3.98e+05 | ops/s | 2.13x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 9.88078e+06 ± 1.726e+04 | 3.17433e+06 ± 3.967e+05 | ops/s | 3.11x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | vector (gemm-routed only) | thrpt | 9.9005e+06 ± 5.037e+04 | 3.17433e+06 ± 3.967e+05 | ops/s | 3.12x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 8.47041e+06 ± 4.8e+04 | 2.76484e+06 ± 9.662e+04 | ops/s | 3.06x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=4 | vector (gemm-routed only) | thrpt | 8.50055e+06 ± 9.503e+04 | 2.76484e+06 ± 9.662e+04 | ops/s | 3.07x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.60706e+07 ± 8.985e+04 | 5.53703e+06 ± 1.254e+05 | ops/s | 2.90x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | vector (gemm-routed only) | thrpt | 1.59526e+07 ± 2.247e+05 | 5.53703e+06 ± 1.254e+05 | ops/s | 2.88x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.10555e+07 ± 4.301e+05 | 4.87142e+06 ± 2.317e+04 | ops/s | 2.27x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | vector (gemm-routed only) | thrpt | 1.09034e+07 ± 1.878e+05 | 4.87142e+06 ± 2.317e+04 | ops/s | 2.24x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.229743 ± 0.0004469 | 1.27436 ± 0.003285 | ms/op | 5.55x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.3637 ± 0.1538 | 122.371 ± 0.277 | ms/op | 5.02x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.84442 ± 0.005116 | 12.2348 ± 0.0243 | ms/op | 6.63x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 220.56 ± 0.5401 | 1229.12 ± 52.47 | ms/op | 5.57x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00576957 ± 4.068e-06 | 0.0404723 ± 0.0002258 | ms/op | 7.01x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.590097 ± 0.001189 | 3.79805 ± 0.009682 | ms/op | 6.44x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0463703 ± 0.0001576 | 0.370923 ± 0.0008543 | ms/op | 8.00x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.16704 ± 0.008465 | 37.7115 ± 0.05668 | ms/op | 6.11x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.14024 ± 0.00034 | 1.27376 ± 0.004152 | ms/op | 9.08x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.7611 ± 0.3825 | 122.153 ± 0.3108 | ms/op | 8.28x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19185 ± 0.005725 | 12.2092 ± 0.02215 | ms/op | 10.24x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 144.639 ± 2.938 | 1224.35 ± 59.24 | ms/op | 8.46x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00495714 ± 3.344e-05 | 0.0402447 ± 4.063e-05 | ms/op | 8.12x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.542477 ± 0.002594 | 3.80444 ± 0.0131 | ms/op | 7.01x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0529522 ± 0.0001558 | 0.37168 ± 0.001576 | ms/op | 7.02x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.63962 ± 0.01391 | 38.6747 ± 1.649 | ms/op | 6.86x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 210572 ± 789.1 | 198026 ± 2188 | ops/s | 1.06x | ahead |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 17765.2 ± 45.79 | 15967.4 ± 356.3 | ops/s | 1.11x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.05285e+06 ± 4.917e+04 | 315640 ± 849.2 | ops/s | 6.50x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 209527 ± 4656 | 32695.3 ± 1290 | ops/s | 6.41x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.13594e+06 ± 5.738e+04 | 325511 ± 362.1 | ops/s | 3.49x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 41223.2 ± 2.506e+04 | 35281.5 ± 18.67 | ops/s | 1.17x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.26899e+06 ± 1997 | 2.25607e+06 ± 6874 | ops/s | 1.01x | ahead |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 222708 ± 2824 | 222978 ± 2747 | ops/s | 1.00x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97829 ± 1028 | 34052 ± 741 | ops/s | 2.87x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 3714.49 ± 10.93 | 1769.73 ± 11.43 | ops/s | 2.10x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 574.748 ± 2.618 | 315.133 ± 3.575 | ops/s | 1.82x | ahead |  |

## Delta against Milestone 1

Milestone 1 is `2026-10-09-breeze-milestone1-trillium-laneB.md` (source ff03eed). It used the same lane, CPU model, JDK, pinning and runner protocol, on a different node of the same type. Milestone 1 was a legacy receipt with no numerical preflight. The `Pinv` fixture changed between the runs (seed 2200 to seed 2100; see `benchmarks/README.md`), so its rows carry no change ratio.

Ratios are gale speed / Breeze speed (>1 means gale is faster). `change` is new ratio / old ratio. Old is milestone 1 ff03eed; new is milestone 2 0ff6440.

| class | op | params | gale backend | old ratio | new ratio | change | old verdict | new verdict |
|---|---|---|---|---:|---:|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | 0.92x | 0.96x | 1.04x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | 0.98x | 0.94x | 0.96x | tie | tie |
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | 0.98x | 1.04x | 1.06x | tie | ahead |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | 0.98x | 0.99x | 1.01x | behind | tie |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | 0.98x | 0.92x | 0.94x | behind | tie |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | 0.95x | 0.96x | 1.01x | behind | behind |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | 11.81x | 11.76x | 1.00x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | 13.62x | 13.69x | 1.01x | ahead | ahead |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | 14.21x | 14.44x | 1.02x | ahead | ahead |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | 0.51x | 0.51x | 1.00x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | 0.99x | 0.99x | 1.00x | tie | tie |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | 0.58x | 0.57x | 0.98x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | 0.97x | 0.97x | 1.00x | behind | tie |
| BlasL2BreezeJmh | Gemv | n=256 | pure | 0.42x | 0.42x | 1.00x | behind | behind |
| BlasL2BreezeJmh | Gemv | n=256 | vector | 0.96x | 1.05x | 1.09x | behind | tie |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | 0.43x | 0.44x | 1.02x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | 0.96x | 1.00x | 1.04x | behind | tie |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | 0.48x | 0.49x | 1.02x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | 0.93x | 0.94x | 1.01x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=256 | pure | 0.34x | 0.34x | 1.00x | behind | behind |
| BlasL2BreezeJmh | GemvT | n=256 | vector | 0.92x | 1.01x | 1.10x | behind | tie |
| BlasL3BreezeJmh | AtA | n=16 | pure | 2.00x | 1.99x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=16 | vector | 1.99x | 1.98x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=256 | pure | 0.60x | 0.61x | 1.02x | behind | behind |
| BlasL3BreezeJmh | AtA | n=256 | vector | 0.62x | 0.61x | 0.98x | behind | behind |
| BlasL3BreezeJmh | AtA | n=64 | pure | 1.06x | 1.07x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | AtA | n=64 | vector | 1.06x | 1.06x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | pure | 1.17x | 1.17x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=16 | vector | 1.17x | 1.17x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=256 | pure | 0.45x | 0.45x | 1.00x | behind | behind |
| BlasL3BreezeJmh | Gemm | n=256 | vector | 1.37x | 1.36x | 0.99x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | pure | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | Gemm | n=64 | vector | 1.12x | 1.12x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | 1.10x | 1.10x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=16 | vector | 1.10x | 1.10x | 1.00x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | 0.44x | 0.44x | 1.00x | behind | behind |
| BlasL3BreezeJmh | GemmTall | n=256 | vector | 1.57x | 1.59x | 1.01x | ahead | ahead |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | 0.62x | 0.62x | 1.00x | behind | behind |
| BlasL3BreezeJmh | GemmTall | n=64 | vector | 0.62x | 0.62x | 1.00x | behind | behind |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | 1.29x | 1.46x | 1.13x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=16 | vector (gemm-routed only) | 1.29x | 1.46x | 1.13x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | 1.11x | 1.93x | 1.74x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=256 | vector (gemm-routed only) | 1.11x | 1.94x | 1.75x | ahead | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | 0.97x | 1.51x | 1.56x | behind | ahead |
| DenseDecompositionBreezeJmh | Det | n=64 | vector (gemm-routed only) | 0.96x | 1.51x | 1.57x | behind | ahead |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | 0.02x | 2.58x | 129.00x | behind | ahead |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | 0.10x | 14.34x | 143.40x | behind | ahead |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | 0.01x | 1.50x | 150.00x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | 0.60x | 1.41x | 2.35x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=16 | vector (gemm-routed only) | 0.60x | 1.35x | 2.25x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | 0.32x | 2.21x | 6.91x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=256 | vector (gemm-routed only) | 0.32x | 2.28x | 7.12x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | 0.42x | 1.78x | 4.24x | behind | ahead |
| DenseDecompositionBreezeJmh | Inv | n=64 | vector (gemm-routed only) | 0.43x | 1.78x | 4.14x | behind | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | 2.28x | 2.95x | fixture changed | ahead | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | 0.93x | 1.06x | fixture changed | behind | ahead |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | 1.07x | 1.27x | fixture changed | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | 2.34x | 3.70x | 1.58x | ahead | ahead |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | 0.95x | 1.17x | 1.23x | behind | ahead |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | 1.08x | 1.30x | 1.20x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | 1.28x | 1.30x | 1.02x | ahead | ahead |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | 1.14x | 1.17x | 1.03x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | 1.29x | 1.30x | 1.01x | ahead | ahead |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | 1.17x | 1.16x | 0.99x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | 1.28x | 1.28x | 1.00x | ahead | ahead |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | 1.16x | 1.19x | 1.03x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=16 | pure | 3.13x | 3.91x | 1.25x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=16 | vector (gemm-routed only) | 3.12x | 3.91x | 1.25x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | pure | 1.40x | 1.63x | 1.16x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=256 | vector (gemm-routed only) | 1.39x | 1.63x | 1.17x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | pure | 2.13x | 2.78x | 1.31x | ahead | ahead |
| FactorizationBreezeJmh | Chol | n=64 | vector (gemm-routed only) | 2.14x | 2.77x | 1.29x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | pure | 1.19x | 1.35x | 1.13x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=16 | vector (gemm-routed only) | 1.18x | 1.36x | 1.15x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | pure | 1.10x | 1.95x | 1.77x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=256 | vector (gemm-routed only) | 1.10x | 1.95x | 1.77x | ahead | ahead |
| FactorizationBreezeJmh | Lu | n=64 | pure | 0.94x | 1.47x | 1.56x | behind | ahead |
| FactorizationBreezeJmh | Lu | n=64 | vector (gemm-routed only) | 0.95x | 1.47x | 1.55x | behind | ahead |
| FactorizationBreezeJmh | Qr | n=16 | pure | 1.30x | 1.54x | 1.18x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=16 | vector (gemm-routed only) | 1.29x | 1.51x | 1.17x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | pure | 1.43x | 1.48x | 1.03x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=256 | vector (gemm-routed only) | 1.45x | 1.49x | 1.03x | ahead | ahead |
| FactorizationBreezeJmh | Qr | n=64 | pure | 0.92x | 1.09x | 1.18x | behind | ahead |
| FactorizationBreezeJmh | Qr | n=64 | vector (gemm-routed only) | 0.92x | 1.09x | 1.18x | behind | ahead |
| FactorizationBreezeJmh | Solve | n=16 | pure | 1.21x | 1.37x | 1.13x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=16 | vector (gemm-routed only) | 1.19x | 1.37x | 1.15x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | pure | 1.10x | 1.68x | 1.53x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=256 | vector (gemm-routed only) | 1.09x | 1.77x | 1.62x | ahead | ahead |
| FactorizationBreezeJmh | Solve | n=64 | pure | 0.92x | 1.37x | 1.49x | behind | ahead |
| FactorizationBreezeJmh | Solve | n=64 | vector (gemm-routed only) | 0.92x | 1.38x | 1.50x | behind | ahead |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | 1.02x | 1.04x | 1.02x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=1024 | vector (gemm-routed only) | 1.02x | 1.04x | 1.02x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | 1.20x | 1.28x | 1.07x | ahead | ahead |
| FactorizationLargeBreezeJmh | Chol | n=512 | vector (gemm-routed only) | 1.21x | 1.28x | 1.06x | ahead | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | 0.69x | 1.38x | 2.00x | behind | ahead |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | 0.73x | 1.39x | 1.90x | behind | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | 1.38x | 1.38x | 1.00x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | vector (gemm-routed only) | 1.97x | 1.89x | 0.96x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | 1.22x | 1.23x | 1.01x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | vector (gemm-routed only) | 1.37x | 1.36x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | 1.07x | 1.85x | 1.73x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=1024 | vector (gemm-routed only) | 1.07x | 1.85x | 1.73x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | 1.08x | 1.85x | 1.71x | ahead | ahead |
| FactorizationLargeBreezeJmh | Lu | n=512 | vector (gemm-routed only) | 1.09x | 1.85x | 1.70x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.50x | 1.45x | 0.97x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=1024 | vector (gemm-routed only) | 2.17x | 2.10x | 0.97x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.41x | 1.39x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Qr | n=512 | vector (gemm-routed only) | 1.44x | 1.42x | 0.99x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | 1.06x | 1.83x | 1.73x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=1024 | vector (gemm-routed only) | 1.06x | 1.83x | 1.73x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | 1.08x | 1.81x | 1.68x | ahead | ahead |
| FactorizationLargeBreezeJmh | Solve | n=512 | vector (gemm-routed only) | 1.08x | 1.81x | 1.68x | ahead | ahead |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | 0.66x | 0.66x | 1.00x | behind | behind |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | 0.91x | 0.90x | 0.99x | behind | behind |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | 0.54x | 0.53x | 0.98x | withheld | withheld |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | 0.73x | 0.73x | 1.00x | withheld | withheld |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | 0.75x | 0.75x | 1.00x | behind | behind |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | 0.85x | 0.82x | 0.96x | withheld | withheld |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.77x | 0.87x | 1.13x | behind | behind |
| LeastSquaresBreezeJmh | Lstsq | n=16 | vector (gemm-routed only) | 0.76x | 0.88x | 1.16x | behind | behind |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.19x | 1.22x | 1.03x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=256 | vector (gemm-routed only) | 1.24x | 1.28x | 1.03x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.04x | 1.23x | 1.18x | ahead | ahead |
| LeastSquaresBreezeJmh | Lstsq | n=64 | vector (gemm-routed only) | 1.04x | 1.23x | 1.18x | ahead | ahead |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | 2.74x | 2.83x | 1.03x | ahead | ahead |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | 0.84x | 1.29x | 1.54x | behind | ahead |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | 14.63x | 24.53x | 1.68x | ahead | ahead |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | 11.80x | 11.81x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | 3.68x | 3.81x | 1.04x | ahead | ahead |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.46x | 3.46x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | 1.59x | 1.59x | 1.00x | ahead | ahead |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | 18.49x | 18.46x | 1.00x | ahead | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | 0.68x | 2.06x | 3.03x | behind | ahead |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | vector (gemm-routed only) | 0.67x | 2.06x | 3.07x | behind | ahead |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | 0.74x | 2.57x | 3.47x | behind | ahead |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | vector (gemm-routed only) | 0.74x | 2.58x | 3.49x | behind | ahead |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | 2.62x | 2.68x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | 1.58x | 1.62x | 1.03x | ahead | ahead |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | 1.65x | 1.66x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | 1.01x | 1.00x | 0.99x | tie | tie |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | 1.49x | 1.51x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | 1.00x | 0.99x | 0.99x | tie | tie |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | 0.98x | 0.99x | 1.01x | behind | tie |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | 1.04x | 1.07x | 1.03x | ahead | ahead |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | 1.03x | 1.05x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | 0.70x | 1.06x | 1.51x | behind | tie |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | 0.74x | 0.76x | 1.03x | behind | behind |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | 0.71x | 0.84x | 1.18x | behind | tie |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | 30.53x | 30.48x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | 30.96x | 30.98x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | 31.22x | 31.20x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | 40.53x | 40.51x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | 40.92x | 40.93x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | 40.95x | 40.96x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | 1.23x | 1.24x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | 1.64x | 1.64x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | 1.63x | 1.63x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | 0.85x | 0.85x | 1.00x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | 0.71x | 0.72x | 1.01x | behind | behind |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | 0.53x | 0.54x | 1.02x | behind | behind |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | 1.78x | 1.83x | 1.03x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | 2.15x | 2.20x | 1.02x | ahead | ahead |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | 1.82x | 1.83x | 1.01x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | 3.27x | 3.27x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | 3.37x | 3.38x | 1.00x | ahead | ahead |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | 3.58x | 3.57x | 1.00x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | pure | 2.96x | 3.52x | 1.19x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=3 | vector (gemm-routed only) | 2.99x | 3.53x | 1.18x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | pure | 2.47x | 2.74x | 1.11x | ahead | ahead |
| SmallDenseBreezeJmh | Det | n=4 | vector (gemm-routed only) | 2.49x | 2.72x | 1.09x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | 1.06x | 1.22x | 1.15x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=3 | vector | 1.05x | 1.21x | 1.15x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | 1.85x | 2.14x | 1.16x | ahead | ahead |
| SmallDenseBreezeJmh | Gemm | n=4 | vector | 1.90x | 2.13x | 1.12x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | pure | 2.94x | 3.11x | 1.06x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=3 | vector (gemm-routed only) | 2.98x | 3.12x | 1.05x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | pure | 2.17x | 3.06x | 1.41x | ahead | ahead |
| SmallDenseBreezeJmh | Inv | n=4 | vector (gemm-routed only) | 2.17x | 3.07x | 1.41x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | pure | 2.60x | 2.90x | 1.12x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=3 | vector (gemm-routed only) | 2.63x | 2.88x | 1.10x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | pure | 2.16x | 2.27x | 1.05x | ahead | ahead |
| SmallDenseBreezeJmh | Solve | n=4 | vector (gemm-routed only) | 2.18x | 2.24x | 1.03x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | 5.57x | 5.55x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | 4.98x | 5.02x | 1.01x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | 6.64x | 6.63x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | 5.56x | 5.57x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | 7.00x | 7.01x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | 6.65x | 6.44x | 0.97x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | 8.08x | 8.00x | 0.99x | ahead | ahead |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | 6.13x | 6.11x | 1.00x | ahead | ahead |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | 9.22x | 9.08x | 0.98x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | 8.44x | 8.28x | 0.98x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | 10.25x | 10.24x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | 8.28x | 8.46x | 1.02x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | 8.10x | 8.12x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | 6.99x | 7.01x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | 7.01x | 7.02x | 1.00x | withheld | withheld |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | 6.72x | 6.86x | 1.02x | withheld | withheld |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | 1.06x | 1.06x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | 1.09x | 1.11x | 1.02x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | 6.61x | 6.50x | 0.98x | ahead | ahead |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | 6.41x | 6.41x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | 3.49x | 3.49x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | 1.16x | 1.17x | 1.01x | tie | tie |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | 1.01x | 1.01x | 1.00x | ahead | ahead |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | 1.00x | 1.00x | 1.00x | tie | tie |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | 1.07x | 1.82x | 1.70x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | 2.25x | 2.87x | 1.28x | ahead | ahead |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | 1.32x | 2.10x | 1.59x | ahead | ahead |

## Supplementary: controlled n=16 eigen recheck

Run separately on debug node tri0035 (Slurm job 2535336, exclusive), same bundle, SHA binding, pinning and lane-B fork arguments, through `supp.sh eig16B` (`docs/verification/w21-simd-spike/x86-avx512/supp.sh`). It used more samples than the sweep: 4 forks × 10 warmup × 10 measurement iterations of 1 s, `-p n=16`, `DenseDecompositionBreezeJmh.(gale|breeze)Eig` and `SymEigenBreezeJmh.(gale|breeze)EigSym`. The scoreboard ran with `--validation --strict` and exited 0 with 4 rows validated. Files: `2026-10-10-breeze-milestone2-trillium-laneB-eig16.{json,netlib.jsonl,validation.jsonl}`.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33766.2 ± 460.2 | 12918.4 ± 43.96 | ops/s | 2.61x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97724.1 ± 208.6 | 34527.7 ± 35.92 | ops/s | 2.83x | ahead |  |

## Provenance

- Date: Slurm job 2535352, start/end/elapsed: 2026-10-10T21:13:03 2026-10-10T23:29:09 02:16:06 
- Benchmarked source: `0ff64405ea981f5a6dbf026a12e2bcc81617f412` (breeze/integration, clean tree). Compiled locally (sbt 1.11.7, JDK 25.0.1, Scala 3.7.4). The exported `benchmarksJVM/Jmh/fullClasspath` was copied to the cluster. No sbt ran on the cluster. Every fork received `-Dgale.bench.sourceSha=0ff64405ea981f5a6dbf026a12e2bcc81617f412` from the bundle's `SHA` file.
- Node: SciNet trillium `tri0468`, partition `compute`, `--exclusive --nodes=1`. CPU: 2 × AMD EPYC 9655 96-core (Zen 5), AVX-512, SMT off, 192 CPUs, NPS4, boost enabled. Load average at job start: 1831.47, 2182.92, 1811.59.
- JDK: Temurin 25+36 (`module load java/25`). `JAVA_TOOL_OPTIONS`, `JDK_JAVA_OPTIONS` and `_JAVA_OPTIONS` were unset.
- Pinning: `numactl --physcpubind=0-7 --membind=0` (the CCD of cpu0 and its NUMA node). It applies to the host JVM and is inherited by every fork.
- Netlib (blas, lapack, vectorModule) over all sidecar records (`2026-10-10-breeze-milestone2-trillium-laneB.netlib.jsonl`): [('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- Numerical preflight (`2026-10-10-breeze-milestone2-trillium-laneB.validation.jsonl`): {'passed': 302} sourceSha: ['0ff64405ea981f5a6dbf026a12e2bcc81617f412'].
- Runner: the reviewed `docs/verification/w21-simd-spike/x86-avx512/run.sh laneB` (unmodified, `GALE_BENCH_ROOT=/scratch/brad/gale-bench/m2`), i.e. `numactl … java --add-modules=jdk.incubator.vector -Dgale.bench.sourceSha=<sha> -cp <classpath> org.openjdk.jmh.Main -foe true -jvmArgsAppend "-Xms4g -Xmx4g -Dgale.bench.sourceSha=<sha> -Dgale.bench.lane=B -Dgale.bench.netlibSidecar=<out>/netlib.jsonl" -p backend=pure,vector -rf json -rff <out>/result.json '.*BreezeJmh.*'`. All other JMH settings (2 forks etc.) come from the annotations.
- Scoreboard: `python3 -I tools/bench/breeze_scoreboard.py --lane B --strict --commit 0ff64405ea981f5a6dbf026a12e2bcc81617f412 --netlib 2026-10-10-breeze-milestone2-trillium-laneB.netlib.jsonl --validation 2026-10-10-breeze-milestone2-trillium-laneB.validation.jsonl 2026-10-10-breeze-milestone2-trillium-laneB.json`. It exited 0, with 0 unpaired rows.
