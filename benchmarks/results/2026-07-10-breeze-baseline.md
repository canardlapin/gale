# Breeze comparison baseline — 2026-07-10

JMH forked sweep (`-f 1 -wi 3 -i 5 -r 500ms -w 500ms -prof gc`), Breeze 2.1.0 on the
pure-Java netlib fallback (no native BLAS). Machine-specific numbers; the ratio column
is the signal. Baseline for the `gale-blas-kernel-tuning` bead.

| class | op | n | gale ops/s | breeze ops/s | gale/breeze |
|---|---|---:|---:|---:|---:|
| BlasL1BreezeJmh | Axpy | 65536 | 33686.0 | 77910.1 | **0.43x** |
| BlasL1BreezeJmh | Axpy | 262144 | 8137.8 | 20572.7 | **0.40x** |
| BlasL1BreezeJmh | Axpy | 1048576 | 1911.2 | 5055.6 | **0.38x** |
| BlasL1BreezeJmh | Dot | 65536 | 19794.7 | 68275.5 | **0.29x** |
| BlasL1BreezeJmh | Dot | 262144 | 4822.1 | 14212.4 | **0.34x** |
| BlasL1BreezeJmh | Dot | 1048576 | 1296.7 | 3448.5 | **0.38x** |
| BlasL1BreezeJmh | Norm | 65536 | 20953.9 | 3847.1 | **5.45x** |
| BlasL1BreezeJmh | Norm | 262144 | 4952.0 | 1426.1 | **3.47x** |
| BlasL1BreezeJmh | Norm | 1048576 | 1250.3 | 401.6 | **3.11x** |
| BlasL2BreezeJmh | Gemv | 256 | 28479.5 | 56966.8 | **0.50x** |
| BlasL2BreezeJmh | Gemv | 1024 | 1328.5 | 4546.3 | **0.29x** |
| BlasL2BreezeJmh | Gemv | 2048 | 321.0 | 1096.3 | **0.29x** |
| BlasL2BreezeJmh | GemvT | 256 | 44675.6 | 53593.4 | **0.83x** |
| BlasL2BreezeJmh | GemvT | 1024 | 2781.6 | 3296.7 | **0.84x** |
| BlasL2BreezeJmh | GemvT | 2048 | 722.0 | 781.1 | **0.92x** |
| BlasL3BreezeJmh | AtA | 16 | 114357.4 | 439987.3 | **0.26x** |
| BlasL3BreezeJmh | AtA | 64 | 1392.7 | 6550.9 | **0.21x** |
| BlasL3BreezeJmh | AtA | 256 | 13.0 | 112.6 | **0.12x** |
| BlasL3BreezeJmh | Gemm | 16 | 673442.4 | 1572783.5 | **0.43x** |
| BlasL3BreezeJmh | Gemm | 64 | 13733.3 | 32336.9 | **0.42x** |
| BlasL3BreezeJmh | Gemm | 256 | 174.7 | 445.1 | **0.39x** |
| BlasL3BreezeJmh | GemmTall | 16 | 170090.4 | 455598.8 | **0.37x** |
| BlasL3BreezeJmh | GemmTall | 64 | 3397.3 | 6718.2 | **0.51x** |
| BlasL3BreezeJmh | GemmTall | 256 | 44.6 | 109.2 | **0.41x** |
| FactorizationBreezeJmh | Chol | 16 | 1433524.3 | 480110.8 | **2.99x** |
| FactorizationBreezeJmh | Chol | 64 | 54951.1 | 22781.0 | **2.41x** |
| FactorizationBreezeJmh | Chol | 256 | 754.4 | 812.2 | **0.93x** |
| FactorizationBreezeJmh | Lu | 16 | 1379580.1 | 1142314.4 | **1.21x** |
| FactorizationBreezeJmh | Lu | 64 | 31984.1 | 37954.0 | **0.84x** |
| FactorizationBreezeJmh | Lu | 256 | 517.1 | 636.5 | **0.81x** |
| FactorizationBreezeJmh | Qr | 16 | 452782.2 | 333859.0 | **1.36x** |
| FactorizationBreezeJmh | Qr | 64 | 10521.9 | 13957.9 | **0.75x** |
| FactorizationBreezeJmh | Qr | 256 | 53.9 | 210.1 | **0.26x** |
| FactorizationBreezeJmh | Solve | 16 | 1022090.0 | 837794.2 | **1.22x** |
| FactorizationBreezeJmh | Solve | 64 | 29305.5 | 35461.4 | **0.83x** |
| FactorizationBreezeJmh | Solve | 256 | 534.6 | 625.9 | **0.85x** |
| LeastSquaresBreezeJmh | Lstsq | 16 | 80661.6 | 103241.5 | **0.78x** |
| LeastSquaresBreezeJmh | Lstsq | 64 | 1461.0 | 1866.7 | **0.78x** |
| LeastSquaresBreezeJmh | Lstsq | 256 | 8.1 | 29.0 | **0.28x** |
| SymEigenBreezeJmh | EigSym | 16 | 62158.4 | 20872.6 | **2.98x** |
| SymEigenBreezeJmh | EigSym | 64 | 1378.2 | 1163.9 | **1.18x** |
| SymEigenBreezeJmh | EigSym | 128 | 157.7 | 231.3 | **0.68x** |
