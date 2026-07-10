# Breeze comparison — after BLAS kernel tuning (2026-07-10)

Same JMH protocol as the baseline (`2026-07-10-breeze-baseline.md`). Kernel changes:
4-lane ddot/daxpy, unrolled row-major gemv, 4x4 register-tiled gemm micro-kernel
(GemmBlock retuned 64→128), dedicated syrk kernel for `a.t * a`. Note: this sweep ran
under background load; treat single-digit percentage moves as noise.

| class | op | n | before g/b | after g/b |
|---|---|---:|---:|---:|
| BlasL1BreezeJmh | Axpy | 65536 | 0.43x | **0.65x** |
| BlasL1BreezeJmh | Axpy | 262144 | 0.40x | **0.63x** |
| BlasL1BreezeJmh | Axpy | 1048576 | 0.38x | **0.68x** |
| BlasL1BreezeJmh | Dot | 65536 | 0.29x | **0.98x** |
| BlasL1BreezeJmh | Dot | 262144 | 0.34x | **0.77x** |
| BlasL1BreezeJmh | Dot | 1048576 | 0.38x | **0.86x** |
| BlasL1BreezeJmh | Norm | 65536 | 5.45x | **4.93x** |
| BlasL1BreezeJmh | Norm | 262144 | 3.47x | **3.71x** |
| BlasL1BreezeJmh | Norm | 1048576 | 3.11x | **3.17x** |
| BlasL2BreezeJmh | Gemv | 256 | 0.50x | **0.95x** |
| BlasL2BreezeJmh | Gemv | 1024 | 0.29x | **0.96x** |
| BlasL2BreezeJmh | Gemv | 2048 | 0.29x | **0.98x** |
| BlasL2BreezeJmh | GemvT | 256 | 0.83x | **0.79x** |
| BlasL2BreezeJmh | GemvT | 1024 | 0.84x | **1.08x** |
| BlasL2BreezeJmh | GemvT | 2048 | 0.92x | **0.90x** |
| BlasL3BreezeJmh | AtA | 16 | 0.26x | **0.36x** |
| BlasL3BreezeJmh | AtA | 64 | 0.21x | **0.74x** |
| BlasL3BreezeJmh | AtA | 256 | 0.12x | **0.77x** |
| BlasL3BreezeJmh | Gemm | 16 | 0.43x | **0.77x** |
| BlasL3BreezeJmh | Gemm | 64 | 0.42x | **0.70x** |
| BlasL3BreezeJmh | Gemm | 256 | 0.39x | **0.75x** |
| BlasL3BreezeJmh | GemmTall | 16 | 0.37x | **0.64x** |
| BlasL3BreezeJmh | GemmTall | 64 | 0.51x | **0.84x** |
| BlasL3BreezeJmh | GemmTall | 256 | 0.41x | **0.80x** |
| FactorizationBreezeJmh | Chol | 16 | 2.99x | **3.35x** |
| FactorizationBreezeJmh | Chol | 64 | 2.41x | **2.23x** |
| FactorizationBreezeJmh | Chol | 256 | 0.93x | **0.70x** |
| FactorizationBreezeJmh | Lu | 16 | 1.21x | **1.19x** |
| FactorizationBreezeJmh | Lu | 64 | 0.84x | **0.82x** |
| FactorizationBreezeJmh | Lu | 256 | 0.81x | **0.86x** |
| FactorizationBreezeJmh | Qr | 16 | 1.36x | **1.31x** |
| FactorizationBreezeJmh | Qr | 64 | 0.75x | **0.83x** |
| FactorizationBreezeJmh | Qr | 256 | 0.26x | **0.23x** |
| FactorizationBreezeJmh | Solve | 16 | 1.22x | **1.27x** |
| FactorizationBreezeJmh | Solve | 64 | 0.83x | **0.82x** |
| FactorizationBreezeJmh | Solve | 256 | 0.85x | **0.62x** |
| LeastSquaresBreezeJmh | Lstsq | 16 | 0.78x | **0.80x** |
| LeastSquaresBreezeJmh | Lstsq | 64 | 0.78x | **0.63x** |
| LeastSquaresBreezeJmh | Lstsq | 256 | 0.28x | **0.27x** |
| SymEigenBreezeJmh | EigSym | 16 | 2.98x | **3.12x** |
| SymEigenBreezeJmh | EigSym | 64 | 1.18x | **1.35x** |
| SymEigenBreezeJmh | EigSym | 128 | 0.68x | **0.81x** |

**Remaining large-n gaps** (candidates for follow-up): QR 0.23x and least-squares 0.27x
at n=256 (gale's QR is unblocked Householder — needs WY/blocked reflectors, a separate
algorithmic change, not a kernel tune); axpy ~0.65x (store-bound); triangular-solve-heavy
paths (Solve/Chol at 256). Dot/gemv are at parity; norm and small-n factorizations stay ahead.
