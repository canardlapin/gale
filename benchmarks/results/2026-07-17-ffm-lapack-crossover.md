# FFM LAPACK crossover — Accelerate / JDK 22 (2026-07-17)

## Environment and protocol

- Apple ARM64, macOS 14.3
- OpenJDK 22
- runtime-discovered macOS Accelerate BLAS/LAPACK
- public Gale factorization/solver facades; native measurements include heap to
  native copies, LAPACK, native to heap copies, and typed-result allocation
- JMH 1.37, average time, one thread, 2 forks, 5 x 500 ms warmups, and
  8 x 500 ms measurements

## Factorization and eigen kernels

Times are microseconds per operation. Speedup is pure Gale divided by FFM.

| Operation | n | Pure Gale | FFM | Speedup |
| --- | ---: | ---: | ---: | ---: |
| LU | 64 | 32.337 | 27.796 | 1.16x |
| LU | 128 | 211.656 | 116.130 | 1.82x |
| LU | 256 | 1,885.326 | 1,147.102 | 1.64x |
| LU | 512 | 14,183.284 | 3,859.524 | 3.68x |
| Cholesky | 64 | 20.455 | 15.411 | 1.33x |
| Cholesky | 128 | 148.977 | 526.177 | 0.28x |
| Cholesky | 256 | 1,359.614 | 1,144.523 | 1.19x |
| Cholesky | 512 | 12,891.207 | 3,322.594 | 3.88x |
| QR (2n x n) | 64 | 133.763 | 106.577 | 1.26x |
| QR (2n x n) | 128 | 1,139.751 | 4,715.555 | 0.24x |
| QR (2n x n) | 256 | 6,700.868 | 17,644.589 | 0.38x |
| QR (2n x n) | 512 | 52,192.100 | 414,599.053 | 0.13x |
| symmetric eigenvalues | 64 | 177.695 | 151.054 | 1.18x |
| symmetric eigenvalues | 128 | 1,136.931 | 550.148 | 2.07x |
| symmetric eigenvalues | 256 | 7,955.065 | 4,878.870 | 1.63x |
| symmetric eigenvalues | 512 | 65,109.310 | 23,509.122 | 2.77x |

LU has a stable measured win, so the Accelerate default is conservatively
`nativeLuMinSize = 128`. Cholesky and QR are not monotone: a single lower-bound
threshold cannot safely describe either route. Both therefore remain disabled
by default and available through explicit `FfmBlasThresholds` overrides. The
symmetric-eigen provider is real and conforming, but Gale does not yet expose a
backend-routed public eigen facade, so no automatic threshold is claimed.

## End-to-end scenarios

| Scenario | n | Pure Gale | Forced FFM | Speedup |
| --- | ---: | ---: | ---: | ---: |
| solve via LU | 64 | 35.627 | 32.596 | 1.09x |
| solve via LU | 128 | 221.298 | 144.594 | 1.53x |
| solve via LU | 256 | 1,925.750 | 1,203.644 | 1.60x |
| least squares via QR | 64 | 154.403 | 114.908 | 1.34x |
| least squares via QR | 128 | 1,203.232 | 4,608.864 | 0.26x |
| least squares via QR | 256 | 7,255.185 | 20,306.885 | 0.36x |

The scenario sweep validates the policy rather than merely timing isolated
calls: default LU routing improves solve, while default QR routing would regress
least squares by 3.8x at n=128 and 2.8x at n=256. The benchmark harness forces
all native thresholds to zero on purpose; these numbers are evidence for the
defaults, not the behavior of the defaults themselves.
