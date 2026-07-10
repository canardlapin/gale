# Breeze comparison — FMA, tiled SYRK, and blocked QR (2026-07-10)

This follow-up measures the two performance changes after
`2026-07-10-breeze-tuned.md`: platform FMA plus assign-oriented dense kernels,
and hybrid blocked/compact-WY Householder QR. Breeze 2.1.0 is using its pure-Java
BLAS/LAPACK fallback, as in the earlier reports.

## Measurement status

The correctness results below are release-grade. The performance results are
**directional only**: an unrelated process held roughly one CPU continuously, so
the machine never reached the quiet state required by the benchmark protocol.
The final focused run used OpenJDK 25.0.1, one fork, seven 500 ms warmups, and
seven 500 ms measurements. A quiet default two-fork refresh remains desirable;
the paired results nevertheless have low within-fork variance and reverse the
large gaps by wide margins.

## Largest-size algorithmic results

| operation | n | Gale ops/s | Breeze ops/s | Gale / Breeze | previous |
|---|---:|---:|---:|---:|---:|
| A^T A | 256 | 217.231 +/- 2.336 | 113.795 +/- 0.627 | **1.91x** | 0.77x |
| QR | 256 | 387.748 +/- 1.519 | 214.768 +/- 0.760 | **1.81x** | 0.23x |
| least squares (m = 4n) | 256 | 65.268 +/- 0.272 | 29.620 +/- 0.097 | **2.20x** | 0.27x |

The small unblocked QR path remains ahead: 1.79x at n=16 and 1.41x at n=64 in
a separate one-fork run with five 500 ms warmups and measurements.

Earlier focused measurements from the same contested session put dot at 0.97x,
GEMM at 1.04x, tall GEMM at 1.05x, and AXPY at 0.64x. AXPY remains the one clear
dense-kernel gap and is store-bandwidth-bound.

## Factorization reassessment

The n=256 dips in the previous background-loaded sweep did not reproduce:

| operation | Gale ops/s | Breeze ops/s | Gale / Breeze |
|---|---:|---:|---:|
| Cholesky | 765.1 | 828.5 | **0.92x** |
| LU | 557.8 | 646.9 | **0.86x** |
| solve | 539.3 | 634.5 | **0.85x** |

These 8–15% gaps do not justify adding blocked Cholesky/LU and pivoting
complexity in this arc. They should be checked again in the eventual quiet
two-fork sweep.

## Correctness gates

- JVM core: 245 passed.
- Scala.js core: 242 passed.
- Laws: 17 passed on JVM and 17 on Scala.js.
- Breeze numerical parity: 25 passed.
- JVM JMH and Scala.js benchmark harnesses compile.
- QR-specific coverage includes tall and wide blocked shapes, partial final
  panels, exact rank deficiency, extreme small/large scales, a near-maximum
  finite Householder overflow case, least squares, and workspace reuse.

## Implementation summary

- JVM `Math.fma` with a Scala.js multiply-add fallback.
- GEMM `beta = 0` assign path and 4x4 tiled assign-only SYRK.
- Unblocked QR below `min(m,n) = 96`; 32-column compact-WY panels above it.
- Stable scaled Householder norms and overflow-safe normalized reflectors.
- Packed `V^T`, compact `T`, and update `W` share one reusable workspace.
- Trailing QR updates use the tuned GEMM path; public compact-reflector and lazy
  `Q` contracts are unchanged.
