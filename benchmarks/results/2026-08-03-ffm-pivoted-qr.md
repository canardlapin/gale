# Copy-inclusive FFM pivoted-QR assessment

Date: 2026-08-03

Decision: no-go. Keep Gale's portable column-pivoted QR authoritative and do
not create a production FFM implementation bead. The isolated `dgeqp3`
prototype fails both the exact downstream performance threshold and Gale's
deterministic permutation contract.

## Scope and environment

The benchmark-only prototype uses LAPACK `dgeqp3` through JDK FFM. It is not
reachable from Gale's production backend. Every timed call includes logical
row-major `DMat` to owned column-major staging, heap-to-native copy, workspace
query and allocation, `dgeqp3`, native-to-heap copy, permutation conversion,
rank decision, and owned Gale reflector, R, tau, and permutation results.

The court ran on a MacBook Pro with an Apple M3 Max (14 cores, 36 GB), macOS
14.3 arm64, JDK 25.0.1, sbt 1.11.7, Scala 3.7.4, JMH 1.37, and the macOS
Accelerate framework version 1.11. `VECLIB_MAXIMUM_THREADS=1` was set before
the sbt/JMH process started. The portable implementation is commit
`bc335cf9212f9c681fd1cbd7fac3ceabb8201635`; the assessment harness was built
on the subsequent evidence commit `9e005e7` without changing core QR code.

The JMH court used one thread, two forks, three 300 ms warmups, five 300 ms
measurements, average time in microseconds, and the GC profiler:

```text
VECLIB_MAXIMUM_THREADS=1 sbt \
  'benchmarksFfm/Jmh/run -f 2 -wi 3 -i 5 -w 300ms -r 300ms -t 1 -prof gc -rf json -rff /absolute/path/result.json gale.bench.FfmPivotedQrJmh.*'
```

Raw receipts: [JMH](./2026-08-03-ffm-pivoted-qr.json) and
[numerical court](./2026-08-03-ffm-pivoted-qr-correctness.txt).

## Copy-inclusive timing

Speedup is portable time divided by FFM time; values above 1 favor FFM.
Normalized allocation is JVM heap only and excludes the native footprint
reported separately below.

| shape | portable us/op | FFM us/op | speedup | portable B/op | FFM B/op |
| --- | ---: | ---: | ---: | ---: | ---: |
| `512 x 3` | 9.752 | 14.859 | 0.66x | 29,096 | 63,392 |
| `512 x 5` | 18.344 | 22.799 | 0.80x | 45,513 | 104,625 |
| `512 x 6` | 23.544 | 26.675 | 0.88x | 53,713 | 125,113 |
| `512 x 8` | 37.160 | 39.382 | 0.94x | 70,129 | 166,081 |
| `512 x 16` | 116.490 | 67.319 | 1.73x | 135,796 | 330,050 |
| `512 x 24` | 231.685 | 123.725 | 1.87x | 201,517 | 494,019 |
| `1024 x 3` | 18.433 | 27.274 | 0.68x | 57,769 | 124,833 |
| `1024 x 5` | 36.076 | 45.706 | 0.79x | 90,570 | 206,785 |
| `1024 x 6` | 47.755 | 53.793 | 0.89x | 106,961 | 247,993 |
| `1024 x 8` | 72.901 | 70.528 | 1.03x | 139,762 | 329,922 |
| `1024 x 16` | 248.196 | 151.113 | 1.64x | 271,024 | 657,732 |
| `1024 x 24` | 492.670 | 251.382 | 1.96x | 402,319 | 985,576 |
| `2048 x 3` | 37.694 | 53.740 | 0.70x | 115,113 | 247,713 |
| `2048 x 5` | 78.331 | 83.224 | 0.94x | 180,682 | 411,586 |
| `2048 x 6` | 100.639 | 95.291 | 1.06x | 213,459 | 493,490 |
| `2048 x 8` | 159.631 | 133.914 | 1.19x | 279,046 | 657,603 |
| `2048 x 16` | 515.222 | 299.841 | 1.72x | 541,457 | 1,313,150 |
| `2048 x 24` | 987.051 | 491.161 | 2.01x | 803,770 | 1,968,698 |
| `4096 x 3` | 75.356 | 94.666 | 0.80x | 229,802 | 493,474 |
| `4096 x 5` | 156.242 | 157.859 | 0.99x | 360,917 | 821,188 |
| `4096 x 6` | 200.646 | 186.897 | 1.07x | 426,492 | 985,016 |
| `4096 x 8` | 319.703 | 258.162 | 1.24x | 557,660 | 1,312,994 |
| `4096 x 16` | 1,032.875 | 789.150 | 1.31x | 1,082,171 | 2,623,950 |
| `4096 x 24` | 2,003.357 | 1,251.551 | 1.60x | 1,606,612 | 3,934,814 |
| `10000 x 3` | 181.815 | 219.733 | 0.83x | 560,448 | 1,201,981 |
| `10000 x 5` | 363.827 | 314.722 | 1.16x | 880,588 | 2,002,098 |
| `10000 x 6` | 479.903 | 585.522 | 0.82x | 1,040,618 | 2,402,285 |
| `10000 x 8` | 785.459 | 826.901 | 0.95x | 1,360,684 | 3,202,443 |
| `10000 x 16` | 2,533.493 | 1,662.172 | 1.52x | 2,640,865 | 6,402,594 |
| `10000 x 24` | 4,874.030 | 2,682.994 | 1.82x | 3,921,052 | 9,602,763 |

The exact regress4s widths are `p=3,5,6`. None reaches the required 1.5x
copy-inclusive speedup at any measured height. Their best point is only 1.16x
at `10000 x 5`; `4096 x 6` is 1.07x, and `10000 x 6` is 0.82x. The point
estimates therefore fail before applying the additional 20% uncertainty safety
margin, and there is no monotone adjacent-size crossover. Wider `p=16/24`
sometimes exceeds 1.5x, but `p=16` falls to 1.31x at `4096`, and neither width
is an exact downstream shape in this epic.

FFM heap allocation is 2.1-2.8x the portable route across the measured court.
This is expected from column-major staging and independent owned Gale result
arrays, but it makes the route unattractive for the M1/M5 allocation-sensitive
work even where vendor computation is faster.

## Native memory and ownership

Each call creates one confined arena and closes it before returning. Peak native
payload includes the `m x p` matrix, `min(m,p)` tau values, `p` pivot integers,
the queried LAPACK work array, the query scalar, and call scalars. It ranges
from 13,424 bytes at `512 x 3` to 1,927,100 bytes at `10000 x 24`. The exact
`4096 x 6` and `10000 x 6` footprints are 198,596 and 481,988 bytes. Arena
allocator bookkeeping and alignment overhead are excluded. These payload
values are absent from JMH's GC allocation metric.

Gale's `NativeDMat` lets a real caller retain native storage for explicit GEMM,
but no current QR factorization consumes that storage or returns a factor whose
ownership can outlive the call. A persistent-native QR comparison would require
a new ownership and factor API rather than reuse by an existing caller, so it
is outside this assessment and cannot justify production routing.

## Numerical contract

The prototype passed reconstruction, thin-probe orthogonality, least-squares
residual, rank-deficient rank, extreme finite scaling, repeated native
determinism, and typed RHS-dimension failure checks. It is confined to the
JDK-only benchmark module, so shared and Scala.js sources cannot reference it.

All 30 ordinary fixtures produced exactly the same permutation as Gale, and an
exact duplicate-column fixture preserved the same tie choice. However, on 256
deterministic near-tied correlated designs, Accelerate and Gale selected
different valid permutations in 45 cases. Gale recomputes exact candidate
norms and owns a deterministic pivot rule; `dgeqp3` uses its vendor algorithm
and cannot promise that same permutation. A numerically valid but different
pivot order cannot silently implement the existing `QROptions` contract.

## Threshold decision

The production-bead threshold is conjunctive. This prototype fails:

- 1.5x copy-inclusive speedup on exact downstream shapes;
- the 20% margin beyond measurement uncertainty on those shapes;
- a monotone crossover across two adjacent tall sizes; and
- deterministic permutation authority on adversarial inputs.

Only Accelerate 1.11 was available and measured. OpenBLAS and MKL remain
unswept, so there is also no cross-family evidence for an explicit supported
policy. Portable QR remains authoritative; no hidden dispatch, native storage
API, or production implementation bead is created.

## Verification

The benchmark-only prototype compiled on JDK 25 and its numerical assessment
main completed with `status=ok`. Gale's full `testAllFull compileAll parityTest
interopBreezeTest benchCompile` gate passed, including 610 JVM and 600
Scala.js core tests, 41 laws per platform, 8 Ravel tests per platform, 45
backend-parity tests, 24 Breeze interop tests, full Scala.js optimization,
Scala-next compilation, and JMH annotation processing. The FFM backend's 21
loader, BLAS, and LAPACK tests also passed. The repository does not define a
`scalafmtCheckAll` task; compilation applied the configured warning policy.
