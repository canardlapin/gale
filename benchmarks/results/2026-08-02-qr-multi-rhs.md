# Portable QR multi-RHS reflector court

Date: 2026-08-02

Decision: retain the original column-first reflector loop for `q < 8` and use
an allocation-free eight-column row-contiguous kernel for `q >= 8`. The shared
portable implementation is selected. No public API, backend dispatch, mutable
ownership rule, pivot/rank decision, or Scala.js dependency changes.

## Court and environment

`QrMultiRhsJmh` prepares a deterministic tall design, response matrix, and
column-pivoted QR factor outside the timed boundary. `applyQtOwned` and
`solveLeastSquaresOwned` include Gale's mandatory RHS copy and owned result;
the latter also includes triangular solve and pivot permutation.
`factorPivotedQr` protects the unchanged factorization path. The full sweep is
`n={512,2048,10000}`, `p={6,24}`, and `q={1,8,16,32,100}`.

The receipts were collected on a MacBook Pro (Apple M3 Max, 14 cores, 36 GB),
macOS 14.3 arm64, JDK 25.0.1, sbt 1.11.7, Scala 3.7.4, and JMH 1.37. JMH used
one thread, two forks, three 300 ms warmups, five 300 ms measurements, average
time in microseconds, and the GC profiler. The implementation is Gale's shared
portable dense kernel; no optional JVM backend is involved. The immutable
baseline revision is `9c0fbc2c6fac124eaabbd5c219500471962d2867`; candidate
receipts use that revision plus the issue-scoped working-tree changes.

Representative commands:

```text
sbt 'benchmarksJVM/Jmh/run -prof gc -rf json -rff benchmarks/results/2026-08-02-qr-multi-rhs-candidate-exact.json -p n=2048 -p p=6 -p q=16 gale.bench.QrMultiRhsJmh.*'
sbt 'benchmarksJVM/Jmh/run -prof gc -rf json -rff benchmarks/results/2026-08-02-qr-multi-rhs-candidate-sweep.json gale.bench.QrMultiRhsJmh.*(applyQtOwned|solveLeastSquaresOwned)'
```

Raw receipts: [baseline exact](./2026-08-02-qr-multi-rhs-baseline-exact.json),
[baseline sweep](./2026-08-02-qr-multi-rhs-baseline-sweep.json),
[selected exact](./2026-08-02-qr-multi-rhs-candidate-exact.json),
[selected sweep](./2026-08-02-qr-multi-rhs-candidate-sweep.json), and
[focused q=1 control](./2026-08-02-qr-multi-rhs-candidate-q1-control.json).

## Exact regress4s shape

The exact `n=2048,p=6,q=16` court improved materially in two independent
selected-candidate runs. Values below are the focused baseline/candidate run;
normalized allocation includes the copy and owned result.

| method | baseline us/op | candidate us/op | speedup | baseline B/op | candidate B/op |
| --- | ---: | ---: | ---: | ---: | ---: |
| applyQT owned | 635.368 | 147.375 | 4.31x | 262,254.6 | 262,219.5 |
| least squares owned | 623.875 | 149.688 | 4.17x | 263,070.6 | 263,059.5 |
| pivoted QR factorization control | 180.088 | 163.858 | 1.10x | 213,514.5 | 213,502.5 |

The full sweep independently measured 144.497 us/op for `applyQT` and 148.358
us/op for least squares at the same shape. Thus the improvement is not a
single focused-run result, allocation does not increase, and the factorization
control remains comfortably inside the five-percent regression limit.

Across all 48 wide-RHS comparisons, every selected candidate improved. The
speedup ranges were 4.55x-5.64x at q=8, 3.82x-5.29x at q=16, 3.83x-4.29x at
q=32, and 2.87x-3.56x at q=100. Normalized allocation did not increase in any
wide case.

## Protected scalar control and crossover

Applying the blocked loop at every width was rejected: its worst q=1 point,
`applyQT` at `n=512,p=6`, regressed from 8.920 to 10.850 us/op (+21.64%). The
raw rejected receipts are [monolithic exact](./2026-08-02-qr-multi-rhs-candidate-monolithic-exact.json)
and [monolithic sweep](./2026-08-02-qr-multi-rhs-candidate-monolithic-sweep.json).

The selected dispatch preserves the original scalar loop below eight columns.
In the canonical full sweep, the worst q=1 regression was +2.03%; the focused
small-shape repeat measured 8.967 us/op for `applyQT` and 8.972 us/op for least
squares. Eight is therefore the measured crossover used by this court; widths
2-7 were not swept and deliberately retain the scalar path. A q=17 test covers
the blocked kernel's partial final block.

Compact-WY QR work was not selected because it accelerates factorization, not
post-factor application to a matrix RHS. A panel/GEMM route would add workspace
and ownership complexity after the direct blocked loop already cleared the
downstream gate. No Vector, FFM, or native candidate was benchmarked or enabled:
the portable shared route wins across the required sweep and avoids an
unnecessary JVM-only crossover policy.

## Correctness and integration evidence

The matrix-RHS tests compare wide `applyQ`/`applyQT` calls with independent
single-column calls and explicit Q, and compare least squares with independent
`DVec` solves. They cover q=1, q=8, a q=17 partial block, pivoted QR, exact and
near rank deficiency, extreme scales, non-contiguous inputs, storage ownership,
and deterministic repeats. `QRSuite` passes 19/19 on both JVM and Scala.js.

The current implementation passed:

```text
sbt testAllFull compileAll parityTest interopBreezeTest benchCompile
```

That run included 596 core JVM tests, 586 core Scala.js tests, 41 laws tests per
platform, 8 Ravel interop tests per platform, 45 parity tests, 24 Breeze interop
tests, full-optimized Scala.js links, Scala-next compilation, and JMH
compilation. After strengthening the single-column differential assertions,
the focused 19-test JVM and Scala.js QR suites were rerun successfully.

The unchanged regress4s JVM/Scala.js compile, test, and benchmark-oracle gates
also pass against this local Gale checkout. Its signed two-fork paired court is:

| method | time ms/op | allocated B/op |
| --- | ---: | ---: |
| M2 complete end to end | 0.704 | 1,514,388 |
| 16 scalar regress4s fits | 5.872 | 9,748,871 |

This is an 8.34x batch-over-scalar speedup, clearing the required 5x gate and
improving the prior 4.74x result without changing regress4s. The retained raw
receipt is [regress4s M2 candidate](./2026-08-02-regress4s-m2-gale-candidate.json).

The implementation changes only private traversal in
`DenseDecompositions.applyReflectorsLeft`; no Gale public API change is
necessary.
