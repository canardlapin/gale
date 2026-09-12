# Packed banded SPD qualification (PHRF-06)

Date: 2026-09-12. Provider base: `288596690a4027157b40d5ffbf36e7d24d98e009`
(Gale main). The independent `release/0.1-stabilization` checkout was not
changed. This note qualifies the banded factor and solve capability, not an
fMRI estimator or a statistical calibration.

## Contract

`BandedCholesky.factorLower` reads row-packed lower bands with
`bands(i, d) = A(i, i-d)` and returns an immutable `ExactSolveFactor`.
Column zero is the diagonal; inactive leading padding is ignored and zeroed.
The general `Banded` adapter reads its lower triangle under the dense
Cholesky symmetry convention. The kernels do not expand the coefficient
matrix to dense square storage.

Pure inputs are preserved. `DMatBuilder.consumeBandedCholesky` transfers
storage and closes the builder on both success and failure. Factor storage
is `n (b+1)` doubles; factorization uses an additional `n` doubles to compute
the symmetric matrix 1-norm. Factorization performs
`sum_i m_i (m_i-1) / 2` inner dot-product steps, where `m_i = min(i,b)`.
This gives `O(n b²)` factorization and `O(n b)` storage for positive bandwidth.
The diagonal case is linear in `n`.

Vector and multiple-RHS solves reuse the factor. Destination solves operate
in place on `MutableDVec` (including internal strided views) or an open
`DMatBuilder`; each RHS has two triangular passes and `O(n b)` work, without
numerical scratch buffers. `solveTranspose` solves the symmetric `Aᵀ` system;
separate `solveLower` and `solveLowerTranspose` expose `L` and `Lᵀ` solves.
Dimensions and finite RHS values are checked before writes; arithmetic
overflow can leave a partially solved destination and returns a typed error.

The log determinant sums logarithms of squared diagonal pivots. Conditioning
uses absolute-value triangular comparison systems to bound the inverse
1-norm, separately from the minimum/maximum pivot ratio. It is a conservative
bound in exact arithmetic, not a directed-rounding certificate or calibrated
admission threshold. The bound can overflow to infinity even when a useful
solve exists; the reciprocal then becomes zero. Near cancellation can make
it pessimistic. Consumers must not use `pivotRatio` as reciprocal condition.

## Independent evidence

- Discrete Laplacian with diagonal 2 and off-diagonal -1: analytic lower
  factor, determinant `n+1`, Green's-function inverse, and quadratic solution
  to an all-ones RHS, including `n=300`.
- Known lower factors at diagonal, narrow and full bandwidth: dense products
  reconstruct the original SPD matrix and independent RHS solutions.
- Algebraic laws across nine size/bandwidth combinations: dense residual,
  superposition, reversal, positive scaling, log-determinant scaling,
  transpose equivalence and conditioning bound versus explicit inverse norm.
- 18 seeded size/bandwidth cases compared with independent Breeze Cholesky
  and general solves, each with four RHS columns. This host uses netlib's
  Java fallback because native JNILAPACK is unavailable.
- Finite-input and pivot failures, absolute pivot tolerance, padding,
  zero-sized systems, zero-column RHS, strided views, destination ownership,
  rejected writes after transfer, and explicit conditioning overflow.
- A 100000-row bandwidth-2 case solves to ones while retaining 300000 factor
  entries (2.4 MB of numerical factor storage); a dense factor would require
  80 GB. This checks the narrow-band execution path, not peak process RSS.

## Validation and measurements

Local host: Apple M3 Max (arm64), Homebrew JDK 25.0.1, Scala 3.7.4,
sbt 1.11.7, Node 26.7.0. Hosted CI uses the repository's JDK 21 / Node 22
lanes. The local runs below were sequential when sharing provider outputs.

| Gate | Result |
| --- | --- |
| Full `coreJVM/test` | 652 passed, including the final 12-test banded suite |
| Full `lawsJVM/test` | 54 passed |
| Full `parity/test` | 84 passed, including 18 banded cases in one test |
| Clean full `coreJS/test` | 642 passed, including the final 12-test banded suite |
| Clean full `lawsJS/test` | 54 passed |
| `coreJS/Test/fullLinkJS`, `lawsJS/Test/fullLinkJS` | Passed |
| `scalafmtCheckAll` | Passed |
| `docsCheck` | Passed; existing Scaladoc and mdoc warnings outside this change remain |
| ScalaFIM `scalafimCompileAll` with explicit local Gale override | Passed, warning-clean |
| ScalaFIM banded provider, profile reduction, condition fit and compact runtime suites | 8 passed on JVM, 8 on JS |

ScalaFIM consumer base was `fd62ea7aabb25ccd415faacfedb7243152861502`.
One concurrent provider/consumer run wrote overlapping compiler outputs and
was discarded after stale TASTy warnings and missing JS IR. The clean,
sequential rerun above replaced that evidence; no source workaround was used.

### JMH scaling and allocations

One thread, one fork, three 500 ms warmups and five 500 ms measurements,
JMH 1.37 average-time mode with `-prof gc`, `-Xms256m -Xmx512m`. The pure
factor benchmark includes the input band copy and matrix-norm calculation.
The conditioning column also evaluates the lazy comparison bound. In-place
solve timings include resetting the RHS values, with no per-invocation setup
that could hide numerical allocations.

Means in microseconds per operation (four-RHS timing is for all four):

| N | b | Factor | Factor + conditioning | One RHS in place | Four RHS in place |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 300 | 4 | 13.7 | 19.1 | 9.0 | 13.6 |
| 300 | 16 | 53.4 | 63.6 | 21.2 | 41.5 |
| 300 | 64 | 404.5 | 436.9 | 74.2 | 147.7 |
| 1200 | 4 | 55.2 | 77.6 | 35.6 | 54.8 |
| 1200 | 16 | 198.5 | 245.1 | 90.1 | 171.3 |
| 1200 | 64 | 1785.5 | 1972.3 | 321.2 | 696.7 |
| 4800 | 4 | 206.0 | 315.8 | 140.9 | 218.8 |
| 4800 | 16 | 817.2 | 985.1 | 369.1 | 740.4 |
| 4800 | 64 | 7293.8 | 8018.7 | 1314.1 | 2971.4 |

The raw [JMH receipt](../benchmarks/results/2026-09-12-banded-cholesky.json)
retains every sample and confidence interval. These are host measurements,
not universal latency guarantees. The largest four-RHS case had a wide
99.9% confidence interval (about ±957 microseconds), retained without trimming.

At fixed bandwidth 64, increasing N from 1200 to 4800 increases factor time
4.09× and one-RHS solve time 4.09×. Factor allocations follow the expected
`8 N (b+2)` bytes plus factor/diagnostic objects; at N=1200, b=64 the measured
pure-factor allocation is 633961 B/op and factor plus conditioning is
643627 B/op. The extra comparison scratch is N doubles. In-place solves
measure roughly 16–57 B/op across the grid, including profiler/background
allocation attributed per operation; they do not allocate N- or Nb-sized
numerical buffers. This is not a zero-object-allocation claim.

A [wider-band run](../benchmarks/results/2026-09-12-banded-cholesky-wide.json)
keeps N=1200 and the same JMH settings, overriding bandwidth to 128 and 256:

| b | Factor (microseconds) | One RHS (microseconds) |
| ---: | ---: | ---: |
| 128 | 6288.4 | 668.6 |
| 256 | 26084.4 | 1257.2 |

Doubling bandwidth from 128 to 256 increases factor time 4.15× and solve
time 1.88×, approaching the quadratic versus linear operation counts.
Finite boundaries and per-row log/pivot/validation work mean wall time need
not follow an exact 4× or 2× ratio. These results support the implemented
banded work/storage bounds; they do not qualify a downstream trial workflow.

Reproduce from the repository root (use an absolute output path because JMH
runs in the benchmark project directory):

```sh
sbt 'benchmarksJVM/Jmh/run -wi 3 -i 5 -w 500ms -r 500ms -f 1 -t 1 -jvmArgsAppend "-Xms256m -Xmx512m" -prof gc -rf json -rff /tmp/banded-cholesky.json .*BandedCholeskyJmh.*'
sbt 'benchmarksJVM/Jmh/run -wi 3 -i 5 -w 500ms -r 500ms -f 1 -t 1 -p n=1200 -p bandwidth=128,256 -jvmArgsAppend "-Xms256m -Xmx512m" -prof gc -rf json -rff /tmp/banded-cholesky-wide.json .*BandedCholeskyJmh.(factorPacked|solveInPlace)'
```
