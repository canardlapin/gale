# Breeze equivalence: the supported replacement slice

Gale is a replacement for a focused part of Breeze, not a source-compatible fork
and not a scientific-computing umbrella. The supported claim is:

> For Scala 3 programs whose numerical core is dense real `Double` linear algebra,
> sparse matrix-vector products, and selected spectral operations, Gale provides
> comparable numerical results, a typed failure model, cross-platform JVM/Scala.js
> execution, and competitive pure-JVM performance.

This document fixes the boundary of that claim. It is both a migration guide and
an acceptance checklist; functionality outside the table is not implied by the
word "equivalence."

## Capability matrix

| Breeze use | Gale equivalent | Evidence | Status / difference |
| --- | --- | --- | --- |
| `DenseVector[Double]` | `DVec` | dense parity and conversion suites | Equivalent numerical values; Gale is immutable-facing. |
| `DenseMatrix[Double]` | `DMat` | dense parity and conversion suites | Equivalent logical values and transpose views; storage layout is deliberately private. |
| vector add/subtract/scale, `dot`, norm, AXPY | `+`, `-`, scalar `*`, `dot`, `norm2`, `axpyInPlace` | paired parity + JMH | Numerically equivalent shared subset. |
| matrix add/subtract, `A * x`, `A * B`, transpose products | same operators | paired rectangular/transpose parity | Numerically equivalent shared subset. |
| `A \ b` square solve | `A.solve(b): Either[...]` | Breeze differential solve tests | Equivalent result; Gale exposes typed failure instead of throwing by default. |
| `A \ B` several right-hand sides | reuse `A.lu`, or `BreezeMigration.solve(A, B)` | interop differential tests | Equivalent result; migration shim factors `A` once. |
| overdetermined `A \ b` / `A \ B` | `A.leastSquares`, or migration shims | unblocked + blocked differential tests | Full-column-rank tall systems; rank deficiency is explicit. |
| `det`, LU, Cholesky, QR | typed `LU`, `Cholesky`, `QR` results | factor reconstruction and differential tests | Equivalent mathematical factors/solutions; signs and pivots may differ. |
| `rank` | `rankEstimate` | full-rank and clearly deficient differential cases | Same result on well-separated rank decisions; algorithms and tolerance rules differ. |
| `cond` | `conditionEstimate` | exact overlap tests on diagonal matrices | **Not the same general quantity:** Breeze uses exact SVD 2-norm; Gale estimates the 1-norm condition number. |
| `eigSym` | `Eigen.eigSymmetric` | eigenvalue, sign-aware vector, and repeated-eigenspace parity | Equivalent dense symmetric result and ascending order. Native Gale reads the lower triangle; the migration shim first applies Breeze's symmetry guard. |
| `svd` extremes | `Svds.svd` partial largest/smallest | tall + wide differential tests | Partial top-k/bottom-k without forming the full decomposition; values to ~1e-7 on the differential fixtures. |
| `svd` full | `A.svd` / `Svds.svd(A, SingularSelection.All)` | full-spectrum differential tests + core reconstruction/orthonormality suites | Economy shapes only (`U` m×k, `Vᵀ` k×n, k = min(m, n)) versus Breeze's full square factors; values agree to ~1e-9, vectors up to sign for well-separated singular values. |
| `pinv` | `A.pinv: Either[...]` | Moore–Penrose condition suite + full-rank differential tests | Elementwise agreement ~1e-8 on full-rank fixtures. Cutoff is the MATLAB/SciPy `pinv` convention `max(m, n)·eps·sigma_max` (NumPy's `pinv` default is a fixed `rcond=1e-15`), which is not Breeze's policy: near the rank boundary the two libraries may truncate differently. |
| `kron` | `A.kron(B)` | rectangular differential tests | Same values to ~1e-12 (plain products; no algorithmic freedom). |
| `CSCMatrix * vector` | `CSR`, `CSC`, `Banded`, `Diagonal` matvec | rectangular transpose differential tests | Equivalent stored-matrix action; formats and canonicalization APIs differ. |
| dense/sparse conversion | `gale.interop.breeze.*` | bit-exact special-value and aliasing tests | Gale to Breeze always copies. Breeze to Gale may be an explicit zero-copy view when strides are positive. |

The parity suite compares mathematical invariants when raw factors are not unique:
LU reconstruction with pivots, `Q R = A`, `Q.t Q = I`, `R.t R = A.t A`,
Cholesky reconstruction, solve residuals, eigenvector signs, and repeated-eigenvalue
subspace projectors. This avoids false failures and vacuous agreement on arbitrary
factor conventions.

## Migration boundary

Add `gale-interop-breeze` while moving a codebase one call site at a time:

```scala
import gale.interop.breeze.*
import gale.interop.breeze.BreezeMigration

val galeView = fromBreezeView(breezeMatrix) // O(1), aliases Breeze storage
val galeCopy = fromBreezeCopy(breezeMatrix) // independent row-major Gale value

val x  = BreezeMigration.solve(a, b)
val xs = BreezeMigration.solve(a, manyRightHandSides)
val ls = BreezeMigration.leastSquares(tallA, manyObservations)
```

The migration shims keep Breeze types at the boundary and throw `LinAlgError` on
failure. They do not reproduce Breeze exception classes. Once a call site is
ported, prefer Gale values and the native `Either[LinAlgError, A]` API so copies
leave the hot path.

`fromBreezeView` is intentionally loud about aliasing. A positive-stride Breeze
matrix or vector can share storage with Gale. A negative-stride vector must use
`fromBreezeCopy`. Gale-to-Breeze conversion is always a copy because Gale never
exposes its platform storage as public API.

## Performance claim

The release-grade baseline is
[`benchmarks/results/2026-07-11-breeze-release-grade.md`](../benchmarks/results/2026-07-11-breeze-release-grade.md):
two JMH forks, identical seeded inputs, and Breeze 2.1's pure-Java netlib fallback.
Gale was faster on 23 of 42 operation/size pairs, within the documented parity
band on 2, and behind on 17. Its strongest results were matrix multiplication,
`A.t * A`, blocked QR, medium/large least squares, and stable 2-norm. Its remaining
pure-kernel gaps were AXPY, medium/large GEMV, LU/solve, and the largest symmetric
eigen case. That score is runtime-conditional: enabling the incubating Vector API
on JDK 22 also enables Breeze's VectorBLAS implementation, so it must not be read
as a universal cross-JDK ranking.

The current same-process JDK 22 control is
[`benchmarks/results/2026-07-17-breeze-jdk22-vector-enabled.md`](../benchmarks/results/2026-07-17-breeze-jdk22-vector-enabled.md).
On this two-double-lane ARM runtime, Gale's optional Vector GEMM is 1.43–1.67x
faster than pure Gale. It beats Breeze VectorBLAS at `n=128` and is within 7.4%
at `n=256` and 4.0% at `n=512`. Vector GEMV is 1.77–1.98x faster than pure Gale,
but remains 17–29% slower than Breeze VectorBLAS. These are competitive results,
not a blanket claim that Gale is faster than Breeze.

The optional JDK 22 FFM backend materially accelerates selected large operations
on the measured Accelerate route, including heap-copy cost: the authoritative
two-fork result is 2.84x pure Gale for GEMM at `n=512`, and native-LU solve is
1.53–1.60x at `n=128–256`. Accelerate was sharply non-monotone at `n=256`, so
automatic GEMM begins at `512^3`; native QR, Cholesky, and GEMV remain default-off.
This is **not** evidence against Breeze configured with native BLAS/LAPACK. That
claim would require a direct same-library comparison.

The JDK Vector backend dispatches contiguous GEMV/GEMM only from the measured
128-square crossover. Transposed or otherwise strided operands retain the pure
implementation. Thresholds are performance policy, not capability detection, and
must be re-measured for materially different JVM/CPU families.

## Explicit non-equivalence

Gale does not claim replacement coverage for:

- Breeze collections, generic scalar/operator machinery, broadcasting, slicing
  syntax, negative-stride Gale views, or source/binary compatibility;
- probability distributions, optimization, statistics, signal processing,
  plotting, machine learning, tensors, or other Breeze modules;
- general complex matrix storage and arithmetic;
- sparse direct factorization or a complete native LAPACK implementation;
- identical exception classes, pivot choices, reflector signs, eigenvector signs,
  or bit-identical results when legal algorithms reassociate floating-point work;
- performance superiority over a native-enabled Breeze installation until a
  controlled same-library benchmark demonstrates it.

## Acceptance gates

The limited equivalence claim is maintained by these gates:

```sh
sbt parityTest                 # Breeze differential/invariant suite
sbt interopBreezeTest          # conversions, migration, aliasing, failures
sbt testAll                    # shared JVM + Scala.js semantics
sbt benchCompile               # paired JMH stays buildable
```

Optional backends add their own conformance suites. A backend is not allowed to
claim a capability merely because it loads: it must pass independent arithmetic
oracles, layout/stride cases, facade routing tests, and a measured crossover that
includes conversion cost.
