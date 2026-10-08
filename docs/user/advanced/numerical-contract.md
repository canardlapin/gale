# Numerical, sparse, and backend contract

This document states what Gale v1 guarantees and where callers must make an
explicit numerical or performance choice. It complements the runnable
[worked examples](../guides/examples.md) and the focused
[Breeze migration guide](../guides/breeze-equivalence.md). The current
portable / backend-only / deferred split is tabulated in
[Shipped vs deferred](shipped-vs-deferred.md).

## Floating-point guarantees

Gale implements IEEE 754 `Double` arithmetic. It guarantees dimensions,
ordering conventions, typed failure modes, and documented mathematical
invariants; it does not guarantee exact real arithmetic or universal
bit-for-bit agreement between legal algorithms.

- The pure single-threaded implementation is deterministic for a fixed Gale
  build and runtime.
- JVM kernels may use `Math.fma`; Scala.js uses the JavaScript number operation.
  Cross-platform results can differ in their final ulps.
- Vector and vendor BLAS/LAPACK backends may reassociate operations. Their
  answers must satisfy conformance tolerances but need not match pure Gale bits.
- Solves and decompositions use scale-aware tests. There is no single absolute
  tolerance valid for every matrix; applications should assess a relative
  residual such as `||A*x-b|| / (||A||*||x|| + ||b||)` against a tolerance
  appropriate to data scale and conditioning.
- Iterative and spectral results carry convergence/residual diagnostics. A
  returned approximation is not a claim that every requested component
  converged; inspect the diagnostics or call the result's `requireConverged`
  helper to enforce residual convergence. For partial spectral extremes,
  `requireExtremeCertified` is the stricter policy: it additionally requires an
  independent membership certificate and returns the typed
  `SpectralExtremeNotCertified` error when residuals pass without one.

Dense LU uses partial pivoting. Cholesky and symmetric eigen read the lower
triangle. QR stores Householder reflectors and reports a numerical rank.
Symmetric eigenvalues are ascending; singular values are descending. Legal
pivot, reflector-sign, eigenvector-sign, and repeated-eigenspace choices are not
part of the identity contract.

### Cross-platform singularity, rank, and backend residuals

LU reports `SingularMatrix` only when a pivot is exactly `0` or `NaN`.
`conditionEstimate` maps that failure to `Right(+∞)` and rejects a rectangular
input with `Left(NonSquareMatrix)`. A matrix that is rank-deficient in exact
arithmetic but whose LU pivots stay nonzero is ill-conditioned, not singular:
`conditionEstimate` stays finite.

Scala.js (and the experimental Wasm lane) lower `PlatformMath.fma` to `a*b+c`,
not IEEE FMA. Reconstructing `UΣVᵀ` with a planted `σ=0` after Householder QR
is therefore only approximately singular on those platforms. Portable
singularity plants are IEEE-exact outer products (for example
`[[1,2,3],[2,4,6],[3,6,9]]`) or exact zero rows. Do not use a QR-reconstructed
reduced SVD as a `SingularMatrix` / `+∞` fixture.

`rankEstimate` is QR numerical rank at `2 · max(m, n) · ε · max|R_ii|`. `pinv`
zeros singular values at or below the MATLAB/SciPy cutoff
`max(m, n) · ε · σ_max`. Prefer exact diagonals and exact outer products so the
keep/drop decision is not lost in reconstruction rounding. The shared
`NumericalPolicySuite` pins these decisions on every core platform.

Caller-selected `SvdCutoff` applies consistently to `truncatedSvd`, `pinv`,
`minimumNormLeastSquares`, row/column spaces and numerical null spaces. Keep
exactly `sigma > cutoff`; the selected rank is independent of the existing
`SVD.rank` metadata, whose default relative tolerance is `1e-10`. Deliberate
truncation defines geometry and solves for the truncated matrix. Ordinary RREF
pivot tolerance has a separate contract. See
[minimum-norm solves and subspaces](../guides/subspaces-and-minimum-norm.md).

An imported Vector or native BLAS/LAPACK backend must agree with
`PureBackend` on the `LinAlgError` class for those IEEE-exact plants and must
keep solve residuals inside the documented conformance tolerance. Factor
entries, pivot indices, and reflector signs may differ. The reusable
`BackendConformanceSuite` enforces that residual and error-class agreement.

Structural failures use `Either[LinAlgError, A]` on total solve/factorization
entry points. Primitive arithmetic operators and `mulInto` methods validate
their preconditions and may throw `LinAlgError`; they are intentionally not
silently totalized.

### Partial generalized symmetric-definite operators

The operator-level partial generalized symmetric-definite solver is
`Eigen.eigSymmetricGeneralized(Aop, Bop, n, Count(...), options,
preconditioner)`. Its portable engine is LOBPCG over typed symmetric `A` and
positive-definite `B` operators. It calls only `A*X`, `B*X`, and the explicit
preconditioner; it does not materialize the operators or silently construct
`B^-1`.

The operator surface accepts
`Count(k, SmallestAlgebraic | LargestAlgebraic)`, `1 <= k < n`, plus an optional
`n × k` initial subspace. Results are ascending, vectors are
`B`-orthonormal, and diagnostics report the true generalized residual and
`B`-Gram error. Non-convergence returns the converged subset in `Right`, never a
false full result. Residual convergence does not by itself certify membership in
the requested global extreme; callers needing that stronger guarantee use
`requireExtremeCertified`. The separately named generalized block-Lanczos
engine is available when the caller supplies an explicit `MetricSolveOperator`;
it never constructs `B^-1` or silently factorizes the metric. The
[operator eigensolver guide](generalized-operator-eigen.md) gives the complete
selection, metric-solve, backend, and work-accounting contract.

### Reductions and elementwise numerics

`DVec`/`DMat` reductions (`sum`, `mean`, `max`, `min`, `argmax`, `argmin`,
per-axis forms, and the norms) and the `Numerics` functions follow these rules
on both JVM and Scala.js:

- **NaN.** Sums, means, and norms follow IEEE arithmetic, so a NaN entry gives
  NaN (and `+Inf + -Inf` gives NaN). `max`/`min` return NaN when any entry is
  NaN; `argmax`/`argmin` return the index of the **first** NaN. A matrix
  `argmax` reports the first NaN in row-major order, independent of layout.
- **Ties.** `argmax`/`argmin` return the first occurrence (row-major for a
  matrix). `-0.0` and `0.0` compare equal.
- **Empty input.** `sum` and `sumExact` of an empty vector or matrix are `0.0`;
  `norm1`, `norm2`, `normInf`, and `normFrobenius` are `0.0`. `mean`, `max`,
  `min`, `argmax`, and `argmin` throw `LinAlgError.EmptyInput`. A per-axis
  reduction throws only when it has a result entry to produce from an empty
  line (for example `max(Axis.Rows)` on a `3×0` matrix); a per-axis `sum` of
  empty lines is `0.0`. `logSumExp` of an empty input is `-Inf`.
- **Determinism of `sum`.** `sum` may reassociate (it keeps several partial
  sums and walks storage in layout order), so it is deterministic for one
  build and platform, but the final bits can differ from a left-to-right loop,
  between a matrix and its transpose view, and between JVM and Scala.js.
  `sumExact` is the exact sum rounded once to the nearest `Double`; it is
  identical on every platform, layout, and permutation of the entries.
  Per-axis reductions and the per-axis log-domain functions likewise
  accumulate in layout order (a row-major matrix reduced per column is read
  row by row), so they can differ in the last bits from the same line reduced
  as a vector.
- **Means.** `mean` is `sum / n` whenever that sum is finite, so ordinary
  inputs give exactly `sum / n`. When the sum overflows or is NaN, the mean is
  recomputed as `Σ(x_i / n)`. Entries near `Double.MaxValue` therefore give a
  finite mean (clamped to `±Double.MaxValue` if rounding at the edge of the
  range overflows), while infinite and NaN entries keep their IEEE results:
  `mean(MaxValue, MaxValue, -Inf)` is `-Inf`, not NaN. Per-axis means apply
  the same rule to each line.
- **Overflow.** `normFrobenius` and `norm2` are scaled, so entries near `1e300`
  give a finite norm and tiny entries are not lost to underflow. An infinite
  entry gives `+Inf` unless a NaN is also present.
- **Log domain.** `logSumExp` shifts by the maximum, so finite inputs never
  overflow. NaN anywhere gives NaN; otherwise any `+Inf` gives `+Inf`; all
  `-Inf` gives `-Inf`. `softmax` and `logSoftmax` return all-NaN for an input
  or line whose maximum is not finite (NaN, `+Inf`, or all `-Inf`).
- **Elementwise functions.** `Numerics.exp`, `log`, `log1p`, and `expm1` apply
  the platform `Math` function to each entry, so they agree bit for bit with
  `a.pointwise.map(math.exp)` on the same platform. `sigmoid` uses the
  overflow-free form `1/(1+t)` or `t/(1+t)` with `t = exp(-|x|)`.

## Storage and allocation boundary

`DVec` and `DMat` are immutable-facing values and views. Their owned platform
storage (`Array[Double]` on JVM, `Float64Array` on Scala.js) is private API.
Callers can export copies and can explicitly create documented views, but no v1
public signature promises `Array[Double]` as the representation.

`Vec[Double]` construction selects primitive `DVec`. Generic `Vec[A]` behavior
is correctness-oriented and is not the primitive throughput path. Ordinary
operators remain allocating and immutable-facing. Profiled repeated pipelines
can opt into the explicit allocation-control tier:

- [dense destinations](dense-destinations.md) for GEMM and fused linear
  combinations;
- [reusable primitive workspaces](workspaces.md) with checked scratch
  requirements; and
- single-owner sparse value destinations used by
  [symbolic sparse plans](sparse-plans.md).

Destinations and workspaces expose operations, not their backing storage. They
are sequential mutable resources and are not safe for concurrent use. Checked
constructors defensively copy caller arrays; only a Gale-owned builder may
transfer its storage to an immutable result.

## Sparse v1 support

Gale v1 supports:

- COO construction with `Sum`, `Last`, or `Error` duplicate policy;
- canonical CSR and CSC storage, transpose views, sparse addition/subtraction,
  scalar and value transforms, diagonal/trace, and sparse matrix-vector action;
- immutable [compressed patterns](sparse-patterns.md), independent numeric
  rebinding, and zero-copy structural transpose between CSR and CSC patterns;
- checked [symbolic union and product plans](sparse-plans.md) for repeated
  numeric replay on one exact pair of compressed patterns;
- banded, diagonal, identity, zero, and permutation structural matrices;
- dense conversion with an explicit maximum-entry guard;
- Matrix Market `coordinate real general` read/write;
- matrix-free `LinearOperator` use with CG, BiCGSTAB, restarted GMRES, CGNR, and
  LSQR, plus Jacobi preconditioning; repeated CG solves can use caller-owned
  `CgWorkspace` storage through `cgWith`.

`rebind`, `mapValues`, scalar multiplication, and symbolic-plan numeric replay
are structure-preserving: they retain every stored position, including explicit
numeric zeros. `pruneZeros`, thresholded `prune`, and `canonicalize` are explicit
structure-changing operations. Symbolic analysis creates a new fixed result
pattern; later replay does not change it. Canonicalization sorts indices,
combines duplicates, and prunes exact zeros; it does not apply a numerical
near-zero threshold unless the caller does so explicitly.

The product plan is a fixed-pattern analyze-once/replay-many facility, not a
general allocating sparse matrix-matrix multiplication facade. Gale v1 does not
claim general sparse matrix-matrix multiplication, an implemented sparse direct
LU/Cholesky/QR provider, complex sparse storage, every Matrix Market
field/symmetry, or a full Breeze sparse-collection replacement. The shared sparse-direct seam defaults to `SparseDirectProvider.none` on
JVM and Scala.js. `import gale.sparse.direct.pure.given` enables portable
sparse Cholesky only (minimum-degree `ProviderDefault`, exact-zero-style
`NotPositiveDefinite`, fill guard). Do not infer sparse LU or QR from the
types or from that import. Sparse QR, pivoting LU, and C Wasm each need
their own spec in a [future version](../../sparse-direct-future.md). The
current plan is [Scala.js sparse-direct](../../sparse-direct-js.md).

## Choosing a backend

| Need | Choice | Runtime and policy |
| --- | --- | --- |
| portability, reproducibility, small/strided work | no import (`PureBackend`) | JVM 21+ and Scala.js; default |
| faster contiguous GEMV/GEMM on JVM | `import gale.backend.jvm.vector.given` | JDK 21+ with Vector incubator module; measured thresholds |
| large native GEMM/LU on a measured library | `import gale.backend.jvm.blas.given` | JDK 22+, native access enabled; runtime discovery and family-specific thresholds |
| copy-free repeated native GEMM | explicit `NativeDMat` | JDK 22+; caller owns lifetime and conversion |
| browser/Node | ordinary Scala.js JavaScript | supported cross-platform route |
| experimental Wasm investigation | `GALE_WASM=1` build | correctness-tested but currently much slower and default-off |

Importing no backend preserves the pure behavior. Importing an accelerator does
not route every primitive: only coarse public seams with measured thresholds are
eligible. A backend that loads but lacks a proven crossover remains
direct-callable and default-off. Unknown OpenBLAS/MKL thresholds are not inferred
from Accelerate measurements.

Native thread count defaults to one. A request for multiple native threads is
accepted only when the selected library exposes a supported thread-control
symbol; otherwise loading fails rather than silently oversubscribing the JVM.

The current thresholds and rejected routes are recorded in the
[backend dashboard](https://github.com/canardlapin/gale/blob/main/benchmarks/dashboard.md).
