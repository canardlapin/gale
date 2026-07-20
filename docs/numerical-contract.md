# Numerical, sparse, and backend contract

This document states what Gale v1 guarantees and where callers must make an
explicit numerical or performance choice. It complements the runnable
[worked examples](examples.md) and the focused
[Breeze migration guide](breeze-equivalence.md).

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

Structural failures use `Either[LinAlgError, A]` on total solve/factorization
entry points. Primitive arithmetic operators and `mulInto` methods validate
their preconditions and may throw `LinAlgError`; they are intentionally not
silently totalized.

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
  LSQR, plus Jacobi preconditioning.

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
field/symmetry, or a full Breeze sparse-collection replacement. The JVM-only
[sparse-direct provider boundary](sparse-direct-provider.md) advertises no
capability until an optional implementation passes its documented gates.

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
[backend dashboard](../benchmarks/dashboard.md).
