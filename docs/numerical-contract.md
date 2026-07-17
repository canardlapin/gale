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
  helper.

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
is correctness-oriented and is not the primitive throughput path. Destination
APIs (`mulInto`, mutable builders, reusable factorization workspaces) exist where
allocation control is part of the public contract.

## Sparse v1 support

Gale v1 supports:

- COO construction with `Sum`, `Last`, or `Error` duplicate policy;
- canonical CSR and CSC storage, transpose views, sparse addition/subtraction,
  scalar and value transforms, diagonal/trace, and sparse matrix-vector action;
- banded, diagonal, identity, zero, and permutation structural matrices;
- dense conversion with an explicit maximum-entry guard;
- Matrix Market `coordinate real general` read/write;
- matrix-free `LinearOperator` use with CG, BiCGSTAB, restarted GMRES, CGNR, and
  LSQR, plus Jacobi preconditioning.

Gale v1 does not claim general sparse matrix-matrix multiplication, sparse direct
LU/Cholesky/QR, complex sparse storage, every Matrix Market field/symmetry, or a
full Breeze sparse-collection replacement. Canonicalization sorts indices,
combines duplicates, and prunes exact zeros; it does not apply a numerical
near-zero threshold unless the caller does so explicitly.

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
