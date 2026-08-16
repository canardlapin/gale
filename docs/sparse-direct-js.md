# Scala.js sparse-direct

Status: **Phase 0 landed** (shared seam, `none` default). No factorization
is shipped. This document is the JS counterpart of
[the provider boundary](sparse-direct-provider.md).

Gale has a staged seam on every core platform:

`provider → workspace → symbolic analysis → numeric factor → solve`

The default provider is `SparseDirectProvider.none` and every entry point
returns `Left(UnsupportedOperation)`. Phase 0 moved the types from
`core/jvm` into `core/shared` so Scala.js resolves the same package.

The Scala.js solution is **the same contract, with a JS-capable provider**.
It is not a second facade, not an async `Promise` API, and not an iterative
solver pretending to be a factorization.

## Decision

1. **Promote the existing types and facade to `core/shared`.** After that move,
   JVM and JS resolve the same `given SparseDirectProvider = none`. Callers
   write the same staged code on both platforms. Absence of a real provider
   remains a typed `UnsupportedOperation`, not a missing classpath.
2. **Do not invent a JS-specific surface.** No `factorAsync`, no
   `js.Promise`, no `WebWorker` workspace type, no routing through
   `Backend` / `Capability.NativeSparse`. Provider selection stays an
   explicit `given SparseDirectProvider`.
3. **The first real JS provider is pure Scala**, not embedded C/Wasm. A
   portable Cholesky, then a static-sparsity LU, can live next to the dense
   factorizations, run on JVM, JavaScript, and the experimental
   `GALE_WASM=1` Scala.js-to-Wasm lane, and still keep `none` as the default
   `given`.
4. **Embedded CSparse / SuiteSparse Wasm is a later optional module**, with
   extra packaging and load-time gates. It is a different Wasm from
   `GALE_WASM=1` and must not be implied by that flag.

Until a provider advertises a capability and passes the gates below, docs
continue to say: Gale ships no sparse direct factorization. Use CG / BiCGSTAB
/ LSQR, or convert a small system to dense.

## Why the API used to be JVM-only

`gale.sparse.direct` started as a single file under `core/jvm`. Portable
sparse (CSR/CSC, patterns, symbolic union/product plans, matvec) and the
iterative solvers already cross-compiled, so JS had the *inputs* a direct
provider needs and none of the *seam*.

That split matches the native-module rule in
[backend architecture](backend-architecture.md): optional native code is
absent from the JS classpath rather than stubbed. It is the wrong split for a
**pure** provider. Phase 0 moved the staged types, facade validation,
diagnostics, and `none` default into `core/shared`.

## Phase 0 — one seam on every platform

`gale.sparse.direct` now lives in `core/shared` without enabling any
factorization.

| Change | Rule |
| --- | --- |
| Types, facade, `SparseDirectProvider.none` | Identical source on JVM and JS |
| Default `given` | Still `none`; empty capability set |
| Error copy | `no provider is installed` (drop the word `JVM`) |
| `AutoCloseable` | Keep. JS has no ARM; callers still `close()` so large factor buffers become collectible |
| Tests | Shared empty-provider lock (today's `SparseDirectParitySuite` first test, plus a core JS suite). JVM keeps the fake dense-oracle provider suite; JS does not run `java.util.concurrent` stress |
| Docs | This file plus a pointer from `sparse-direct-provider.md` |

Phase 0 is API hygiene. It does not advertise LU, Cholesky, or QR. It makes
the later provider an import, not a new package name.

```scala
import gale.sparse.direct.*

SparseDirect.capabilities // Set.empty on JVM and JS
SparseDirect.newWorkspace() // Left(UnsupportedOperation(...))
```

A future caller enables work with an explicit given, the same way Vector and
FFM backends are imported:

```scala
import gale.sparse.direct.pure.given

val workspace = SparseDirect.newWorkspace().orThrow
val symbolic = SparseDirect.analyze(a.pattern, SparseDirectFactorization.Cholesky, workspace).orThrow
val factor = SparseDirect.factor(symbolic, a, workspace).orThrow
val solved = SparseDirect.solve(factor, rhs, workspace).orThrow
```

`aWithNewValues` must still `sharesPatternStorage` with the analyzed
`CSRPattern`. That rule is already the portable sparse-plan rule; JS does not
relax it.

## What a JS provider can actually be

Three implementations are compatible with the existing trait. Only the first
is the recommended Scala.js solution.

### A. Pure Scala provider (recommended)

A left-looking sparse Cholesky, then a static-sparsity sparse LU, written
against `CSR` / `CSRPattern` and platform `DoubleArray` / `IndexArray`
(`Float64Array` / `Int32Array` on JS). No JNI, no FFM, no Emscripten.

**Why this is the JS solution**

- The public API is already staged for exactly this split: symbolic fill
  estimate and column permutation are independent of numeric values; numeric
  factor binds a value vector on the same pattern; solve is a separate
  lifetime.
- Gale's dense LU already reports `SingularMatrix` only on an exact-zero or
  NaN pivot. A static-sparsity sparse LU can use the same policy, so JS and
  JVM agree on the error class for IEEE-exact plants.
- The same source runs on JVM (useful as a portable reference oracle),
  ordinary Scala.js, and `GALE_WASM=1`. One implementation, three runtimes.
- Browser and Node have no honest SuiteSparse story. A pure provider does not
  pretend otherwise.

**Ship order**

| Step | Family | Capabilities to advertise only after gates pass |
| --- | --- | --- |
| 1 | Sparse Cholesky | `Cholesky`, `UserOrdering`, `TransposeSolve`, `MultipleRhs` |
| 2 | Static-sparsity LU | `LU` plus the same optional features |
| 3 | Sparse QR | Deferred. SPQR-class rectangular QR is a different project |

Cholesky first: SPD systems are the common browser/Node workload (FEM,
graphs, kernels, covariance), the symbolic phase is exact (no numeric
pivoting), and `NotPositiveDefinite` already exists.

**Cholesky sketch**

Symbolic, square canonical `CSRPattern`:

1. Use the existing zero-copy pattern transpose so column access is CSC.
2. Compute the elimination tree (Liu).
3. Count nonzeros per column of `L` and allocate a compressed factor pattern.
4. `Natural` is identity. `ProviderDefault` is AMD (or a documented
   minimum-degree stand-in) and must report the permutation actually used.
   `User` applies the caller permutation and requires `UserOrdering`.
5. Write `predictedFactorNnz`. If it exceeds the fill guard (below), return
   `Left(InvalidArgument)` before allocating.

Numeric:

1. Left-looking: for column `j`, assemble `A(:, j)` into work, solve
   `L(0:j, 0:j) x = A(0:j, j)`, set `L_jj = sqrt(a_jj - x·x)`.
2. Non-positive or NaN pivot → `Left(NotPositiveDefinite(j))`.
3. Factor lifetime is independent of the symbolic handle, as the JVM contract
   already requires.

Solve:

- Forward `L y = P b`, back `Lᵀ x = y` (and the transpose operation, which
  for SPD is the opposite triangle order).
- Multiple RHS is repeated column solves; no extra capability machinery.
- Residual in `SparseSolveDiagnostics` is optional; if computed, use
  `‖b - A x‖` with `A` applied through portable CSR matvec, not a dense
  conversion.

**LU sketch (after Cholesky)**

Static sparsity: the symbolic pattern of `L` and `U` is computed from the
graph of `A` (or `A` with the chosen column permutation) and does **not**
change during numeric factorization. No threshold pivoting, no delayed
pivots, no SuperLU-style numeric restructuring.

Numeric left-looking LU: exact-zero or NaN pivot → `Left(SingularMatrix(j))`.
That matches [the numerical contract](user/advanced/numerical-contract.md#cross-platform-singularity-rank-and-backend-residuals).
A reconstructed near-singular sparse matrix is ill-conditioned, not
`SingularMatrix`.

Row permutation is identity until a later pivoting provider exists. The
numeric factor still returns a Gale-owned `Permutation` copy, never an
aliased `Int32Array`.

This will fail more often than UMFPACK on poorly ordered nonsymmetric
problems. That is acceptable and must be documented: JS LU is a portable
direct solver, not a SuiteSparse replacement. SciPy `splu` remains the
*eventual* JVM-native differential target
([`SparseDirectParitySuite`](../parity/src/test/scala/gale/parity/SparseDirectParitySuite.scala)),
not the JS correctness oracle. JS LU is scored against Gale dense LU on the
same numbers, plus residual tests, plus IEEE-exact singular plants.

**QR** stays `UnsupportedOperation` until a separate spec. Do not advertise
`SparseDirectCapability.QR` from a thin “convert to dense and Householder”
fallback.

### B. Dense-fallback provider (oracle only)

The JVM test suite already factors a CSR by `toDense()` and dense LU. That
is a conformance oracle, not a sparse-direct implementation.

A shipped `SparseDirectProvider.denseFallback(maxOrder: Int)` may exist as a
named testing/debug tool if it:

- advertises `LU` / `Cholesky` only for `n ≤ maxOrder`;
- rejects larger patterns with `UnsupportedOperation` or `InvalidArgument`
  before allocating `n²`;
- is **not** the `pure` given and is **not** documented as the JS solution.

It must not be how Gale claims “sparse direct on Scala.js.”

### C. Embedded C Wasm (CSparse / SuiteSparse) — later, optional, gated

This is **not** `GALE_WASM=1`. That flag compiles Gale's Scala to Wasm. An
embedded solver is C compiled to Wasm, instantiated from JavaScript, with
its own linear memory.

It can implement `SparseDirectProvider` only if it looks like every other
provider after load:

| Constraint | Rule |
| --- | --- |
| Sync public API | `analyze` / `factor` / `solve` stay `Either[LinAlgError, _]`. Instantiation is `load(): Either[LinAlgError, SparseDirectProvider]` (or `Future` at the module edge only). No `given` is published until `load` succeeds |
| Copies | CSR values/indices are copied into the Wasm heap; permutations are copied out to Gale `Permutation`. No `Float64Array` alias escapes |
| Index width | `Int32` only, matching JS `IndexArray` |
| Threads | `BackendConfig.singleThreaded`. SharedArrayBuffer workers are out of v1 |
| Licence / ABI | Recorded before the module is advertised. SuiteSparse LGPL vs CSparse BSD is a product choice, not an implementation detail |
| Failure | Missing module, wrong ABI, or OOM → no `given`, so the default `none` remains. Never a mid-solve linkage exception |
| Packaging | Separate optional JS artifact. No Wasm binary in `gale-core` |

Until those gates have evidence, this option is **no-go**. It is also the
wrong first step: it does not run on the Scala.js-to-Wasm lane, it makes CI
depend on a binary blob, and it does not help the portable numerical-policy
story.

## JS contract amendments

The JVM document still applies. These are the JS-specific tightenings.

**Memory.** Symbolic analysis must estimate fill and refuse to allocate an
unbounded factor. A documented guard such as
`predictedFactorNnz > max(inputNnz * fillFactor, absoluteNnzCap)` returns
`Left(InvalidArgument)` with the two counts in the message. Browser and Node
heaps are the reason; the same guard is useful on JVM.

**Storage.** Factor and workspace buffers are platform `DoubleArray` /
`IndexArray`. Public signatures stay `DVec`, `MutableDVec`, `DMat`,
`DMatBuilder`, `Permutation`. No `Float64Array` in the provider trait.

**Threads.** `BackendConfig.singleThreaded` only. The JVM rule “concurrent
solves need distinct workspaces” remains, but JS tests interleave
sequentially rather than using `Future` thread pools. A workspace is still
not safe for overlapping calls (including overlapping JS microtasks).

**Lifecycle.** `close()` drops references to large typed arrays so the JS GC
can reclaim them. Idempotent. Closed workspace / analysis / factor is
rejected before provider code, same as JVM. There are no file descriptors to
leak; the resource-safety gate is “no retained buffers after close” plus
idempotence.

**Determinism.** A pure provider sets `SparseSymbolicDiagnostics.deterministic = true`.
No wall-clock fields (already forbidden).

**Ordering.** `ProviderDefault` must name its algorithm in docs (AMD or
minimum degree). It must not silently change between releases without a
capability or version note. `Natural` is bit-stable given the same pattern.

## What not to do

- Do not teach JS callers a different sequence than
  `newWorkspace → analyze → factor → solve`.
- Do not route sparse direct through `given Backend` or
  `Capability.NativeSparse`.
- Do not implement “direct” as one CG call with a tight tolerance.
- Do not plant JS singularity via a reconstructed reduced SVD; use IEEE-exact
  outer products, as in `NumericalPolicySuite`.
- Do not put Emscripten, npm `.wasm` blobs, or SuiteSparse into `gale-core`.
- Do not treat `GALE_WASM=1` as a sparse-direct accelerator. It is a Scala
  backend, and a pure provider already runs there.
- Do not advertise `QR` or pivoting LU to look complete.

## Go / no-go gates

Phase 0 (shared seam) may land without a provider. A JS-visible capability
(`Cholesky` or `LU` on `SparseDirectProvider.pure`) may be advertised only
when every applicable row passes.

| Gate | Required evidence |
| --- | --- |
| Seam | Types and facade live in `core/shared`; JS and JVM empty-provider locks pass; error text is platform-neutral |
| Capability honesty | `pure` advertises only implemented families. QR stays off until it exists. `none` remains the default `given` |
| Contract conformance | Shared provider suite: symbolic reuse on `rebind`, pattern-change rejection, closed-resource rejection, normal/transpose/multi-RHS solves, zero-sized matrices, malformed handles closed before error return |
| Numeric oracle | Cholesky reconstructs `L Lᵀ ≈ A` and solves with a small residual vs portable CSR matvec. LU matches Gale dense LU on the same numbers for well-ordered fixtures; IEEE-exact singular / non-SPD plants return the same error *class* as dense |
| Fill guard | A deliberately dense fill pattern is rejected with `InvalidArgument` rather than exhausting the JS heap |
| Resource safety | Repeated analyze/factor/solve/close on JS does not retain factor buffers after `close()`; close is idempotent |
| Cross-platform | `coreJVM` and `coreJS` run the same pure-provider tests. Optional `GALE_WASM=1` is a follow-up, not a blocker |
| Performance | Documented as a portable correctness path, not a SuiteSparse replacement. A reuse crossover (analyze once, factor many) is measured on Node; one favourable matrix is insufficient to claim acceleration |
| Packaging | Pure provider is ordinary Scala in core (or a Scala-only module). No native / Wasm binary on the default JS artifact |

Embedded C Wasm, if ever attempted, additionally inherits every packaging,
licence, ABI, and “no `given` until `load` succeeds” row from the JVM
SuiteSparse table in `sparse-direct-provider.md`.

## Recommended landing sequence

1. **Phase 0 (landed)** — `SparseDirect.scala` in `core/shared`;
   `SparseDirectSeamSuite` empty-provider lock on JVM and JS; platform-neutral
   error copy.
2. **Phase 1** — `SparseDirectProvider.pure` with Cholesky, `Natural` +
   `ProviderDefault` AMD, `UserOrdering`, transpose and multi-RHS, fill
   guard. Explicit import, not the default `given`.
3. **Phase 2** — static-sparsity LU on the same provider, exact-zero pivot
   policy, dense-LU residual oracle, IEEE-exact `SingularMatrix` plants.
4. **Not scheduled** — sparse QR; threshold / delayed pivoting LU; embedded
   CSparse Wasm; any implicit routing from `Backend`.

The portable iterative solvers remain the default sparse solve on Scala.js
until a caller imports a capable provider. Direct factorization is an opt-in
for repeated same-pattern SPD (then square nonsymmetric) systems, not a
replacement for CG on matrix-free operators.
