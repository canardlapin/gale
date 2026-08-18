# Scala.js sparse-direct

Status: **Phase 0 and Phase 1 landed** (shared seam; portable Cholesky
behind `import gale.sparse.direct.pure.given`). The default `given` is
still `none`. LU, QR, and C Wasm are not shipped. This document is the JS
counterpart of [the provider boundary](sparse-direct-provider.md).

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
| 3 | Sparse QR | Deferred. Not a third family on the same engine; see [Why later items are different projects](#why-later-items-are-different-projects) |

Cholesky first: SPD systems are the common browser/Node workload (FEM,
graphs, kernels, covariance), the symbolic phase is exact (no numeric
pivoting), and `NotPositiveDefinite` already exists.

**Cholesky sketch**

Symbolic, square canonical `CSRPattern`:

1. Use the existing zero-copy pattern transpose so column access is CSC.
2. Compute the elimination tree (Liu).
3. Count nonzeros per column of `L` and allocate a compressed factor pattern.
4. `Natural` is identity. `ProviderDefault` is minimum degree with
   lowest-index tie-break (a documented AMD stand-in) and reports the
   permutation actually used. `User` applies the caller permutation and
   requires `UserOrdering`.
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

**QR** stays `UnsupportedOperation` until a separate spec. It is not a
third family on the Cholesky/LU engine (see below). Do not advertise
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
2. **Phase 1 (landed)** — `SparseDirectProvider.pure` with Cholesky,
   `Natural` + `ProviderDefault` minimum degree (lowest-index tie-break),
   `UserOrdering`, transpose and multi-RHS, fill guard
   `max(inputNnz * 64, 8_000_000)`. Explicit import, not the default
   `given`. Evidence: `PureSparseCholeskySuite`.
3. **Phase 2** — static-sparsity LU on the same provider, exact-zero pivot
   policy, dense-LU residual oracle, IEEE-exact `SingularMatrix` plants.
4. **Future version (spec first)** — sparse QR, threshold / delayed
   pivoting LU, and embedded C Wasm. Each needs its own spec in a later
   version; they are not leftover families on the Phase 1–2 engine. See
   [the future-version note](sparse-direct-future.md). Implicit `Backend`
   routing stays out.

The portable iterative solvers remain the default sparse solve on Scala.js
until a caller imports a capable provider. Direct factorization is an opt-in
for repeated same-pattern SPD (then square nonsymmetric) systems, not a
replacement for CG on matrix-free operators.

## Why later items are different projects

Phase 1–2 share one engine and one contract:

- square `CSR` / `CSRPattern`;
- symbolic analysis computes a **fixed** factor sparsity from the graph;
- numeric factorization only writes values into that pattern (and may fail
  with `NotPositiveDefinite` or `SingularMatrix`);
- `factor(analysis, rebound)` reuses the symbolic handle because
  `sharesPatternStorage` means the graph did not change;
- Gale's exact-zero / non-positive pivot policy;
- ordinary Scala in `gale-core`, no extra artifact.

Sparse QR, pivoting LU, and C Wasm each drop a different one of those
assumptions. Putting them on the same checklist as “also implement QR”
would advertise a capability the Cholesky/LU engine cannot honestly
provide.

### Sparse QR is a different factorization, not a third enum case

The facade already allows rectangular QR and asks the factor for
`rhsRows` / `solutionRows`. That is the *seam*. The *algorithm* is not
“Cholesky of `AᵀA`” and not “dense Householder on `toDense()`.”

- **Different graph.** Sparse Cholesky/LU symbolics are the elimination
  tree of a square pattern. Sparse QR (SPQR / multifrontal Householder)
  symbolics are the column intersection graph of `A`, equivalently the
  pattern of `AᵀA`, plus a frontal assembly order. Fill, ordering (COLAMD
  vs AMD), and `predictedFactorNnz` are a different analysis.
- **Different factors.** The portable result is Householder/Givens
  reflectors plus `R`, or a Q-less least-squares apply. It is not `L`/`U`
  triangles and a triangular solve. Dense Gale already stores reflectors
  and applies `Qᵀ` without forming `Q`; a sparse provider has to invent
  the compressed analogue and new solve diagnostics.
- **Different failure policy.** Dense `leastSquares` is tall-only; wide
  systems are `Left(UnsupportedOperation("underdetermined least squares"))`.
  Rank drop is `RankDeficient` at the QR cutoff
  `2 · max(m,n) · ε · max|R_ii|`. Shipping sparse QR reopens both
  decisions (underdetermined, numerical rank) for compressed inputs. That
  is product work, not a leftover `SparseDirectCapability.QR` bit.
- **Normal equations are already shipped.** `cgnr` / `lsqr` solve sparse
  least squares without a QR factor. Forming `AᵀA` and calling Phase 1
  Cholesky would square the condition number and must not be advertised
  as `QR`.
- **Size.** SuiteSparseQR is a library. A honest pure SPQR is a
  standalone implementation project (frontals, rank-revealing, apply-Q).
  A `toDense()` Householder behind `SparseDirectCapability.QR` fails the
  capability-honesty gate.

So QR waits for its own spec: rectangular contract, rank policy, and an
engine. It does not ride along once LU exists.

### Threshold / delayed pivoting is a different LU contract

Phase 2 LU is **static sparsity**: the `L`/`U` pattern is fixed at
analyze time; numeric left-looking only fills it; a zero/NaN pivot is
`SingularMatrix`. That matches dense Gale's exact-zero policy and the
`analyze` → `rebind` → `factor` reuse rule.

Threshold pivoting (UMFPACK / SuperLU / SciPy `splu`) and delayed
pivoting change that contract:

- **The pattern is no longer a function of the graph.** A row swap during
  numeric factorization changes the sparsity of `L`. The symbolic handle
  cannot be applied to new values without either (a) accepting a worse
  pivot and possible `SingularMatrix`, or (b) re-analyzing. (b) voids
  “changing values on one exact pattern reuses symbolic state.”
- **The pivot policy is no longer exact-zero.** Threshold pivoting treats
  a small-but-nonzero `u_ii` as a pivot failure or a swap candidate.
  That is a different `LinAlgError` story from
  `NumericalPolicySuite` / dense LU. A JS provider that used a threshold
  would disagree with JVM dense LU on the IEEE-exact plants the trust
  work just locked.
- **The data structure is different.** Delayed pivoting and supernodal /
  multifrontal assembly are the body of SuperLU/UMFPACK, not a flag on
  left-looking LU. SciPy `splu` is the *eventual JVM-native* differential
  target, which implies a SuiteSparse/SuperLU wrapper project, not a
  tweak to `SparseDirectProvider.pure`.

A later pivoting provider is welcome. It should be a named capability or
a distinct provider (`pure` vs `umfpack`), with its own reuse and error
rules written down first. It is not “Phase 2 plus a threshold.”

### C Wasm is a packaging project, not a solver family

Embedded CSparse/SuiteSparse Wasm does not add a factorization to the
pure engine. It is a **load, ABI, memory, and licence** project that
happens to implement `SparseDirectProvider` after `load()` succeeds.

- It is a different Wasm from `GALE_WASM=1` (C heap vs Scala compiled to
  Wasm). The pure provider already runs on the Scala.js-to-Wasm lane; a
  `.wasm` blob does not.
- The public Gale API is synchronous `Either`. Instantiation is async.
  The work is `load(): Either[LinAlgError, SparseDirectProvider]` plus
  copy-in/copy-out of CSR and permutations, not Cholesky math.
- CI, artefact size, LGPL vs BSD, `Int32` index width, and “no `given`
  until load succeeds” are the go/no-go rows. Those do not appear in
  Phase 1–2.
- Wrapping C could later supply QR or pivoting LU. That still does not
  make QR or pivoting “part of the pure JS provider.” It makes them
  features of a second, optional module.

Treat C Wasm as `gale-backend-js-csparse` (or similar) if it ever
clears the gates. Do not block portable Cholesky on it, and do not
pretend it is how Scala.js grows `SparseDirectCapability.QR`.
