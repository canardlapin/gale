# Backend Architecture (gale v0.5 acceleration)

## Purpose

This document is the design phase for four v0.5 acceleration beads, one per
section. Each section is self-contained enough to close its bead, and all four
compose into a single acceleration story: a general kernel/factorization backend
(§A), the threading policy every backend obeys (§B), the off-heap dense storage
the FFM backend needs (§C), and the staged plan for gale's own pure `gemm` up to
the point native takes over (§D).

It **composes with, and must not contradict, `docs/spectral-backend-boundary.md`**
(committed; skeleton implemented at `SpectralBackend.scala`). That doc solved the
*coarse, rare-call* boundary — a whole eigendecomposition behind a
`given SpectralBackend`. This doc solves the *general* boundary, whose hard new
problem is the **hot path**: `dot`/`axpy`/`gemm` are called from tight loops and
28 static call sites, so the seam that lets a backend accelerate them must add
**zero overhead when no backend is imported**. Where this doc reuses the spectral
doc's reviewed mechanisms (`given` resolution with an always-present fallback, a
`compose` combinator, the thread-safety requirement, conformance-by-laws), it
cross-references rather than re-deciding.

Design only; no code, one document.

## Status and provenance

| § | Bead | Scope |
|---|---|---|
| A | `gale-v0-5-backend-contract` | The general `Backend` capability contract + kernel dispatch seam |
| B | `gale-v0-5-threading-policy` | `BackendConfig`, pool ownership, oversubscription rules |
| C | `gale-v0-5-native-dmat` | `NativeDMat` over FFM `MemorySegment` |
| D | `gale-v0-5-blocked-gemm-roadmap` | Pure `gemm` stages → native crossover (feeds `gale-v0-5-native-thresholds`) |

Downstream/sibling beads referenced: `gale-v0-5-blas-ffm` (the FFM BLAS/LAPACK
backend that implements §A/§C), `gale-v0-5-vector-backend` (the JDK Vector API
backend, the key §D pure lever), `gale-v0-5-native-thresholds` (fed by §D),
`gale-epic-v0-5-acceleration` (the epic), and the scenario-benchmark bead (§D
measurement). Verified against the tree at commit `8333bde`: **none** of
`Backend`, `Capability`, `DenseDoubleKernel`, `BackendConfig`, `NativeDMat`, or a
`Layout` type exist yet — the PRD § Backend Requirements / § Backend Performance
Strategy blocks are sketches this document turns into a contract.

## Shared constraints (inherited; not relitigated)

1. **`gale-core` stays pure** — no native, no FFM, no JVM-only dep in core (PRD
   Doctrine 4). Every native implementation is a separate optional JVM module
   under the `gale.backend.*` convention.
2. **Backend selection is explicit through `given` imports** (PRD Doctrine 8).
3. **Pure JVM/JS is a serious backend, not a fallback** (PRD § Backend
   Performance Strategy). Native wins only where "architecture-specific blocking,
   packing, SIMD microkernels, prefetching, and threading dominate" — large dense
   `gemm`, blocked factorization, SVD/eigen. Level-1/2 memory-bound kernels stay
   pure.
4. **Thresholded, benchmark-derived dispatch** — "Calling native code for vectors
   of length three is a bug, not an optimization." No hard-coded guesses.

---

## §A. The general backend contract (`gale-v0-5-backend-contract`)

### A.1 The problem the spectral boundary did not have

The spectral boundary attaches to *coarse* operations: one `given SpectralBackend`
call decides a whole decomposition, so per-call `given` resolution is free
relative to the O(n³) work. The general kernel boundary is different in kind.
`DoubleKernels` is a `private[gale] object` invoked from **28 static call sites**
(Matrix.scala, Factorizations.scala, Vec.scala, MutableVec.scala,
TriangularSolve.scala), many in tight numeric code. Those sites are static,
monomorphic, and JIT-inlined today — that is *why* the pure L1/L2 kernels are
competitive. Any seam that threads a `given Backend` indirection through
`ddot`/`daxpy` would put a virtual call and a capability check on the hottest
path in the library, for an operation (native L1) the PRD explicitly says pure
JVM should win. **The seam must therefore be placed where native can help and
kept out of where it cannot.**

### A.2 Decision — dispatch by call granularity: coarse entry points seam, hot inner loops stay static-pure

The 28 call sites split on **call granularity, not BLAS level** — a single call
whose own work is ≥ O(n²) can absorb one `given` resolution; a kernel invoked from
a tight loop cannot. The boundary is drawn there:

| Class | Sites | Backend seam? |
|---|---|---|
| **Coarse entry** — `dgemm`, `dsyrkRowMajor`, standalone `dgemv` (`DMat.*(DVec)`), and the factorization drivers `lu`/`qr`/`cholesky`/`solve`/`leastSquares` | `DMat.*`, `Factorizations.*` entry points | **Yes** — few, high-level, ≥ O(n²) operations where native/SIMD can dominate; one resolution amortizes over the whole op. |
| **Hot inner loop** — `ddot`, `daxpy`, `dscal`, `dnrm2`, `dcopy`, `dadd`, `dsub`, in-loop `dtrsv`, and gemv/dot *as microkernel steps* | Vec, MutableVec, factorization scalar inner loops | **No** — static `DoubleKernels.*`, forever. Memory-bound and loop-called; a `given` resolution (or a per-element native call) is pure overhead where pure JVM already wins. |

The public L1 `Vec` operations (`dot`, `+`, `-`, `*`, `norm2`) also stay static:
each is a single call but memory-bound, and holding the "L1 stays pure" line
crisply is worth more than a bandwidth-limited SIMD sliver on them.

The seam is at the **facade**, not inside `DoubleKernels`. `DoubleKernels` stays
the pure, portable, always-present kernel set (and the JS path). A coarse facade
method resolves a `given Backend` and dispatches:

```scala
// DMat.*(that) after the seam lands — illustrative:
def *(that: DMat)(using backend: Backend): DMat =
  if backend.acceleratesGemm                   // {NativeBlas ∪ Vectorized} — NOT NativeBlas alone
     && flops(rows, that.cols, cols) >= backend.thresholds.nativeGemmMinFlops
  then backend.denseDouble.gemm(/* … */)     // native OR vectorized path
  else DoubleKernels.dgemm(/* … */)          // today's tuned pure path, byte-identical
```

The predicate is `acceleratesGemm` (`capabilities` contains `NativeBlas` **or**
`Vectorized`), not `contains(NativeBlas)` — the `Vectorized` backend's whole lever
is a SIMD `gemm` (A.2.1, §D.2), so a native-only test would silently drop it to the
pure path. With the default `given Backend = Backend.pure` (§A.4),
`acceleratesGemm` is `false` (`pure`'s capabilities are `{WasmFriendly,
Deterministic}` — no gemm accelerator), so `A * B` with no import calls
`DoubleKernels.dgemm` exactly as today. The guarantee that this costs nothing does
**not** rest on constant-folding a Set lookup — a shared `Backend` call site can be
bi-/megamorphic once a real backend coexists with `pure`. It rests on three
things: the check runs only on coarse ops (a few ns against a µs–ms gemm), the
argument is an enum singleton carrying a primitive `Long` threshold (no
allocation/boxing on the path), and **A-R1 is a measured REQUIREMENT** — the
no-backend dispatch MUST be shown allocation-free and within JMH noise of the
direct static call (the PRD's "backend dispatch overhead" benchmark). If a
`using`-parameter seam cannot meet that bar, the fallback is to specialize
`Backend.pure` as a `final` singleton the JIT devirtualizes; the L1/L2 non-seam
backstops it regardless.

Why per-call `(using Backend)` on coarse ops rather than a module-level install:
a mutable global `var backend` (install-on-import) is the Registry the spectral
doc already rejected — order-dependent, hostile to test isolation and to Doctrine
8. Per-call resolution on an O(n³) op is free; on L1/L2 it does not happen at all.

### A.2.1 Why the Vector-API backend does not seam L1/L2 (a design question, resolved)

A natural objection: the JDK Vector API's canonical use is exactly SIMD
`dot`/`axpy`, and it is *pure JVM* (no FFM), so it escapes the "native L1 call is
overhead" argument. Should a `Vectorized` backend therefore accelerate L1/L2?
Reaching a kernel called from a tight inner loop *transparently* would require an
**install-once kernel table inside `DoubleKernels`** (a `var activeKernels:
DenseDoubleKernel` swapped on import) — the only mechanism that reaches an
inner-loop kernel without threading a `given` into the loop.

That mechanism is **rejected**, for two reasons that outrank the upside:

1. **Doctrine 8 / A.8.2.** An install-once global is order-dependent, hostile to
   test isolation, and is precisely the mutable Registry the spectral doc and A-2
   already rejected. The only acceptable form would be a `given`-derived,
   lexically-scoped kernel set — never a process-global `var`.
2. **The upside is bandwidth-bound.** L1/L2 are memory-bound, so SIMD's win there
   is a small, bandwidth-limited fraction (the measured axpy gap is ~0.63x against
   *pure* Breeze — a store-bandwidth ceiling, not a compute one), not the 4–8x it
   delivers on compute-bound `gemm`. It does not justify touching the hot path or
   introducing the global.

So the `Vectorized` backend's lever is **`gemm` (and `syrk`) at the coarse seam**.
The factorization drivers reach that same gemm through their `(using Backend)`
(A.7) — but **only where a driver is itself blocked** (its trailing-submatrix
update is a `gemm`). Today that is **QR** (`factorBlockedQR` calls
`DoubleKernels.dgemm`); **LU and Cholesky are currently unblocked** (right-looking,
L1/L2 panel updates — `dtrsv`/`dcopy`/`dnrm2`, no internal `gemm`), which is exactly
why the release-grade sweep shows QR ahead of pure Breeze while Solve/LU trail. So a
`Vectorized` backend accelerates standalone `A*B` and blocked-QR's internal gemm
today; it accelerates LU/Cholesky internals only once those drivers are blocked (a
§D roadmap item), not by importing a backend. Standalone `gemv` is coarse and *may*
be routed (the Vector path copies nothing), guarded by a measured `nativeGemvMinWork`.
Loop-called L1/L2 stays pure-scalar. If a scenario benchmark later shows a compelling
L1/L2 SIMD win, that reopens as a separate, measured decision — with a scoped kernel
set, not a global.

**The two-level factorization model (how a `(using Backend)` driver dispatches).**
A driver has two independent ways to use a backend:

- **Whole-op replacement** — a `NativeLapack` backend's `denseDouble.luFactor` /
  `choleskyFactor` / `qrFactor` replaces the entire factorization (LAPACK
  `getrf`/`potrf`/`geqrf`). The pure hooks default to `PureFallback.*`, so a backend
  that does *not* advertise `NativeLapack` leaves them pure.
- **Internal-gemm routing** — a *pure but blocked* driver keeps its own algorithm and
  routes only its trailing-update `gemm` through the `(using Backend)` facade (the
  same `acceleratesGemm` test as `DMat.*`), so a `Vectorized` backend with no
  `luFactor` at all still accelerates the blocked driver's dominant flops.

The driver takes whole-op replacement when the backend advertises the matching
`Native*` factor capability and the size clears `nativeFactorizationMinSize`;
otherwise it runs the pure blocked path with backend-routed internal gemm. A
`Vectorized`-only backend never takes the first branch (it has no `luFactor`), only
the second, and only for a blocked driver.

### A.3 The contract types

```scala
package gale.backend   // the trait + enums live in gale-core (pure); impls in gale.backend.jvm.*

enum Capability:
  case NativeBlas, NativeLapack, NativeSparse     // native providers
  case Vectorized                                  // JDK Vector API SIMD kernels
  case Multithreaded                               // uses a compute pool (see §B)
  case WasmFriendly, Deterministic                 // pure-path properties

/** The accelerable dense-Double kernel surface (PRD § Backend Requirements). The
  * FULL surface is declared so a native-BLAS backend can implement all of it, but
  * the DISPATCH POLICY (A.2, thresholds) decides what is actually ROUTED: only the
  * coarse ops — `gemm`, `syrk`, standalone `gemv`, and the factorization drivers —
  * dispatch here. Loop-called L1/L2 (`dot`/`axpy`/… as inner-loop steps) never
  * route through a backend (A.2, A.2.1); those methods exist for completeness and
  * for the coarse standalone paths, not to seam the hot loops. All methods assume
  * validated dimensions (unsafe).
  */
trait DenseDoubleKernel:
  def dot(n: Int, x: DoubleArray, xOff: Int, xStr: Int, y: DoubleArray, yOff: Int, yStr: Int): Double
  def axpy(/* … */): Unit
  def scal(/* … */): Unit
  def gemv(/* … */): Unit
  def gemm(/* … */): Unit
  def syrk(/* … */): Unit                            // a.t*a; pure keeps the half-flops dsyrkRowMajor
  def trsm(/* … */): Either[LinAlgError, Unit]
  // WHOLE-OP factorization hooks (two-level model, A.2.1): only a NativeLapack backend
  // fills these; every other backend — incl. Vectorized — leaves the default, so the
  // pure blocked driver runs and routes only its internal gemm through the seam.
  def luFactor(/* … */): Either[LinAlgError, Unit]   = PureFallback.luFactor(/* … */)
  def choleskyFactor(/* … */): Either[LinAlgError, Unit] = PureFallback.choleskyFactor(/* … */)
  def qrFactor(/* … */): Either[LinAlgError, Unit]   = PureFallback.qrFactor(/* … */)

trait BackendThresholds:                            // measured, per backend family (§D, PRD)
  def nativeGemmMinFlops: Long                       // the primary crossover (§D.4); two variants, heap vs NativeDMat
  def nativeGemvMinWork: Long                        // standalone `A*x` only; conservative/∞ when copy cost dominates
  def nativeFactorizationMinSize: Int                // per-routine (LU/Cholesky/QR differ — see Open questions)
  // No `nativeDotMinLength`: loop-called L1 is never routed (A.2.1), so it has no threshold.

trait Backend:
  def name: String
  def capabilities: Set[Capability]
  def denseDouble: DenseDoubleKernel
  def thresholds: BackendThresholds
  def config: BackendConfig                          // threading (§B); read at construction
  def spectral: Option[SpectralBackend] = None       // the bridge (A.5)
  /** Does this backend provide an accelerated `gemm`/`syrk`? The coarse gemm seam
    * (A.2) routes on THIS, not on `NativeBlas` alone, so a `Vectorized` backend
    * triggers the gemm path (finding that a native-only test silently dropped it). */
  final def acceleratesGemm: Boolean =
    capabilities.contains(Capability.NativeBlas) || capabilities.contains(Capability.Vectorized)
```

**Capability discovery is a `Set`, consistent with `SpectralCapability`** (spectral
doc §1.1): the same "Set for the yes/no question, total methods for the work"
idiom, so the two contracts read as one family. `Capability` classifies
*acceleration class* (does this backend make `gemm` fast?); `SpectralCapability`
classifies *operations* (can it do QZ?). They stay distinct enums — a point the
spectral doc already fixed.

### A.4 Resolution: `given` + `Backend.pure` fallback + `compose` — mirroring spectral

Reused verbatim from the spectral doc's reviewed solution (§2.1–§2.2):

```scala
object Backend:
  /** Always-in-scope, lowest priority. Unlike SpectralBackend.none (which returns
    * Left(Unsupported)), `pure` COMPUTES — its denseDouble forwards to DoubleKernels
    * and its capabilities are {WasmFriendly, Deterministic}, no native. With no
    * acceleration import this is the only `given Backend`, so every coarse facade
    * call takes the pure branch — today's behaviour, byte-identical. */
  given pure: Backend = PureBackend

  /** Union capabilities; dispatch each op to the first part that has the needed
    * capability (earlier parts win). Same shape and rationale as
    * SpectralBackend.compose — two `given Backend`s in scope are a compile-time
    * ambiguity, so a multi-backend user imports the VALUES and declares one
    * composite given. */
  def compose(primary: Backend, rest: Backend*): Backend
```

The naming difference is deliberate and load-bearing: the spectral fallback is
**`none`** (its job is to say "unsupported"); the kernel fallback is **`pure`**
(its job is to *compute* with the portable kernels). Both are the always-present,
import-overridable default; both use `compose` for the multi-provider case; both
carry the **thread-safety REQUIREMENT (A-R2 = spectral G1): a `given Backend` is a
shared singleton and MUST be safe for concurrent invocation.**

### A.5 The `Backend.spectral` bridge — how a general backend exposes a `SpectralBackend`

The spectral doc (§2.3) promised composition "without nesting": a spectral op
resolves `summon[SpectralBackend]`, and `Backend` *may* expose
`def spectral: Option[SpectralBackend]` so one module and one import satisfy both.
This doc fixes that bridge:

- The FFM module (`gale-backend-jvm-blas-ffm`, or a `-lapack` sibling) provides
  **both** a `given Backend` (kernels/factorizations) and a `given SpectralBackend`
  (QZ/GSVD/dense-accel), and its `Backend.spectral` returns `Some(thatSpectral)`.
- Resolution stays independent: kernel facades `summon[Backend]`, spectral facades
  `summon[SpectralBackend]`. The bridge is a *convenience* for code holding a
  `Backend` that wants its spectral half, and a *self-imposed advertising contract*:
  **a backend MUST NOT advertise `NativeLapack` unless it also provides `spectral`
  (A-R3)**. This is a rule on the module author, not an automatic implication —
  LAPACK-the-library has eig/SVD, but a module could bind only `getrf`/`potrf`/
  `geqrf` without `geev`/`gesdd`/`ggev`; A-R3 forbids advertising the capability in
  that case, so a `NativeLapack` claim and the `spectral` half never disagree.
  Verified by the conformance suite.
- `compose` composes the bridge too: the composite's `spectral` is the `compose`
  of the parts' spectral providers (or `None` if none present).

### A.6 Conformance — reuse the laws, extend to kernels

Consistent with the spectral doc (§2.6), conformance is **by law, not by
bit-identity**. A `Backend`'s `denseDouble` conforms iff, over the shared
ScalaCheck matrix generators, its results satisfy the existing `MatrixLaws` /
`SolverLaws` bundles (associativity of `gemm` against the reference, factorization
reconstruction `‖PA − LU‖`, `‖A − QR‖`, `‖A − LLᵀ‖`, triangular-solve residuals)
within a documented tolerance. Kernels **may reassociate** (the FMA/vectorized
paths already do vs older gale), so equality is law-equivalence, and
**cross-platform determinism is per-platform** (the stance already committed in
`110e041` after the FMA kernels). A `BackendConformanceSuite` in `gale-laws`
(the reusable testkit) drives every advertised `Capability` and asserts the laws;
this is the same kit the spectral `SpectralBackendConformanceSuite` extends.

### A.7 Migration sketch

| Change | Class | No-import behaviour |
|---|---|---|
| Add `gale.backend` package to core: `Capability`, `Backend`(+`pure`,`compose`), `DenseDoubleKernel`, `BackendThresholds`, `BackendConfig` | additive | n/a |
| `DMat.*(that)` and `a.t*a` (syrk) gain `(using Backend)` | source-evolve | `Backend.pure` → `DoubleKernels.dgemm`/`dsyrkRowMajor`, byte-identical (A-R1) |
| Factorization drivers `lu`/`qr`/`cholesky`/`solve`/`leastSquares` gain `(using Backend)` | source-evolve | `Backend.pure` → today's pure factorization |
| **All L1/L2 sites** | untouched | static `DoubleKernels.*`, zero change |
| `gale-backend-jvm-blas-ffm`: `given Backend` + `given SpectralBackend`, thresholds, FFM kernels | new module | absent unless imported (and absent on JS) |

Pre-1.0, adding `(using Backend = pure)` to a shipped coarse method is acceptable
source evolution (same stance as the spectral migration); the invariant is
behavioural identity of the no-import path (A-R1).

### A.8 Non-goals (§A)

1. No backend seam inside L1/L2 kernels or inside `DoubleKernels` (A.2).
2. No mutable global backend install / service-loader / classpath scan (Doctrine 8).
3. No lazy expression graph — fusion stays BLAS-shaped (`axpy`/`gemv`/`gemm`/one
   pointwise pass), per PRD.
4. No backend-specific tuning knobs on the pure facade signatures (thresholds live
   on the backend, not on `DMat.*`).

---

## §B. Threading policy (`gale-v0-5-threading-policy`)

### B.1 `BackendConfig` and who owns the pool

```scala
final case class BackendConfig(
    jvmThreads: Int,                     // size of gale's ONE fixed compute pool; 1 = single-threaded
    nativeThreads: Int,                  // handed to the native library (openblas/MKL/OMP)
    allowNestedParallelism: Boolean = false
)
object BackendConfig:
  val singleThreaded: BackendConfig = BackendConfig(1, 1)          // the DEFAULT
  def dedicated(cores: Int): BackendConfig = BackendConfig(cores, cores)
```

Decisions:

- **gale owns exactly one fixed-size compute pool**, sized `jvmThreads`, created
  at backend construction, **not virtual threads** (PRD: "CPU-bound kernels should
  not use virtual threads"). The pool is owned by the `Backend` instance (or a
  process-global pool the backend references) and shut down when the backend is
  closed. Concurrent facade calls submit to it; this is safe under the A-R2
  thread-safety requirement.
- **Config reaches backends at construction, never per call** — consistent with
  the spectral doc's G2 deferral (spectral facades take no per-call threading
  arg). A `given Backend`/`given SpectralBackend` is built *with* a `BackendConfig`.
- **`nativeThreads` is set explicitly on the native library at construction**
  (`openblas_set_num_threads` / `MKL_NUM_THREADS` / `OMP_NUM_THREADS`). **B-R1
  (REQUIREMENT): gale MUST NOT let a native library default to "all cores."** That
  default, combined with a JVM pool, is the classic oversubscription bug; gale
  pins native threads to `nativeThreads`.

### B.2 Oversubscription rules

- `allowNestedParallelism = false` (default): within a single operation, **at most
  one layer parallelizes.** When gale dispatches an op to the native library
  (which uses `nativeThreads`), the JVM pool MUST NOT also parallelize that op;
  and gale's own parallel pure paths (if enabled) do not nest native calls inside
  parallel regions.
- **Embedded in a parallel host (Spark, ZIO, Akka, a thread-per-request server):**
  the safe default is `BackendConfig.singleThreaded` — **gale adds no parallelism;
  the host owns the cores.** This is why the default is `(1, 1)`, not
  `(availableProcessors, …)`: N host tasks each running gale must not each spawn an
  M-thread pool. Standalone batch users opt into `BackendConfig.dedicated(cores)`.
- **B-R2 (REQUIREMENT): gale MUST NOT spawn a pool larger than `jvmThreads`**, and
  the total intended concurrency of any single op is `max(jvmThreads-region,
  nativeThreads)`, never their product. **Scope:** under the leaning
  process-global-pool decision (Open questions), "`jvmThreads`" here means the
  process pool's *configured* size under a **first-config-wins** rule — a later
  backend constructed with a smaller `jvmThreads` cannot shrink an already-sized
  shared pool, and likewise a single `nativeThreads` pin (B-R1) is process-global.
  B-R2/B-R1 are therefore guaranteed **per process**, not independently per backend;
  the Open question must settle the first-config-wins/reconfigure rule that makes
  this precise.

### B.3 What the pure kernels do, and the JS story

- **Pure kernels are single-threaded by default in v0.5.** The portable-correctness
  bar (deterministic, cross-platform) is single-threaded. A **parallel-outer-loop
  pure `gemm`** (partition the row blocks across the pool) is **opt-in** behind
  `jvmThreads > 1`, and is a §D roadmap item — *not* the default path and *not*
  required for v0.5 correctness. (PRD: "Large JVM kernels **may** use a fixed
  compute pool" — a permission, not a requirement.)
- **JS / Wasm: single-threaded, always.** `BackendConfig` is a JVM concept; on JS
  there is no pool and no native library. `jvmThreads` is ignored (effectively 1),
  and the parallel pure path does not exist on JS. This falls out of §C's "no
  native module on JS" and needs no separate machinery.

### B.4 Non-goals (§B)

1. No virtual threads for CPU-bound kernels.
2. No per-call threading configuration (construction-time only, G2).
3. No auto-detection of the host's parallelism — the safe default is
   single-threaded; the user declares `dedicated(cores)` for standalone use.
4. No nested parallelism by default.

---

## §C. Native `DMat` design (`gale-v0-5-native-dmat`)

### C.1 Decision — a separate `NativeDMat` type, not a native-backed `DoubleArray`

The PRD names `NativeDMat` as a distinct type, and this doc follows it decisively
over the alternative (making `DoubleArray` a heap-or-native union behind its
opaque type). The alternative is rejected on three verified grounds:

1. **The opaque-type + inline-extension invariant.** `DoubleArray` is
   `opaque type DoubleArray = Array[Double]` with `private[gale] inline`
   `apply`/`update`/`length` extensions that compile to bare array access. A
   native-or-heap union would force every one of those inline accessors — the
   hottest code in the library — to branch heap-vs-native per element. That is
   exactly the pure-path tax §A refuses.
2. **JS portability.** `DoubleArray` is `Float64Array` on JS. A `MemorySegment`
   variant cannot exist on JS; overloading `DoubleArray` would fracture the
   cross-platform storage abstraction that the whole core rests on.
3. **The interop shadowing contract.** `core/{jvm,js}/interop.scala` carries a
   documented invariant: the public `toArray`/`toArrayRowMajor` (JVM) /
   `toFloat64Array` (JS) extensions deliberately shadow same-named `private[gale]`
   members, with `userland.InteropSuite` as the canary. **`NativeDMat` is
   additive** — new type, new `toNative`/`toHeap` methods, none of which share a
   name with those members — so **that invariant is untouched (verified: the
   interop files define only `fromArrayUnsafe`/`toArray`/`toArrayRowMajor`/
   `toFloat64Array`).** A `DoubleArray` union would instead put native handling
   right through that shadowed surface.

### C.2 Shape, and where it lives

```scala
// gale.backend.jvm.native (JVM-only module) — NOT gale-core.
enum Layout:
  case RowMajor, ColMajor

final class NativeDMat private[backend] (
    private[backend] val memory: java.lang.foreign.MemorySegment,
    val rows: Int,
    val cols: Int,
    val leadingDimension: Int,      // stride between columns (ColMajor) / rows (RowMajor)
    val layout: Layout              // default ColMajor — the BLAS/LAPACK native convention
) extends AutoCloseable
```

- **`NativeDMat` and `Layout` live in the JVM native module, not core** — core
  stays pure/FFM-free (Doctrine 4). This resolves a PRD **Open Decision**
  ("Whether `NativeDMat` lives in the BLAS backend module or a small JVM-native
  storage module shared by BLAS and future sparse direct solvers") in favour of a
  small JVM-native module that `gale-backend-jvm-blas-ffm` and a future
  `gale-backend-jvm-sparse` both depend on. Cross-referenced, not re-opened.
- **`Layout` is introduced here** — it does not exist in core today (verified: no
  `enum Layout`; `DMat` uses raw strides + `isContiguousRowMajor` predicates).
  Default `ColMajor` because LAPACK is column-major; a heap `DMat` (row-major by
  construction) is transposed on copy-in when the native routine wants column-major.

### C.3 Ownership and lifetime — Arena-scoped, explicit, checked

- **`NativeDMat` is `AutoCloseable` and Arena-backed.** Its `MemorySegment` is
  allocated from a caller-supplied `java.lang.foreign.Arena`. **Whoever owns the
  Arena closes it**; `NativeDMat.close()` is a convenience that frees its segment.
- **C-R1 (REQUIREMENT): use-after-free is a checked `IllegalStateException`, not
  undefined behaviour.** FFM already enforces this — a `MemorySegment` from a
  closed `Arena` throws on access — so the requirement is "rely on the Arena scope,
  never on raw addresses," and gale must not cache a segment past its Arena.
- Two usage patterns, both explicit:
  - **Confined/scoped** (`Arena.ofConfined`, try-with-resources) for a single
    native computation — the common case; the segment is freed at scope exit.
  - **Shared** (`Arena.ofShared` / a long-lived Arena) for a `NativeDMat` reused
    across many native calls (e.g. an iterative solver's operator), freed when the
    owner disposes the Arena.

### C.4 Copy boundaries

```scala
extension (a: DMat)  def toNative(using Arena): NativeDMat     // heap → off-heap (+ maybe transpose to ColMajor)
extension (n: NativeDMat) def toHeap: DMat                     // off-heap → heap (row-major DMat)
```

- **Conversions are explicit** (PRD: "Conversions must be explicit or performed
  only above documented thresholds"). No implicit heap↔native coercion.
- **Copy cost is counted in the threshold decision (§A, §D).** A native op on a
  heap `DMat` above threshold does copy-in → compute → copy-out; the crossover
  size must include both copies (PRD: "A native call that copies heap-backed
  inputs into off-heap memory must include that cost in its threshold decision").
  When the inputs are already `NativeDMat`, the copy terms drop and the crossover
  is lower — §D distinguishes the two thresholds.
- **No zero-copy view between `DMat` and `NativeDMat`.** A heap `Array[Double]`
  cannot alias off-heap `MemorySegment` storage (nor the reverse), so every
  `toNative`/`toHeap` is a genuine copy — there is deliberately no `fromNativeView`.
  This is the off-heap mirror of the `gale → Breeze` copy-only rule shipped in
  `gale-interop-breeze`: an encapsulated store never aliases out, in either
  direction across the FFM boundary.

### C.5 The JS non-story

There is no FFM, no `MemorySegment`, no `Arena`, and no `NativeDMat` on JS. The
native module is JVM-only and simply absent from the JS classpath; JS code uses
heap `DMat`/`Float64Array` throughout. Nothing in core or the JS build references
native storage, so there is nothing to stub.

### C.6 Non-goals (§C)

1. No native-backed `DoubleArray` (C.1); native storage is a distinct type.
2. No `NativeDMat` in core; no FFM dependency in core (Doctrine 4).
3. No implicit heap↔native conversion; no hidden copies.
4. No manual `free`/pointer arithmetic exposed — lifetime is Arena scope only (C-R1).
5. No change to `core/{jvm,js}/interop.scala` or the shadowing invariant (C.1.3).

---

## §D. Blocked `gemm` roadmap (`gale-v0-5-blocked-gemm-roadmap`)

### D.1 Where the pure kernel is today (verified at `8333bde`)

`DoubleKernels.dgemm` branches once on layout:

- **Row-major** (unit column stride on A, B, C): `dgemmRowMajor`, an i-k-j loop
  built on a **4×4 register-tiled microkernel** (`gemmPanel`); above `64³ =
  262144` elements it switches to **`dgemmBlockedRowMajor`** with cache blocking
  at **`GemmBlock = 128`** (comment records 128 measured best with the 4×4 panel).
- **Any strided / transposed layout**: `dgemmStrided`, a general i-j-k
  dot-product loop with FMA, honouring arbitrary strides.
- `a.t * a` is routed to a dedicated assign-only **`dsyrkRowMajor`**.

Measured (`benchmarks/results/2026-07-11-breeze-release-grade.md`, vs Breeze 2.1
**pure-Java** netlib fallback): **gale is already ahead of pure Breeze on the
whole gemm family** — Gemm 1.11–1.19x, GemmTall 1.03–1.47x, AtA 2.08–2.98x. So in
the PRD's stage ladder `TinyUnrolledGemm → PureBlockedGemm → JvmVectorGemm →
NativeBlasGemm`, gale sits solidly at **PureBlockedGemm**. The open question is
which of the remaining stages earn their complexity, and where native takes over.

### D.2 Further pure stages, ranked by expected value

| Stage | What | Expected value | Cost | Gate |
|---|---|---|---|---|
| **D-1 Panel packing / copy-blocking** | Pack A and B blocks into contiguous L1/L2-resident scratch before the microkernel (GotoBLAS/BLIS structure), removing stride + TLB penalties at large `n` | Real at `n ≳ 256`, where B currently streams from L2/L3 (the bench note flags this); keeps the win growing instead of flattening | Medium; a scratch `Workspace` (PRD already sketches `Workspace`) | measure vs current `GemmBlock=128` |
| **D-2 Vector API microkernel (`JvmVectorGemm`)** | Rewrite the 4×4 (→ 8×4) panel with explicit `jdk.incubator.vector` FMA lanes | The single biggest pure lever — explicit SIMD is what F2J/native exploit; routed at the coarse gemm seam, it also lifts standalone **gemv** (0.67–0.90x, coarse-routable via `nativeGemvMinWork`). It does **not** touch standalone **axpy** (0.63x): that is a public L1 op held pure per A.2.1, so its store-bandwidth gap is *accepted*, not closed (unless the A.2.1 scoped-kernel-set escape hatch is later invoked) | Medium; **JVM-only, optional module** (`gale-v0-5-vector-backend`); incubator-API availability caveat | after D-1; behind `Capability.Vectorized` |
| **D-3 Wider register panels (8×4 / 6×8)** | Larger microkernel for more register reuse | Marginal **without** the Vector API (scalar JVM register pressure caps it) — fold into D-2, not worthwhile standalone | Low | with D-2 |
| **D-4 Parallel outer loop** | Partition row blocks across gale's compute pool | Largest lever at very large `n`, but **gated on §B** (single-threaded default) and opt-in `jvmThreads>1` | Medium | after §B; opt-in, not default |

Recommendation: **D-1 then D-2 are the worthwhile pure stages**; D-3 folds into
D-2; D-4 is opt-in and policy-gated. Anything beyond D-2 competes with native and
should defer to the crossover rule (D.4) rather than chase OpenBLAS in pure Scala
(explicit PRD non-goal: "Beating OpenBLAS/BLIS/MKL/Accelerate on large GEMM is not
a v1 goal").

### D.3 Measurement protocol

- **Existing JMH Breeze pairs** (`benchmarksJVM/Jmh/run .*Breeze.*`), the harness
  behind `benchmarks/results/*`, with `@Param` size sweeps, ≥2 forks, warmup +
  steady-state windows, and `-prof gc` allocation checks (all already in place).
  Pure stages (D-1..D-3) are validated here against pure-Breeze and against the
  prior gale baseline (reassociation allowed; `parityTest` + `testAll` must stay
  green; JS must not regress — re-run `benchSmokeJS`).
- **Scenario benchmarks** (the scenario-benchmark bead): real solver/factorization
  workloads, so gemm changes are judged in JIT-realistic context, not only
  isolated throughput (PRD: "JMH microbenchmarks are necessary but not
  sufficient").
- **Native / copy-inclusive benchmarks** (for D.4 crossover): FFM `dgemm` timed
  **including heap→native copy-in and copy-out** (PRD "heap-to-native conversion
  cost"), and separately with inputs already `NativeDMat`. Native comparison is vs
  **native BLAS** (OpenBLAS/MKL/Accelerate) — outside the current pure-Breeze
  sweeps, which is why it is a distinct protocol.

### D.4 The decision rule for native takeover (feeds `gale-v0-5-native-thresholds`)

Native FFM `dgemm` should take over exactly when it is faster end-to-end:

```
route to native  ⇔  t_native(n) + t_copyIn(n) + t_copyOut(n)  <  t_pure(n)
```

Because the copies are `O(n²)` and the multiply is `O(n³)`, the copy term
amortizes with size, so there is a crossover `n*` above which native wins **even
including copies**. Consequences fixed here:

1. **Two thresholds, both measured** — `nativeGemmMinFlops` for heap inputs
   (copy-inclusive) and a lower one for `NativeDMat` inputs (copy-free). The
   backend exposes them via `BackendThresholds` (§A.3).
2. **D-R1 (REQUIREMENT): thresholds are measured per platform and backend family,
   conservative by default** (PRD: "Measure native crossover thresholds instead of
   guessing them"; "Calling native code for vectors of length three is a bug").
   The default `nativeGemmMinFlops` is set *above* the measured crossover so
   borderline sizes stay on the already-competitive pure path and avoid native
   thread/copy complexity.
3. Below threshold, and on JS always, the tuned pure kernel (D.1, plus D-1/D-2
   when they land) runs. The pure path is never removed — it is the correctness
   reference (§A.6) and the whole of the JS product.

### D.5 Non-goals (§D)

1. Not beating native BLAS in pure Scala (PRD).
2. No parallel pure gemm by default (§B; opt-in only).
3. No native routing below the measured, conservative crossover (D-R1).
4. No new gemm layout/threshold knobs on the public API — tuning stays internal
   or on `BackendThresholds`.

---

## Appendix — decision log and cross-references

| ID | Decision | Trace |
|---|---|---|
| A-1 | Kernel seam at coarse L3/driver facades only; L1/L2 stay static-pure | §A.2 — pure JVM wins on memory-bound L1/L2 (PRD); keeps the hot path indirection-free |
| A-2 | Per-call `(using Backend)` on coarse ops; **no** mutable global install | §A.2 — Registry rejected (Doctrine 8, spectral doc precedent) |
| A-2b | Vector-API backend does **not** seam L1/L2; install-once kernel table in `DoubleKernels` rejected (Doctrine 8) — bandwidth-bound upside doesn't justify it; lever is `gemm`/`syrk` at the coarse seam | §A.2.1 |
| A-3 | `given Backend = Backend.pure` fallback (computes) + `compose`; mirrors `SpectralBackend.none`/`compose` | §A.4 — `pure` vs `none` naming is intentional (compute vs unsupported) |
| A-4 | Full `DenseDoubleKernel` surface declared; **dispatch policy** routes only coarse ops (gemm/syrk/standalone gemv/factorization drivers); loop-called L1/L2 never routes | §A.3, §A.2.1 — Vector-API lever is `gemm`, not memory-bound L1; no `nativeDotMinLength` |
| A-5 | `Backend.spectral: Option[SpectralBackend]` bridge; `NativeLapack ⇒ spectral.isDefined` | §A.5 — composes with spectral doc §2.3 without nesting |
| A-R1 | **REQUIREMENT**: no-import dispatch measured allocation-free and within JMH noise of the direct call | §A.2 — the PRD "backend dispatch overhead" benchmark |
| A-R2 | **REQUIREMENT**: `given Backend` is a shared singleton, concurrency-safe | §A.4 = spectral G1 |
| B-1 | One fixed compute pool (not virtual threads), sized `jvmThreads`; owned by the backend | §B.1 (PRD) |
| B-2 | Config at construction, never per call | §B.1 = spectral G2 |
| B-R1 | **REQUIREMENT**: native thread count pinned to `nativeThreads`; never let the native lib grab all cores | §B.1 |
| B-R2 | **REQUIREMENT**: never spawn a pool > `jvmThreads`; single-op concurrency is `max`, not product | §B.2 |
| B-3 | Default `BackendConfig.singleThreaded`; pure kernels single-threaded; parallel pure gemm opt-in | §B.2–B.3 — safe under a parallel host |
| C-1 | Separate `NativeDMat` type, not a native-backed `DoubleArray` | §C.1 — opaque-type/inline invariant, JS portability, interop shadowing |
| C-2 | `NativeDMat` + `Layout` live in a JVM-native module, not core; default `ColMajor` | §C.2 — resolves a PRD Open Decision |
| C-R1 | **REQUIREMENT**: Arena-scoped lifetime; use-after-free is a checked exception, not UB | §C.3 |
| C-3 | Explicit `toNative`/`toHeap`; copy cost counted in thresholds; interop.scala untouched | §C.4, §C.1.3 |
| D-1 | Worthwhile pure stages are packing (D-1) then Vector API (D-2); wider panels fold into D-2; parallel is opt-in | §D.2 |
| D-R1 | **REQUIREMENT**: native crossover thresholds measured per platform/family, conservative, copy-inclusive; two thresholds (heap vs NativeDMat) | §D.4 (PRD) |

**Cross-references to `docs/spectral-backend-boundary.md`:** resolution/`none`
fallback (spectral §2.1 ↔ A.4), `compose` combinator (spectral §2.2 ↔ A.4),
thread-safety G1 (spectral §1.1 ↔ A-R2), threading deferral G2 (spectral §1.1 ↔
B-2), conformance-by-laws (spectral §2.6 ↔ A.6), the `Backend.spectral` bridge
(spectral §2.3 ↔ A.5), and `SpectralCapability` vs `Capability` staying distinct
(spectral §1.1 ↔ A.3). This document does not restate or alter any spectral
decision; §C's constructor-visibility note is consistent with the spectral doc's
`private[spectral]` step 0 (both keep native/backend modules out of the sealed
core types).

### Open questions (for the implementation beads)

- **Whether the compute pool is per-`Backend` or a process-global gale pool the
  backend references.** Leaning process-global (one pool, sized by the first
  constructed backend's `jvmThreads`) to make B-R2 trivial across multiple
  backends; needs a documented "first config wins / reconfigure" rule.
- **`nativeFactorizationMinSize` per routine.** LU, Cholesky, and QR have
  different crossovers (the release-grade sweep shows QR pure is *ahead* of pure
  Breeze at all sizes while Solve/LU trail at 64–256); the native factorization
  threshold is likely per-routine, not one number. Feeds `gale-v0-5-native-thresholds`.
- **Vector API module vs incubator stability.** `jdk.incubator.vector` is still an
  incubator; the D-2 module must pin a JDK and document the fallback to the scalar
  panel when the module/JDK is absent.
