# Spectral Backend Boundary (gale v0.3.5 → v0.5)

## Purpose

`docs/spectral-parity.md` gated the **types**: every spectral result value, ADT,
and diagnostic traces to a MATLAB/SciPy capability row, and every dropped
capability traces to an explicit "out". This document gates the **backend**: it
defines the boundary through which an optional production engine (native
LAPACK/ARPACK/SLEPc-style, or any future accelerated provider) supplies the
computations that the v0.3.5 pure core deliberately left as `Left(...)` seams —
**without leaking a single native symbol, option, or dependency into the pure
`gale-core` API**.

The parity matrix marked a family of rows **backend**: § 5 (generalized
nonsymmetric / QZ), § 9 (rank-deficient GSVD), and the shift-invert targets of
§ 6–7. For those rows v0.3.5 shipped the *result types* and returns
`Left(UnsupportedOperation)` / `Left(RankDeficient)` where a backend would
compute. This document decides exactly what fills those seams, how it is
resolved into scope, how it is tested, and what it must never do.

This is a design, not an implementation. It specifies the contract the eventual
implementation must satisfy; like the parity doc, its decisions are meant to be
traceable and enumerable.

## Status and provenance

- Tracker bead: `gale-v0-3-5-spectral-backend-boundary` — "Define optional
  backend boundary for production spectral engines such as native
  LAPACK/ARPACK/SLEPc-style providers without leaking them into the pure core
  API."
- **Must compose with** the planned general acceleration contract, not rival it:
  - `gale-v0-5-backend-contract` — the general `Capability` / `Backend` /
    `DenseDoubleKernel` kernel-acceleration contract (PRD § Backend
    Requirements), and the owner of the shared `BackendConfig` threading policy.
  - `gale-epic-v0-5-acceleration` — the v0.5 epic that owns `gale-v0-5-blas-ffm`
    (FFM BLAS/LAPACK), `gale-v0-5-vector-backend` (JDK Vector API),
    `gale-v0-5-threading-policy` (`BackendConfig`), and the
    `NativeDMat`/threshold work.
- Related deferral bead: `gale-left-eigenvectors`. **Its dense half already
  landed** — dense nonsymmetric left and left-and-right eigenvectors ship in pure
  core (`Eigen.nonsymmetricVectorFlags` accepts all four `EigenVectors` cases;
  left vectors are the conjugated rows of `V⁻¹`, and a defective matrix surfaces
  as `Left(SingularMatrix)`). Only the *iterative* left-vector case remains a seam
  (§ 3, S4).

## Fixed constraints (inherited from the parity doc; not relitigated)

1. **`gale-core` stays pure.** No native, no Breeze, no Cats/ZIO/Spire, no
   JVM-only dependency in core (PRD Design Doctrine 4; Acceptance Criteria). The
   backend boundary is an *interface* in core; every *implementation* lives in a
   separate optional JVM module.
2. **Backend selection is explicit through `given` imports** (PRD Design
   Doctrine 8; `import gale.backend.jvm.blas.given`). No service-loader, no
   classpath scanning, no runtime download.
3. **Real inputs; `Complex` is a boundary value type only.** No complex storage
   or kernel tier — even a complex-capable native library must hand results back
   through gale's real structure-of-arrays packing.
4. **Either-first failure model.** Structural/precondition violations are
   `Left(LinAlgError)`; non-convergence and degeneracy are `Right(result +
   SpectralDiagnostics)` with an explicit `requireConverged` opt-in
   (parity doc § Convergence & failure semantics).
5. **Canonical ordering/layout is gale's, not the engine's** (parity doc
   § "Ordering guarantees gale commits to"). The output order is fixed by gale so
   results are reproducible across backends and platforms; a backend supplies
   spectra, never the order.

---

## 0. What v0.3.5 already froze — and the gaps this boundary must span

The shipped `gale.spectral` surface (`core/shared/src/main/scala/gale/spectral/`)
is the fixed ground this boundary attaches to. Four properties of it shape every
decision below.

**0.1 The result types have `private[gale]` constructors — which is *not* a wall,
and must be tightened.** `EigenDecomposition`, `NonsymmetricEigenDecomposition`,
`SVD`, and `GeneralizedSVD` are all constructed `private[gale]` (some with
construction-time `require` invariants — the nonsymmetric one validates
real-Schur conjugate-pair packing in its body). It is tempting to claim a backend
module therefore cannot build them, but **Scala `private[gale]` is visible to
*every* subpackage of `gale`, including the PRD's `gale.backend.*` modules** — so
a backend can call these constructors today. The wall does not exist yet.

The boundary makes it exist, as **migration step 0** (§ 5.0): tighten the four
constructors from `private[gale]` to `private[spectral]` (the enclosing
`gale.spectral` package). That genuinely excludes `gale.backend.*` while keeping
the facades — `Eigen` and `Svds`, both in package `gale.spectral` — able to
construct. Nothing outside `gale.spectral` constructs these types today, so the
tightening is a non-breaking, mechanical change. **Only after step 0** is the
design's central discipline — *a backend returns raw numeric factors; the facade
in `gale.spectral` canonicalizes them and assembles the sealed result* —
enforced by the compiler rather than by convention. That discipline is what keeps
ordering, sign, packing, residual measurement, and shape validation with the
facade (§ 1.5, § 2.5), where a backend cannot bypass them.

**0.2 `GeneralizedEigenvalue` is a producer-less, container-less type.**
`GeneralizedEigenvalue(alpha: Complex, beta: Double)` (with `value`, `isInfinite`,
`isFinitePair`) exists and is unit-tested, but nothing produces it: there is **no
`GeneralizedEigenDecomposition` container, no `eigGeneralizedNonsymmetric`
facade, and no `qz` entry point anywhere in the source.** The parity doc § 5
defined only the scalar `(α, β)` value type in v0.3.5. So the QZ boundary is not
"wire a backend into an existing method" — it must **introduce a new result
container and a new facade** alongside the backend (§ 1.3, § 5). This is the
largest single piece of new public surface the boundary implies, and defining its
shape *and its ordering machinery* (§ 1.3) now is exactly what keeps the eventual
backend from inventing an ad-hoc one.

**0.3 The pure dense kernels are a `private[gale] object`, not an abstraction.**
`DenseSpectralKernels.{symmetricEigen, nonsymmetricEigen, ...}` are hard object
calls inside the facades (`Eigen.eigSymmetric` calls
`DenseSpectralKernels.symmetricEigen(a, wantVectors)` directly). Nothing today is
polymorphic over "which engine computes the dense decomposition." An
*accelerated* dense path (a backend replacing `syev`/`geev`/`gesdd`) therefore
requires the facade to gain a resolution point *before* the hard object call —
not a rewrite of the kernels (§ 5). The pure kernels remain the default and the
JS/no-backend path.

**0.4 `UnsupportedOperation` shipped with one field, not the PRD's two.**
The frozen error is `LinAlgError.UnsupportedOperation(operation: String)`. The
PRD § Error Model sketch had `UnsupportedOperation(operation: String, backend:
String)`. **Decision: keep the one-arg form.** Every existing seam already uses
it; adding a field to a frozen case class is a breaking change for no benefit.
Backend context (which capability, whether a backend was registered) is carried
in the `operation` string, e.g. `"generalized nonsymmetric eigen (QZ): no
spectral backend registered"`. This is the single traceable discrepancy between
the PRD error sketch and the shipped ADT, and it is resolved in favour of the
shipped ADT.

---

## 1. The contract — `SpectralBackend`

### 1.1 Trait shape, capability discovery, and the trait-level requirements

```scala
package gale.spectral

/** An optional provider of spectral computations the pure core does not ship.
  * Lives in gale-core as an interface; every implementation is in a separate
  * optional JVM module and is brought into scope by an explicit `given` import
  * (§ 2). Methods take already-validated dense inputs and return *raw numeric
  * carriers* (§ 1.2) — never the sealed result types, which the facade builds.
  *
  * Thread-safety (a REQUIREMENT of the trait, G1): a given `SpectralBackend` is
  * a shared singleton resolved once and reused across all call sites and
  * threads. An implementation MUST be safe for concurrent invocation — either
  * stateless, or internally synchronized. A backend wrapping a non-reentrant
  * native library MUST serialize internally.
  */
trait SpectralBackend:
  def name: String
  def capabilities: Set[SpectralCapability]

  // Every operation has a default `Left(UnsupportedOperation)` body; a backend
  // overrides only the ones it supports and lists exactly those in `capabilities`.
  def denseSymmetricEigen(a: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
    unsupported("dense symmetric eigen")
  def denseNonsymmetricEigen(a: DMat, vectors: EigenVectors): Either[LinAlgError, RawNonsymmetricEigen] =
    unsupported("dense nonsymmetric eigen")
  def denseSvd(a: DMat, wantVectors: Boolean): Either[LinAlgError, RawSvd] =
    unsupported("dense SVD")
  def generalizedNonsymmetricEigen(a: DMat, b: DMat, vectors: EigenVectors): Either[LinAlgError, RawGeneralizedEigen] =
    unsupported("generalized nonsymmetric eigen (QZ)")
  def rankDeficientGsvd(a: DMat, b: DMat, wantVectors: Boolean): Either[LinAlgError, RawGsvd] =
    unsupported("rank-deficient GSVD")
  def shiftInvertOperator(a: DMat, b: Option[DMat], sigma: Double): Either[LinAlgError, DoubleLinearOperator] =
    unsupported("shift-invert operator")

  protected final def unsupported(op: String): Either[LinAlgError, Nothing] =
    Left(LinAlgError.UnsupportedOperation(s"$op: $name backend does not provide it"))

enum SpectralCapability:
  case DenseSymmetricEigen, DenseNonsymmetricEigen, DenseSvd   // acceleration of shipped pure ops
  case GeneralizedNonsymmetricEigen                            // QZ / ggev — the §5 seam
  case GeneralizedSchur                                        // optional (Q,Z,AA,BB) factors
  case RankDeficientGsvd                                       // ggsvd3 hard cases — the §9 seam
  case ShiftInvertSolve                                        // (A−σB) providers for §6–7 targets
  case IterativeGeneralized                                    // B-inner-product / ARPACK-class generalized
  // Deliberately absent: ComplexShift (needs a complex tier; out — see §4).

object SpectralBackend:
  /** The always-in-scope, lowest-priority fallback: no capabilities, every
    * operation `Left(UnsupportedOperation)`. Reproduces today's seam behaviour
    * exactly when no acceleration module is imported (§ 2.1). */
  given none: SpectralBackend with
    def name = "none"
    def capabilities = Set.empty
```

**Discovery decision — a capability `Set`, not `Option`-returning provider
methods.** The task poses "a capabilities set vs `Option`-returning methods." We
choose the **Set for discovery + total `Either`-returning methods for
invocation**, related by the conformance obligation:

> A backend that lists a `SpectralCapability` in `capabilities` MUST implement
> the corresponding method (return non-`UnsupportedOperation` for structurally
> valid input), and MUST NOT list a capability it does not honour.

This is stated as a **conformance obligation, not a type-level guarantee** (D-d):
at the type level a claimed-but-unimplemented capability is indistinguishable
from one that always returns `Left(UnsupportedOperation)`. The conformance suite
(§ 2.6) is what verifies the biconditional — it drives every method whose
capability is advertised and fails the backend if any returns
`UnsupportedOperation` on valid input.

Rationale for Set + Either over `Option`-returning accessors, in order of weight:

1. **One failure channel.** The whole library is `Either[LinAlgError, _]`-first.
   `Option`-returning provider accessors (`def qz: Option[QzProvider]`) fork the
   idiom: callers would juggle `Option` at discovery and `Either` at invocation.
   A `Set` for the yes/no question plus `Either` for the result keeps exactly one
   error channel — the same one the seams already speak.
2. **Composes with the general contract verbatim.** `gale-v0-5-backend-contract`
   already defines `Backend.capabilities: Set[Capability]`. Mirroring that shape
   (a `Set` of a spectral-specific enum) means the two contracts read as one
   family, not two conventions.
3. **The facade's current behaviour is already "call, handle `Left`."** Each
   shipped seam *is* an `Either`-returning method that returns
   `Left(UnsupportedOperation)`. Routing that call through
   `summon[SpectralBackend].method(...)` with the same return type is a drop-in.

**Why `SpectralCapability` is separate from the general `Capability` enum.** The
general `Capability` (`NativeBlas, Vectorized, Multithreaded, WasmFriendly,
Deterministic, …`) classifies *kernel acceleration* ("does this backend make
`gemm` fast?"). Spectral capabilities are *operation-level* ("can this backend do
QZ at all?") — a different granularity that does not belong in the kernel enum.
The two remain orthogonal and are both discovered through `given` imports (§ 2).

**Threading config (G2).** Threading configuration — the PRD's `BackendConfig`
from `gale-v0-5-threading-policy` — does **not** flow through the spectral
boundary in v1. Spectral facade methods take no per-call threading argument and
`SpectralOptions` gains no threading field (§ 4, non-goal 3). A backend reads its
thread policy from the shared `BackendConfig` owned by `gale-v0-5-backend-contract`
at construction time; per-operation threading control is deliberately deferred so
oversubscription policy never leaks into the pure option types.

### 1.2 Raw carriers and the operation signatures, against the frozen result types

Backends return **public, invariant-free numeric carriers** in `gale.spectral`;
the facade turns each into the corresponding `private[spectral]` sealed result
(after the step-0 tightening) following canonicalization (§ 2.5) and diagnostics
(§ 1.5). Carriers hold only what a native routine naturally produces.

```scala
final case class RawSymmetricEigen(values: DVec, vectors: DMat)          // ascending not required; facade sorts
final case class RawNonsymmetricEigen(                                    // real-Schur SoA, as geev/dgeev emit
    re: DVec, im: DVec, rightPacked: DMat, leftPacked: Option[DMat])
final case class RawSvd(sigma: DVec, u: DMat, vt: DMat)                   // any order; facade fixes descending
final case class RawGeneralizedEigen(                                     // the QZ payload; (α,β) SoA
    alphaRe: DVec, alphaIm: DVec, beta: DVec,
    rightPacked: DMat, leftPacked: Option[DMat],
    schur: Option[QzSchur])                                              // present iff GeneralizedSchur capable
final case class RawGsvd(u: DMat, v: DMat, x: DMat, c: DVec, s: DVec)     // c,s as direct column norms
final case class QzSchur(aa: DMat, bb: DMat, q: DMat, z: DMat)           // optional; not required by the boundary
```

Signature-to-result-type map (every method's output is assembled by the facade
into the named frozen type):

| Backend method | Raw carrier | Facade assembles → | Fills seam |
|---|---|---|---|
| `denseSymmetricEigen` | `RawSymmetricEigen` | `EigenDecomposition` | acceleration of dense `Eigen.eigSymmetric` |
| `denseNonsymmetricEigen` | `RawNonsymmetricEigen` | `NonsymmetricEigenDecomposition` | acceleration of dense `Eigen.eigNonsymmetric` (all four vector modes, including left, already pure — § 3 S3) |
| `denseSvd` | `RawSvd` | `SVD` (full dense — pure Golub–Kahan–Reinsch kernel ships as of v0.5) | acceleration of full dense SVD (parity § 3) |
| `generalizedNonsymmetricEigen` | `RawGeneralizedEigen` | **`GeneralizedEigenDecomposition`** (new, § 1.3) | parity § 5 — QZ |
| `rankDeficientGsvd` | `RawGsvd` | `GeneralizedSVD` | parity § 9 — rank-deficient GSVD |
| `shiftInvertOperator` | `DoubleLinearOperator` | consumed by the existing Lanczos/Arnoldi loop | parity § 6–7 — targeted selection |

Note `shiftInvertOperator` returns gale's own `DoubleLinearOperator`: a backend
factorizes `A − σB` (or `A − σI`) once and exposes `(A − σB)⁻¹` as an operator,
which drops straight into the **already-shipped** iterative
`Eigen.eigSymmetric(op, n, …)` / `eigNonsymmetric(op, n, …)` paths. The backend
supplies the fast solve; gale's pure Krylov iteration and canonicalization are
unchanged.

### 1.3 The new result container for QZ, and its canonical-order machinery

Because § 0.2 shows no container exists, the boundary must define one now, *with
its ordering rules*, since the existing nonsymmetric sorter cannot be reused
(D-b). It mirrors `NonsymmetricEigenDecomposition` (SoA real storage, real-Schur
packed vectors, typed accessors, `private[spectral]` constructor) but carries the
**projective `(α, β)`** spectrum so infinite eigenvalues from a singular/
rank-deficient `B` are representable — the whole reason `GeneralizedEigenvalue`
exists.

```scala
final class GeneralizedEigenDecomposition private[spectral] (
    private[spectral] val alphaRe: DVec,
    private[spectral] val alphaIm: DVec,
    private[spectral] val beta: DVec,                  // exactly 0.0 marks an infinite eigenvalue
    private[spectral] val rightVectorsPacked: DMat,    // real-Schur packing, as NonsymmetricEigenDecomposition
    private[spectral] val leftVectorsPacked: Option[DMat],
    val diagnostics: SpectralDiagnostics
):
  def size: Int
  def eigenvalue(i: Int): GeneralizedEigenvalue        // = GeneralizedEigenvalue(Complex(alphaRe,alphaIm), beta)
  def isInfinite(i: Int): Boolean                      // beta(i) == 0.0
  def isRealPair(i: Int): Boolean                      // alphaIm(i) == 0.0 (finite real λ, or infinite)
  def eigenvector(i: Int): (DVec, DVec)                // decode packing exactly like the nonsymmetric type
  def leftEigenvector(i: Int): (DVec, DVec)
  def requireConverged: Either[LinAlgError, GeneralizedEigenDecomposition]
  def requireExtremeCertified: Either[LinAlgError, GeneralizedEigenDecomposition]
```

**Canonical-order machinery this container needs (D-b).** The parity doc's
generalized-nonsymmetric row says "as nonsymmetric, comparing `α/β` with
infinities (`β = 0`) placed by criterion." `nonsymDenseIndices` only sorts
`Complex`, so a **projective sorter** must be added:

- `generalizedIndices(alphaRe, alphaIm, beta, order)` sorts by the criterion on
  `α/β` **without dividing by `β`**, via cross-multiplication, so `β = 0` never
  produces `∞`/`NaN`. For a magnitude order, compare pair *i* vs *j* by
  `|α_i|·|β_j|` against `|α_j|·|β_i|`; for a real-part order, by
  `α_i.re·β_j` against `α_j.re·β_i` (LAPACK returns `β ≥ 0`, so the sign of the
  cross-product is meaningful). An infinite eigenvalue (`β_i = 0`) then compares
  as the largest magnitude automatically: `|α_i|·β_j > |α_j|·0`.
  **`β ≥ 0` is a hard producer contract**, load-bearing for the comparator's
  totality (a negative `β` inverts cross-multiplied comparisons and breaks
  transitivity): it is documented on `RawGeneralizedEigen` and the facade
  '''guards''' it before sorting — any negative `β` is rejected as
  `Left(InvalidArgument)` naming the backend contract violation, never a
  sort-time crash.
- **Multiple infinite eigenvalues are mutually tied** (all "= ∞"). Break the tie
  deterministically — by descending `|α|`, then by original producer index — and
  document it as an arbitrary-but-fixed order, exactly as the parity doc fixes an
  arbitrary layout for otherwise-incomparable spectra.
- **Complex-infinite (`alphaIm ≠ 0` with `β = 0`) is not representable and not
  produced.** The real generalized Schur form emits infinite eigenvalues as
  `1×1` blocks (a zero on the `BB` diagonal with a real `AA` entry); complex
  conjugate pairs are always *finite* `2×2` blocks. So the packing invariant is
  **`β(i) = 0 ⟹ alphaIm(i) = 0`**: an infinite eigenvalue owns a single real
  column. The constructor validates this (a `require`, like the nonsymmetric
  packing check); a backend that violates it is a conformance failure, not a
  supported input.
- **Generalized pairing predicate.** The nonsymmetric packing check tests
  `re(i+1) = re(i) ∧ im(i+1) = −im(i)`. The generalized check adds **`β`
  equality**: a conjugate pair satisfies
  `alphaRe(i+1) = alphaRe(i) ∧ alphaIm(i+1) = −alphaIm(i) ∧ beta(i+1) = beta(i)`,
  with the positive-imaginary member first. Without the `β` clause two pairs that
  share an `α` but differ in `β` could be mispacked.

As with the pure paths, this layout is imposed by the facade on the backend's raw
output; the backend supplies `(α, β)` and vectors, gale supplies the order. **The
container and its facade ship as part of the boundary work but every public path
that yields one requires a `given SpectralBackend`** (§ 5), so no producer-less
type reaches the pure surface.

### 1.4 Failure contract — how backend outcomes map to `LinAlgError`

The mapping is identical to the pure paths (parity doc § Convergence & failure
semantics); the backend inherits the rule, it does not get to invent a new one.

| Situation | Result | Notes |
|---|---|---|
| Capability absent / no backend registered | `Left(UnsupportedOperation(op))` | Produced by `SpectralBackend.none` or a backend's default method. The **only** new-vs-today difference is that with a capable backend imported it becomes a `Right`. |
| Structural precondition (non-square, shape disagreement, empty) | same `Left` variant as pure (`NonSquareMatrix`, `DimensionMismatch`, `InvalidArgument`) | **The facade validates shape *before* calling the backend** — a backend never sees malformed input, so these never originate inside a backend. |
| `B` not positive-definite (generalized symmetric-definite acceleration) | `Left(NotPositiveDefinite)` | Facade's Cholesky gate runs first, exactly as today. |
| Rank-deficient pencil, now *computed* (QZ / GSVD backend) | `Right` with `Infinite`/`Zero` typed values (GSVD) or `β = 0` infinite eigenvalues (QZ) + `diagnostics.rank` | The `Left(RankDeficient)` seam becomes a `Right`: the backend's whole job is to *handle* the case the pure path refused. Facade returns `Left(RankDeficient)` only if **no** rank-deficient-capable backend is present. |
| Backend ran, did not fully converge (iterative / ARPACK-class) | `Right(partial result)`, `allConverged = false` | Same as `SolverResult.NotConverged` / pure Lanczos. Never a `Left`. `requireConverged` is the caller's residual-policy throw; `requireExtremeCertified` additionally enforces certified spectral membership. |
| LAPACK `info < 0` (illegal argument to the native routine) | `Left(InvalidArgument("<backend> reported illegal argument N to <routine>"))` | Indicates a gale↔backend marshalling bug; surfaces as `InvalidArgument`, and is a conformance failure (§ 2.6), not a user-facing normal path. |
| Native library unavailable / fails to load | **cannot occur at call time** | Availability is a *registration* responsibility (§ 2.2): if the native library cannot load, the module does not provide its `given`, so resolution falls back to `none` and the call is a clean `Left(UnsupportedOperation)`. A registered backend must not throw for environmental reasons mid-call. |

Backends **must not** throw for any of the above; `LinAlgError` extends
`RuntimeException`, and the total-API contract is that these surface as values.
A native segfault-class fault is out of the failure model entirely (it is a bug
in the native provider, not a modelled outcome).

### 1.5 Diagnostics obligations — the facade re-derives, the backend counts

`SpectralDiagnostics(requested, converged, residuals, orthogonalityError,
iterations, rank, extremalityCertified)` must be **honest**. To make honesty
structural rather than trusted, the boundary splits responsibility:

- **The facade computes `residuals` and `orthogonalityError`** from the vectors
  the backend returned, using **the same mathematics the `SpectralLaws`
  residual/orthogonality checks use** — but implementing that arithmetic itself.
- **The facade sets `extremalityCertified`** only when the backend result comes
  from a full-space reduction (or a future backend supplies a separately
  validated extremality certificate). Ritz residuals alone never set it.
  The facade **cannot import `SpectralLaws`**: `gale-laws` depends on
  `gale-core`, not the reverse (a core→laws import would be circular), and the
  law functions are throwing `Unit` assertions, not residual-returning
  functions. So the small residual/orthogonality arithmetic (`‖A v − λ v‖`,
  `‖β A x − α B x‖`, `‖MᵀM − I‖`, …) is duplicated deliberately in the facade;
  the conformance suite (§ 2.6, which *is* in `gale-laws`) then independently
  re-checks the assembled results via the laws. A backend therefore **cannot
  report a small residual it did not achieve** — gale re-derives it, and the
  laws re-check it.
- **The backend supplies `iterations`, `requested`/`converged`, and `rank`** —
  facts only it knows (an iterative/ARPACK engine's iteration count and which
  requested pairs met tolerance; a rank-revealing factorization's numerical
  rank). These arrive as a small `BackendConvergence(requested, converged,
  iterations, rank: Option[Int])` alongside the raw carrier for the iterative
  methods; dense one-shot methods report `iterations = 0`, `converged =
  requested`, and let the facade derive `rank` where meaningful (as the pure SVD
  path already does).

Consequence: the diagnostics a backend "must populate honestly" reduce to the
counts it alone owns; everything checkable is checked by gale.

---

## 2. Registration and resolution

### 2.1 Resolution: a `given SpectralBackend` with an always-present `none` fallback

**Decision: `given`-based resolution, with `SpectralBackend.none` as a
companion-object fallback.** Rejected alternatives:

- *A mutable `Registry` object with `register(backend)`* — global mutable state,
  order-dependent, hostile to the "explicit `given` imports" doctrine and to
  test isolation. Rejected.
- *An explicit `backend =` value argument on every call* — verbose, and it
  cannot be defaulted to "whatever the user imported" (a default argument is
  fixed at the definition site). Rejected as the *primary* mechanism, though the
  `using` clause below is exactly the type-directed version of it.

The mechanism:

```scala
// gale-core, companion of SpectralBackend — lowest priority, always in scope:
object SpectralBackend:
  given none: SpectralBackend = /* no capabilities, all Left(UnsupportedOperation) */

// gale-backend-jvm-lapack (optional module) — higher priority when imported:
package gale.backend.jvm.lapack
val backend: SpectralBackend = new LapackSpectralBackend(/* FFM handles */)
given SpectralBackend = backend   // convenience for the single-backend case
```

A backend-aware facade method takes `(using backend: SpectralBackend)`. With **no
import**, the only candidate is the companion `none`, so the method behaves
byte-for-byte as it does today (the seam's `Left`, or the pure kernel). With
`import gale.backend.jvm.lapack.given` in scope, that given outranks `none`
(import-scope beats companion-scope in Scala 3 given prioritization) and wins. No
scanning, no registry, no global state — resolution is lexical and explicit,
satisfying PRD Design Doctrine 8.

### 2.2 The multi-provider story — an explicit composite, never two givens

§ 2.6 and § 3 envision two capability-*disjoint* backend modules (an FFM-LAPACK
dense provider and a later ARPACK-class iterative provider). **Importing both
`given`s is a compile-time ambiguity** — Scala cannot pick between two
`given SpectralBackend` of equal priority. The boundary resolves this explicitly
(D-a):

> **At most one `given SpectralBackend` is ever in scope.** To use several
> capability-disjoint engines together, the user constructs a *composite*
> explicitly and provides it as the single given:

```scala
import gale.backend.jvm.lapack       // brings the `backend` VALUE, not a given
import gale.backend.jvm.arpack       // brings the `backend` VALUE, not a given
given SpectralBackend =
  SpectralBackend.compose(lapack.backend, arpack.backend)   // one given, no ambiguity
```

`compose` is part of the core contract:

```scala
object SpectralBackend:
  /** A composite whose `capabilities` is the UNION of the parts, dispatching each
    * operation to the FIRST part whose `capabilities` contains the needed
    * capability. Earlier parts win on overlap (documented precedence). Its own
    * thread-safety follows from the parts' (§ 1.1). */
  def compose(primary: SpectralBackend, rest: SpectralBackend*): SpectralBackend
```

Each acceleration module therefore exports its backend as a **plain `val`** (for
composition) *and* a convenience `given` (for the common single-backend case);
the multi-backend user imports the `val`s and declares one composite given. This
keeps "explicit `given` imports" intact while making the multi-provider story
realizable rather than an ambiguity.

### 2.3 Module layout and composition with `gale-v0-5-backend-contract`

```text
gale-core                    JVM + JS   SpectralBackend trait + compose, SpectralCapability,
                                        Raw* carriers, GeneralizedEigenDecomposition,
                                        SpectralBackend.none — pure, no native.
gale-backend-jvm-blas-ffm    JVM only   given Backend (DenseDoubleKernel via FFM BLAS/LAPACK)
                                        AND SpectralBackend value+given (QZ, GSVD, dense accel,
                                        shift-invert) — the FFM LAPACK provider.
gale-backend-jvm-arpack      JVM only   (later) SpectralBackend value+given with the iterative /
                                        IterativeGeneralized / large-scale ShiftInvertSolve
                                        capabilities — ARPACK/SLEPc-class.
```

Composition rule: `SpectralBackend` is a **peer** of `Backend`, resolved by the
same `given` mechanism, not nested inside it — a spectral op needs
`summon[SpectralBackend]`, not a `Backend` member walk. For the convenience of a
caller holding a `Backend`, `Backend` **may** expose
`def spectral: Option[SpectralBackend] = None` as an optional bridge (a single
native module overrides it), but the facade's resolution path is
`summon[SpectralBackend]`. This keeps the two contracts independently importable
and independently testable while letting one module satisfy both. It invents no
rival capability system: kernel acceleration stays in `Capability`, spectral
operations in `SpectralCapability`, threading in the shared `BackendConfig`.

### 2.4 The JS / Scala.js story

Scala.js (and the Wasm profile) have **no native backends** — none of the
acceleration modules cross-compile to JS. On JS the only `given SpectralBackend`
in scope is `none`. Therefore:

- Every backend-scoped operation on JS returns `Left(UnsupportedOperation)` (QZ,
  rank-deficient GSVD, shift-invert) — the *identical* value the JVM returns with
  no import. The pure paths (symmetric/nonsymmetric dense eig including left
  vectors, partial eig/SVD, generalized symmetric-definite, full-rank GSVD) run
  unchanged on JS.
- This is not a JS deficiency to apologize for; it is the boundary working as
  designed. The pure core is the whole product on JS, and the backend is a JVM
  performance/coverage option. Cross-platform determinism holds trivially on the
  pure surface (§ 2.5).

### 2.5 Determinism: the facade canonicalizes; equality is by law, not by bits

A native engine legitimately produces results that differ from the pure kernels
and from each other: eigenvector **signs/phases** are free, **ordering within
ties** is unspecified, the **partial-convergence subset** an iterative engine
returns can differ, and floating-point reassociation shifts low bits. The
boundary's determinism story has two layers:

1. **Layout/order is re-imposed by the facade, not trusted from the backend.**
   The facade applies gale's canonical ordering to every backend result using the
   *same* helpers the pure path uses (`denseSelectionIndices`/`assembleDense` for
   symmetric; `nonsymDenseIndices`/`assembleNonsymDense` for nonsymmetric,
   preserving adjacent conjugate pairs; the new `generalizedIndices` of § 1.3 for
   QZ; descending sort for singular values; `Infinite`-first/`Zero`-last for
   GSVD). So the parity doc's ordering table holds regardless of which engine
   computed the spectrum. A backend that emits a different order is *silently
   corrected*, not passed through.
2. **Equality across engines is "passes the laws," not "bit-identical."** Signs
   and phases are *not* canonicalized to a bit-pattern; instead the conformance
   contract is stated over sign-agnostic invariants — residuals, orthogonality,
   ordering, membership, the CS identity — which every valid decomposition
   satisfies. Within one backend on one platform, results are deterministic (same
   input → same output); *across* backends only law-equivalence is promised, and
   this is documented as the guarantee. (This mirrors the existing gale-vs-gale
   reassociation stance in the kernel-tuning beads.)

Not canonicalizing sign is deliberate: a global sign rule is undefined at exact
ties, does not extend to a *phase* convention for complex conjugate-pair columns,
and buys nothing the residual/orthogonality laws don't already cover.

### 2.6 Conformance: `gale-laws` / `SpectralLaws` *is* the conformance kit

`gale-laws` is a reusable cross-built module (`laws` in `build.sbt`) that already
encodes every guarantee the parity doc commits to. A `SpectralBackend` conforms
iff, routed through the facade, its results pass the `SpectralLaws` family for
each advertised capability. A new `SpectralBackendConformanceSuite` (in
`gale-laws`) drives this: parameterized by a `SpectralBackend`, it reads
`capabilities`, runs exactly the laws each capability implies over the shared
ScalaCheck matrix generators, and — closing the § 1.1 biconditional — asserts
that every advertised capability's method returns a `Right` on valid input.

| `SpectralCapability` | Laws (from `SpectralLaws`) that must hold |
|---|---|
| `DenseSymmetricEigen` | `symmetricResidual`, `orthonormalColumns`, `ascending`, `symmetricMembership` |
| `DenseNonsymmetricEigen` | `nonsymmetricResidual`, `nonsymmetricOrdering`, `nonsymmetricMembership` (+ a left-vector residual `‖wᴴA − λwᴴ‖` when left vectors are returned) |
| `DenseSvd` | `svdResidual`, `orthonormalColumns` (U and V), `descending`, `singularMembership`, `singularRankConsistent` |
| `GeneralizedNonsymmetricEigen` | **`generalizedNonsymmetricResidual`**, **`generalizedEigenOrdering`** — *new laws, below* |
| `GeneralizedSchur` (optional) | reconstruction `‖Q AA Zᵀ − A‖`, `‖Q BB Zᵀ − B‖`; orthonormality of `Q`, `Z` — *new law* |
| `RankDeficientGsvd` | `gsvdReconstruction`, `csIdentity`, `gsvdWellDeterminedOrthonormal`, `gsvdDescendingRatio` (existing — written to handle `Zero`/`Infinite`, not yet exercised by a rank-deficient producer) |
| `ShiftInvertSolve` | `symmetricResidual`/`nonsymmetricResidual` on the targeted subset, `subsetOfSpectrum` |
| `IterativeGeneralized` | `generalizedResidual`, `bOrthonormal`, `subsetOfSpectrum` |

**The QZ law family the kit must add (traceable to § 0.2).** `SpectralLaws` has
**no law exercising `GeneralizedEigenvalue` inside a decomposition** — the
existing `generalizedResidual` is for the real-λ symmetric-definite
`EigenDecomposition`, not the projective `(α, β)` QZ result. The boundary work
must add, in `gale-laws`:

- `generalizedNonsymmetricResidual(a, b, d: GeneralizedEigenDecomposition, tol)` —
  the **homogeneous** residual `‖β·A x − α·B x‖ ≤ tol·max(1, ‖A‖, ‖B‖)`, in real
  arithmetic. The homogeneous form (not `‖A x − (α/β)·B x‖`) is what makes it
  finite for every case:
  - **Finite** eigenvalue `λ = α/β` (`β ≠ 0`): `β·A x − α·B x = 0` is equivalent
    to `A x = λ·B x` — the correct pencil residual, with no division.
  - **Infinite** eigenvalue (`β = 0`): the expression reduces to `‖−α·B x‖ =
    |α|·‖B x‖`, i.e. it checks that the eigenvector lies in `null(B)` — the
    correct test for an infinite generalized eigenvalue. (It is **`‖B x‖`, not
    `‖A x‖`** — this is the fix to the earlier draft, which had the reduction
    backwards.)
- `generalizedEigenOrdering(d, order, tol)` — the criterion is monotonic in
  `α/β` under `generalizedIndices` (compared by cross-multiplication, § 1.3),
  infinities are placed at the criterion's extreme with the documented
  `|α|`-descending tiebreak among themselves, and every conjugate pair is
  adjacent with the positive-imaginary member first and equal `β`.

These are additions to `gale-laws`, not to `gale-core`; the frozen v0.3.5 core
surface is untouched by the conformance work.

---

## 3. Seam inventory — every shipped `Left(...)` mapped to a backend operation

Enumerating every place v0.3.5 refuses, and what fills it. "Fill class" states
whether v0.5's **FFM LAPACK** backend (dense, `gale-v0-5-blas-ffm`) suffices, an
**ARPACK-class** iterative provider is required, it is **pure-deferrable** (no
backend needed), or it is **out**. Loci are given by method, not line number
(line numbers drift).

| # | Seam (method / branch) | Today's return | Backend op that fills it | Fill class |
|---|---|---|---|---|
| S1 | `Eigen.eigSymmetric(op, n, sel, opts, target)`, `target = Some(_)` branch | `Left(UnsupportedOperation("shift-invert / targeted selection"))` | `shiftInvertOperator(a, None, σ)` → feeds existing Lanczos | **Pure-deferrable** (wire gale's own LU/iterative solve — the scaladoc says "until the `LinearSolvePlan` wiring lands"); *optionally* FFM-accelerated factorization |
| S2 | `Eigen.eigNonsymmetric(op, n, sel, opts, target)`, `target = Some(_)` branch | `Left(UnsupportedOperation(...))` | `shiftInvertOperator(a, None, σ)` → feeds existing Arnoldi | **Pure-deferrable**; optionally FFM-accelerated |
| S3 | dense nonsymmetric left / left-and-right eigenvectors | **ships in pure core** (`nonsymmetricVectorFlags` accepts all four modes; left = conjugated rows of `V⁻¹`; defective ⇒ `Left(SingularMatrix)`) | `denseNonsymmetricEigen(a, LeftAndRight)` — **acceleration only** | **FFM LAPACK** (`dgeev` left+right), a speed option; *no seam remains* — the `gale-left-eigenvectors` dense half is done |
| S4 | `Eigen.eigNonsymmetric(op, …)` iterative, left vectors: `validateArnoldiVectors` rejects `Left`/`LeftAndRight` | `Left(UnsupportedOperation("… a Krylov basis for A gives no left vectors; use the dense eigNonsymmetric"))` | two-sided Krylov, or `Aᵀ`-Arnoldi re-pairing — **the genuine left-vector seam** | **ARPACK-class** / **pure-deferrable** (not a dense-LAPACK seam) |
| S5 | **QZ / generalized nonsymmetric eigen — no facade exists** (parity § 5; type present, § 0.2) | *no entry point*; parity doc specifies `Left(UnsupportedOperation)` | `generalizedNonsymmetricEigen(a, b, vectors)` → new `GeneralizedEigenDecomposition` | **FFM LAPACK** (`ggev`; `gges`+`tgsen` if `GeneralizedSchur`) |
| S6 | `Svds.gsvd(a, b, vectors)` rank-deficient pencil (`m+p < n`, or measured `rank < n`) | `Left(RankDeficient(rank, n))` | `rankDeficientGsvd(a, b, wantVectors)` → `GeneralizedSVD` with `Infinite`/`Zero` | **FFM LAPACK** (`ggsvd3`: `ggsvp3` + `tgsja`) |
| S7 | Full **dense** SVD (parity § 3) — CLOSED for coverage in v0.5: the pure Golub–Kahan–Reinsch kernel ships behind `Svds.svd(a, All)` / `DMat.svd` | pure dense kernel computes (economy factors) | `denseSvd(a, wantVectors)` → `SVD` | **FFM LAPACK** (`gesdd`/`gesvd`) — acceleration only |
| S8 | Accelerated dense symmetric/nonsymmetric eig for production scale (PRD Backend Performance Strategy) — **symmetric half WIRED (v0.5)**: `Eigen.eigSymmetric(a, sel, vectors)` takes `(using SpectralBackend)` and routes through `SpectralBackend.routesDenseSymmetricEigen(n)` (capability ∧ `n ≥ denseSymmetricEigenMinSize`, the measured spectral threshold — Accelerate 128 per `benchmarks/results/2026-07-17-ffm-lapack-crossover.md`, unswept families `Int.MaxValue`); the facade validates before the gate, re-imposes ascending order, re-derives residuals, falls back to the pure kernel on a provider `Left`, and throws `InvalidArgument` on malformed factors. Nonsymmetric half remains pure-only | pure kernels run below threshold / with no import (byte-identical) | `denseSymmetricEigen` / `denseNonsymmetricEigen` | **FFM LAPACK** (`dsyev` shipped; `geev` pending) |
| S9 | Iterative **generalized** eigen (`B`-inner-product Lanczos, large/sparse; parity § 6 "in-b", no operator facade shipped) | *no operator facade* | `IterativeGeneralized` provider, or pure B-Lanczos | **Pure-deferrable** (small) / **ARPACK-class** (large) |
| S10 | Complex shift σ (off the real axis; parity § Explicitly OUT) | n/a — `SpectralTarget.sigma` is `Double` | would need complex solves / complex tier | **Out** (§ 4) |

**Which the v0.5 FFM LAPACK backend fills directly:** S5 (QZ), S6 (rank-deficient
GSVD), S7 (full dense SVD), S8 (accelerated dense eig), the factorization behind
S1/S2, and S3 as a pure acceleration. This is the bulk of the boundary's value
and the first backend to build.

**Which need an ARPACK/SLEPc-class iterative provider (a later, separate
module):** the *large-scale* forms of S1/S2 (shift-invert at scale), **S4 (the
one genuine left-vector seam — iterative)**, and S9 (iterative generalized).
These are matrix-free/sparse and are not what a dense LAPACK backend addresses.

**Which need no backend at all (pure-deferrable, do in pure gale first):** the
*small-dense* forms of S1/S2 (wire the existing `LinearSolvePlan` into the shipped
Lanczos/Arnoldi) and, if desired, S4/S9's small cases. The boundary explicitly
does **not** force these behind a native backend; it just *also* lets a backend
serve them faster. **S3 is already done in pure core** and needs no backend at
all except for speed.

---

## 4. Non-goals — what the boundary deliberately does not cover

1. **No plugin/classloader/ServiceLoader magic.** Resolution is `given` imports
   only (§ 2.1). No classpath scanning, no `META-INF/services`, no reflection.
2. **No runtime downloads or auto-provisioning of native libraries.** Native
   availability is a build/deploy concern; a module either loads at registration
   or does not provide its `given`/`val` (§ 1.4, § 2.2).
3. **No backend-specific options in the pure option types.** `SpectralOptions`,
   `SpectralTarget`, `EigenSelection`, `SingularSelection`, `LinearSolvePlan`
   stay engine-neutral. No `lapackDriver = "evr"`, no `arpackNcv`, no
   `nobalance`, no `'buckling'`/`'cayley'` modes, **and no threading field**
   (threading lives in the shared `BackendConfig`, G2). Backend tuning that has no
   engine-neutral meaning lives inside the backend, never in `gale-core`.
4. **No new complex storage or kernel tier.** Even a complex-capable native
   library returns results through gale's real SoA packing; `Complex` stays a
   boundary value type. Complex-shift targeting (S10) is therefore out.
5. **No behavioral change to the no-backend path.** With no acceleration import,
   every method must be byte-for-byte what v0.3.5 ships — same `Left`, same pure
   kernel output (§ 5). The boundary is *additive*.
6. **No mandatory public Schur/QZ factor API.** `GeneralizedSchur` (the
   `(Q, Z, AA, BB)` factors) is an *optional* capability a backend may expose;
   the boundary does not require a public generalized-Schur result type in this
   round (parity § 5 marks it "backend / deferred"; public Schur is PRD "Later").
7. **The backend does not own ordering, sign, canonicalization, residual
   measurement, shape validation, or threading policy.** Those are the facade's
   (or the shared config's), by construction (§ 1.5, § 2.5, G2). A backend is a
   numeric engine, not a policy source.
8. **No eigenvalue-cluster reordering (`ordqz`/`tgsen`) surface** beyond what a
   `GeneralizedSchur`-capable backend may do internally (parity § 5).

---

## 5. Migration sketch — the minimal change-set when the first backend lands

The cost of this boundary must be visible now. Here is the entire change-set for
the first real backend (`gale-backend-jvm-lapack`, filling S5/S6/S7/S8 and
accelerating S3), split into "byte-identical" and "additive."

**Pre-1.0 note.** gale's binary-compatibility policy begins at v1.0 (roadmap).
Until then, tightening a constructor (step 0) or adding a
`(using SpectralBackend = ...)` clause to a shipped method is acceptable *source*
evolution; the invariant the boundary guarantees is **behavioral** identity of
the no-import path, which holds because `none` reproduces today's exact returns.

### 5.0 Step 0 — tighten the result-type constructors (prerequisite)

Change the constructors of `EigenDecomposition`, `NonsymmetricEigenDecomposition`,
`SVD`, and `GeneralizedSVD` from `private[gale]` to `private[spectral]` (§ 0.1).
This walls `gale.backend.*` out of them, so the raw-carrier / facade-assembly
discipline becomes compiler-enforced. Non-breaking: nothing outside
`gale.spectral` constructs these types today. Everything below assumes step 0.

### 5.1 What stays byte-identical (no acceleration import)

- All shipped pure kernels (`DenseSpectralKernels.*`) — untouched.
- Every facade call **without** a backend import: `SpectralBackend.none` is
  summoned, so `eigSymmetric`/`eigNonsymmetric`/`eigSymmetricGeneralized`/`svd`/
  `gsvd` return exactly today's values — the same `Right(pure result)` and the
  same `Left(UnsupportedOperation)` / `Left(RankDeficient)` at the seams. This
  includes dense nonsymmetric **left** vectors, which already ship (S3).
- All frozen result types except the *additive* `GeneralizedEigenDecomposition`.
- The entire JS surface (§ 2.4).

### 5.2 What grows a backend-aware overload / entry point

| Facade | Change | Behaviour with `none` |
|---|---|---|
| `Eigen.eigSymmetric(a, sel, vectors)` | **LANDED (S8)** — gains `(using b: SpectralBackend)`; consults `b.denseSymmetricEigen` when `b.routesDenseSymmetricEigen(n)` (capability ∧ measured size threshold), else the pure call; provider `Left` falls back to pure, malformed factors throw | pure kernel — identical |
| `Eigen.eigNonsymmetric(a, sel, vectors)` | gains `(using b)`; routes **acceleration** (S3/S8, all four vector modes) to `b.denseNonsymmetricEigen` | pure kernel (incl. left vectors) — identical |
| `Eigen.eigSymmetric/eigNonsymmetric(op, n, sel, opts, target)` | when `target = Some(ShiftInvert(σ, plan))`, obtain the solve op (S1/S2) — from `b.shiftInvertOperator` if `ShiftInvertSolve`-capable, else the pure `LinearSolvePlan` wiring — and run the existing Krylov loop | `target` still `Left(UnsupportedOperation)` until the pure wiring lands; identical |
| `Svds.svd(a, sel, vectors)` | gains `(using b)` for `DenseSvd` (S7) | pure — identical |
| `Svds.gsvd(a, b, vectors)` | gains `(using bk)`; on the rank-deficient branch (S6) calls `bk.rankDeficientGsvd` if capable, else the shipped `Left(RankDeficient)` | `Left(RankDeficient)` — identical |
| **new** `Eigen.eigGeneralizedNonsymmetric(a, b, vectors)(using bk)` | brand-new facade (S5) returning the new `GeneralizedEigenDecomposition`; `Left(UnsupportedOperation)` with no capable backend | `Left(UnsupportedOperation("… QZ: no spectral backend registered"))` |

The facade dispatch pattern is uniform and small:

```scala
def eigGeneralizedNonsymmetric(a: DMat, b: DMat, vectors: EigenVectors)(using
    backend: SpectralBackend
): Either[LinAlgError, GeneralizedEigenDecomposition] =
  validateSquarePencil(a, b).flatMap { n =>                    // 1. facade validates shape (Left as usual)
    backend.generalizedNonsymmetricEigen(a, b, vectors).map {  // 2. raw compute, or clean Left from `none`
      raw => assembleGeneralized(a, b, raw, order, vectors)    // 3. facade sorts (generalizedIndices),
    }                                                          //    packs, measures residuals, builds result
  }
```

### 5.3 New code the backend module itself adds (outside `gale-core`)

- `gale-backend-jvm-lapack`: `val backend: SpectralBackend` + convenience `given`,
  the FFM bindings, and the raw-carrier producers. **Zero access to
  `private[spectral]`** — it only builds `Raw*` carriers.
- `gale-laws`: the QZ law family (§ 2.6) and `SpectralBackendConformanceSuite`.

### 5.4 The one-time additive core surface (ships with the boundary, gated behind a backend)

- `SpectralBackend` (+ `compose`), `SpectralCapability`, the `Raw*` carriers,
  `BackendConvergence`.
- `GeneralizedEigenDecomposition` (§ 1.3), its `generalizedIndices` sorter, and
  the `eigGeneralizedNonsymmetric` facade.
- Optional `Backend.spectral: Option[SpectralBackend]` bridge (§ 2.3).

Total footprint: step-0 constructor tightening, a trait + enum + combinator + six
carriers + one new result type + its sorter + one new facade in core,
`(using SpectralBackend)` threaded through five existing facades with
`none`-identical behaviour, and the conformance additions in `gale-laws`. Nothing
native enters `gale-core`; nothing changes for a caller who does not import a
backend.

---

## Appendix — decision log (traceable)

| ID | Decision | Rationale / trace |
|---|---|---|
| D1 | Backends return `Raw*` carriers; facade builds the sealed results — enforced by tightening constructors to `private[spectral]` (step 0) | § 0.1, § 5.0 — `private[gale]` is visible to `gale.backend.*`, so it is not a wall until tightened; then facade-owned canonicalization/diagnostics is compiler-enforced |
| D2 | Discovery = `Set[SpectralCapability]`; invocation = total `Either` methods; the biconditional is a **conformance obligation**, not a type guarantee | § 1.1, § 2.6 (D-d) — claimed-but-unimplemented is indistinguishable at the type level; the suite verifies it |
| D3 | `SpectralCapability` separate from general `Capability` | § 1.1 — operation-level vs kernel-acceleration granularity |
| D4 | Keep one-arg `UnsupportedOperation(operation)`; encode backend context in the string | § 0.4 — shipped ADT is frozen; PRD's two-arg sketch not adopted |
| D5 | `given`-based resolution + companion `SpectralBackend.none` fallback | § 2.1 — PRD Doctrine 8; no registry/global state; no-import = today's behaviour |
| D6 | Multi-provider via explicit `SpectralBackend.compose`; **never two givens** | § 2.2 (D-a) — two givens are a compile-time ambiguity; modules export a `val` for composition |
| D7 | Facade re-imposes canonical order/packing on backend output, incl. a new projective `generalizedIndices` sorter for QZ | § 1.3, § 2.5 (D-b) — parity ordering table holds across engines |
| D8 | Cross-engine equality = "passes the laws," not bit-identical; sign/phase not canonicalized | § 2.5 — residual/orthogonality laws are sign-agnostic |
| D9 | `SpectralLaws` (in `gale-laws`) is the conformance kit; facade **re-implements** the same arithmetic (cannot import laws — core↛laws, and laws are throwing assertions) | § 1.5 (D-c), § 2.6 |
| D10 | New `GeneralizedEigenDecomposition` + its packing invariants (`β=0 ⟹ alphaIm=0`; pair predicate adds `β` equality) defined now | § 0.2, § 1.3 (D-b) |
| D11 | Shift-invert and *iterative* left vectors (S4) are the real seams; dense left vectors (S3) already ship; QZ and rank-deficient GSVD are backend-required | § 3 (C2) — separates shipped, pure-deferrable, and backend-required |
| D12 | Backend `given`s are shared singletons; the trait REQUIRES concurrency-safety | § 1.1 (G1) |
| D13 | Threading config does not flow through the boundary in v1; backends read the shared `BackendConfig` at construction | § 1.1, § 4 (G2) |
| D14 | QZ homogeneous residual `‖βAx − αBx‖`: finite ⇒ pencil residual; infinite (β=0) ⇒ `|α|·‖Bx‖` (x ∈ null(B)), **not** `‖Ax‖` | § 2.6 (C3) — the flagship new conformance law, corrected |

### Open questions (for the implementation bead, not blocking this design)

- **Module name for the first spectral backend.** `gale-backend-jvm-lapack` vs
  folding the spectral provider into `gale-backend-jvm-blas-ffm` (one FFM module
  for both kernels and LAPACK spectral). Leaning to the latter; named `-lapack`
  here for exposition.
- **`BackendConvergence` vs widening the carriers.** Whether the iterative
  count/rank travels as a separate small record or as fields on each `Raw*`
  carrier. Cosmetic; does not affect the boundary.
- **`denseSvd` (S7) as "acceleration" or "new capability."** Full dense SVD is
  pure-deferred (needs a bidiagonal-SVD kernel), so a backend providing it fills
  a *coverage* gap, and may warrant its own pure implementation later regardless.
