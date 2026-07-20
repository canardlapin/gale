# Spectral Capability Parity Matrix (MATLAB / SciPy → gale v0.3.5)

## Purpose

This document is a capability-parity matrix between MATLAB's and SciPy's spectral
APIs and the spectral surface gale intends to ship in **v0.3.5**. It exists to
**gate the design of the v0.3.5 spectral result types** (`EigenDecomposition`,
`SVD`, generalized spectral results, `SpectralDiagnostics`, and the selection /
ordering / options ADTs). Every type, enum, and diagnostic the spectral-types
bead defines should be traceable back to a row in one of the tables below, and
every capability we deliberately drop should be traceable to an explicit "out"
decision in the closing section.

Gale targets MATLAB/SciPy's spectral **capability space**, not their syntax.
Stringly-typed mode flags (`'LM'`, `'SM'`, `'LA'`, `'largestabs'`, …) become
typed ADTs; return tuples become named result values; convergence flags and
exceptions become gale's Either-first / diagnostics idiom.

## Fixed constraints (decided; not relitigated here)

These are inputs to the matrix, not open questions:

1. **Real inputs, complex only in outputs.** v0.3.5 spectral inputs are real
   (`DMat`/`DVec`). A boxed boundary value type `Complex(re: Double, im: Double)`
   exists only at the API boundary — there is no complex storage tier and no
   complex kernel tier.
2. **Nonsymmetric eigen packing is caller-invisible.** Nonsymmetric eigen results
   store structure-of-arrays real data (`re: DVec`, `im: DVec`) plus a real `DMat`
   of eigenvectors in LAPACK real-Schur packed convention, exposed only through
   typed accessors (`eigenvalue(i): Complex`, `isRealPair(i)`, `eigenvector(i)`).
   The packing never appears in a public signature.
3. **Symmetric and SVD APIs are complex-free.** Real eigenvalues, real singular
   values, real vectors — no `Complex` anywhere in those signatures.
4. **Phase split.** Phase **a**: dense + partial **symmetric eigen** and
   **partial SVD** (full dense SVD stays deferred — see § 3).
   Phase **b**: **nonsymmetric eigen**, **generalized eigen**, and **GSVD**.
5. **Internal dense kernels.** Householder tridiagonalization + implicit-shift
   tridiagonal QL/QR (symmetric); Hessenberg reduction + Francis double-shift QR
   (nonsymmetric). Partial methods (Lanczos/Arnoldi) solve their small projected
   problems with these same kernels.

## Error-model idiom this matrix maps onto

gale's existing precedents, which the spectral types must be consistent with:

- **`LinAlgError`** — sealed ADT: `NonSquareMatrix`, `DimensionMismatch`,
  `SingularMatrix`, `NotPositiveDefinite`, `RankDeficient`, `InvalidArgument`,
  `UnsupportedOperation`, `DidNotConverge`, … Total APIs return
  `Either[LinAlgError, A]`; `.orThrow` is the throwing convenience.
- **`FactorizationDiagnostics(info: Int, rank: Option[Int])`** with
  `isSuccess = (info == 0)` — the LAPACK-style "non-fatal problem recorded as a
  number, result still returned" pattern.
- **`SolverResult.{Converged, NotConverged}`** — both carry `x`, `iterations`,
  `residual`; **`NotConverged` returns the best iterate, it is not a `Left`**.
  `orThrow` raises `DidNotConverge`.
- **`conditionEstimate`** returns **`Right(Double.PositiveInfinity)`** for a
  singular matrix — a degenerate-but-meaningful answer is a `Right`, not a `Left`.

These three precedents fix the convergence-mapping rule used throughout
(§ Convergence & failure semantics): **structural/precondition violations are
`Left`; non-convergence and degeneracy are `Right` + diagnostics.**

Legend for the "gale v0.3.5 plan" column: **in-a** / **in-b** = shipped in phase
a / b; **deferred** = out of v0.3.5, on the roadmap (PRD "Later"); **backend** =
API + result type defined in v0.3.5, pure implementation behind an optional
backend boundary; **out** = not planned (rationale in the closing section).

---

## 1. Dense symmetric / Hermitian eigen (`eigh`)

LAPACK: `syev`/`syevd`/`syevr`/`syevx` (standard), `sygv`/`sygvd`/`sygvx`
(generalized). Real symmetric ⇒ real eigenvalues, orthonormal eigenvectors.

| Capability | MATLAB | SciPy | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| All eigenvalues (values only) | `eig(A)` (symmetric detected; ascending) | `eigh(A, eigvals_only=True)` | **in-a** `eigSymmetric(A, EigenSelection.All, ValuesOnly)` | Both ecosystems return **ascending**; gale guarantees ascending-algebraic. |
| All eigenpairs | `[V,D]=eig(A)` | `eigh(A)` → `(w, V)` | **in-a** | `A V = V diag(λ)`, `V` orthonormal. Complex-free. |
| Top/bottom `k` by algebraic value | not in dense `eig` (use `eigs`) | not by `k`; use `subset_by_index` | **in-a** `EigenSelection.Count(k, LargestAlgebraic\|SmallestAlgebraic)` | MATLAB has **no dense subset** — forces users onto iterative `eigs` even for dense. gale offers it directly (slice of the full spectrum in phase a). |
| Subset by index range | — | `eigh(A, subset_by_index=[il, iu])` (0-based, inclusive) | **in-a** `EigenSelection.IndexRange(from, to)` | Maps to `syevr`/`syevx`. Ascending-rank indices. Symmetric-only. |
| Subset by value interval | — | `eigh(A, subset_by_value=(vl, vu))` half-open `(vl, vu]` | **in-a** `EigenSelection.ValueInterval(lower, upper)` (`(lower, upper]`) | Count of results is data-dependent; diagnostics report how many fell in range. Symmetric-only. |
| Eigenvector layout (values vs matrix) | `eig(A,'vector')` / `'matrix'` | always `(w, V)` | **in-a** — always `(DVec, DMat)` | gale returns values as `DVec`, vectors as `DMat` columns; no `'matrix'` diagonal form. |
| Driver / algorithm choice | (internal) | `driver='ev'\|'evd'\|'evr'\|'evx'` | **out** (v0.3.5) | Single kernel (tridiag QL/QR). Driver selection is an internal optimization, not public surface. |
| Lower-vs-upper stored triangle | — | `lower=True\|False` | **in-a** (implicit) | gale reads one triangle by assumption, mirroring the `Cholesky` precedent (lower triangle only). Expose which triangle via option if needed; default lower. |

**Ordering guarantee.** gale returns symmetric eigenvalues **ascending-algebraic
always**, regardless of the selection order. The `order` in
`EigenSelection.Count` decides *which* `k` are returned, never the output order —
identical to `eigh`/`eigsh`, and the least-surprising choice. `LargestAlgebraic`
and `SmallestAlgebraic` mean the ordinary order on the real line;
`LargestMagnitude`/`SmallestMagnitude` are also valid selectors for symmetric
inputs (they differ from algebraic only when the spectrum straddles zero).

**Phase-a scope decision.** The planned tridiagonal QL/QR kernel *is* a full
symmetric eigensolver, and the partial-symmetric path depends on it for its
projected (tridiagonal) problem. We therefore **expose full dense
`eigSymmetric` in phase a** even though the PRD lists "Full symmetric
eigendecomposition" under *Later*; dense subset selection then reduces to slicing
the full ascending spectrum. A true subset driver (bisection + inverse iteration,
`syevr`-style) is a later performance optimization, not a phase-a blocker.

---

## 2. Dense nonsymmetric eigen (`eig`)

LAPACK: `geev` (values + optional left/right vectors), preceded by Hessenberg
reduction (`gehrd`) and Francis QR (`hseqr`). Real inputs can have complex
eigenvalues in conjugate pairs.

| Capability | MATLAB | SciPy | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| Eigenvalues only | `eig(A)` (order **not** guaranteed) | `eig(A, right=False)` → `w` (unordered) | **in-b** `eig(A, EigenVectors.ValuesOnly)` | Complex outputs via SoA `(re, im)`; typed `eigenvalue(i): Complex`. |
| Right eigenvectors | `[V,D]=eig(A)`, `A V = V D` | `eig(A)` → `(w, vr)` | **in-b** `EigenVectors.Right` | Real `DMat` in real-Schur packing; `eigenvector(i)` returns a complex pair view for complex eigenvalues. |
| Left eigenvectors | `[V,D,W]=eig(A)`, `Wᴴ A = D Wᴴ` | `eig(A, left=True)` → `(w, vl, vr)` | **shipped (dense)** — conjugated rows of `V⁻¹`, unit 2-norm, `wᴴA = λwᴴ`; defective/near-defective `a` → `Left(SingularMatrix)` (deliberate divergence from `dgeev`, which returns near-parallel vectors). Iterative path still rejects (`UnsupportedOperation`) — a one-sided Krylov basis gives no left vectors. | gale's convention: unit 2-norm, not biorthonormal (`wᴴv = 1` is NOT forced). |
| Skip balancing | `eig(A,'nobalance')` | (balancing internal to `geev`) | **out** (v0.3.5) | Balancing is on by default internally; `'nobalance'` is a niche numerical knob, deferred. |
| Ordering of results | not documented | not documented | **in-b** — gale-canonical | gale imposes a deterministic order the ecosystems lack (see below). |

**Ordering guarantee.** Complex spectra have no natural total order, so gale
imposes one: sort by the selection criterion (default `LargestMagnitude` ⇒
descending `|λ|`), tie-broken by descending real part then descending imaginary
part. **Conjugate pairs are kept adjacent with the positive-imaginary member
first**, which is mandatory given the real-Schur SoA packing — `isRealPair(i)` /
the pair accessors depend on that adjacency.

---

## 3. Dense SVD (`svd`)

LAPACK: `gesvd` / `gesdd`. `A = U Σ Vᵀ`, singular values descending.

| Capability | MATLAB | SciPy | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| Singular values only | `svd(A)` (descending) | `svdvals(A)` / `svd(A, compute_uv=False)` | **in** (v0.5) | Dense Golub–Kahan–Reinsch kernel, values-only assembly. |
| Full SVD | `[U,S,V]=svd(A)` (`U` m×m, `V` n×n) | `svd(A, full_matrices=True)` → `(U, s, Vh)` | **in — economy factors** (v0.5) | `Svds.svd(a, All)` / `DMat.svd`; gale returns economy `U` m×k, `Vᵀ` k×n (k = min(m,n)), not full square factors. |
| Economy / thin SVD | `svd(A,'econ')`, `svd(A,0)` | `svd(A, full_matrices=False)` | **in** (v0.5) | The shipped shape. |
| LAPACK driver choice | (internal) | `lapack_driver='gesdd'\|'gesvd'` | **out** | Not public surface. |
| Partial SVD (top/bottom `k`) | `svds` | `svds` | **in-a** | See § 6 — this is the SVD gale actually ships in v0.3.5. |

**Scope note.** v0.3.5's SVD deliverable is **partial** (§ 6). Full dense SVD is
deferred because it needs Householder bidiagonalization + an implicit-QR
**bidiagonal SVD** kernel, which is *not* in the phase-a kernel plan. Partial SVD
deliberately avoids that kernel: Lanczos bidiagonalization produces a small
bidiagonal `B_k`, and its SVD is obtained by running the **existing tridiagonal
QL/QR kernel** on the **Jordan–Wielandt augmented tridiagonal** — the
perfect-shuffle permutation of `[[0, B_kᵀ], [B_k, 0]]`, whose eigenvalues are
`±σ_i`. The augmented form is essential: the normal-equation alternative
(`T = B_kᵀ B_k`, eigenvalues `σ²`) squares the condition number and destroys the
smallest singular values, which § 8 ships as `SingularOrder.Smallest`. So partial
SVD reuses the symmetric kernel and needs no new dense-SVD machinery — with the
caveat that a generic symmetric QL/QR on the augmented form delivers absolute
accuracy `~ε‖A‖`, not the high *relative* accuracy of a dedicated zero-shift/dqds
bidiagonal kernel (`bdsqr`/`bdsdc`); acceptable for the PRD's portable-correctness
scope, and documented as such.

---

## 4. Generalized symmetric-definite eigen (`A x = λ B x`, `B` SPD)

LAPACK: `sygv`/`sygvd`/`sygvx`. Reduces to standard symmetric via Cholesky of `B`
(`B = LLᵀ`, solve `L⁻¹ A L⁻ᵀ y = λ y`, `x = L⁻ᵀ y`). Real eigenvalues, `B`-orthonormal vectors.

| Capability | MATLAB | SciPy | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| Sym-definite generalized eigenpairs | `eig(A,B,'chol')` (auto when `A` sym, `B` SPD) | `eigh(A, B)` (type 1) | **shipped** `eigSymmetricGeneralized(A, B, selection)` | Cholesky-of-`B` reduction onto the phase-a symmetric kernel; `B`-orthonormal vectors; error amplifies ~`κ(B)` (documented). |
| Problem-type variants | — | `eigh(A, B, type=1\|2\|3)` (`ABx=λx`, `BAx=λx`) | **out** (v0.3.5) | PRD names only `A x = λ B x` (type 1). Types 2/3 deferred. |
| Subset by index / value | — | `subset_by_index` / `subset_by_value` with `gvx` driver | **in-b** (via `EigenSelection`) | Same selection ADT as § 1; ascending-algebraic output. |
| `B` not positive-definite | error | `LinAlgError`-equivalent | **`Left(NotPositiveDefinite)`** | Cholesky of `B` fails ⇒ same `Left` the dense `Cholesky` path already returns. |

---

## 5. Generalized nonsymmetric eigen / QZ (`A x = λ B x`, general)

LAPACK: `ggev` (values + vectors), `gges` (generalized Schur / QZ), `tgsen`
(reorder). Eigenvalues are **projective** `(α:β)`, `λ = α/β`; `β = 0` ⇒ infinite
eigenvalue (singular pencil / rank-deficient `B`).

| Capability | MATLAB | SciPy | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| Generalized eigenvalues `(α, β)` | `eig(A,B)` / `[AA,BB]=qz` diag | `eig(A,B)`; `eig(...,homogeneous_eigvals=True)` → `(α,β)` | **backend** — `GeneralizedEigenvalue(alpha: Complex, beta: Double)` | Result type defined in v0.3.5; **must** carry `(α,β)` to represent `∞`. |
| Left/right generalized vectors | `[V,D,W]=eig(A,B)` | `eig(A, B, left=True)` | **backend** | Same packing story as § 2. |
| Generalized Schur (QZ) | `[AA,BB,Q,Z]=qz(A,B)`, `'real'`/`'complex'` | `qz(A, B, output='real'\|'complex')` | **backend / deferred** | Public Schur form is PRD "Later"; the QZ backend can expose it opportunistically. |
| Reorder / cluster selection | `ordqz(...)` incl. `'lhp'/'rhp'/'udi'/'udo'` | `ordqz(..., sort='lhp'/'rhp'/'iuc'/'ouc'\|callable)` | **out** (v0.3.5) | Eigenvalue-cluster reordering (`tgsen`) is out; note the region-keyword parity for a future design. |
| Infinite eigenvalue handling | `β=0` in `BB` diagonal | `homogeneous_eigvals` | **backend** — `isInfinite`, `isFinitePair` accessors | Central design consequence — see finding #1. |

**Backend boundary.** Pure QZ (Hessenberg-triangular reduction + generalized
Francis) is heavy and out of v0.3.5's pure scope. gale defines the **API and
result type** in v0.3.5 and routes the computation to an optional backend; with
no backend present the call returns `Left(UnsupportedOperation("generalized nonsymmetric eigen"))`,
mirroring the existing `UnsupportedOperation` precedent.

---

## 6. Partial symmetric eigen (`eigsh`)

Block Krylov with soft locking, thick restarting, and full reorthogonalization.
SciPy wraps ARPACK (IRLM); MATLAB `eigs` replaced ARPACK with Krylov–Schur in
R2017b. Constraint: `k < n`. The initial `ncv` satisfies `k < ncv ≤ n`, defaults
to `min(n, max(2k+1, 20))`, and grows toward `n` across non-converged restarts.

| Capability | MATLAB `eigs` | SciPy `eigsh` | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| `k` largest magnitude (default) | `eigs(A,k)` = `'largestabs'` (default) | `which='LM'` (default) | **in-a** `Count(k, LargestMagnitude)` | Defaults agree: LM. |
| `k` smallest magnitude | `'smallestabs'` | `which='SM'` | **in-a** `Count(k, SmallestMagnitude)` | `SM` without shift-invert is slow/fragile in both — gale documents preferring shift-invert around 0. |
| `k` largest/smallest algebraic | `'largestreal'`/`'smallestreal'` (≡ algebraic for symmetric) | `which='LA'`/`'SA'` | **in-a** `Count(k, LargestAlgebraic\|SmallestAlgebraic)` | For real symmetric, real-part = algebraic. |
| Both ends | `'bothendsreal'` | `which='BE'` (⌈k/2⌉ high, ⌊k/2⌋ low) | **in-a** `Count(k, BothEnds)` | **New enum case** — not in the PRD's `EigenOrder`. |
| Shift-invert near σ | `sigma` scalar (auto-factorizes `A-σB`) | `sigma=σ, mode='normal'` (+`'buckling'`/`'cayley'`) | **in-a** `SpectralTarget.ShiftInvert(σ, plan)` (σ **real**) | gale requires an **explicit** `LinearSolvePlan`; it never auto-factorizes. `'buckling'`/`'cayley'` modes **deferred**. |
| Generalized `A x = λ B x` | `eigs(A,B,k)` | `eigsh(A, k, M=B)` (`B` SPD) | **in-b** | Needs `B`-inner-product Lanczos; grouped with generalized eigen (phase b). |
| Matrix-free operator | `eigs(Afun,n,...)` | `LinearOperator` | **in-a** | gale already has `DoubleLinearOperator` with `applyTo`. |
| Start vector / subspace / tol / maxiter | `StartVector`, `SubspaceDimension`, `Tolerance`, `MaxIterations` | `v0`, `ncv`, `tol`, `maxiter` | **in-a** `SpectralOptions` | The caller vector is the first block column; deterministic orthogonal probes fill a block at least `k` wide. Default initial `ncv = min(n, max(2k+1, 20))`. |
| Convergence reporting | `[V,D,flag]` (`flag` 0/1) | raises `ArpackNoConvergence` (carries partial results) | **in-a** `Right(result + SpectralDiagnostics)` | See § Convergence. Never a `Left` for non-convergence. |

**Ordering guarantee.** As in § 1, partial symmetric results are returned
**ascending-algebraic**, matching `eigsh` (which sorts ascending; MATLAB `eigs`
does not rigorously document its within-set order). This is deliberately stronger
than either ecosystem.

**Multiplicity guarantee.** `Count(k, order)` counts algebraic multiplicity.
Repeated Ritz roots retain independent orthonormal directions from their
eigenspace; convergence and tests are defined by per-pair residuals and the
returned subspace projector, never by a particular basis inside a repeated
eigenspace. A result cannot report `allConverged` unless all `k` requested pairs
meet the residual criterion. For a matrix-free partial solve, that residual
criterion establishes convergence only within the explored Krylov space; it does
not prove that an invariant starting subspace contains the requested global
extreme. `diagnostics.convergenceStatus` distinguishes `ResidualConverged` from
`ExtremeCertified`; the latter currently requires a full-space reduction. Exact
Krylov breakdown replenishes the orthogonal complement rather than silently
substituting another distinct root.

---

## 7. Partial nonsymmetric eigen (`eigs`)

Arnoldi with implicit restarting — SciPy wraps ARPACK (IRAM); MATLAB uses
Krylov–Schur since R2017b. Constraint: `k < n-1`.

| Capability | MATLAB `eigs` | SciPy `eigs` | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| `k` largest/smallest magnitude | `'largestabs'`(default)/`'smallestabs'` | `which='LM'`(default)/`'SM'` | **in-b** `Count(k, LargestMagnitude\|SmallestMagnitude)` | Complex outputs via SoA packing. |
| `k` largest/smallest real part | `'largestreal'`/`'smallestreal'` | `which='LR'`/`'SR'` | **in-b** `Count(k, LargestRealPart\|SmallestRealPart)` | — |
| `k` largest/smallest imag part | `'largestimag'`/`'smallestimag'` | `which='LI'`/`'SI'` | **in-b** (`LargestImaginary`/`SmallestImaginary`) or **deferred** | **New enum cases** absent from PRD; low priority (see finding #2). |
| Both ends (real / imag) | `'bothendsimag'` (nonsymmetric); `'bothendsreal'` availability differs by problem type | (no direct `BE` for `eigs`) | **deferred** | Asymmetry: `eigsh` has `BE`, `eigs` does not. |
| Shift-invert near σ | `sigma` scalar | `sigma=σ` (+ `Minv`/`OPinv`) | **in-b** `ShiftInvert(σ, plan)`, σ **real** | Complex σ (targeting off the real axis) would need complex solves — **out** (finding #5). |
| Convergence reporting | `[V,D,flag]` | `ArpackNoConvergence` (partial) | **in-b** `Right + diagnostics` | Same rule as § 6. |

**Ordering guarantee.** Same canonical nonsymmetric order as § 2 (by selection
criterion; conjugate pairs adjacent, positive-imag first). SciPy `eigs` leaves
order unspecified — gale is stricter.

---

## 8. Partial SVD (`svds`)

| Capability | MATLAB `svds` | SciPy `svds` | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| `k` largest singular values (default) | `svds(A,k)` = `'largest'` | `which='LM'` (default), `solver='arpack'/'lobpcg'/'propack'` | **in-a** `SingularSelection.Count(k, Largest)` | Default agrees: largest. |
| `k` smallest singular values | `'smallest'` | `which='SM'` | **in-a** `Count(k, Smallest)` | Fragile with normal-equation methods; gale documents accuracy limits, prefers Lanczos bidiag. |
| Smallest **nonzero** | `'smallestnz'` | — (via `which='SM'` + rank handling) | **out**; handled via rank reporting | gale reports exact/near-zero singular values through a `rank` diagnostic rather than a separate order mode. |
| Near σ | `sigma` scalar | — (`svds` has no `sigma`) | **deferred** | MATLAB-only; low demand. |
| Start vector / subspace / tol / maxiter | `LeftStartVector`/`RightStartVector`, `SubspaceDimension`, `Tolerance`, `MaxIterations` | `v0`, `ncv`, `tol`, `maxiter` | **in-a** `SpectralOptions` | Single `startVector` (right); left derived. |
| Vectors optional | `[U,S,V]` vs `S` | `return_singular_vectors=True\|'u'\|'vh'\|False` | **in-a** — values-only vs full | gale: values-only or full `(U, s, Vᵀ)`; no one-sided `'u'`/`'vh'` mode in v0.3.5. |
| Convergence reporting | `[U,S,V,flag]` | `ArpackNoConvergence` | **in-a** `Right + diagnostics` | Same rule. |

**Ordering guarantee (critical inconsistency).** **SciPy `svds` returns singular
values ascending; MATLAB `svds` returns them descending; both dense `svd`s return
descending.** gale returns **descending always** — consistent with dense SVD
everywhere and with MATLAB. `Smallest` selects the `k` smallest and still reports
them descending. This single choice resolves a real cross-ecosystem trap.

---

## 9. Generalized SVD (`gsvd`)

Shape: `A = U C Xᵀ`, `B = V S Xᵀ`, `Cᵀ C + Sᵀ S = I`, generalized singular values
`c_i / s_i`. `X` is square and **not** orthogonal. LAPACK: `ggsvd3`
(Kogbetliantz-style: `ggsvp3` preprocessing + `tgsja`; the cosine–sine
decomposition is the mathematical relative, not the LAPACK algorithm).

| Capability | MATLAB | SciPy | gale v0.3.5 plan | Notes |
|---|---|---|---|---|
| Generalized singular values | `gsvd(A,B)` (vector; **ascending**) | **none** (no high-level `gsvd`) | **in-b** (values), **descending** | SciPy asymmetry — see below. gale normalizes to descending like § 8. |
| Full GSVD factors | `[U,V,X,C,S]=gsvd(A,B)` | `cossin` (CS decomposition) / `lapack.?ggsvd3` (low-level only) | **in-b** (full-rank), **deferred/backend** (rank-deficient) | Pure impl targets full-column-rank `[A;B]`; hard rank-deficient cases behind backend. |
| Economy form | `gsvd(A,B,0)` | — | **deferred** | — |
| Zero / infinite / rank-deficient generalized values | implicit in `C`,`S` diagonals | — | **in-b** — typed `GeneralizedSingularValue` (`Finite(c/s)`, `Infinite` when `s=0`, `Zero` when `c=0`) | PRD explicitly requires "explicit zero, infinite, and rank-deficient cases". |

**SciPy asymmetry (accuracy note).** **SciPy ships no first-class `gsvd`.** The CS
decomposition is available as `scipy.linalg.cossin`, and `?ggsvd3` only as a
low-level `scipy.linalg.lapack` wrapper. MATLAB has first-class `gsvd`. gale's
GSVD therefore validates primarily against **MATLAB and LAPACK `ggsvd3`
directly**, not against a SciPy high-level reference. MATLAB's vector-form
`gsvd(A,B)` is **ascending** (inconsistent with its own `svd`, which is
descending) — one more reason gale fixes its own canonical order.

---

## Selection & ordering semantics

### The selection modes across ecosystems

Collapsing every ecosystem knob into their underlying intents:

- **By count + criterion.** "Give me `k` extreme values." Criterion ∈ {largest
  magnitude, smallest magnitude, largest algebraic, smallest algebraic, largest
  real, smallest real, largest imag, smallest imag, both ends}. Source: MATLAB
  `eigs`/`svds` `sigma` words; SciPy `eigs`/`eigsh`/`svds` `which`.
- **By index range.** "Give me the eigenvalues ranked `il..iu` in ascending
  order." Symmetric only. Source: SciPy `eigh(subset_by_index=…)`; LAPACK
  `syevr`/`syevx`. **No MATLAB dense equivalent.**
- **By value interval.** "Give me every eigenvalue in `(vl, vu]`." Count is
  data-dependent. Symmetric only. Source: SciPy `eigh(subset_by_value=…)`;
  LAPACK `syevr`/`syevx`.
- **By target / shift-invert.** "Give me the values nearest σ." Source: MATLAB
  scalar `sigma`; SciPy `sigma=…`.
- **All.** Dense full decomposition.

### Proposed gale typed equivalents

```scala
enum EigenOrder:
  case LargestMagnitude, SmallestMagnitude     // legal: symmetric & nonsymmetric
  case LargestAlgebraic, SmallestAlgebraic     // legal: symmetric only (real-line order)
  case LargestRealPart,  SmallestRealPart      // legal: nonsymmetric only
  case BothEnds                                 // NEW (eigsh 'BE'); legal: symmetric only
  // deferred (nonsymmetric, low demand): LargestImaginary, SmallestImaginary

enum SingularOrder:
  case Largest, Smallest

// Unifies count/index/interval/all. Range & Interval are symmetric-only
// (validated at the boundary → Left(InvalidArgument) if used on a
// nonsymmetric problem).
enum EigenSelection:
  case All
  case Count(k: Int, order: EigenOrder)
  case IndexRange(from: Int, to: Int)            // ascending-algebraic rank, inclusive
  case ValueInterval(lower: Double, upper: Double) // half-open (lower, upper]

enum SingularSelection:
  case All
  case Count(k: Int, order: SingularOrder)

enum SpectralTarget:
  case Around(value: Double)                     // real only in v0.3.5
  case ShiftInvert(sigma: Double, plan: LinearSolvePlan)  // real σ; explicit solve
                                                 // (PRD's `solver` field renamed `plan`)

enum EigenVectors:
  case ValuesOnly, Right, Left, LeftAndRight     // symmetric API restricts to ValuesOnly | Vectors
```

The PRD's count-only `SpectralSelection(count, order, target)` cannot express
index-range or value-interval selection; **`EigenSelection` supersedes it**. The
`target` moves into `SpectralOptions`/`SpectralTarget` so that shift-invert
composes with any selection.

**Order legality is validated at the boundary**, like `IndexRange`/`ValueInterval`:
an `EigenOrder` used on the wrong problem type (`BothEnds` or the algebraic
orders on a nonsymmetric problem; the real-part orders on a symmetric problem)
returns `Left(InvalidArgument)` rather than silently reinterpreting the intent.
The magnitude orders are legal everywhere.

### Ordering guarantees gale commits to (documented, deterministic)

MATLAB and SciPy are mutually inconsistent (svds ascending vs descending; eigsh
ascending vs eigs unspecified; gsvd ascending vs svd descending). gale therefore
publishes its own guarantees:

| Result | Canonical output order |
|---|---|
| Symmetric eigen (dense & partial) | **ascending algebraic**, always |
| SVD (partial) | **descending** singular value, always |
| GSVD | **descending** generalized singular value; `Infinite` (s=0) first, `Zero` (c=0) last |
| Nonsymmetric eigen | by **selection criterion**; ties → desc. real, then desc. imag; **conjugate pairs adjacent, +imag first** |
| Generalized nonsym eigen | as nonsymmetric, comparing `α/β` with infinities (β=0) placed by criterion |

The selection `order` chooses *membership*; the table above fixes the *layout*.
This decoupling is the single most important ordering decision — it makes results
reproducible across backends and platforms (JVM / Scala.js / Wasm).

---

## Convergence & failure semantics

### How the ecosystems behave

- **MATLAB** returns a **flag** as an extra output: `[V,D,flag]=eigs(...)`,
  `[U,S,V,flag]=svds(...)`; `flag = 0` iff all requested pairs converged,
  otherwise `1`. Whatever converged is still returned. `FailureTreatment`
  (`'replacenan'|'keep'|'drop'`) controls how unconverged slots are filled.
- **SciPy sparse** raises **`ArpackNoConvergence`** — but the exception object
  carries `.eigenvalues` / `.eigenvectors` (or singular triplets) for the subset
  that *did* converge. So "failure" still ships partial results, just through the
  exception instead of the return value.
- **SciPy dense** (`eig`/`eigh`/`svd`) essentially always converges; a LAPACK
  `info > 0` surfaces as `LinAlgError`.

### gale's Either-first / diagnostics mapping

gale already has three precedents that dictate the rule:
`FactorizationDiagnostics.info` (non-fatal problem → number, result still
returned), `SolverResult.NotConverged` (best iterate returned, **not** a `Left`),
and `conditionEstimate` (singular → `Right(∞)`). Applying them:

**`Left(LinAlgError)` — structural / precondition violations (no meaningful
result):**

| Condition | Error |
|---|---|
| Non-square eigen input | `NonSquareMatrix` |
| `A`/`B` shape disagree (generalized) | `DimensionMismatch` |
| `k ≤ 0`, or `k ≥ n` (`eigsh`) / `k ≥ n-1` (`eigs`) | `InvalidArgument` (message suggests the dense API — no silent dense fallback, per PRD non-goal) |
| `IndexRange`/`ValueInterval` on a nonsymmetric problem | `InvalidArgument` |
| `EigenOrder` illegal for the problem type (`BothEnds`/algebraic on nonsymmetric; real-part on symmetric) | `InvalidArgument` |
| Generalized sym-definite with `B` not SPD | `NotPositiveDefinite` (Cholesky of `B` fails) |
| Shift-invert with `A-σB` (or `A-σI`) singular (σ at an eigenvalue) | `SingularMatrix` (from the inner solve) |
| Generalized nonsymmetric / QZ with no backend | `UnsupportedOperation` |

**`Right(result + SpectralDiagnostics)` — convergence & degeneracy (result, or
partial result, is meaningful):**

- **Residual convergence:** `Right(result)`, `diagnostics.allConverged = true`,
  every per-pair residual below tolerance. For iterative partial solvers this is
  explicitly convergence in the explored subspace, with
  `convergenceStatus = ResidualConverged`, not a universal extremality proof.
- **Extreme-certified convergence:** `convergenceStatus = ExtremeCertified`
  additionally proves membership in the requested spectral extreme. Gale
  currently issues this status only for full-space reductions.
- **Partial convergence:** `Right(result containing only the converged pairs)`,
  `diagnostics.allConverged = false`, `requested`/`converged` counts recorded,
  per-pair residuals recorded. This is the typed analogue of
  `ArpackNoConvergence`'s partial results and of `SolverResult.NotConverged` —
  **never a `Left`.**
- **Zero convergence:** `Right(empty result)`, `converged = 0`,
  `allConverged = false`. Consistent with always returning a value; the caller
  decides severity.
- **Rank deficiency (SVD/GSVD):** reported in `diagnostics.rank` and via typed
  `Zero`/`Infinite` generalized-value cases — a `Right`, mirroring
  `conditionEstimate`'s `Right(∞)`.

**Typed enforcement helpers.** `requireConverged` preserves the historical
residual policy: it returns `DidNotConverge(iterations, worstResidual)` when
`!allConverged` but accepts `ResidualConverged`. `requireExtremeCertified` is the
strict policy: it accepts only `ExtremeCertified`, returns `DidNotConverge` for
`NotConverged`, and returns the distinct typed error
`SpectralExtremeNotCertified` when residuals converged without a global
extremality certificate. The strict helper enforces a certificate already present
in the diagnostics; it does not create one.

**One-line rule:** *structural violations are `Left`; non-convergence and
degeneracy are `Right` + diagnostics, with explicit residual-only or
extreme-certified enforcement selected by the caller.*

---

## Implications for spectral-types (recommended definitions)

Each item traces to matrix rows; **bold** items are additions or corrections to
the PRD's sketch.

1. **`Complex(re: Double, im: Double)`** — boundary value type only (fixed
   constraint 1). Used by nonsymmetric eigen accessors (§ 2, § 7).

2. **`GeneralizedEigenvalue(alpha: Complex, beta: Double)`** with
   `value: Complex` (= α/β), `isInfinite: Boolean` (β≈0), `isFinitePair`. A plain
   `Complex` **cannot** represent the infinite eigenvalues produced by a singular
   / rank-deficient `B`; the projective `(α:β)` form is required (§ 5). *This is
   the highest-impact type decision.*

3. **`EigenSelection`** ADT — `All | Count(k, EigenOrder) | IndexRange(from, to) |
   ValueInterval(lower, upper]`. **Supersedes the PRD's count-only
   `SpectralSelection`**, which cannot express `eigh` subset-by-index /
   subset-by-value (§ 1, § 4). `IndexRange`/`ValueInterval` are symmetric-only.

4. **`EigenOrder`** — PRD's six cases **plus `BothEnds`** (§ 6, `eigsh 'BE'`).
   `LargestImaginary`/`SmallestImaginary` (§ 7) are deferred, low-priority.

5. **`SingularSelection`** = `All | Count(k, SingularOrder)`; **`SingularOrder` =
   `Largest | Smallest`.** No value/index subset (no ecosystem precedent worth
   cloning) (§ 3, § 8).

6. **`SpectralTarget`** = `Around(Double) | ShiftInvert(sigma: Double, plan:
   LinearSolvePlan)` — **σ is real only** in v0.3.5; shift-invert takes an
   **explicit** solve plan (gale never auto-factorizes, unlike MATLAB `eigs`)
   (§ 6, § 7).

7. **`SpectralOptions`** = `tolerance`, `maxIterations`, `subspaceDimension`
   (ncv/p, default `min(n, max(2k+1, 20))`), `startVector`, `returnVectors:
   EigenVectors` (§ 6–8).

8. **`EigenVectors`** = `ValuesOnly | Right | Left | LeftAndRight` (§ 2). The
   **symmetric API restricts to `ValuesOnly` | vectors** (left = right for
   symmetric — do not expose the 4-case enum there).

9. **`EigenDecomposition`** (symmetric): `eigenvalues: DVec` (ascending),
   `eigenvectors: DMat`, `diagnostics: SpectralDiagnostics` — complex-free.

10. **`NonsymmetricEigenDecomposition`**: SoA `(re: DVec, im: DVec)` + real `DMat`
    (real-Schur packing) + optional left vectors; accessors `eigenvalue(i):
    Complex`, `isRealPair(i)`, `eigenvector(i)`, `leftEigenvector(i)` (§ 2, § 7).

11. **`SVD`**: `singularValues: DVec` (descending), `u: DMat`, `vt: DMat`,
    **`rank: Int`** (numerical rank at tolerance — addition to the PRD's `SVD`
    sketch), `diagnostics` (§ 8).

12. **`GeneralizedSingularValue`** = `Finite(ratio: Double) | Infinite | Zero`;
    `GSVD` result carries the `(U, V, X, C, S)` factors and ordered generalized
    values (§ 9).

13. **`SpectralDiagnostics`** — the diagnostics contract that makes non-convergence
    a `Right`: `requested: Int`, `converged: Int`, `allConverged: Boolean`,
    `residuals: DVec` (per-pair), `orthogonalityError: Double`, `iterations: Int`,
    `rank: Option[Int]`, `extremalityCertified: Boolean`, and derived
    `convergenceStatus: NotConverged | ResidualConverged | ExtremeCertified`.
    Plus residual-based `requireConverged` and strict
    `requireExtremeCertified: Either[LinAlgError, Self]` on each result. Directly
    mirrors `FactorizationDiagnostics` + `SolverResult` while making the stronger
    spectral-membership claim explicit.

### Explicitly OUT of v0.3.5 (with rationale)

- **Complex inputs** — fixed constraint; no complex storage/kernel tier. Complex
  appears only in outputs.
- **Complex shift σ** (targeting off the real axis, MATLAB/SciPy allow it) —
  needs complex linear solves; no complex kernels ⇒ out.
- **Full dense SVD / economy SVD** — needs a bidiagonal-SVD kernel absent from the
  phase-a plan; PRD "Later". (Partial SVD ships via the Jordan–Wielandt augmented
  tridiagonal — absolute-accuracy tradeoff noted in § 3.)
- **`eigh` problem types 2 & 3** (`ABx=λx`, `BAx=λx`) — PRD names only type 1.
- **QZ reordering / `ordqz`** (`tgsen`) — eigenvalue-cluster reordering; out.
- **`'buckling'` / `'cayley'` shift-invert modes** (`eigsh`) — specialized
  structural-mechanics modes; plain shift-invert only.
- **`'nobalance'`, driver selection, one-sided singular-vector modes** — internal
  or niche knobs, not public surface.
- **LOBPCG / PROPACK solver families** — SciPy `svds`/`eigsh` alternative engines;
  gale uses block Krylov/Arnoldi. Reference only.
- **Polynomial / quadratic eigenproblems (`polyeig`)** — out of the linear
  spectral scope entirely.
- **Nonsymmetric Krylov–Schur sophistication** — the symmetric path now uses a
  thick-restarted block method; a corresponding Arnoldi/Krylov–Schur upgrade is
  still a robustness goal, not part of the public type surface.
- **`smallestnz` singular selection** — folded into rank reporting rather than a
  distinct order mode.
- **Pure generalized nonsymmetric (QZ) and hard rank-deficient GSVD** — API +
  result types defined; computation behind an optional backend
  (`UnsupportedOperation` with no backend).
