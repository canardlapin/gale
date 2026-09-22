# Sparse-direct: future-version note

Status: **desire only.** Do not implement these surfaces from this note.
Each item needs its own spec in a **future version** before any capability
is advertised.

Portable sparse Cholesky is already an explicit import
(`gale.sparse.direct.pure.given`). That engine is square, fixed-pattern,
and exact-zero / non-positive on pivots. The three items below each drop a
different one of those invariants. They are not leftover enum cases on the
same checklist. The technical split is in
[Scala.js sparse-direct](sparse-direct-js.md#why-later-items-are-different-projects).

| Future item | Why a new spec | Must not be treated as |
| --- | --- | --- |
| Sparse QR | Different graph (`AᵀA` / column intersection), Householder/Givens factors, `RankDeficient` and underdetermined-LS product decisions | `SparseDirectCapability.QR` on the Cholesky provider, Cholesky of `AᵀA`, or `toDense()` Householder |
| Threshold / delayed pivoting LU | Numeric row swaps change the `L`/`U` pattern and the exact-zero `SingularMatrix` policy; SuperLU/UMFPACK-class data structures | “Phase 2 static-sparsity LU plus a threshold” |
| Embedded C Wasm (CSparse / SuiteSparse) | Load/ABI/licence/memory project; different Wasm from `GALE_WASM=1`; sync `Either` only after `load()` | A `.wasm` blob in `gale-core`, or `GALE_WASM=1` as a sparse-direct accelerator |

When a future version wants one of these, write the spec first (contract,
errors, reuse, packaging, go/no-go gates), then implement. Do not advertise
the capability from a thin fallback in the meantime.
