# Shipped vs deferred surfaces

This table describes the **current tree**, not the MATLAB/SciPy capability plan
in `docs/spectral-parity.md`. Use it to see which public calls have a portable
implementation, which require an imported backend, and which are locked to
`Left(UnsupportedOperation)` (or a documented stand-in) until work lands.

Lock tests live in `UnsupportedSurfaceSuite`. Gale-owned numerical definitions
that Breeze cannot referee (`rankEstimate`, `pinv` cutoff, `conditionEstimate`)
are pinned by `GaleNumericalContractSuite`.

| Surface | Status | Caller sees | Evidence |
| --- | --- | --- | --- |
| Dense `+`/`−`/`*`/`solve`/`lu`/`cholesky`/`qr`/`det` | portable | factors or `Either` | core factorization suites |
| Tall least squares | portable | `Either` (`RankDeficient` when rank drops) | `QRSuite` |
| Underdetermined least squares | deferred | `Left(UnsupportedOperation)` | `UnsupportedSurfaceSuite` |
| `rankEstimate` / `conditionEstimate` / `pinv` | portable (Gale contract) | QR rank; Hager `κ₁` (singular → `Right(+∞)`, rectangular → `Left(NonSquareMatrix)`); MATLAB/SciPy `pinv` cutoff | `GaleNumericalContractSuite` |
| Dense symmetric / nonsymmetric eigen | portable | `Right` + diagnostics; dense path is `ExtremeCertified` | dense eigen suites |
| Dense left eigenvectors | portable | `wᴴA = λwᴴ`; defective → `Left(SingularMatrix)` | `EigNonsymmetricLeftVectorSuite` |
| Iterative Lanczos / Arnoldi (no `target`) | portable | `Right` of residual-converged pairs; `requireExtremeCertified` is stricter | Lanczos / Arnoldi suites |
| Iterative `ShiftInvert` / `Around` | deferred | `Left(UnsupportedOperation)` | `UnsupportedSurfaceSuite` |
| Arnoldi `Left` / `LeftAndRight` | deferred | `Left(UnsupportedOperation)` | `UnsupportedSurfaceSuite` |
| Dense generalized symmetric-definite `eigSymmetricGeneralized(A, B)` | portable | `B`-orthonormal, `ExtremeCertified` | `EigSymmetricGeneralizedSuite`, `GeneralizedEigenTrustSuite` |
| Operator LOBPCG / generalized Lanczos | portable | residuals + `B`-Gram; extreme membership only when certified | `GeneralizedEigenTrustSuite` |
| Generalized nonsymmetric / QZ | backend-only | `Left(UnsupportedOperation)` with no QZ backend | `UnsupportedSurfaceSuite` |
| Full-column-rank GSVD | portable | `Right` | `GeneralizedSvdSuite` |
| Rank-deficient GSVD | backend-only | `Left(RankDeficient)` unless a capable backend is imported | `UnsupportedSurfaceSuite` |
| Sparse CSR/CSC arithmetic and matvec | portable | values / views | sparse suites, `parity/` |
| Sparse direct LU / Cholesky / QR | backend-only | `Left(UnsupportedOperation)` (no provider in this build) | `docs/sparse-direct-provider.md` |
| Vector / FFM BLAS backends | backend-only | same answers within conformance; not bit-identical | backend suites |
| Complex matrix storage | out | not an API | `docs/user/guides/breeze-equivalence.md` |

Status words:

- **portable** — implemented in `gale-core` with no acceleration import.
- **backend-only** — the public method exists; the default `given` declines.
- **deferred** — the public method exists and is locked to
  `UnsupportedOperation` until the named wiring lands.
- **out** — not part of the real-`Double` v1 surface.

When a deferred row starts succeeding, update this table and the lock suite in
the same change. Do not treat a green residual as a shipped extreme, and do not
treat a Breeze-overlapping fixture as Gale's cutoff policy.
