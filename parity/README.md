# Breeze parity tests

The parity module runs Gale and Breeze 2.1.0 on the same numerical inputs. It
compares results that have the same mathematical meaning. For factorizations,
where signs or bases may differ, the tests compare reconstructions, residuals,
or subspace projectors instead of comparing raw factors.

Run the suite from the repository root:

```sh
sbt parityTest
```

`EverydayOpsParitySuite` uses ScalaCheck to vary matrix shapes and data seeds.
The factorization and spectral suites use fixed adversarial and
well-conditioned fixtures. A parity test should state the shared mathematical
contract and use a tolerance that accounts for the algorithms being compared.

## Coverage checklist

Status values: **covered** (differential test present), **out** (no honest Breeze
reference or deliberate non-goal).

| Operation | Breeze API | Gale API | Suite | Status |
| --- | --- | --- | --- | --- |
| Dense ± / * / axpy / dot / scale | `+`, `-`, `*`, `dot` | same, including `A * α` | `DenseOpsParitySuite` | covered |
| Slice / strided-view products | `A(i until j, …) * x` | `slice` then `*` / `col` as `x` | `DenseOpsParitySuite` | covered |
| Construct / slice / gather / update / pointwise | indexing, `:*`, etc. | `slice`, `gather*`, `updated`, `pointwise`, `zipMapExact` | `EverydayOpsParitySuite` | covered |
| Vector zeros / fill / tabulate | `DenseVector.zeros/fill/tabulate` | `Vec.zeros/fill/tabulate` | `EverydayOpsParitySuite` | covered |
| det / solve / LU / Chol / QR / lstsq / inv | `det`, `\`, `lu`, `cholesky`, `qr`, `inv` | `det`, `solve`, `lu`, `cholesky`, `qr`, `leastSquares`, `solve(I)` | `FactorizationParitySuite` | covered |
| Cholesky / reused-QR solve paths | `\` on SPD / tall | `cholesky.solve`, `qr.solveLeastSquares` | `FactorizationParitySuite` | covered |
| rank / cond (clear cases) | `rank`, `cond` | `rankEstimate`, `conditionEstimate` | `FactorizationParitySuite` | covered (overlap only) |
| Sparse matvec / transpose-matvec | `CSCMatrix *` | Banded/CSR/CSC/COO/Diagonal `*` | `BandedSparseParitySuite` | covered |
| Identity / zero / permutation matvec | `CSCMatrix` structural | `Sparse.identity/zero/permutation` | `BandedSparseParitySuite` | covered |
| Sparse + / − / scale / sparse-dense `*` | `CSCMatrix` arithmetic | CSR/CSC `+`, `-`, `*`, CSR `* DMat` | `BandedSparseParitySuite` | covered |
| `pinv(A)*b` vs `A \\ b` (full-rank) | `pinv`, `\` | `pinv` then `*` | `FullSvdParitySuite` | covered |
| Symmetric eigen (dense + Lanczos) | `eigSym` | `Eigen.eigSymmetric` | `SpectralParitySuite` | covered |
| Nonsymmetric eigen | `eig` | `Eigen.eigNonsymmetric` | `NonsymmetricEigenParitySuite` | covered |
| Partial / full SVD, `pinv`, `kron` | `svd`, `pinv`, `kron` | `Svds.svd`, `pinv`, `kron` | `SvdQrParitySuite`, `FullSvdParitySuite` | covered |
| Blocked QR / lstsq | `qr`, `\` | `qr`, `leastSquares` | `SvdQrParitySuite` | covered |
| Iterative solve (solution equivalence) | dense `\` | `cg` / `bicgstab` / `gmres` / `lsqr` / `cgnr` | `IterativeSolveParitySuite` | covered (workload replaceability, not algorithm parity) |
| Generalized symmetric eigen | — | `Eigen.eigSymmetricGeneralized` | — | out (no Breeze public `eigh(A,B)`) |
| Sparse direct factorization | SuiteSparse / native | — | — | out (Gale non-goal) |
| Near-cutoff rank / `pinv` / `cond` | policy-dependent | policy-dependent | — | out (deliberate; pinned in core) |

Published conversion and migration shims live in `interop-breeze` (`sbt
interopBreezeTest`), not in this differential harness.

## Migration pain points

This table records cases found while writing parity tests where the Breeze
expression is materially easier to write or read. It is an API-design input,
not a claim that Gale should copy Breeze.

| Operation | Breeze | Gale | Why the Gale expression is harder | Status |
| --- | --- | --- | --- | --- |
| Matrix literal | `DenseMatrix((1.0, 2.0), (3.0, 4.0))` | `Matrix(2, 2)(1.0, 2.0, 3.0, 4.0)` | Gale requires the dimensions and a flattened row-major value list. The compiler cannot check that visual rows have equal lengths because the rows are not present in the expression. | Consider a row-based constructor that keeps `Matrix(rows, cols)` for generated or flattened input. |
| Contiguous matrix slice | `a(1 until 4, 2 until 5)` | `a.slice(1, 4, 2, 5)` | Two range values show which endpoints belong together. Gale's four adjacent integers are easier to transpose or misread. | Consider an overload that accepts row and column ranges. |
| Exception-style solve | `a \ b` | `a.solve(b).orThrow` | Gale returns `Either[LinAlgError, A]`, which preserves the failure type but adds ceremony in programs that deliberately use exceptions. | Keep the typed default. Evaluate a clearly named throwing convenience only if migrations repeatedly add local wrappers. |

Add a row only when the shorter Breeze form also makes the operation clearer.
Do not list deliberate differences merely because their spelling differs.
