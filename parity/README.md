# Breeze parity tests

The parity module runs Gale and Breeze 2.1.0 on the same numerical inputs. It
compares results that have the same mathematical meaning. For factorizations,
where signs or bases may differ, the tests compare reconstructions, residuals,
or subspace projectors instead of comparing raw factors.

Where Breeze has no honest public reference, the same module compares Gale to
checked-in NumPy / SciPy fixtures (the R counterparts are `geigen`,
`kappa` / `rcond`, `pracma::pinv`, `Matrix::lu`, and the Krylov solvers in
`Matrix`). CI does not need Python or R; regenerate fixtures from the
repository root with:

```sh
python3 parity/scripts/generate_numpy_references.py
```

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
reference or deliberate non-goal). **SciPy** means a NumPy/SciPy fixture is
the reference.

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
| Sparse + / − / scale | `CSCMatrix` arithmetic | CSR/CSC `+`, `-`, `*` | `BandedSparseParitySuite` | covered |
| Sparse identity / zero / permutation | eye / zeros / perm CSC | `Sparse.identity` / `zero` / `permutation` | `BandedSparseParitySuite` | covered |
| Sparse × dense product | `CSCMatrix * DenseMatrix` | `CSR * DMat` | `BandedSparseParitySuite` | covered |
| Sparse inspect (`apply` / row / col / `t` / trace / `toDense`) | CSC accessors | CSR/CSC accessors | `BandedSparseParitySuite` | covered |
| `pinv(A)*b` vs `A \\ b` (full-rank) | `pinv`, `\` | `pinv` then `*` | `FullSvdParitySuite` | covered |
| Symmetric eigen (dense + Lanczos) | `eigSym` | `Eigen.eigSymmetric` | `SpectralParitySuite` | covered |
| Nonsymmetric eigen | `eig` | `Eigen.eigNonsymmetric` | `NonsymmetricEigenParitySuite` | covered |
| Partial / full SVD, `pinv`, `kron` | `svd`, `pinv`, `kron` | `Svds.svd`, `pinv`, `kron` | `SvdQrParitySuite`, `FullSvdParitySuite` | covered |
| Blocked QR / lstsq | `qr`, `\` | `qr`, `leastSquares` | `SvdQrParitySuite` | covered |
| Iterative solve (solution equivalence) | dense `\` | `cg` / `bicgstab` / `gmres` / `lsqr` / `cgnr` | `IterativeSolveParitySuite` | covered (workload replaceability vs Breeze `\\`) |
| Iterative algorithm (Krylov diagnostics) | — | `cg` / `bicgstab` / `gmres` / `lsqr` | `IterativeAlgorithmParitySuite` | SciPy (`sparse.linalg`; solution + residual band + iteration band) |
| Generalized symmetric eigen | — | `Eigen.eigSymmetricGeneralized` | `GeneralizedSpectralParitySuite` | SciPy (`eigh(A, B)` type 1) |
| GSVD (full-column-rank) | — | `Svds.gsvd` | `GeneralizedSpectralParitySuite` | SciPy (Gram-pencil `eigh(AᵀA, BᵀB)`; no high-level `gsvd`) |
| QZ / generalized nonsymmetric | — | `Eigen.eigGeneralizedNonsymmetric` | `GeneralizedSpectralParitySuite` | covered (unsupported-contract lock; SciPy `qz` / `eig(A,B)` is the future target) |
| Sparse direct factorization | SuiteSparse / native | `SparseDirect` seam | `SparseDirectParitySuite` | SciPy (`splu` vs dense LU) + empty-provider lock |
| Near-cutoff rank / `pinv` / `cond` | policy-dependent | `rankEstimate`, `pinv`, `conditionEstimate` | `NearCutoffParitySuite` | SciPy (`pinv` MATLAB `rtol`, SVD rank, `cond(A, 1)` lower bound) |

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
