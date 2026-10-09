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

CI runs `parityTest` and `interopBreezeTest` as separate timed steps on the
`breeze-interop` job (`timeout-minutes: 20`). Both modules fork the test JVM
to isolate Breeze initialization from sbt. The timeouts bound stalled runs;
netlib still probes native implementations and falls back to Java when they
are unavailable. These settings do not force a particular BLAS/LAPACK backend.

`EverydayOpsParitySuite` uses ScalaCheck to vary matrix shapes and data seeds.
`ReductionsNumericsParitySuite` and `SparseVectorParitySuite` do the same for
sizes 1 to 64, plus fixed 1K and 64K cases. To replay a printed failing seed,
run the failing suite alone:
`sbt -Dgale.parity.seed=<seed> "parity/testOnly gale.parity.<Suite>"`. Where Gale's documented
semantics differ from Breeze, these suites check both results. The
[Breeze migration guide](../docs/user/guides/breeze-equivalence.md) lists those
cases.
The factorization and spectral suites use fixed adversarial and
well-conditioned fixtures. A parity test should state the shared mathematical
contract and use a tolerance that accounts for the algorithms being compared.

## Cross-platform Breeze goldens

Breeze runs only on the JVM, so the parity suites above say nothing about
Gale's Scala.js numerics. The golden corpus closes that gap. A JVM generator
runs Breeze over a fixed seeded corpus and writes the inputs and Breeze's
results into Gale's shared test sources. `BreezeGoldenSuite` then replays them
through Gale's public API on both the JVM and Scala.js, with no Breeze
dependency.

| File | Role |
| --- | --- |
| `parity/src/test/scala/gale/parity/GenerateBreezeGoldens.scala` | Generator, and the freshness check |
| `core/shared/src/test/scala/gale/golden/BreezeGoldens.scala` | Generated corpus (do not edit) |
| `core/shared/src/test/scala/gale/golden/GoldenData.scala` | Case types and the hex decoder |
| `core/shared/src/test/scala/gale/golden/GoldenTolerance.scala` | Error units shared by the replay and the freshness check |
| `core/shared/src/test/scala/gale/golden/BreezeGoldenSuite.scala` | Replay and tolerance policy |

```sh
sbt breezeGoldens        # regenerate the corpus after changing the generator or Breeze
sbt breezeGoldensCheck   # regenerate in memory; fail if the checked-in corpus is stale
sbt "coreJVM/testOnly gale.golden.BreezeGoldenSuite" "coreJS/testOnly gale.golden.BreezeGoldenSuite"
```

**Coverage.** The corpus has 69 cases (about 93,000 doubles, about 1.5 MB of
source). Sizes run from 1 to 64, with a few near 100. It covers:

- dense `A·x`, `Aᵀ·x`, `A·B`, `Aᵀ·D`, `dot` and `axpy`;
- LU solve, `det` and the inverse, on random, diagonally dominant and
  prescribed-κ matrices (κ = 1e2, 1e6, 1e10);
- the Cholesky factor `L` and solve, including Hilbert(8) and κ = 1e8;
- the sign-normalized QR `R` and least squares, including κ = 1e6;
- `eigSym` eigenvalues and eigenvectors (random, a repeated-eigenvalue
  spectrum, and a graded spectrum with κ = 1e10), plus values only at n = 96;
- singular values, the full-rank `pinv`, and `kron`;
- vector and matrix reductions: `sum`, `mean`, `max`, `min`, `argmax`,
  `argmin`, the 1-, 2- and ∞-norms, and the Frobenius norm;
- `exp`, `log`, `sigmoid`, `logSumExp` and `softmax`. Breeze's `softmax(v)` is
  log-sum-exp, and the generator builds the true softmax as
  `exp(v − softmax(v))`.
- `SparseVector` sparse and dense `dot`, `+`, `−`, `sum`, `max`, `min` and the
  norms.

**Encoding.** Every double is stored as its 16-hex-digit IEEE-754 bit pattern,
so both platforms decode exactly the same inputs and references. The generator
builds inputs with plain IEEE arithmetic, `sqrt` and `StrictMath` only, so the
inputs are bit-identical on every host. The data is split into string literals
and part objects so that no JVM constant or class initializer exceeds its size
limit.

**Tolerance policy.** Results are never compared bit-for-bit. Gale and Breeze
use different algorithms and summation orders, and on Scala.js
`PlatformMath.fma` is `a*b+c`. Every check asserts `error ≤ C · unit`, with a
single `C = 8`. There are two kinds of check.

- **Forward checks** compare gale with Breeze. The unit has the shape of the
  standard forward-error bound, `n · κ · ε · scale`, and is computed from the
  case for each reference array (`GoldenTolerance.forwardUnit`). Here `κ` is
  Breeze's 2-norm condition number for solves, inverses, `det`, Cholesky, QR
  and `pinv`, and 1 for products and reductions. Least squares adds the
  `κ² ‖r‖ / (‖A‖ ‖x‖)` term.
- **Backward checks** test gale's own result, without κ. They are
  `‖b − A x̂‖∞ ≤ C n ε (‖A‖∞ ‖x̂‖∞ + ‖b‖∞)` for LU and Cholesky solves,
  `‖A X − I‖∞ ≤ C n ε ‖A‖∞ ‖X‖∞` for the inverse, and
  `‖L Lᵀ − A‖max ≤ C n ε max aᵢᵢ` for Cholesky. Least squares uses the
  normal-equations residual
  `‖Aᵀ r̂‖∞ ≤ C m ε ‖A‖₁ (‖A‖∞ ‖x̂‖∞ + ‖b‖∞ + ‖r̂‖∞)`. On ill-conditioned
  inputs the forward bounds are loose. These checks catch a solver that loses
  digits there.

Some operations are exact or correctly rounded in both libraries: `max`, `min`,
`argmax`, `argmin`, vector `normInf`, `kron`, and sparse `+` and `−`. These
have a unit of 0 and must match exactly. Eigenvectors are compared through the
projector of each eigenvalue cluster, with unit `n ε ‖A‖ / gap`. A cluster
joins two neighbouring eigenvalues only when their gap is at most
`64 n ε ‖A‖` or at most `1e-8 |λ|`, so a graded spectrum is not merged into one
projector. The suite prints the worst observed `error / unit` for each family
and output on each platform. It lists the conditioned cases (κ, Hilbert,
graded) separately. Do not loosen `C` to hide a failure; a failure that occurs
only on JS is a numerical bug to report.

**Freshness.** Breeze's references depend on the host. netlib chooses native,
SIMD or scalar Java BLAS for the machine, and libm intrinsics differ by
architecture. A regenerated corpus is therefore byte-identical only on a host
like the one that generated it. `breezeGoldensCheck` therefore passes in
either of two cases:

- the regenerated source is byte-identical to the checked-in file; or
- the case list, the structure and every input bit are identical, and each
  reference moved by no more than 64 of the replay's own forward units (the
  replay allows 8). Eigenvectors may move by 64 projector units. The exact
  operations have a band of 0.

Any other change fails, and the failure message names the case. The
`breeze-interop` CI job runs the check after `parityTest`. To share the unit
definitions, the `parity` project depends on the core test classes
(`test->test`).

## Coverage checklist

Status values: **covered** (differential test present), **out** (no honest Breeze
reference or deliberate non-goal). **SciPy** means a NumPy/SciPy fixture is
the reference.

| Operation | Breeze API | Gale API | Suite | Status |
| --- | --- | --- | --- | --- |
| Dense ± / * / axpy / dot / scale | `+`, `-`, `*`, `dot` | same, including `A * α` | `DenseOpsParitySuite` | covered |
| Slice / strided-view products | `A(i until j, …) * x` | `slice` then `*` / `col` as `x` | `DenseOpsParitySuite` | covered |
| Reductions, per-axis reductions, norms | `sum`, `mean`, `max`, `argmax`, `sum(A(::, *))`, `norm(v, p)` | `sum`, `sumExact`, `mean`, `max`, `argmax`, `sum(Axis.Cols)`, `norm1`/`norm2`/`normInf`/`normFrobenius` | `ReductionsNumericsParitySuite` | covered (divergences pinned) |
| Elementwise and log-domain numerics | `breeze.numerics.*`, `softmax` (= log-sum-exp) | `Numerics.exp`/`log`/`log1p`/`expm1`/`sigmoid`/`logSumExp`/`softmax`/`logSoftmax` | `ReductionsNumericsParitySuite` | covered (divergences pinned) |
| Sparse vectors | `SparseVector`, `VectorBuilder`, `CSCMatrix * SparseVector` | `SparseVector`, `rowSparse`/`colSparse`, `CSR * SparseVector` | `SparseVectorParitySuite` | covered (divergences pinned) |
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
| L-BFGS / L-BFGS-B (x*, f*, residual, active set) | `LBFGS`, `LBFGSB` with an equivalent stopping rule | `LBFGS.minimize`, `LBFGSB.minimize` | `OptimizerParitySuite` | covered (timings and callback counts informational, in `parity/target/optimizer-parity.md`) |
| Generalized symmetric eigen | — | `Eigen.eigSymmetricGeneralized` | `GeneralizedSpectralParitySuite` | SciPy (`eigh(A, B)` type 1) |
| GSVD (full-column-rank) | — | `Svds.gsvd` | `GeneralizedSpectralParitySuite` | SciPy (Gram-pencil `eigh(AᵀA, BᵀB)`; no high-level `gsvd`) |
| QZ / generalized nonsymmetric | — | `Eigen.eigGeneralizedNonsymmetric` | `GeneralizedSpectralParitySuite` | covered (unsupported-contract lock; SciPy `qz` / `eig(A,B)` is the future target) |
| Sparse direct factorization | SuiteSparse / native | `SparseDirect` seam | `SparseDirectParitySuite` | SciPy (`splu` vs dense LU) + empty-provider lock |
| Near-cutoff rank / `pinv` / `cond` | policy-dependent | `rankEstimate`, `pinv`, `conditionEstimate` | `NearCutoffParitySuite` | SciPy (`pinv` MATLAB `rtol`, SVD rank, `cond(A, 1)` lower bound); Gale definitions in core `GaleNumericalContractSuite` |
| Cross-platform replay (JVM + Scala.js) | dense ops, factorizations, `eigSym`, `svd`, `pinv`, `kron`, reductions, numerics, `SparseVector` | same | `BreezeGoldenSuite` (core shared tests) | covered (generated Breeze goldens) |

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
