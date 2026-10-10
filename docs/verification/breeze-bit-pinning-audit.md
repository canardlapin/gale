# Bit-pinning audit for W1 kernel reassociation (plan item W1.7)

- **Date:** 2026-10-08
- **Source revision:** `f4ff88d` (worktree `gale-breeze-w53`)
- **Scope:** test sources in `core/{shared,jvm,js}`, `laws/shared` (tests and the reusable
  `gale.laws` kits), `parity`, `backend-jvm-{vector,native,blas-ffm}`, `interop-breeze`, and
  `interop-ravel`. The audit covers assertions on the outputs of LU (pivots, det, solve,
  `SingularMatrix(k)`), Cholesky, triangular solve, gemv/matvec, and symmetric
  eigen/tridiagonalization.
- **Method:** grep over all 168 test files for `assertEquals` on doubles, sequences and arrays,
  `== `, `sameElements`, `doubleTo(Raw)LongBits`, `.toSeq` and `valuesRowMajor` comparisons,
  zero-tolerance `assertEqualsDouble(…, 0.0)`, typed-error indices, and exact iteration
  counts. Each hit was read in context. The kernel routes were traced in
  `linalg/Factorizations.scala`, `linalg/DenseCholeskyWorkspace.scala`,
  `kernel/DoubleKernels.scala`, `linalg/Matrix.scala:310-365` and
  `spectral/DenseSpectralKernels.scala`. No build or test run was made.

## Contract lines that permit reassociation

- `docs/user/advanced/numerical-contract.md:14-15`: "it does not guarantee exact real
  arithmetic or universal bit-for-bit agreement between legal algorithms."
- `numerical-contract.md:17-18`: "The pure single-threaded implementation is deterministic
  for a fixed Gale build and runtime." Determinism holds within a build, so a new build may
  change bits.
- `numerical-contract.md:37-39`: "Legal pivot, reflector-sign, eigenvector-sign, and
  repeated-eigenspace choices are not part of the identity contract."
- `numerical-contract.md:43`: "LU reports `SingularMatrix` only when a pivot is exactly `0`
  or `NaN`." A blocked LU must keep this rule.
- `docs/backend-architecture.md:318-319`: "conformance is **by law, not by bit-identity**".
  Lines 325-327: "Kernels **may reassociate** (the FMA/vectorized paths already do vs older
  gale), so equality is law-equivalence, and **cross-platform determinism is per-platform**".
- `backend-architecture.md:598-600`: the pure stages are validated "against the prior gale
  baseline (reassociation allowed; `parityTest` + `testAll` must stay green …)".

## Route facts that decide the handling

| Route | Kernel today | Relevance |
|---|---|---|
| `Factorizations.lu` (`:213-266`) | Inline right-looking elimination without fma. No other caller exists, and there is no LU workspace route. | W1.3 replaces this route above n of about 96. |
| `cholesky` (`:271-291`) and `DenseCholeskyWorkspace.factor`/`testPositiveDefinite` (`DenseCholeskyWorkspace.scala:40,128`) | Both call `DoubleKernels.dpotrfLower` (`DoubleKernels.scala:341`). | W1.4 must enter blocking through one function that both call. |
| LU `solve(DVec)` / `solve(DMat)` (`:1183`, `:1214-1271`) | `dtrsv` on each column, so the matrix and vector results are bit-identical today. | W1.2 moves the matrix route to `dtrsmLeft`. |
| Cholesky `solve(DVec)` (`:1273`) and the workspace solve (`DenseCholeskyWorkspace.scala:97-98`) | `dtrsv` | W1.2 does not touch these routes. |
| Cholesky `solve(DMat)` (`:1297-1339`) | Inline loops with `value -= l*x` and no fma. | W1.2 changes these bits at every n if `dtrsm` uses fma. |
| `TriangularSolve.lower/upper` | `dtrsv` | W1.2 does not touch this route. |
| `symmetricEigen`, `symmetricEigenWith` and `tridiagonalize` (`DenseSpectralKernels.scala:98,174,193`) | All three call `tred2` (`:226`). | W1.5 must keep the single shared call. |
| `DMat * DVec` (`Matrix.scala:323,339,355`) and `PureBackend.gemv` (`Backend.scala:150`) | `dgemvRowMajor`, `dgemvColMajor` or `dgemv` | Since W1.1 (`breeze/hc-gemv`), `PureBackend.gemv` selects its kernel by layout in the same order as `DMat.*`. The two routes are therefore bit-identical, which `BackendContractSuite` pins exactly. Before W1.1 the pure backend always called the generic `dgemv`. |

## Findings

| File:line | Op | What is pinned | Broken by | Handling | Justification |
|---|---|---|---|---|---|
| `core/shared/.../linalg/LUSuite.scala:28` | LU pivots | Integer pivots `Seq(1,0)`, n=2 | none (n < 96) | KEEP | Integer pivot result below the W1.3 threshold. W1.3 also requires the `pivots(i)` = original-row convention to survive. |
| `LUSuite.scala:49` | LU | `Left(SingularMatrix(1))`, 2×2 IEEE-exact rank-1 | none | KEEP | An exact zero pivot in any order: partial pivoting picks row `(2, 4)`, the multiplier is `0.5`, and `2 - 0.5·4 = 0` exactly. The index is an integer below threshold. |
| `linalg/NumericalPolicySuite.scala:52-54, 58, 62-65` | LU | Error class on the 3×3 `exactRank1` input, `SingularMatrix(0)` for `zeros(2,2)`, and the exact-zero-row case | none | KEEP | IEEE-exact plants (`:20-22` documents this). Small n. |
| `NumericalPolicySuite.scala:68-75` | LU | A 2×2 input with `4+1e-12` must factor (`Right`) | none | KEEP | A structural decision with a large margin. Small n. |
| `NumericalPolicySuite.scala:82`, `linalg/GaleNumericalContractSuite.scala:178` | LU → conditionEstimate | `+∞` for the 3×3 exact outer product | none | KEEP | Depends on an exact-zero pivot (contract `:43`). 3×3 is below threshold. |
| `GaleNumericalContractSuite.scala:153,155,164` | LU matrix solve (`exactCond1`, `:63`) and Hager | `est <= exact·(1+1e-8)` and `\|est-4\|<1e-9` | W1.2 (matrix solve ulps) | KEEP (already tolerance) | The relative slack of 1e-8 is far above the ulp shifts. No change needed. |
| `laws/shared/src/main/.../BackendConformanceSuite.scala:264-280, 322-324` | LU, Cholesky | Error *class* parity with PureBackend on exact plants | none | KEEP | The class does not depend on summation order for exact plants. |
| `linalg/CholeskySuite.scala:38` | Cholesky | `NotPositiveDefinite(1)`, 2×2 indefinite | none | KEEP | Integer index. The pivot `1-4=-3` is far from 0. |
| `CholeskySuite.scala:66-69` | Cholesky | `NotPositiveDefinite(1)` with tol 1e-12 and pivot 1e-14 | none | KEEP | Integer index with a 100× margin. Diagonal input, so the result is exact. |
| `CholeskySuite.scala:77-83` | Cholesky solve (vector and strided matrix RHS) | Typed error, input bits unchanged (`doubleToRawLongBits` `:79,:82`), factor unchanged `:83` | W1.2 (matrix RHS path) | KEEP | Input non-mutation and finiteness pre-check (plan W1.2: "Keep the transactional finiteness checks"). |
| `CholeskySuite.scala:84` | Cholesky vector solve | `Seq(2.0, 3.0)` exactly, identity factor | none | KEEP | Exactly representable (identity: `x/1`, `0·x`). The vector route is not rerouted. |
| `CholeskySuite.scala:91-94` | Cholesky solve | Overflow becomes a typed Left, and inputs are unchanged | W1.2 | KEEP | The post-solve finiteness check and non-mutation must survive `dtrsm`. |
| `linalg/DenseCholeskyWorkspaceSuite.scala:26-28, 33-35` | workspace factor/solve | NaN in the upper triangle and the 99/88/77/66 sentinels are untouched (tol 0.0) | W1.4 | KEEP | Input non-mutation, and only the lower triangle is written. |
| `DenseCholeskyWorkspaceSuite.scala:36` | workspace solve | Factor-array bits unchanged after three solves | none | KEEP | Input non-mutation (the solve copies the factor into scratch). |
| `DenseCholeskyWorkspaceSuite.scala:45` | workspace factor | Bits of a failed input unchanged (0, -1, NaN, +∞, indefinite 2×2) | W1.4 if it folds the finiteness check into the panel pass | KEEP | Transactional failure (the `factorLowerInPlace` doc). Blocking must still factor in scratch and copy back only on success. |
| `DenseCholeskyWorkspaceSuite.scala:50-51` | workspace solve | RHS unchanged on NaN | none | KEEP | Input non-mutation. |
| `DenseCholeskyWorkspaceSuite.scala:72` | workspace solve | `0.5` exactly (1×1, `2/2/2`) | none | KEEP | Exactly representable. |
| `DenseCholeskyWorkspaceSuite.scala:82-86` | `testPositiveDefinite` | Enum class plus input `a` unchanged | W1.4 | KEEP | Classification and non-mutation. `NonFiniteInput` comes from the input scan (`DenseCholeskyWorkspace.scala:35`), not from the kernel. |
| `DenseCholeskyWorkspaceSuite.scala:93` | workspace solve | `Double.MaxValue` RHS unchanged after a non-finite solution | none | KEEP | Non-mutation. |
| `linalg/MatrixPropertiesSuite.scala:78` | Cholesky | `NotPositiveDefinite(1)`, 2×2 | none | KEEP | Integer index, small n. |
| `optim/ConstrainedRayleighSuite.scala:78-81` | Cholesky | `DenominatorNotPositiveDefinite(1)` for `diag(1,0)` | none | KEEP | Exact zero pivot. Integer index. |
| `linalg/TriangularSolveSuite.scala:42, 51` | `dtrsv` | `SingularMatrix(1)` and `SingularMatrix(2)` on exact zero diagonals | none | KEEP | Integer index. `TriangularSolve` stays on `dtrsv`. |
| `linalg/KernelRegressionSuite.scala:31` | gemv (row-major, col-major, strided) | `Seq(3.0, 7.0)` exactly, with NaN in `y` | W1.1 | KEEP | Integer-valued: every order gives the same result. The NaN-ignoring `beta==0` assignment must survive the 4-row tiles. |
| `linalg/DenseCoreSuite.scala:57` | gemv | `Seq(8.0, 20.0)` exactly | W1.1 | KEEP | Integer-valued 2×3. |
| `DenseCoreSuite.scala:92, 97` | gemv `mulInto` (row-major, transpose) | Exact values plus -7 sentinels outside the output view | W1.1 | KEEP | Integer-valued, and the `y` stride/offset must be respected. |
| `linalg/LinearOperatorSuite.scala:12, 16` | gemv `applyTo` / `transposeApplyTo` | `Seq(8,20)` and `Seq(-3,-3,-3)` exactly | W1.1 | KEEP | Integer-valued. |
| `backend/BackendSeamSuite.scala:157-228` | gemm/syrk/gemv routing witnesses | Doubled versus pure at **1e-12** tolerance | W1.1 | KEEP (already tolerance) | No bit pin. The class doc (`:11`) calls the whole suite the "byte-identical witness"; that is wording only. |
| `BackendSeamSuite.scala:248-250, 257, 265, 277-278, 285` | LU/Cholesky/QR routing | Integer call counts and rank | none | KEEP | Integer routing results. |
| `backend/BackendContractSuite.scala:270-272` | `PureBackend.gemv` versus `DMat.*` | Exact equality (tightened from 1e-12 in W1.1), plus a bitwise test over row-major and transposed views with tails | W1.1 | TIGHTENED | Since W1.1 both routes run the same layout-selected kernel (`Backend.scala:150` and `Matrix.scala:323`). |
| `spectral/DenseSymmetricWorkspaceSuite.scala:61` | eigSym values-only | **Workspace equals ordinary route, bit for bit**, n=9 | none if both stay on one `tred2` | KEEP | Route agreement through the shared kernel (plan W1.5). n=9 is below 64, so this test does **not** check sharing above the threshold. |
| `DenseSymmetricWorkspaceSuite.scala:67-68, 85-86` | eigSym | The same route twice gives the same bits | none | KEEP | Determinism for a fixed build (contract `:17-18`). |
| `DenseSymmetricWorkspaceSuite.scala:40, 44, 46, 59, 76` | eigSym workspace | Exact scratch sizes: 56 = n²+n, 7 = n, 0, and the measured capacity | W1.5 if the latrd panel W (n×nb) is drawn from the workspace | CONDITIONAL | This is an API contract (`symmetricEigenRequirement`, `DenseSpectralKernels.scala:156`), not floating point. Keep it if W1.5 fits the existing `(d, e, eOffset, workspace)` layout. Otherwise the requirement change is a reviewed API decision and must not be loosened silently. |
| `spectral/EigSymmetricBackendSeamSuite.scala:36,43` via `:107, :140, :184, :210, :217-224, :244` | eigSym facade | Exact eigenvalues and eigenvectors, provider route versus pure | none | KEEP | Both sides call `DenseSpectralKernels.symmetricEigen`. The test providers (`:58,:74,:234`) call the same kernel, so equality holds at any n while the kernel is shared. |
| `spectral/EigSymmetricGeneralizedSuite.scala:185` | generalized eigSym (Cholesky + tred2) | Upper-triangle garbage gives bit-identical eigenvalues, n=6 | none | KEEP | Both runs make identical lower-triangle copies (`Factorizations.scala:275-282`, `symmetrizedLowerRowMajor`). Blocked kernels must read only those copies. |
| `spectral/EigSymmetricGeneralizedBackendSuite.scala:81, 126-127` | LOBPCG (Rayleigh-Ritz eigSym and Cholesky) | Routed equals pure, bit for bit | none | KEEP | Both sides run the same pure LOBPCG code. The projected problems are small. |
| `spectral/LobpcgSuite.scala:151-152` | LOBPCG | Two identical runs are bit-identical | none | KEEP | Determinism (contract `:17-18`). |
| `spectral/EigSymmetricGeneralizedOperatorSuite.scala:137-138` | generalized | Result unchanged after the operators are corrupted | none | KEEP | Result ownership and non-mutation. |
| `spectral/EigSymmetricDenseSuite.scala:128-129` | eigSym values-only | Diagnostics `0.0` | none | KEEP | Structural: no vectors means no residual arithmetic. |
| `spectral/DenseSpectralKernelsSymmetricSuite.scala:229` | tridiagonal QL | `DidNotConverge(0)` at a zero sweep budget | none | KEEP | Budget count from `solveTridiagonal`, not `tred2`. |
| `optim/CMAESSuite.scala:17, 22-23` | CMA-ES (eigSym n=2) | Repeat runs are bit-identical | none | KEEP | Determinism for a fixed build. |

**Named leads with no W1 bit pin:**

- `linalg/GemmTuningSuite.scala:46` pins exact `gram(i, j) == gram(j, i)` for `a.t * a`.
  That symmetry is exact because the `dsyrk` path mirrors one triangle and does not
  recompute it. The values themselves are compared at 1e-12 (`:40`). W1.4's new
  `dsyrkLowerUpdate` is an internal Cholesky trailing update and does not replace the
  public gram route. KEEP. If W1.4 reroutes `a.t * a` through it, the new kernel must still
  mirror the triangle.
- `linalg/PerfDoctrineSuite.scala:18-19` pins exact elementwise add/sub against a reference
  loop. Those operations are correctly rounded single operations with no W1 route. KEEP.
- The plan (`docs/breeze-competitiveness-plan.md`, W1.5 and W1.7) cites
  `DenseSymmetricWorkspaceSuite:52`; that line is the test name. The exact-agreement
  assertion is at `:61`.
- `GaleNumericalContractSuite.scala:172`, the line the plan cites, is a comment. The
  assertion it introduces is at `:178` (see the table).

**Excluded after reading:** parity suites, `laws` property suites, `VectorGemmSuite`,
`FfmLapackSuite` and `BreezeMigrationSuite` all use tolerances. FFM `:123` pins a LAPACK
result, not a pure kernel. `TinyKernelSuite` det and matvec use `Tiny.scala` cofactor and
unrolled code, not LU or gemv. The QR, banded Cholesky, sparse, and LM/L-BFGS-B zero-tolerance
pins have no W1 route. `MatrixLaws.assertExact` is used only for `A.t.t`.

## Summary

| Handling | Count |
|---|---|
| KEEP | 39 rows: integer indices and pivots, input non-mutation, exactly representable values, route agreement or determinism, and 3 rows that already use tolerances |
| CONVERT | **0** |
| CONDITIONAL | 1 row: the eigSym scratch-requirement sizes |

No test compares a W1-affected output against a hard-coded value that is not exactly
representable, so no test needs to move to a tolerance. The real risk is the opposite: there
is **no coverage above any blocking threshold**. Every dense LU, Cholesky and symmetric-eigen
test uses n ≤ 25. Parity uses n ∈ {4..25}, `laws` uses n ≤ 12 and core uses n ≤ 14. The only
n ≥ 96 tests are QR (`QRSuite:103-148`, `BackendSeamSuite:234,284`) and matrix-free Lanczos.
Without new tests, W1.3 to W1.5 would ship without exercising the blocked path at all.

W5.3 (same branch, `breeze/w53-parity-hardening`) closes this gap at the Breeze-parity
level only. `FactorizationHardeningParitySuite` covers LU, Cholesky, inverse, QR and least
squares at n = 96, 128, 200 and 300. `SpectralHardeningParitySuite` covers eigSym at the
same sizes. All comparisons use c·n·ε·κ, Weyl or Davis–Kahan tolerances, so they stay valid
under reassociation. The core-suite exact twins recommended below are still needed.

W5.3 also found that the unblocked tridiagonal QL kernel
(`DenseSpectralKernels.solveTridiagonal`) used to exhaust its 30-sweeps-per-value cap for
n ≥ about 96 when the spectrum had a highly repeated eigenvalue at 0. Commit `5237c48` on
this branch fixed it with a norm-scaled deflation test. The fix changes no existing
exact-equality result: `coreJVM`, `coreJS`, `lawsJVM` and `parity` all pass unchanged. W1.5
should be measured against the post-fix kernel.

## Recommendations by W1 item

- **W1.1 gemv.** No existing test breaks, because every exact gemv pin is integer-valued.
  - Bit-identical output is achievable and preferred. Keep each row's four-accumulator
    sequence in `dgemvRowMajor` (`:262-282`). In `dgemvColMajor`, nest the four column
    updates as `fma(s3,a3,fma(s2,a2,fma(s1,a1,fma(s0,a0,y))))`.
  - If the implementation instead reassociates, contract `:14-15` allows it.
  - **As implemented (branch `breeze/hc-gemv`):** `dgemvColMajor` follows the nesting
    above and is bit-identical to the one-column sweep. `dgemvRowMajor` deviates from the
    four-accumulator recommendation. A four-row tile that kept four accumulators per row
    (sixteen chains) spilled registers and ran 0.87× at n=256. The kernel therefore uses
    two accumulators per row over column pairs, `(p0 + p1)`, followed by an `fma` for an
    odd last column. The leftover rows (`rows % 4`) use exactly the same per-row
    arithmetic as the tile, so each output row depends only on that row and `x`, not on
    its position: `(A * x)(i)` is bit-identical to `A.slice(i, i + 1, 0, cols) * x` on
    JVM and JS (`KernelRegressionSuite`, "row-major gemv rows are bit-identical to their
    1-row slice products"). Compared with the pre-W1.1 kernel, row-major gemv bits change
    for non-integer inputs. This reassociation is within contract `:14-15`. No exact pin
    breaks, and the Breeze goldens still pass.
  - Keep the `beta==0` assignment so the `KernelRegressionSuite:10-32` NaN test holds.
  - Add one tolerance test with rows ≥ 9 and a remainder (`rows % 4 ≠ 0`, cols ≥ 8) against a
    naive or `BigDecimal` reference, so the tile and remainder paths run on `PureBackend`.
    The `BackendConformanceSuite` gemv test (`:109`, 5×7) runs only for the vector and FFM
    backends.
- **W1.2 trsm.**
  - Preserve `CholeskySuite:72-97` (raw-bit non-mutation and pre/post finiteness Lefts) and
    `GaleNumericalContractSuite:63` behaviour.
  - Cholesky `solve(DMat)` bits will change at every n, because the current loop does not
    use fma. No test pins those bits.
  - LU matrix-solve columns will stop matching the vector solve bit for bit. No test pins
    that either; do not add one.
  - Add tolerance tests with n > block size (for example 100×7 RHS) for the LU and Cholesky
    matrix solves, both against the vector solve column by column at a relative residual
    level. Add a strided or transposed RHS non-mutation case at that size.
  - Parity `FactorizationParitySuite:77,157` uses n ≤ 25. Add n ≈ 100.
  - **As implemented (branch `breeze/hc-solve`):** `DoubleKernels.dtrsmLeft` is a
    left-looking blocked solve with 8-row diagonal blocks; the off-diagonal update is one
    `dgemm(-1, 1)` per block. Both the LU and Cholesky matrix solves use it (Cholesky `Lᵀ`
    through swapped strides, packed row-major above one block). For n ≤ 8 each column is
    bit-identical to `dtrsv`, so the LU matrix solve keeps its old bits there. The Cholesky
    matrix solve changes bits at every n, because its old loop did not use fma. Above n = 8
    both routes change bits. All of this is within contract `:14-15`. The diagonal is scanned before any write, so `SingularMatrix(k)` and
    `NotPositiveDefinite(k)` report the same index as before. `TrsmKernelSuite` and the
    n = 100 × 7 tests in `LUSuite` and `CholeskySuite` cover the blocked path. The n = 96 to 300
    multi-RHS cases in `FactorizationHardeningParitySuite` cover it against Breeze.
- **W1.3 blocked LU.**
  - `LUSuite:28,49`, `NumericalPolicySuite:51-75` and `GaleNumericalContractSuite:178` stay
    exact because they are below the threshold.
  - Add tests above the threshold (n ≈ 128) that stay exact:
    1. A planted exact zero column at k=100 must give `SingularMatrix(100)` in any order
       (contract `:43`).
    2. A strictly diagonally dominant matrix with an integer pivot vector known by
       construction, such as a row permutation of a dominant matrix, must return the exact
       pivots and `pivots(i)` = original row.
  - Add tolerance tests for the `‖PA−LU‖` reconstruction and for `det` relative error.
    Extend parity LU and det to n ≈ 128.
- **W1.4 blocked Cholesky.**
  - Route both `Factorizations.cholesky` (`:283`) and `DenseCholeskyWorkspace`
    (`:40`, `:128`) through the same blocked entry point, for example by making
    `dpotrfLower` dispatch.
  - Add an exact route-agreement test above the threshold: `factorLowerInPlace` lower
    triangle == `cholesky(A).lower`, as raw bits. This is the analogue of
    `DenseSymmetricWorkspaceSuite:61`, and no such test exists today.
  - Keep `DenseCholeskyWorkspaceSuite:36,45` and `:97-103` (overflow fails closed) when you
    fold the finiteness check. A non-finite value must still produce Left with no partial
    write.
  - Add `NotPositiveDefinite(k)` at k > nb with an exact plant, for example a diagonal
    matrix with a zero at index 100, and a tolerance reconstruction test.
- **W1.5 blocked tridiagonalization.**
  - Keep `DenseSymmetricWorkspaceSuite:61` exact, and add an n ≈ 96 twin for values-only
    and for vectors.
  - Keep `EigSymmetricBackendSeamSuite`, which calls the shared kernel.
  - Decide the panel scratch explicitly. If it does not fit the existing layout, update
    `symmetricEigenRequirement` and `DenseSymmetricWorkspaceSuite:40-46,59,76` as an API
    change.
  - Add eigen residual and orthogonality tolerance tests at n ≈ 96 and 200, and extend
    `SpectralParitySuite` (n ≤ 24 today).
  - `DenseSpectralKernelsSymmetricSuite:232` (values-only versus full at 1e-12, n=14) needs
    an above-threshold twin, because the `wantQ` and values-only blocked paths may diverge
    by ulps.
