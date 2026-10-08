# Breeze competitiveness plan

Review date: 2026-10-08. Source revision: `f4ff88d943227cc318c415ad9b1393b05310c8d6`.
Status: approved 2026-10-08; in execution.

Mote epic: `bd-01M4EJ8JTE95TM3HQ7ZYC358EN`; one task per workstream item (W0–W6).


## Context

Breeze 2.1.0 is gale's main competitor: it is mostly retired and JVM-only, but battle-tested through Spark. The last full head-to-head (`benchmarks/results/2026-07-11-breeze-release-grade.md`) had gale ahead on 23/42 pairs. Since then none of the gap kernels has improved. An audit on 2026-10-08 found five problems:

- **The baseline is mislabelled.** Since 07-17, JMH forks receive `--add-modules=jdk.incubator.vector` (`build.sbt:415`, added in df62a78), so any default Breeze run since then gets netlib **VectorBLAS**, while the README and build comment still say F2J. The 07-11 release-grade sweep (80611e3) predates both `.jvmopts` and that option, so it most likely ran scalar `Java11BLAS` against scalar gale: its gaps (axpy, gemv, LU, solve) are genuine scalar-vs-scalar gaps. No receipt recorded the class until W0 (verified 2026-10-08 by w0-bench).
- **Kernel gaps remain.** LU and Cholesky are unblocked. There is no trsm. Row-major gemv is untiled and the transposed gemv is not unrolled. `tred2` is unblocked. Cholesky got an extra copy plus finiteness pass on 10-03 (2c794d8) that has not been re-measured.
- **The everyday API is missing.** Gale has no reductions, no exp/log/softmax/logsumexp, no public norm1/normInf/Frobenius and no SparseVector. Benchmarks can't be won on operations gale doesn't have.
- **Parity tests use easy inputs.** Inputs are well-conditioned and n ≤ 25, so future blocked paths would never be exercised. There are no edge cases (empty, 1×1, NaN/Inf, singular, ill-conditioned, views). LBFGS/LBFGSB have no Breeze parity. Nothing on Scala.js is checked against an external reference.
- **Benchmark coverage is thin.** There is no Breeze benchmark for SVD, eig, inv, det, pinv, sparse, optimizers, elementwise, small fixed-size, multi-RHS, or n ≥ 512.

**Decisions taken (user, 2026-10-08):**

- Reopen ADR A-2b so SIMD can serve L1 and reductions, behind a benchmark gate.
- `sum` is fast by default, with `sumExact` as the exact variant.
- Empty `max`/`argmax`/`mean` throw `LinAlgError`.
- In scope: reductions, elementwise numerics, the three norms, SparseVector.
- Out of scope: FFT, distributions, stats, broadcasting. Broadcasting is a candidate for a later version.

**Outcome sought:** a reproducible two-lane scoreboard in which gale ≥ Breeze on every L2, L3, factorization and eigen pair at n ≥ 64. Any exception must be documented with its cause. The new op families must be ≥ 1.0× Breeze. Every new or changed path must be parity-tested against Breeze, including edge cases and sizes above the blocking thresholds.

## Ground rules

- **Repo state.** The main tree has another session's uncommitted optim work (CMAES, ALM, and LBFGS edits). Each workstream runs in its own git worktree and branch, with one sbt writer per tree. Never stage unrelated files.
- **Disk.** Run `df -h /System/Volumes/Data` before full builds and stop if less than 20 GB is free. Delete each worktree's `target/` dirs on exit.
- **Tracking.** Create one mote epic, "Breeze competitiveness", with one task per numbered item below. Coordinate with the open task `bd-01M427PJCGSXP2VCCE97EEMT5T` (dense Cholesky workspace, currently *doing*), which also calls `dpotrfLower`.
- **Plan copy.** Save a copy of this plan in the repo as `docs/breeze-competitiveness-plan.md`, following the existing `docs/*-plan.md` header convention.
- **Gates for every merge.**
  - `sbt scalafmtCheckAll compileAll testAll parityTest interopBreezeTest docs/mdoc`
  - Plus `vectorBackendTest` for W2 and anything that touches `Backend`.
  - Plus `docsCheck` for documentation changes.
  - A fresh-context review subagent pass before merge.

---

## W0 — Measurement first (do before any kernel work)

1. **Log Breeze's BLAS.** Add `dev.ludovic.netlib.blas.BLAS.getInstance().getClass.getName` (plus the LAPACK equivalent) to a `@Setup(Level.Trial)` in `BreezeBenchData.scala`, written to stderr and to a sidecar file. Every result file records it.
2. **Run two lanes.** Use one pinned JDK, 25 LTS (Homebrew `openjdk@25`, the July receipt JDK; no JDK 21 is installed locally), recorded in each receipt.
   - **Lane A, out-of-box.** Breeze runs whatever non-SIMD BLAS netlib picks, and gale runs pure.
     - JMH forks inherit the host's `--add-modules`, and `jvmArgsAppend` cannot remove it. Pass the CLI `-jvmArgs` (which replaces the inherited args) through a `breezeLaneA` sbt alias. Run with `-p backend=pure`, because the vector option would throw `NoClassDefFoundError` without the module.
     - The receipt **rejects** a VectorBLAS class and also a native class: netlib tries native first, so a Linux host with `libblas.so.3` would otherwise quietly get native BLAS.
   - **Lane B, SIMD.** Fork with `--add-modules`, so Breeze runs VectorBLAS and gale runs with `import gale.backend.jvm.vector.given`. Add a `@Param("pure","vector")` backend switch to the gale side of the Breeze benches, chosen in `@Setup`.
3. **Scoreboard tool.** `tools/bench/breeze_scoreboard.py` reads JMH `-rf json`, pairs the gale and Breeze benchmarks by (op, size), and emits a markdown table: ratio, 99.9% CI overlap → ahead/tie/behind, lane, JDK and netlib class. This replaces the hand-written tables.
4. **Fix the docs.** Correct the "pure-Java F2J" claims (`build.sbt:418-421`, `benchmarks/README.md:136`) and the July receipt's caveat.
5. **Expand Breeze bench coverage** (new JMH classes in `benchmarks/jvm/src/main/scala/gale/bench/`):
   - Factorizations and eigSym at n = 512 and 1024.
   - Multi-RHS solve (n=256, k=64).
   - inv, det, pinv, full SVD, nonsymmetric eig.
   - 3×3 and 4×4 gemm/solve/inv.
   - CSC/CSR matvec and matmul vs `CSCMatrix`.
   - Elementwise add/sub/hadamard and the W3 reductions and numerics as they land.
   - SparseVector dot/axpy/+ (from W4).
   - LBFGS on Rosenbrock-100 and logistic regression (time to tolerance plus evaluation count).
6. **Re-baseline.** Run both lanes, 2 forks, on a quiet machine, and save as `benchmarks/results/2026-10-xx-breeze-two-lane-baseline.{json,md}`. Every later item is judged against this file.

## W1 — Pure-kernel performance (shared JVM/JS, in `core/shared/.../gale/kernel/DoubleKernels.scala` and `linalg/Factorizations.scala`)

Order by expected gain against risk:

1. **gemv.**
   - Row-major `dgemvRowMajor` (`DoubleKernels.scala:245`): use 4-row tiles that share each loaded `x[j]`, each row keeping fma accumulators, and handle the remainder rows.
   - `dgemvColMajor` (`:293`): unroll 4 columns per pass over `y`.
   - Keep the strided fallbacks.
2. **trsm.** Add `dtrsmLeft(uplo, unitDiag, n, nrhs, alpha, A(off, ld), B(off, ld))` for row-major sub-blocks, blocked at 32 to 64. The diagonal block uses a dtrsv-style loop; off-diagonal blocks use `dgemm(alpha=-1, beta=1)`.
   - Route the LU matrix solve (`Factorizations.scala:1214-1271`) and the Cholesky matrix solve (`:1297-1339`) through it.
   - Keep the transactional finiteness checks (`finiteCholeskyValues`, `:1341`) and the guarantee that inputs are never mutated (`CholeskySuite:72-97`).
3. **Blocked LU (getrf-style)** above a threshold of about n ≥ 96; tune it, and keep small n on the existing path so `LUSuite`'s 2×2 pivot and `SingularMatrix(1)` cases stay unchanged.
   - Factor a panel of width nb with partial pivoting, and use fma in the panel update.
   - Apply each row swap across the full row, so the `pivots(i)` = original-row convention (`:17`) and `parity` are preserved.
   - Compute U12 with `dtrsmLeft(lower, unit)` and A22 with `dgemm(-1, 1)` through `dispatchGemm` (`Factorizations.scala:982`).
   - `SingularMatrix(k)` must report the same global k.
4. **Blocked Cholesky (lower, right-looking)** with `dpotrfLower` as the diagonal-block kernel.
   - Add `dpotrfLowerAt(n, a, off, ld, tol)`, an offset and leading-dimension variant.
   - Compute L21 with `dtrsm` (right side, lower transposed). Add this variant, or transpose-pack and reuse the left variant.
   - Add a new `dsyrkLowerUpdate(C -= A·Aᵀ, lower triangle only, alpha/beta)` kernel with a 4×4 register tile. This avoids the slow `dgemmStrided` NT path (`:474`).
   - Re-measure, and if needed remove, the 2c794d8 copy and finiteness overhead by folding the check into the panel pass.
   - Do this after the workspace task (`bd-01M427…`) lands, or coordinate with its owner.
5. **Blocked tridiagonalization** (dsytrd: latrd panel plus a `dsyr2kLowerUpdate`) for n ≥ about 64 in `spectral/DenseSpectralKernels.scala:226`.
   - Ordinary and workspace eigen routes must keep calling one shared kernel (`DenseSymmetricWorkspaceSuite:52` asserts they agree exactly).
   - Fit the `(d, e, eOffset, workspace)` layout.
   - This is the largest item; do it last.
6. **Accept-or-revert rule.** Each item ships only if the scoreboard shows ≥ 1.0× (lane A) at its target sizes and no regression greater than 3% elsewhere. JS gets the same code and has its own gate: no more than 5% regression at n ≤ 256 in `benchSmokeJSFull`. If JS fails that gate, make the blocking threshold platform-specific through `core/{jvm,js}` platform constants.
7. **Bit-pinning audit (before W1.3 to W1.5).** List the tests that assert exact equality on LU, Cholesky or eigen outputs: the `LUSuite` pivots, `DenseCholeskyWorkspaceSuite:36,45`, `DenseSymmetricWorkspaceSuite:52`, `GaleNumericalContractSuite:172` and the `BackendSeamSuite` witnesses. Workspace and ordinary routes must keep sharing one kernel, so exact agreement between them still holds. Tests that compare against hard-coded values switch to tolerance only where the numerical contract already allows reassociation. Never loosen them to make a failing result pass.

## W2 — Reopen ADR A-2b: SIMD for L1 and reductions (gated)

1. **Spike** (in `backend-jvm-vector`, measured only). `VectorBackend` currently forwards `dot`, `axpy` and `nrm2` to the pure kernels, so the spike must first write real SIMD versions of them. Then add SIMD `dsum`, `dmax/argmax` and a `dexp` map (Vector API `lanewise(EXP)`). Measure them with JMH against the scalar kernels and against Breeze VectorBLAS at n = 1K, 64K and 1M.
   - **Gate:** the SIMD kernels must be ≥ 1.25× scalar at n ≥ 4K and ≥ 1.0× Breeze VectorBLAS. The gate applies to each op separately.
   - Vector `EXP` is often slower than scalar `Math.exp`, so expect `exp` to fail. That is an acceptable result.
   - Also measure dispatch overhead under requirement A-R1: the cost of a `using Backend` call at small n (16 to 256) against a direct static call.
   - If every op fails the gate, record the results and close the ADR as "keep A-2b".
2. **ADR**, written as an amendment in `docs/backend-architecture.md` with a new decision id. It defines:
   - Public **facade** L1 and reduction ops gain `(using backend: Backend)`. The companion `given pure` keeps call sites source-compatible (this is pre-M1 with no MiMa, so the binary break is acceptable but must be recorded in `api-stability.md`).
   - **Loop-called internal L1 stays static** `DoubleKernels.*`. Factorization inner loops get no dispatch overhead, which keeps A-2b's original rationale.
   - The explicit-given selection (no ServiceLoader, no auto-default) is **unchanged**.
   - SIMD `exp` must meet a documented accuracy bound and law-equivalence to the pure kernel, not bit identity.
3. **Implement**, if the gate passed:
   - Extend `DenseDoubleKernel` (`backend/Backend.scala:30`) with `sum`, `maxIndex/minIndex`, and `expInto`/`logInto` maps, with pure defaults pointing at the W3 static kernels.
   - Override them in `VectorBackend`, along with the existing forwarded `dot`, `axpy` and `nrm2` (`VectorBackend.scala:31-53`).
   - Add the vector backend to `BackendConformanceSuite`.

## W3 — New public API: reductions, elementwise numerics, norms (`core/shared`)

**Kernels** (static, `DoubleKernels`, each with a contiguous fast path and a strided path):

- `dsum`: 4 accumulators, unrolled, the `ddot` style (`:17`).
- `dmaxIndex/dminIndex`: first-occurrence index. A NaN propagates: the first NaN wins.
- `dasum`: abs sum, for norm1.
- `damax`: max abs, for normInf.
- `dexpInto/dlogInto/dsigmoidInto`.
- `dlogSumExp`: two-pass max-shift. All −∞ gives −∞; any +∞ gives +∞; NaN gives NaN.

**DVec** (`linalg/Vec.scala`, class members like `dot`/`norm2`):

- `sum`, `sumExact` (via `gale.numeric.ExactSum`; map `CapacityExceeded` to a `LinAlgError`), `mean`, `max`, `min`, `argmax`, `argmin`, `norm1`, `normInf`.
- Empty input: `max/min/argmax/argmin/mean` throw a new `LinAlgError.EmptyInput(op)`. Pre-M1 it is fine to add a case to the sealed class. `sum(empty)=0`, `norm*(empty)=0`.

**DMat** (`linalg/Matrix.scala`):

- `sum`, `mean`, `max`, `min`, `argmax` (returns `(row, col)`).
- Per-axis versions via a new `enum Axis { case Rows, Cols }`, as in `A.sum(Axis.Cols): DVec`. Use the `isContiguousRowMajor` fast path, otherwise one strided call per row (copy `addSub` at `:623`).
- Matrix `norm1` (max column abs sum), `normInf` (max row abs sum), `normFrobenius` (scaled with the `dnrm2` recurrence, overflow-safe).

**Elementwise numerics.** PRD.md:614 requires elementwise ops to be explicit. The existing `A.pointwise.map(f)` (`syntax/Syntax.scala:30-60`, DMat only, built on `tabulate`) stays as the general mechanism for arbitrary functions. The named, kernel-backed fast paths go in a `gale.linalg.Numerics` object. Its guide and scaladoc must say that `Numerics.exp(A)` gives the same result as `A.pointwise.map(math.exp)` but runs faster, so the two don't look like competing APIs.

- `exp`, `log`, `log1p`, `expm1`, `sigmoid` (stable: branch on sign, use `exp(-|x|)`, as in `DenseObjectives.scala:87-104`), `logSumExp`, `softmax`, `logSoftmax`.
- Each has DVec and DMat overloads; softmax and logSoftmax also take an axis.
- `sum` stays fast and `sumExact` stays exact.

**Migrate private duplicates** to the public implementations:

- `Factorizations.scala:1957` `norm1`.
- `Eigen.scala:1452` naive `frobeniusNorm`, which has an overflow risk. Check the bit-sensitive tests.
- Leave `OptimNumerics.normInf` (max-abs, a different meaning) alone; at most rename it.

**Naming trap:** Breeze's `softmax(v)` returns log-sum-exp. Document this in the migration table.

## W4 — SparseVector (`core/shared/.../gale/sparse/SparseVector.scala`)

- **Representation:** `length`, sorted unique `IndexArray` indices, `DoubleArray` values.
  - Builder from pairs: sort, then resolve duplicates with the existing `DuplicatePolicy` (`Sum`/`Last`/`Error`, `Sparse.scala:11`).
  - Non-finite values follow the existing `SparseValuePolicy` (`Sparse.scala:14-23`).
  - Explicit zeros need a new, documented rule: keep them by default, which matches Breeze, and add `compact()` to drop them.
  - Constructors from DVec (drop zeros) and empty.
- **Ops:**
  - Access: `apply`, `activeSize`, `foreachActive`, `toDense`.
  - Arithmetic: `dot(SparseVector)` (two-pointer merge), `dot(DVec)`, `+`/`-` (merge, the pattern of `CSR.zipValues` at `:654`), `* scalar`, `axpyInto(MutableDVec)`.
  - Norms: `norm1/norm2/normInf`.
  - Reductions: `sum`, and `max/min`, which must account for implicit zeros when `activeSize < length`.
  - `mapActive`, explicit that it applies to nonzeros only.
- **Interop:** `CSR * SparseVector → DVec`, `CSC.colSparse(j)`, `CSR.rowSparse(i)`. These are additive; the existing dense `row`/`col` stay.
- **Docs:** add a `gale.sparse` row note in `api-stability.md`.

## W5 — Correctness against Breeze (`parity/`, plus cross-platform goldens)

1. **New parity suites, ScalaCheck-driven with seeds logged, sizes 1 to 64.**
   - `ReductionsNumericsParitySuite`: sum, mean, max/min/argmax, norm(v,1|2|Inf), A Frobenius/1/Inf, exp/log/sigmoid, logsumexp (vs Breeze `softmax`), softmax.
   - `SparseVectorParitySuite`, against Breeze `SparseVector`.
   - Edge cases are explicit: length 1, ±Inf, NaN, all −∞, ties for argmax. Where gale and Breeze disagree, assert gale's documented semantics and list the divergence in the suite and the migration doc.
2. **`OptimizerParitySuite`.** Compare gale LBFGS and LBFGSB with Breeze `LBFGS` and `LBFGSB` on quadratics, Rosenbrock-n and bounded logistic regression.
   - Compare the optimum x* and f* within a conditioning-scaled tolerance, and the projected-gradient norm at exit. Do not compare iterates.
   - Set both libraries' convergence tolerances explicitly to equivalent criteria. Report time-to-tolerance and evaluation count as information, not as a parity assertion.
   - Start only after the uncommitted optim work in the main tree (LBFGS, LBFGSB and StrongWolfe edits) has landed or its owner has been consulted.
3. **Factorization and spectral hardening.**
   - Add sizes **above every blocking threshold** (n = 96, 128, 200, 300). W1 items cannot merge without this.
   - Add 1×1 and empty inputs.
   - Add singular error paths: both libraries must fail. Compare failure, not index.
   - Add ill-conditioned inputs (Hilbert 6–12, `withSpectrum` with κ up to 1e12) with tolerance scaled by κ·ε.
   - Pass transposed and strided views into lu/chol/qr/solve/eigSym.
   - *Optional (beyond the core ask; do if time allows):* wide (m < n) lstsq, rank-deficient pinv and lstsq against Breeze (min-norm), and CSR/CSC with explicit zeros and duplicate entries.
   - Convert Factorization/SvdQr/Spectral fixed cases to ScalaCheck generators where they can be (`ParitySupport` already has `spd/symmetric/withSpectrum`).
4. **Cross-platform Breeze goldens.**
   - Add an sbt task in `parity` (JVM) that runs Breeze over a fixed seeded corpus covering dense ops, factorizations, eigen, the new reductions and numerics, and SparseVector. It emits generated Scala source, the same way `NumpyScipyFixtures.scala` is generated, into `core/shared/src/test/scala/gale/golden/BreezeGoldens.scala`.
   - Add a `BreezeGoldenSuite` in core shared tests that replays the corpus on **both JVM and JS**. This is the first external-reference check of JS numerics.
   - Tolerances scale with condition number times ε (c·κ·ε). They are never bit-exact, because `PlatformMath.fma` on JS is `a*b+c` (`core/js/.../PlatformMath.scala:8`), so results differ from the JVM in the last bits.
   - CI checks that the goldens are fresh: regenerate them and diff.

## W6 — Docs and evidence

- `docs/user/advanced/numerical-contract.md`: the NaN policy for reductions, the empty-input policy, the determinism of `sum` (per-platform, reassociation allowed) vs `sumExact` (correctly rounded), and the SIMD `exp` accuracy bound if W2 lands.
- New guide `docs/user/guides/reductions-and-numerics.md` (mdoc), a SparseVector section in the sparse guide, and `directory.conf` and `index.md` entries.
- Extend the existing `docs/user/guides/breeze-equivalence.md` (do not create a parallel migration doc): a Breeze → gale table covering names, semantics and divergences, including the softmax naming trap and the empty-input behaviour.
- `docs/api-stability.md`: an "Everyday-ops additions before M1" section, and an ADR note if W2 lands.
- `benchmarks/dashboard.md`: generated from the scoreboard tool.
- Final receipt `benchmarks/results/2026-1x-xx-breeze-two-lane-final.{json,md}`.

## Sequencing

```
W0 (baseline) ──┬─> W1.1 gemv ─> W1.2 trsm ─> W1.3 LU ─> W1.4 Chol ─> W1.5 tridiag
                │        (each W1 merge gated by W5.3 large-n parity)
                ├─> W3 kernels+API ─┬─> W2 spike ─> ADR ─> W2 impl
                │                   └─> W5.1 reductions parity
                ├─> W4 SparseVector ─> W5.1 sparse parity
                └─> W5.2 optimizer parity, W5.3 hardening, W5.4 goldens (independent)
W6 continuous; final two-lane scoreboard closes the epic.
```

Parallelism: W1, W3/W4 and W5.2–5.4 touch disjoint files and can run in separate worktrees at the same time. W2 waits for the W3 kernels so that SIMD overrides cover the reductions.

## Verification

- **Per item:** the full gate list from the ground rules, plus `sbt "testOnly"` on the touched suites (`LUSuite`, `CholeskySuite`, `TriangularSolveSuite`, `GemmTuningSuite`, `KernelRegressionSuite`, `DenseSymmetricWorkspaceSuite`, `BackendConformanceSuite`), plus `coreJS/test`.
- **Perf:** the scoreboard tool on both lanes against the W0 baseline. A ratio is claimed only when the CIs don't overlap. Each receipt records the JDK, the netlib class, the commit SHA and the machine.
- **Correctness:** `sbt parityTest` (Breeze, JVM) and `BreezeGoldenSuite` on JVM and JS. Test seeds are printed on failure.
- **Epic done when:**
  - Lane A has no L2, L3, factorization or eigen pair at n ≥ 64 below 1.0× without a documented cause.
  - New op families (reductions, numerics, norms, SparseVector) are ≥ 1.0× in lane A at n ≥ 1K.
  - Small fixed-size work (3×3, 4×4) and n < 64 are report-only. Per-call overhead dominates there, so they are tracked but are not a pass/fail target.
  - W5 suites are green on CI.
  - The docs gate (`docsCheck`) passes.
  - A fresh-context reviewer signs off on the final receipt.
