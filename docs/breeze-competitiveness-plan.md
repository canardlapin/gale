# Breeze competitiveness plan

Review date: 2026-10-08. Source revision: `f4ff88d943227cc318c415ad9b1393b05310c8d6`.
Status: approved 2026-10-08; in execution.

Mote epic: `bd-01M4EJ8JTE95TM3HQ7ZYC358EN`; one task per workstream item (W0–W6).

The original context below describes the October 8 baseline. The
[October 10 follow-up](#october-10-follow-up--correctness-and-performance-comparisons)
sets the next execution order after the first hill-climbs.


## Status 2026-10-11

Milestone 3 receipts (`90c74a3`) measure source `96c95e5` on trillium (Zen 5,
JDK 25), with every scoreboard run under `--validation --strict`:

- Lane A, gale pure vs Breeze `Java11BLAS`: 124 ahead, 13 tie, 4 behind,
  10 withheld of 151 pairs (Milestone 2: 122/12/7/10).
- Lane B headline, gale vector vs Breeze `VectorBLAS`: 118 ahead, 14 tie,
  9 behind, 10 withheld (Milestone 2 re-scored: 114/16/11/10). Secondary, gale
  pure vs `VectorBLAS`: 46/0/11/1 of 58.

Merged on 2026-10-11: hc-gemm-vec (`c3ceeb5`), hc-lbfgs (`c9e2b14`),
hc-lstsq16 (`51d74c9`), W2.2 docs (`96c95e5`, ADR A-2c), all measured by
Milestone 3; and hc-sigmoid2 (`6066cdb`), which merged after the receipts but
changes only tests, benchmarks and a comment in the kernel. In flight:
hc-max2, hc-gemv-vec and hc-lbfgs-step. The remaining gaps and their owners are
in [`benchmarks/dashboard.md`](../benchmarks/dashboard.md).

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
   - Ordinary and workspace eigen routes must keep calling one shared kernel (`DenseSymmetricWorkspaceSuite:61` asserts they agree exactly).
   - Fit the `(d, e, eOffset, workspace)` layout.
   - This is the largest item; do it last.
6. **Accept-or-revert rule.** Each item ships only if the scoreboard shows ≥ 1.0× (lane A) at its target sizes and no regression greater than 3% elsewhere. JS gets the same code and has its own gate: no more than 5% regression at n ≤ 256 in `benchSmokeJSFull`. If JS fails that gate, make the blocking threshold platform-specific through `core/{jvm,js}` platform constants.
7. **Bit-pinning audit (before W1.3 to W1.5).** List the tests that assert exact equality on LU, Cholesky or eigen outputs: the `LUSuite` pivots, `DenseCholeskyWorkspaceSuite:36,45`, `DenseSymmetricWorkspaceSuite:61`, `GaleNumericalContractSuite:172` and the `BackendSeamSuite` witnesses. Workspace and ordinary routes must keep sharing one kernel, so exact agreement between them still holds. Tests that compare against hard-coded values switch to tolerance only where the numerical contract already allows reassociation. Never loosen them to make a failing result pass.

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

## October 10 follow-up — correctness and performance comparisons

Planning reference: `breeze/integration` at `e8f1f8732acc79df0b91d0c0dc4d4e39477ed978`.
Status: next-stage plan requested 2026-10-10; implementation has not started under
this follow-up. Keep the existing epic and workstream ownership.

**Outcome.** Every reported speed ratio measures a successfully completed,
numerically validated operation. Publish an updated two-lane comparison, then
reduce its remaining losses through separately measured changes.

**Current evidence.** The complete Milestone 1 receipts measure `ff03eed`, before
the max, solve/inverse and D&C follow-ups. Their old losses are not current-head
verdicts. Nonsymmetric eig remains owned by HC-1. Direct nonsymmetric Breeze
parity stops at n=12; full-SVD factor parity is mostly small; JS goldens omit
nonsymmetric eig and large spectral factors. Current receipts cover Zen 5/JDK 25,
Java11BLAS or VectorBLAS with F2jLAPACK, not native-backed Breeze.

### 1. Validate the exact workloads before timing

Extend the existing Breeze benchmark setup helpers and decomposition/factorization
JMH classes. Before any timed iteration, run both implementations on the exact
fixture and selected backend. Require successful completion, finite outputs where
appropriate, the requested result shape/count, and operation-specific accuracy:
scaled solve/inverse residuals, factor reconstruction, eigenpair residuals and
SVD reconstruction/orthogonality. Respect sign and cluster-basis ambiguity.

Record the preflight result beside the existing SHA, seed, backend and netlib
identity. A `Left`, nonconvergence, or failed accuracy check aborts the trial; it
must never produce a competitive timing. Keep validation outside timed methods
and preserve the existing lazy-result consumption contract. Diagnostics timing,
if measured later, is a separately named workload.

Acceptance: an intentionally rejected result demonstrably fails validation;
all benchmark fixtures pass preflight on the selected lane. Run `benchCompile`,
focused validation checks, and short JMH smoke runs. Smoke timings are not receipts.
If current eig fails, retain the reproducer and return it to HC-1 before timing.

### 2. Close the spectral correctness gaps and finish HC-1

Extend the existing nonsymmetric, full-SVD and spectral-hardening parity suites,
and the Breeze golden generator/replay. Inventory existing cases first; add only
the missing discriminating cases.

- Nonsymmetric eig: cover the benchmark fixtures at n=16/64/256, planted real and
  complex spectra, strong nonnormality, near-defective cases and scale extremes.
  Check normalized residuals for both libraries; use forward eigenvalue/vector
  comparisons only when conditioning makes them informative. Preserve packing,
  ordering and values-only contracts. Coordinate with the existing HC-1 owner.
- Full SVD/pinv: extend to n=128/256 and selected tall/wide shapes, repeated or
  clustered singular values, graded scales and rank deficiency. Use residuals,
  orthogonality, cluster projectors and Moore–Penrose identities; align cutoff
  policy explicitly before comparing pseudoinverses.
- Symmetric eigen: validate the timed n=512/1024 fixtures; retain the existing
  threshold, tiny/huge-scale, workspace-ownership and bit-pin regressions.
- JVM/JS goldens: add a compact nonsymmetric corpus, another D&C vector case
  above n=48, and full-SVD factor evidence. Store vectors or projectors only where
  their comparison is well-defined; constrain generated-source and CI costs.

Acceptance: independent analytic or Breeze evidence plus numerical invariants,
with tolerances justified by scale, precision, conditioning and spectral gaps.
No tolerance widening to conceal failures. Keep a small deterministic CI corpus;
run expensive benchmark-size/adversarial cases in a separately invocable stress
lane if their measured runtime would burden normal CI. Re-run `parityTest`,
`breezeGoldensCheck` and JVM/JS core tests after the changes. HC-1 must pass this
evidence and the existing merge gates before inclusion in Milestone 2.

### 3. Freeze and publish Milestone 2 evidence

After steps 1–2 and the HC-1 merge, freeze one source SHA. Run both existing JMH
lanes on exclusive trillium resources using the established JDK, CPU/memory
pinning, two forks and per-class warmup/measurement settings. Preserve complete
raw JSON, fork logs, preflight evidence and BLAS/LAPACK sidecars. Use the existing
runner and strict scoreboard; record exact commands and scheduler job IDs.

Report deltas against Milestone 1, especially eig, eigSym, SVD/pinv, inverse,
multi-RHS solves, max and n=64 factorizations. Keep pure-versus-vector control
rows, but summarize the comparable vector lane without double-counting them.
Continue withholding comparisons with mismatched work or stopping rules. Capture
allocation evidence for the spectral and solve workloads where workspace or
lazy diagnostics changed, using a separately identified profiling run if needed.

Acceptance: every timed fixture passes preflight; strict pairing/backend checks
pass; each remaining loss has an owner and an evidence-based explanation or next
experiment. Do not transplant branch-local ratios into this receipt. A failure
blocks the affected claim and is reported explicitly.

### 4. Optimize the remaining measured losses, one change at a time

Use Milestone 2 to select sizes and confirm priorities. First investigate SIMD
`AᵀA` (M1 n=256: 0.62×), tall GEMM (n=64: 0.62×), and sigmoid (0.53–0.85×).
Profile layout/packing and dispatch for the products; profile scalar branching
and exponential evaluation for sigmoid. Preserve sigmoid's extreme-input and
NaN semantics. Treat GEMV's smaller 3–8% deficits as a later tuning target.

For each candidate: state the bottleneck hypothesis, retain a same-machine
baseline, change one kernel or threshold, run focused correctness checks, and
measure both the target and neighboring sizes/shapes. Apply the existing W1
accept-or-revert rule: achieve the target Breeze ratio, no reproducible >3%
regression elsewhere, and no >5% shared-code JS regression at n≤256. Reject or
revise a failed candidate; do not start an unbounded tuning loop. Large blocked
LU/Cholesky/tridiagonalization remain later work unless M2 makes them a priority.

### 5. Complete the gated SIMD routing decision

Continue W2 rather than opening another SIMD effort. The retained AVX-512 spike
passes for norm, sum and argmax, but dispatch fails A-R1 at small sizes. Resolve
that cost while preserving explicit backend selection and custom-backend
contracts, then refresh the gate against current scalar kernels. Measure
monomorphic and mixed-backend call sites and allocation, including small inputs.

Acceptance: satisfy both the existing per-operation speed gate and A-R1 before
routing a passing operation. A failed gate leaves that operation unrouted and
documented. Do not relax the gate implicitly. SIMD exp remains excluded from
default routing because its output is JIT-tier/position dependent. Scope every
result to the measured architecture/species width; AVX-512 does not qualify AVX2.

### 6. Compare optimizer time to equivalent accuracy

Reuse the shared stopping criterion and independent residual recomputation from
`OptimizerParitySuite` in the JMH optimizer workload. Retain fixed-iteration
measurements as a separate cost comparison. Start with Rosenbrock and regularized
logistic regression; add bounded logistic/LBFGSB with its projected-gradient
criterion once the shared harness is established.

Record final objective/residual, status, iterations, function/gradient calls,
allocations and elapsed time. Both solvers must reach the same requested accuracy
before a time-to-solution winner is declared. Preserve early termination/failure
information; no success ratio for a failed solve. Only then profile and optimize
callback copies, line-search work or vector kernels implicated by the evidence.

### Integration, qualification and next action

Use one primary writer and the existing independently reviewed merge process;
do not duplicate the active HC-1 work. Preserve public numerical, ownership and
workspace contracts. Apply the ground-rule gates above, plus `benchCompile` for
benchmark changes, `breezeGoldensCheck` for corpus changes, `vectorBackendTest`
for backend changes and `benchSmokeJSFull` for shared-kernel performance changes.
Run affected gates after the final edit, then retain their exact-source receipts.

Refresh the final scoreboard after accepted optimizations. Close relevant Motes
only when their acceptance evidence is present; documented exceptions remain
visible. Native-backed Breeze and additional CPU families need separate explicit
receipts before broadening the claim; provisioning new native backends is outside
this follow-up. Preserve unrelated main-worktree files and stage only owned paths.

Next action: implement step 1 on an isolated branch from current integration,
while coordinating benchmark-sized eig validation with the existing HC-1 owner.

## Before the final two-lane receipt

Do not cut `benchmarks/results/2026-1x-xx-breeze-two-lane-final.{json,md}`
until all of the following hold:

1. **In-flight work has merged or been dropped:** hc-max2, hc-gemv-vec and
   hc-lbfgs-step.
2. **Every Milestone 3 gap is closed or documented here with its cause**
   (see the dashboard's "Still behind" list):
   - lane A `max` at n = 1M, 0.76× (HC-3, `bd-01M4GY1QTAE17DV1A1990V88B8`).
     This one fails an epic criterion: a new op family below 1.0× at n ≥ 1K;
   - lane B vector `gemv`/`gemvT` (0.95–0.96×) and the gemv-bound L-BFGS
     logistic objective at a fixed budget (0.96×) (HC-4,
     `bd-01M4PAPJGS4EF4DXMEAS8NVY01`);
   - `dot` at n = 64K, 0.97× in both lanes: the paired rerun
     (`bd-01M4PAPKQTBWWQSVC5PJ5R9CFT`) decides whether this is a kernel gap;
   - lane B sparse `dot` with a dense vector, 0.96× (no tracker item yet);
   - `sigmoid` at n = 1K and 64K is already documented below for lane A;
     lane B measures 0.95–0.96× for the same backend-insensitive kernel.
3. **The new op families are settled at ≥ 1.0×.** Several lane-A rows are
   ties at 1.00× (`exp`, `logSumExp`, sparse `dot` with a dense vector). The
   reviewer must decide whether a CI-overlapping tie meets the bar.
4. **A fresh sweep of the final source SHA:** both lanes plus the eig16 and
   lane-A hugemethods supplements, on the reviewed runner, with validation
   sidecars bound to that SHA and scoreboards run with `--validation --strict`
   and zero unpaired rows.
5. **W5 suites are green on CI** for that SHA, and `docsCheck` passes.
6. **A fresh-context reviewer signs off** on the final receipt; then
   regenerate the dashboard from it and close W6 and the epic.

## Documented gaps

Accepted exceptions to the outcome above, each with its measured cause.

- **`sigmoid`, n = 1K and 64K: 0.96–0.97× Breeze (lane A).** The hybrid kernel
  (70fb262) uses Breeze's form `1/(1 + exp(-x))` for `x >= -36` and `t/(1 + t)`
  with `t = exp(x)` below. The per-element range test is the whole gap: on
  trillium, the kernel without it (`rawDirect`) is 4.4% faster at 1K, 3.7% at
  64K and 1.9% at 1M, and gale's facade is within 0.4% of the bare kernel.
  Both bit-identical prescans that would skip the test were slower (3.0–9.4%).
  The test is kept because Breeze's form overflows `exp(-x)` below about -709.78
  and returns 0 there, where the true sigmoid is the (subnormal) `exp(x)`; it
  also loses up to 2.4 ulps below -36. At n = 1M gale is 1.42× Breeze, whose map
  slows to 5.10 ns per element there against gale's 3.59. Receipts:
  `benchmarks/results/2026-10-11-hc-sigmoid-trillium-r2.md` and
  `benchmarks/results/2026-10-11-hc-sigmoid2-path-trillium.md`;
  `SigmoidDigestSuite` pins the kernel's bits.
