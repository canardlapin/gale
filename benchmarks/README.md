# Gale Benchmarks

The current cross-backend summary and evidence links live in the
[backend dashboard](dashboard.md).

JMH (JVM) benchmarks covering the dense kernels plus sparse and solver scenarios:

- `DenseKernelJmh` — `dot`, `axpy`, `gemv` over `n` in {256, 4096, 65536}.
- `GemmJmh` — dense `n x n` matrix product over `n` in {32, 128, 512}
  (`n = 512` hits the blocked row-major path).
- `SpmvJmh` — CSR sparse matrix-vector product at 1% density over `n` in
  {1024, 16384}.
- `SolverJmh` — conjugate gradient on a 2D Laplacian (64x64 grid, 4096x4096 SPD).
- `BlockSymmetricEigenJmh` — matrix-free block Krylov with an eight-dimensional
  repeated top eigenspace over `n` in {128, 512, 2048}.
- `GeneralizedLobpcgJmh` — matrix-free generalized symmetric LOBPCG over
  `n` in {128, 512, 2048}, `k` in {4, 8, 16}, clustered diagonal and
  stiffness/mass pencils, and identity/Jacobi/block-Jacobi preconditioning.
  `GeneralizedLobpcgWorkReceipt` records exact iteration and A/B/preconditioner
  work separately from JMH timing.
- `DenseTransformJmh` — boxed matrix/vector exports versus direct primitive
  copies from transposed and strided views over 1,024 and 16,384 values; use
  `-prof gc` for bytes/op.
- `SparseInteropJmh` — boxed `COO.entries` versus reusable primitive COO/CSR
  traversal over 1,024 and 16,384 stored entries; use `-prof gc` for bytes/op.
- `AllocationArchitectureJmh` — dense destination, factorization workspace,
  spectral scratch, and sparse structure/value allocation baselines over `n` in
  {64, 128}; use `-prof gc` and compare the paired allocating/reuse scenarios.
- `QrMultiRhsJmh` — copy-inclusive pivoted-QR factor application and
  least-squares matrix RHS solves over `n` in {512, 2048, 10000}, `p` in
  {6, 24}, and `q` in {1, 8, 16, 32, 100}, plus a protected factorization
  control. Owned and reusable-workspace solve routes are measured for both
  matrix and vector RHS inputs. The exact regress4s fixture is
  `n=2048,p=6,q=16`.
- `TallPivotedQrJmh` — copy-inclusive portable column-pivoted factorization over
  `n` in {512, 1024, 2048, 4096, 10000} and `p` in {3, 5, 6, 8, 16, 24}; the
  exact regress4s factor shapes include `4096x6`, `2048x6`, and `1024x5`.
- `ScaledQrConstructionJmh` — transient builder consumption, algebraic
  row-scaled QR, and scaled matrix-RHS workspace solves at the exact regress4s
  M5 shape `n=1024,p=5,q=8`.

Most established kernel benchmarks use 2 forks with 5x500ms warmup and
5x500ms measurement. `GeneralizedLobpcgJmh` uses one fork, 1x200ms warmup, and
2x200ms measurement because its 54-scenario product is an end-to-end solver
matrix; its receipt distinguishes the short development sweep from a
trustworthy default-settings run.

Run everything:

```bash
sbt "benchmarksJVM/Jmh/run"
```

Run a single benchmark (regex) at one size, e.g. a quick smoke of `dot`:

```bash
sbt "benchmarksJVM/Jmh/run -f 0 -wi 1 -i 1 -p n=256 gale.bench.DenseKernelJmh.dot"
```

## Allocation profiling

JMH ships a GC profiler that reports bytes allocated per operation. Add
`-prof gc` to any run to surface allocation regressions (the dense kernels and
SpMV should report near-zero `gc.alloc.rate.norm`):

```bash
sbt "benchmarksJVM/Jmh/run -prof gc -p n=4096 gale.bench.DenseKernelJmh.dot"
sbt "benchmarksJVM/Jmh/run -prof gc .*DenseTransformJmh.*"
sbt "benchmarksJVM/Jmh/run -prof gc .*SparseInteropJmh.*"
sbt "benchmarksJVM/Jmh/run -prof gc -p n=128 gale.bench.AllocationArchitectureJmh.*"
sbt "benchmarksJVM/Jmh/run -prof gc -p n=2048 -p p=6 -p q=16 gale.bench.QrMultiRhsJmh.*Owned"
sbt "benchmarksJVM/Jmh/run -prof gc -p n=2048 -p p=6 -p q=16 gale.bench.QrMultiRhsJmh.solveLeastSquares.*"
sbt "benchmarksJVM/Jmh/run -prof gc -p n=4096 -p p=6 gale.bench.TallPivotedQrJmh.factorPivotedQr"
sbt "benchmarksJVM/Jmh/run -prof gc gale.bench.ScaledQrConstructionJmh.*"
sbt "benchmarksJVM/Jmh/run -prof gc gale.bench.GeneralizedLobpcgJmh.solve"
sbt "benchmarksJVM/Jmh/runMain gale.bench.GeneralizedLobpcgWorkReceipt"
```

Verify the whole suite compiles (annotation processing included) without running
it:

```bash
sbt benchCompile
```

## JDK 22 FFM benchmarks

The separate `benchmarksFfm` project keeps JDK 22 FFM sources out of the JDK 21
benchmark build:

- `FfmGemmJmh` — heap-copy-inclusive GEMM plus a copy-free `NativeDMat` control.
- `FfmGemvJmh` — heap-copy-inclusive standalone matrix-vector multiply.
- `FfmLapackJmh` — LU, Cholesky, tall QR, and symmetric eigenvalues.
- `FfmSolverScenarioJmh` — end-to-end dense solve and tall least squares.
- `FfmDispatchJmh` — warmed cost of selecting and declining a native provider.

```bash
sbt "benchmarksFfm/Jmh/run .*FfmGemvJmh.*"
sbt "benchmarksFfm/Jmh/run .*FfmLapackJmh.*"
sbt "benchmarksFfm/Jmh/run -prof gc -p n=128 .*FfmSolverScenarioJmh.*"
```

All native crossover receipts include heap/native copies unless explicitly
labelled `NativeDMat`. A backend is enabled by default only when a complete
two-fork sweep supports a conservative threshold on that library family.

## Breeze comparison (paired gale-vs-Breeze)

`benchmarksJVM` also carries Breeze 2.1.0 in **compile scope for this module only**
(it is `publish`-skipped and never a dependency of `core`/`laws`, so gale-core
stays 100% Breeze-free). Each operation has one `@Benchmark` per library over
identical seeded `@Setup` data at matching `@Param` sizes:

- `BlasL1BreezeJmh` — `dot`, in-place `axpy`, 2-`norm` (`n` in {65536, 262144, 1048576}).
- `BlasL2BreezeJmh` — `gemv` and transpose `gemvT`, both allocating (`n` in {256, 1024, 2048}).
- `BlasL3BreezeJmh` — square `gemm`, tall `gemmTall` (`4n x n`), transpose-product `AtA` (`n` in {16, 64, 256}).
- `FactorizationBreezeJmh` — `solve`, `lu` (factorization only: gale `lu` vs breeze `LU.primitive`/`dgetrf`), `chol`, `qr` (no `Q` materialised) (`n` in {16, 64, 256}).
- `LeastSquaresBreezeJmh` — overdetermined `m = 4n` least-squares: gale `leastSquares` vs breeze backslash (`n` in {16, 64, 256}).
- `SymEigenBreezeJmh` — symmetric eigen with vectors: gale `Eigen.eigSymmetric(All)` vs breeze `eigSym` (`n` in {16, 64, 128}).

Backend-sensitive `gale*` methods also take a `GaleBackendState` whose
`@Param backend` is `pure` or `vector` (the `backend-jvm-vector` SIMD backend);
the Breeze twins do not, so each Breeze row runs once per size. Coverage:

| gale methods | backend param | scoreboard label |
|---|---|---|
| `gemv`, `gemvT`, `gemm`, `gemmTall`, `AtA` | `pure`, `vector` | `pure` / `vector` |
| `lu`, `chol`, `solve`, `qr`, `lstsq` | `pure`, `vector` | `vector (gemm-routed only)`: the backend is reached only where a product routes through its gemm |
| `dot`, `axpy`, `norm`, `eigSym` | none | `backend-insensitive`: the operation takes no `Backend` (eigen resolves a `SpectralBackend` the Vector backend does not supply), so these always run pure gale; W2 brings them back when L1 routing lands |

### Two lanes

Breeze's `dev.ludovic.netlib` picks native JNI BLAS first, then `VectorBLAS`
whenever `jdk.incubator.vector` is resolvable, then scalar Java BLAS
(`Java11BLAS`); LAPACK is `F2jLAPACK` unless native LAPACK loads. The JMH config
adds `--add-modules=jdk.incubator.vector`, and forks inherit it, so a plain
`Jmh/run` gives Breeze **SIMD VectorBLAS**. Run one of the two lanes instead:

| Lane | sbt alias | Fork JVM args | Breeze BLAS | gale backend |
|---|---|---|---|---|
| A, out-of-box | `breezeLaneA` | `-jvmArgs -Dgale.bench.lane=A` (replaces **all** inherited fork args) | scalar Java BLAS | `pure` |
| B, SIMD | `breezeLaneB` | inherited `--add-modules` plus `-Dgale.bench.lane=B` | `VectorBLAS` | `pure` and `vector` |

Both aliases set `-rf json`; append a result file, JMH options and a benchmark
regex. Use one pinned JDK (25 LTS) for both lanes:

```bash
sbt "breezeLaneA -rff target/breeze-laneA.json .*BreezeJmh.*"
sbt "breezeLaneB -rff target/breeze-laneB.json .*BreezeJmh.*"
```

`-p backend=vector` in lane A fails that trial with an explicit
`IllegalStateException` (the Vector API module is absent); the other benchmarks
still run. Run single pairs with the same aliases, for example
`sbt "breezeLaneA -prof gc -p n=256 .*BlasL3BreezeJmh.*Gemm$"`.

Each Breeze bench records the netlib classes it resolved, once per trial, to
stderr (`[breeze-netlib] ...`) and as one JSON line appended to the sidecar
`benchmarks/jvm/target/breeze-netlib.jsonl` (the fork's working directory is
`benchmarks/jvm`; override with `-Dgale.bench.netlibSidecar=<path>`). Each record
also carries the JDK version and whether `jdk.incubator.vector` was resolved. Both
lane aliases delete the sidecar first (`benchmarksJVM/breezeNetlibReset`), so every
run writes a fresh one: copy it next to the JSON receipt before the next run.
Running the two lanes concurrently in one checkout is unsupported, because the
reset and both lanes share that fixed path: use a separate worktree per lane, or
give each run its own `-Dgale.bench.netlibSidecar` (which the reset does not delete).

### Scoreboard

`tools/bench/breeze_scoreboard.py` turns a lane's JMH JSON and sidecar into the
markdown scoreboard, replacing hand-written tables (`--strict` fails on unpaired
rows):

```bash
python3 tools/bench/breeze_scoreboard.py --lane A \
  --netlib benchmarks/jvm/target/breeze-netlib.jsonl target/breeze-laneA.json \
  -o benchmarks/results/<date>-breeze-laneA.md
python3 -I tools/bench/test_breeze_scoreboard.py   # self-test on synthetic fixtures
```

Pairing convention: a paired benchmark is `<Class>.gale<Op>` / `<Class>.breeze<Op>`
(`<Op>` upper-case first); rows pair on class, op, and params with gale's
`backend` removed. New Breeze benches must follow it. The ratio is gale speed over
Breeze speed (>1 means gale is faster; time modes are inverted), and the verdict
is `ahead`/`behind` only when the 99.9% confidence intervals do not overlap. The
header records the lane, JDK, fork JVM args, netlib BLAS/LAPACK classes, and
commit. Any receipt is rejected when the sidecar and results disagree (JDK, a
result without a sidecar record, or a stale record without a result) or a score
is not positive. A lane-A receipt is **rejected** if Breeze resolved `VectorBLAS`
or a native BLAS/LAPACK (a Linux host with `libblas.so.3` would otherwise quietly
get native BLAS), if any fork resolved the Vector module (via the JMH args or, for
example, `JDK_JAVA_OPTIONS`), or if gale ran a non-`pure` backend; a lane-B receipt
is rejected unless every fork had the module and Breeze resolved `VectorBLAS`.

### What each lane measures

Lane A compares portable pure-JVM gale with out-of-box scalar Breeze. Lane B
compares both gale backends with Breeze on its SIMD `VectorBLAS`. **Native-BLAS
Breeze (system OpenBLAS/MKL via JNI) is a separate, much faster target and is out
of scope here.** Note that netlib's Java BLAS is itself a mature, heavily
loop-unrolled BLAS, so beating it is a real bar. Receipts before the two lanes
existed did not record Breeze's BLAS class; see the correction in
`results/2026-07-11-breeze-release-grade.md`.

### Sessions only smoke-compile

CI / agent sessions never run a full sweep. They prove the harness builds
(`sbt benchCompile`, which also runs JMH annotation processing) and, at most, run
a single tiny unforked pair as a liveness check. Trustworthy numbers require a
quiet machine, the default 2 forks, and full warmup — an unforked `-f 0` run (what
a liveness check uses) is explicitly *not* a measurement.

## Scala.js smoke runner

The Scala.js runner is a wall-clock smoke harness (not JMH):

```bash
sbt benchSmokeJS        # fastLinkJS
sbt benchSmokeJSFull    # fullLinkJS (optimised) run
```

The runner reports explicit owned-result, destination/workspace reuse, immutable
pattern reuse, and symbolic-plan analysis/reuse counts per operation. Those are
public-contract construction counters, not a claim about every allocation
performed by the JavaScript engine. The benchmark harness is intentionally
small: it keeps performance and construction regressions visible while the
kernel layer is still taking shape.
