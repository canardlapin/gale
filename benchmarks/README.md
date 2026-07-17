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

Each benchmark uses 2 forks with 5x500ms warmup and 5x500ms measurement.

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

Run the whole Breeze sweep (all pairs, all sizes, default 2 forks):

```bash
sbt "benchmarksJVM/Jmh/run .*Breeze.*"
```

Run one class, or a single pair at one size with the allocation profiler:

```bash
sbt "benchmarksJVM/Jmh/run .*FactorizationBreezeJmh.*"
sbt "benchmarksJVM/Jmh/run -prof gc -p n=256 .*BlasL3BreezeJmh.*Gemm$"
```

### Pure-JVM baseline caveat

Breeze here runs on its **pure-Java netlib fallback** (`dev.ludovic.netlib`'s
F2J/Java BLAS — it logs "native BLAS not found … return java instance" at
startup). That is the deliberate comparison target: gale is a pure-JVM / Scala.js
library, so the honest baseline is pure-JVM Breeze. **Native-BLAS Breeze (system
OpenBLAS/MKL via JNI) is a separate, much faster target and is out of scope
here** — closing that gap is deferred to gale's v0.5 acceleration work. Read every
number strictly against the pure-JVM baseline, and note that pure-Java netlib is
itself a mature, heavily loop-unrolled BLAS — beating it is a real bar, not a
formality.

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

The benchmark harness is intentionally small. It exists to keep performance
regressions visible while the kernel layer is still taking shape.
