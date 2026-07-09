# Gale Benchmarks

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

## Scala.js smoke runner

The Scala.js runner is a wall-clock smoke harness (not JMH):

```bash
sbt benchSmokeJS        # fastLinkJS
sbt benchSmokeJSFull    # fullLinkJS (optimised) run
```

The benchmark harness is intentionally small. It exists to keep performance
regressions visible while the kernel layer is still taking shape.
