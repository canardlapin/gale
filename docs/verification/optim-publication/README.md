# Optimization publication integration, 2026-10-07

The optimization portfolio was rebased onto upstream `2fada93`, which adds
`BoxQuasiNewton` and ExactSum resource estimates. Upstream already owns the
finite vector `BoxBounds` API. The new matrix-shaped bounds are therefore named
`MatrixBoxBounds`; L-BFGS-B, its shared line search, tests, benchmark, and public
guide use that name. The existing vector API and its callers are preserved.
No solver algorithm, numerical tolerance, fixture, or benchmark timing boundary
changed during integration.

The [performance report](../optim-performance/README.md) remains the historical
measured snapshot. Its source hashes precede this integration; its 98 artifact
hashes still match byte-for-byte. Raw log whitespace and machine-generated JSON
and CSV formatting are preserved for that reason. This directory records checks
of the integrated source using [source and artifact hashes](source-manifest.json).

All required portable gates passed with exit status zero: `scalafmtCheckAll`,
`compileAll`, `testAllFull`, `parityTest`, `interopBreezeTest`, `docsCheck`,
`benchCompile`, and the optimization tests under `FullOptStage`. Counts are
883 core JVM tests, 870 core JavaScript tests, 54 law tests per platform,
85 Breeze parity tests, 29 interop tests, and 96 optimized JavaScript optimization
tests. Existing documentation warnings remain. Full logs and command/exit
metadata are in [logs/](logs/integration-gates.log).

Both comparison runners were rerun using the same pinned Python environment:
all 19 JVM and all 19 optimized Node endpoints pass independent objective,
feasibility, and stationarity checks. The ten extension cases also pass their
reported certificate checks. SciPy's large diagonal and rotated bounded cases
still fail the common endpoint quality gate. The new standalone samples are
retained as integration evidence; the matched forked JMH results in the
performance report were not rerun or relabeled.

- [JVM first-order checks](firstorder/gale-jvm-comparison-v1.md)
- [Node first-order checks](firstorder/gale-node-comparison-v1.md)
- [JVM LM and L-BFGS-B checks](extensions/jvm-comparison.md)
- [Node LM and L-BFGS-B checks](extensions/node-comparison.md)

Hosted CI, downstream adapter migration, native backend, and Wasm qualification
remain outside these local results.
