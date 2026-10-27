# Axis margins for Gale 0.2

Status: first implementation slice on `feature/0.2-axis-margins`. These APIs are
not part of the published 0.1.0-M1 artifact. No release/tag/publication is part of
this change; admission to the 0.2 release workflow is a separate release task.

## Contract

`Margin[A <: Axis]` is invariant and distinct from `DVec`. It stores an axis and
an immutable vector. `v.along(axis)` returns `Margin[axis.type]` without copying
values; `toDVec` explicitly erases orientation, again without copying. Mapping,
scalar arithmetic, and same-axis addition/subtraction/zipMap retain orientation.
Known opposite axes fail to compile; abstract/runtime axes are also checked
before combination. Axis types do not establish semantic identity or ordering
of features from unrelated data sets.

`Margin.ByRow` and `Margin.ByCol` name the two fixed-axis types. Fused operations
require a shared axis type: independent `Margin[?]` fields deliberately cannot
be combined. Store fixed-axis aliases or retain both margins in a generic
parameter object to preserve their relationship. Runtime choices work when
both operands are tagged using the same stable axis value. A wildcard fallback
would admit statically opposite axes, so none is provided.

`byRow`, `byCol`, and `byAxis(axis)` expose tagged sum/mean/min/max reductions.
The original `DMat` reduction signatures and algorithms remain unchanged.
`std`/variance are not introduced here: their statistical and numerical contracts
need a separate slice. This implementation accepts caller-supplied scales.
The guide leads with centering and weights, and identifies its supplied-scale
example as calibration rather than estimated z-scores.

Matrix-margin addition/subtraction are direct. Multiplication/division require
`pointwise`, preserving matrix-vector multiplication. Scalar pointwise arithmetic
is explicit too. There are no matrix-singleton-dimension expansion rules, no
plain-vector broadcasting, no implicit outer products, and no sparse broadcasts.

Scalar pointwise multiplication delegates to existing `DMat.*`; scalar addition,
subtraction, and literal division use `DMat`'s unary kernel path with an owned
builder result. Division cannot delegate to multiplication by a reciprocal:
that changes rounding and can overflow for subnormal divisors.

A result has exactly the input matrix shape. Every operation checks the margin
length against its named axis, including zero-size matrices. Immutable results
own one new matrix buffer. Tagging owns one small wrapper, not repeated elements.
The public positive-stride contract is unchanged; zero read increments exist
only inside broadcast traversal. No general lazy expression graph is introduced.

## Numerical and ownership rules

`affine(scale, shift)` evaluates `(x * scale) + shift` in that order.
`standardize(center, scale)` evaluates `(x - center) / scale` in that order.
Both make one traversal when supplied their parameters; computing those parameters
has separate costs. Neither rewrites division as reciprocal multiplication,
contracts FMA, or suppresses reads for a zero coefficient. NaN, infinities, signed
zeros, and zero/negative scales follow the written IEEE operations. This is not
a claim of identical arbitrary reductions/backends or NaN payload bits.

`addInto`, `subInto`, pointwise `mulInto`/`divInto`, `affineInto`, and
`standardizeInto` replace an open `DMatBuilder`. They validate axes, dimensions,
closed state, and storage aliasing before any write. A failed validation leaves
an open destination unchanged. After `result()`, the builder is closed. There is
no new mutable matrix type or externally writable broadcast view.

## Evidence and deferred work

Shared JVM/JS suites cover compiler boundaries, square/rectangular matrices,
empty dimensions, strided vectors, transposed/sliced matrices, destination errors,
alias rejection, independence, IEEE edge cases, and order-sensitive examples.
The subtraction/division regression `(10000000000000002 - 10000000000000000) / 3`
must yield `2/3`; reciprocal-affine arithmetic yields `0.5` instead.

`MarginJmh` compares chaining, fused allocating operations, destination reuse,
and an explicit loop for both axes and row-major/transposed layouts. Timing and
allocation results must be reported separately from compilation/test success.

Deferred: variance/std algorithms, sized margins carrying both axis and length,
semantic feature axes, in-place builder arithmetic, enhanced existing Diagonal
operator paths, and matrix-free sparse centering. These should reuse the same
contracts without weakening the established vector/matrix operation meanings.

## Initial local qualification, 27 October 2026

Implementation base: `80dd605` (the owner's documentation commit). Feature branch:
`feature/0.2-axis-margins`. This is a local implementation, not a published 0.2
artifact. The release/version configuration remains unchanged.

- JDK 21.0.13, Node 22.23.4, macOS ARM64.
- Full core/law tests: 1,298 JVM core, 1,272 JS core, and 54 laws on each platform
  (2,678 total), passing again after the row-read optimization.
- Twelve new test groups cover numerical/ownership behavior and compiler
  boundaries. They also pass under full-optimized Scala.js, alongside the four
  existing syntax tests.
- Breeze parity: 182 passed; Breeze interop: 29 passed. These ran before the
  final row-read optimization; the changed margin paths were rechecked afterward.
- JVM/JS Scaladoc, guide-site generation, 37 executable guide files, benchmark
  compilation, and formatting checks passed. New margin source/test/benchmark
  files are explicitly included in the repository's formatting gate.

Bounded JMH measurements used 256×256 matrices, both axes, row-major/transposed
layouts, one fork, two 500 ms warmups and three 500 ms measurements, with the GC
profiler and a 512 MiB heap. Final measured bytes/op across those four cases:

| Method | Bytes per operation |
| --- | ---: |
| Two-operation chain | 1,048,752–1,048,754 |
| Fused standardize | 524,377–524,378 |
| Explicit allocating loop | 524,377–524,378 |
| Reused destination | 0.38–0.93 |

These measurements support two result buffers versus one, and no matrix-buffer
allocation in the destination path. Sub-byte figures are profiler overhead,
not a claim that JVM execution literally allocates fractional objects.
Latency changed substantially between the two short runs, including for the
unchanged explicit-loop control. They do not support a general speedup claim.
Row-wise paths now load constant margin values once per row, reusing the existing
map kernel; column-wise paths retain strided traversal. Both implementations'
receipts are retained rather than discarding the adverse/noisy results.

Raw evidence: `.agent-work/axis-margins-20261027/`, including `court.log`,
`final-court.log`, `hoisted-checks.log`, `format-final.log`, `benchmark.json`, and
`benchmark-hoisted.json`. `source-manifest.json` binds implementation/test hashes.

## Review revisions, 27 October 2026

- Added `Margin.ByRow`/`Margin.ByCol` aliases and standardized the column facade
  spelling to `byCol` throughout source, examples, tests, and benchmarks.
- Made the guide's missing variance/std estimator explicit, separated centering
  from supplied calibration, and marked the guide and index as an unreleased
  preview requiring a source checkout.
- Documented the intentional rejection of independently erased `Margin[?]`
  operands in fused calls. Added executable examples and compiler tests for
  retaining the shared axis in stored generic parameter objects.
- Removed scalar traversal from `MarginOps`. Multiplication reuses `DMat.*`;
  addition, subtraction, and literal division use `dmapInto` with a builder.
  Scalar checks now include transposes, slices, empty matrices, and the finite
  result of dividing a subnormal by itself.
- Clarified strided reads from transposed inputs and limited performance claims
  to result-buffer allocation. The broadcast/fused traversal did not change.

Reverification used one sbt process at a time: 17 focused tests passed on each of
JVM, fast-optimized Scala.js, and full-optimized Scala.js (13 margin groups plus
4 existing syntax groups). Formatting, benchmark compilation, all 37 executable
guide files, JVM/JS Scaladoc, and guide-site generation passed. The earlier full
2,678-test run, parity runs, and JMH results above were not repeated for this
revision; their original receipts remain separate.

Review evidence: `review-checks.log`, `review-optimized-docs.log`, their exit-status
metadata, and `review-source-manifest.json` under the same local evidence directory.
