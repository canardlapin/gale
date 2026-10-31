# Gale

[Getting started](docs/user/getting-started.md) · [Guides](docs/user/guides/examples.md) ·
[Reference](docs/user/reference/index.md) · [Benchmarks](benchmarks/dashboard.md)

Gale is a Scala 3 linear algebra library for building numerical applications
on the JVM and Scala.js. It provides dense and sparse matrices, linear solvers,
spectral methods, and optimization in a portable core. Numerical solves return
typed failures and iterative methods expose convergence diagnostics.

> **Status:** `0.1.0-M1` is available from Maven Central for JVM and Scala.js.
> This is the first milestone release, not a 1.0 compatibility freeze.
> See the [release manifest](docs/release-manifest.md) for the milestone's scope
> and [API stability policy](docs/api-stability.md) for compatibility boundaries.

## Quick start

Use Scala **3.7.4** and JDK **21** for the baseline build. Add the JVM core to
`build.sbt`:

```scala
scalaVersion := "3.7.4"
libraryDependencies += "io.github.canardlapin" %% "gale-core" % "0.1.0-M1"
```

For Scala.js or an sbt cross-project, use `%%%` instead of `%%` with the
appropriate Scala.js/cross-project plugin enabled. The portable core needs no
Vector API import or native library. For source integration, see
[the pinned source setup](docs/user/getting-started.md#add-gale-from-source).

This example fits an intercept and slope to observations following `1 + 2x`:

```scala
import gale.linalg.*

val design = Matrix(4, 2)(
  1.0, 0.0,
  1.0, 1.0,
  1.0, 2.0,
  1.0, 3.0
)
val observations = Vec(1.0, 3.0, 5.0, 7.0)

val coefficients = design.leastSquares(observations).orThrow
coefficients.toSeq // Seq(1.0, 2.0), within floating-point error
```

`leastSquares` returns `Either[LinAlgError, DVec]`; `.orThrow` is convenient in
a first example or test. Application and library code can keep the `Either` and
handle rank deficiency or a shape mismatch explicitly.

The same workflow is compiled and executed in the
[getting-started guide](docs/user/getting-started.md).

## What Gale covers

- **Dense and structured algebra:** vectors, matrices, strided views, sparse
  COO/CSR/CSC formats, banded and specialized matrix representations.
- **Linear systems:** square and least-squares solves, reusable LU/Cholesky/QR
  factors, and matrix-free CG, BiCGSTAB, GMRES, CGNR, and LSQR.
- **Spectral methods:** dense eigenvalue and singular-value decompositions,
  partial iterative methods, and generalized symmetric eigenproblems.
- **Optimization:** first-order and proximal methods, nonlinear least squares,
  and box-constrained local minimization. Start with the
  [optimization guide](docs/user/guides/nonlinear-and-bounded-optimization.md).
- **Numerical control:** exact accumulation, paired residual comparisons,
  owned results, reusable workspaces, and explicit destination buffers. See
  [fit comparisons](docs/user/guides/numerical-workspaces.md) and
  [allocation control](docs/user/advanced/index.md).
- **Optional integrations:** Breeze conversions, Ravel copy conversions,
  JVM Vector API acceleration, and FFM BLAS/LAPACK backends.

The portable core is the default. Acceleration requires an explicit backend
import and runtime setup; see [Enabling SIMD](docs/user/advanced/numerical-contract.md#enabling-simd).
Breeze interop and JVM backends are provisional, outside the first milestone
bundle. The Ravel adapter is source-only and resolves the exact Ravel revision
in `build.sbt`; local provider development requires an explicit
`-Dgale.ravel.build=/absolute/path/to/ravel` override.

## Fit and boundaries

Gale targets the real-`Double` linear algebra slice used by scientific and data
libraries, plus a bounded owned complex Hermitian tier in `gale.complex`.
The [Hermitian guide](docs/user/guides/complex-hermitian.md) defines its portable
PSD factor, rank-aware solve and capacity contracts. Gale is not source-compatible
with Breeze and does not replace its statistics, probability, signal-processing,
plotting, machine learning, tensor, or general complex-number modules.

Portable sparse Cholesky is available through the explicit
`gale.sparse.direct.pure.given` import or `SparseDirectProvider.pure`; the default
provider advertises no sparse-direct capability. Sparse direct LU and QR remain
unimplemented. Dense SVD returns the complete economy decomposition.

The shared code targets Scala 3.7.4. Core JVM use requires JDK 21; finalized
FFM modules require JDK 22 and explicit native access. CI defines separate JVM
and Scala.js/Node test lanes. Experimental WebAssembly is default-off and is
not part of the current compatibility or performance promise.

Floating-point algorithms preserve documented shapes, ordering, failure modes,
and numerical invariants—not bit identity across legal implementations.
Backends may reassociate arithmetic, and partial iterative results must be
interpreted through their convergence diagnostics.

## Documentation

- [Getting started](docs/user/getting-started.md) — setup, a regression fit,
  error handling, and reusable factors.
- [Core concepts](docs/user/core-concepts.md) — values, failures, operators,
  diagnostics, and backend behavior.
- [Worked examples](docs/user/guides/examples.md) — dense, sparse, spectral,
  optimization, and migration workflows.
- [Breeze migration](docs/user/guides/breeze-equivalence.md) — API equivalents
  and differences.
- [Advanced topics](docs/user/advanced/index.md) — numerical guarantees,
  ownership, allocation control, and matrix-free generalized eigensolving.
- [Modules and platforms](docs/user/reference/modules-and-platforms.md) —
  build projects, runtimes, and artifact boundaries.
- [Troubleshooting](docs/user/troubleshooting.md) — shape, rank, convergence,
  backend, and runtime failures.

Run `sbt docsCheck` to build JVM/Scala.js Scaladoc and render the executable
mdoc/Laika guides. The same compiled guides are published as a Starlight site
from [`website/`](website/README.md) to GitHub Pages; there is no stable hosted
Scaladoc URL yet.

## Development

Use JDK 21, sbt 1.11.7 (pinned in `project/build.properties`), and Node 22 for
Scala.js tests. Keep the repository's `.jvmopts` in place.

| Command | What it checks |
| --- | --- |
| `sbt compileAll` | Core and laws on JVM/JS, then local core publication and the Scala 3.8.4 consumer compile probe |
| `sbt scalafmtCheckAll` | Formatting within the configured release scope |
| `sbt testAll` | Core and laws tests on JVM and Scala.js |
| `sbt testAllFull` | `testAll`, then optimized linking of JS test bundles (not optimized test execution) |
| `sbt parityTest interopBreezeTest` | Breeze differential parity and conversion tests |
| `sbt interopRavelTest` | Ravel copy conversions on JVM and Scala.js |
| `sbt vectorBackendTest` | Optional Vector API backend |
| `sbt benchCompile` | JVM and JS benchmark compilation; no timing claim |
| `sbt docs/mdoc` | Compile and execute guide examples |
| `sbt docsCheck` | Scaladoc plus executable, rendered guides |

`compileAll` deliberately scopes compilation; plain root `sbt compile` also
aggregates optional native modules. For those, select JDK 22 and install a
supported BLAS/LAPACK runtime (CI uses OpenBLAS), then run:

```sh
sbt nativeBackendTest blasFfmBackendTest benchFfmCompile
```

The test configuration enables native access for the forked native tests.
Applications using FFM must also enable it explicitly. Experimental WebAssembly
uses `GALE_WASM=1` and Node 25+; it is an advisory CI lane, outside the current
compatibility promise. See the [CI workflow](.github/workflows/ci.yml),
[release policy](docs/release-policy.md), and
[benchmark dashboard](benchmarks/dashboard.md) for the full checks and evidence.

## License

[Apache-2.0](LICENSE)
