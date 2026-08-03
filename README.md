# Gale

[Guide](docs/user/index.md) · [Getting started](docs/user/getting-started.md) ·
[API reference](docs/user/reference/index.md) ·
[Breeze migration](docs/user/guides/breeze-equivalence.md) ·
[Benchmarks](benchmarks/dashboard.md)

Gale is a linear algebra library for Scala 3 that runs the same dense, sparse,
solver, and spectral API on the JVM and Scala.js. Use it when a numerical
library needs real-`Double` matrices, typed failure and convergence evidence,
or matrix-free algorithms without taking a JVM-only dependency.

> **Status:** Pre-release and source-only. The build is `1.0.0-SNAPSHOT`; Gale
> is not available from Maven Central yet.

## Quick start

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

## Use Gale before the first release

Pin a commit as an sbt source dependency. Choose `coreJVM` or `coreJS` for a
single-platform build:

```scala
lazy val galeRevision = "<commit sha>"
lazy val galeBuild = uri(
  s"https://github.com/canardlapin/gale.git#$galeRevision"
)
lazy val galeCoreJVM = ProjectRef(galeBuild, "coreJVM")

lazy val app = project.dependsOn(galeCoreJVM)
```

For a cross-project, use the matching `coreJVM` and `coreJS` project references
on its platform projections. A local checkout can instead run
`coreJVM/publishLocal` and `coreJS/publishLocal`.

The intended public dependency, once a release exists, is:

```scala
libraryDependencies +=
  "io.github.canardlapin" %%% "gale-core" % "<published-version>"
```

Do not substitute `1.0.0-SNAPSHOT` unless that snapshot was actually published
to the resolver used by your build.

## What Gale covers

- Build and transform immutable-facing dense vectors and matrices, including
  strided transpose, row, column, and slice views.
- Solve square and least-squares systems with reusable LU, Cholesky, or QR
  factors and vector or matrix right-hand sides.
- Store sparse matrices as COO, CSR, CSC, banded, diagonal, identity,
  permutation, or zero structures, with explicit canonicalization rules.
- Run CG, BiCGSTAB, GMRES, CGNR, and LSQR against dense matrices, sparse
  matrices, or custom `DoubleLinearOperator`s.
- Compute dense and partial eigenvalue and singular-value decompositions,
  including matrix-free generalized symmetric problems through LOBPCG and
  generalized block Lanczos.
- Opt into allocation-controlled builders, destinations, and workspaces while
  keeping ordinary returned results owned.
- Add JVM Vector API or FFM BLAS/LAPACK acceleration explicitly; the portable
  core remains the default and Scala.js stays independent.

## Fit and boundaries

Gale targets the real-`Double` linear algebra slice used by scientific and data
libraries. It is not source-compatible with Breeze and does not replace
Breeze's statistics, probability, signal-processing, plotting, machine
learning, tensor, or general complex-number modules. Sparse direct
factorization is not implemented in the current core.

The shared code targets Scala 3.7.4. Core JVM use requires JDK 21; finalized
FFM modules require JDK 22 and explicit native access. CI executes the JVM and
Scala.js/Node lanes independently. Experimental WebAssembly is default-off and
is not part of the current compatibility or performance promise.

Floating-point algorithms preserve documented shapes, ordering, failure modes,
and numerical invariants—not bit identity across legal implementations.
Backends may reassociate arithmetic, and partial iterative results must be
interpreted through their convergence diagnostics.

## Documentation

- [Getting started](docs/user/getting-started.md) — install from source and
  complete one dense workflow.
- [Core concepts](docs/user/core-concepts.md) — values, failures, operators,
  diagnostics, and backend behavior.
- [Worked-example map](docs/user/guides/examples.md) — choose a dense, sparse,
  operator, spectral, optimization, or migration task.
- [Advanced topics](docs/user/advanced/index.md) — numerical guarantees,
  allocation control, ownership, sparse structure reuse, and matrix-free
  generalized eigensolving.
- [Modules and platforms](docs/user/reference/modules-and-platforms.md) —
  artifacts, runtimes, and current publication status.
- [Troubleshooting](docs/user/troubleshooting.md) — diagnose shape, rank,
  convergence, backend, and runtime failures.

Build JVM and Scala.js Scaladoc and execute/render the complete mdoc/Laika guide:

```sh
sbt docsCheck
```

The generated site is a local and CI artifact. GitHub Pages publication and a
stable hosted Scaladoc URL are not configured yet.

## Development

```sh
sbt compileAll
sbt testAllFull
sbt parityTest interopBreezeTest
sbt benchCompile
sbt docsCheck
```

JDK 22 native gates are separate:

```sh
sbt nativeBackendTest blasFfmBackendTest benchFfmCompile
```

The CI workflow also exercises Vector API on JDK 21 and 22, FFM/OpenBLAS on JDK
22, and an allow-failure experimental WebAssembly lane. See the
[compatibility and artifact policy](docs/release-policy.md) for the intended v1
boundary and the [benchmark dashboard](benchmarks/dashboard.md) for qualified
backend evidence.

## License

[Apache-2.0](LICENSE)
