# Getting started

Gale's public API is small enough to begin with two imports:

```scala mdoc:silent
import gale.linalg.*
import gale.spectral.*
```

The examples on this page are compiled and executed as part of `docsCheck`.

## Use Gale from this checkout

Gale is not yet publicly released. For now, use it as a source dependency in a
multi-project build or run `publishLocal` from this checkout. The planned core
coordinates are:

```scala
libraryDependencies +=
  "io.github.canardlapin" %% "gale-core" % "<published-version>"
```

Use `%%%` instead of `%%` in a Scala.js or cross-project dependency. Do not
substitute the site's snapshot version into a build unless that artifact was
actually published locally.

## Build a matrix and multiply

Dense literals are row-major. Ordinary arithmetic returns owned immutable
values:

```scala mdoc
val a = Matrix(2, 2)(
  1.0, 2.0,
  3.0, 4.0
)
val x = Vec(1.0, 1.0)

(a * x).toSeq
```

Construct larger values with `tabulate` rather than assembling intermediate
collections:

```scala mdoc
val diagonal = Matrix.tabulate(4, 4): (row, col) =>
  if row == col then row.toDouble + 1.0 else 0.0

(diagonal * Vec.fill(4)(2.0)).toSeq
```

## Handle numerical failure explicitly

Operations that can fail return `Either[LinAlgError, A]`. Keep the `Either` in
library code and decide at the application boundary how to report failure:

```scala mdoc
val system = Matrix(2, 2)(
  3.0, 1.0,
  1.0, 2.0
)
val rhs = Vec(9.0, 8.0)

system.solve(rhs).map(_.toSeq)
```

For a tutorial or a test where failure should abort immediately, importing
`gale.linalg.*` also provides `.orThrow`:

```scala mdoc
val solution = system.solve(rhs).orThrow
solution.toSeq
```

## Compute a symmetric eigendecomposition

The dense facade returns eigenvalues in ascending algebraic order. Property
checks and iterative operator routes are separate APIs; this small dense example
asks for all right eigenvectors:

```scala mdoc
val symmetric = Matrix(2, 2)(
  2.0, 1.0,
  1.0, 2.0
)

val eig = Eigen
  .eigSymmetric(
    symmetric,
    EigenSelection.All,
    EigenVectors.Right
  )
  .orThrow

eig.eigenvalues.toSeq
```

Diagnostics are part of the result rather than inferred from successful return:

```scala mdoc
(
  eig.diagnostics.converged,
  eig.diagnostics.requested,
  eig.diagnostics.allConverged
)
```

## Choose the next guide

- Read [core concepts](core-concepts.md) for the value, error, operator, and
  backend model used throughout Gale.
- Use [worked examples](guides/examples.md) for dense solves, PCA, graph Laplacians,
  sparse formats, sized values, and migration examples.
- Use the
  [matrix-free generalized eigensolver guide](advanced/generalized-operator-eigen.md)
  when `A` and `B` should be operators rather than dense matrices.
- Read the [numerical contract](advanced/numerical-contract.md) before depending on
  determinism, tolerances, sparse canonicalization, or optional backends.
- Read the [Breeze migration guide](guides/breeze-equivalence.md) before treating Gale as a
  replacement for an existing Breeze workload.
