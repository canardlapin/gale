# Get your first result

This page takes Gale from an sbt dependency to one interpreted numerical
result. The example fits a straight line, so the expected coefficients are
visible without knowing any Gale internals.

## Add Gale before the first release

Gale is currently source-only. Pin an exact revision and choose the platform
projection your build uses:

```scala
lazy val galeRevision = "<commit sha>"
lazy val galeBuild = uri(
  s"https://github.com/canardlapin/gale.git#$galeRevision"
)
lazy val galeCoreJVM = ProjectRef(galeBuild, "coreJVM")

lazy val app = project.dependsOn(galeCoreJVM)
```

Use `coreJS` for a Scala.js project. For a cross-project, depend on the matching
`coreJVM` and `coreJS` references from its platform projections.

The intended public coordinate is
`"io.github.canardlapin" %%% "gale-core" % "<published-version>"`, but it is
not installable from Maven Central until a release is actually published. Do
not use the site's snapshot version as though it were a public artifact.

## Fit a small regression

Most dense work starts with one import:

```scala mdoc:silent
import gale.linalg.*
```

Matrix literals are row-major. Here the first column is an intercept and the
second is a predictor:

```scala mdoc
val design = Matrix(4, 2)(
  1.0, 0.0,
  1.0, 1.0,
  1.0, 2.0,
  1.0, 3.0
)
val observations = Vec(1.0, 3.0, 5.0, 7.0)

val coefficients = design.leastSquares(observations).orThrow
coefficients.toSeq
```

The result is an owned `DVec`. The design and observations remain unchanged,
and later Gale operations cannot mutate the coefficients through another
ordinary public handle.

## Keep numerical failure explicit

`leastSquares` returns `Either[LinAlgError, DVec]`. The previous call used
`.orThrow` because failure would make the tutorial invalid. Library code can
retain the error value:

```scala mdoc
val dependent = Matrix(3, 2)(
  1.0, 2.0,
  2.0, 4.0,
  3.0, 6.0
)

dependent.leastSquares(Vec(1.0, 2.0, 3.0)).isLeft
```

A rank-deficient design, singular square matrix, non-positive Cholesky pivot,
or mismatched right-hand side has a distinct `LinAlgError` case.

Primitive arithmetic such as `a + b`, `a * x`, or `x.dot(y)` is intentionally
different: it validates shape preconditions and throws `LinAlgError` when they
are violated. Total numerical entry points such as solves, factorizations, and
spectral facades use `Either` where a structural numerical failure is part of
ordinary control flow.

## Inspect a factor when you need it

The one-line `leastSquares` call factors the design internally. Retain a QR
factor when several right-hand sides share the same design or when rank and
pivot information matter:

```scala mdoc
val qr = design.qr(QROptions(pivoting = QRPivoting.Column))

(
  qr.diagnostics.rank,
  qr.columnPermutation.toIndexSeq,
  qr.solveLeastSquares(observations).orThrow.toSeq
)
```

Explicit QR options use Gale's portable deterministic pivot and rank policy.
The default no-options `qr` may use an imported backend for an eligible shape.

## Choose the next page

- Read [Core concepts](core-concepts.md) for values, errors, operators,
  diagnostics, and backends.
- Continue the same model-fitting workflow in
  [Dense systems and least squares](guides/dense-systems.md).
- Choose another task from the [worked-example map](guides/examples.md).
- Use [Troubleshooting](troubleshooting.md) when a shape, rank, convergence, or
  runtime problem blocks the first result.
