# Get your first result

This page takes Gale from an sbt dependency to one interpreted numerical
result. The example fits a straight line, so the expected coefficients are
visible without knowing any Gale internals.

## Add the released core

Gale `0.1.0-M1` is available from Maven Central. Use Scala 3.7.4 and JDK 21
for the baseline build. In a JVM project's `build.sbt`:

```scala
scalaVersion := "3.7.4"
libraryDependencies += "io.github.canardlapin" %% "gale-core" % "0.1.0-M1"
```

For Scala.js or a cross-project, use `%%%` in place of `%%` with the matching
build plugin enabled. No Vector API or native-library setup is needed for the
portable core artifact. Continue to [the regression example](#fit-a-small-regression).

## Add Gale from source

For source integration instead of an artifact dependency, the following pins
the exact commit behind `v0.1.0-M1`:

```scala
scalaVersion := "3.7.4"

lazy val galeRevision = "62e4998c03a1eae1bed55cefa2d131f7d9e33568"
lazy val galeBuild = uri(
  s"https://github.com/canardlapin/gale.git#$galeRevision"
)
lazy val galeCoreJVM = ProjectRef(galeBuild, "coreJVM")

lazy val app = project.in(file(".")).dependsOn(galeCoreJVM)
```

For an existing project, add `.dependsOn(galeCoreJVM)` to its definition.
Add `--add-modules=jdk.incubator.vector` to `.jvmopts` when consuming the
source build. Use `coreJS` for a Scala.js project. For a cross-project, depend on the matching
`coreJVM` and `coreJS` references from its platform projections.

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
assert(math.abs(coefficients(0) - 1.0) < 1e-12)
assert(math.abs(coefficients(1) - 2.0) < 1e-12)
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
