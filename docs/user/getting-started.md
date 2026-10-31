# Get your first result

This page takes Gale from an sbt dependency to a first program: make a
matrix, multiply it by a vector, and solve a linear system.

## Add the released core

Gale `0.1.0-M1` is available from Maven Central. Use Scala 3.7.4 and JDK 21
for the baseline build. In a JVM project's `build.sbt`:

```scala
scalaVersion := "3.7.4"
libraryDependencies += "io.github.canardlapin" %% "gale-core" % "0.1.0-M1"
```

For Scala.js or a cross-project, use `%%%` in place of `%%` with the matching
build plugin enabled. No Vector API or native-library setup is needed for the
portable core artifact. Continue to [the first program](#make-a-matrix-and-use-it).

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

## Make a matrix and use it

Most dense work starts with one import:

```scala mdoc:silent
import gale.linalg.*
```

`Matrix(rows, cols)` takes its values in row-major order, so the code reads
like the matrix on paper. Indices start at zero.

```scala mdoc
val a = Matrix(3, 3)(
  2.0, 1.0, 0.0,
  1.0, 3.0, 1.0,
  0.0, 1.0, 4.0
)
val x = Vec(1.0, 2.0, 3.0)

a(1, 2)
a.row(1).toSeq
(a * x).toSeq
(a.t * a + Matrix.eye(3)).valuesRowMajor
```

`*` is the matrix product, `a.t` is the transpose, and `+` adds entry by
entry. `DVec` and `DMat` have no readable `toString`, so these examples print
them with `toSeq` and `valuesRowMajor`.

## Solve a system

`a.solve(b)` finds the `x` for which `a * x` equals `b`:

```scala mdoc
val b = a * x
val solved = a.solve(b).orThrow
solved.toSeq
assert(solved.toSeq.zip(x.toSeq).forall((s, e) => math.abs(s - e) < 1e-12))
```

The solution recovers `x`.

## Keep numerical failure explicit

`solve` returns `Either[LinAlgError, DVec]`, not a bare vector. The call
above used `.orThrow`, which suits a script or a test. Library code can keep
the error and decide what to do with it:

```scala mdoc
val singular = Matrix(2, 2)(
  1.0, 2.0,
  2.0, 4.0
)

singular.solve(Vec(1.0, 1.0)).isLeft
```

A singular matrix, a non-positive Cholesky pivot, a rank-deficient
least-squares design, and a mismatched right-hand side each have their own
`LinAlgError` case.

Arithmetic is different. `a + b`, `a * x`, and `x.dot(y)` return plain
values. They check shapes and throw a `LinAlgError` if the shapes do not
match, since a mismatch there is a bug in the calling code. Solves,
factorizations, and spectral methods return `Either` because a singular or
indefinite matrix is an ordinary outcome that correct code must handle.

## Choose the next page

- [Vectors and matrices](guides/vectors-and-matrices.md) covers
  construction, indexing, views, and arithmetic in more detail.
- [Core concepts](core-concepts.md) explains values, errors, operators,
  diagnostics, and backends.
- [Dense systems and least squares](guides/dense-systems.md) fits a
  regression, inspects a pivoted QR factor, and reuses factors.
- [Find a worked example](guides/examples.md) routes by problem.
- [Troubleshooting](troubleshooting.md) helps when a shape, rank,
  convergence, or runtime problem blocks the first result.
