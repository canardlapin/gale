# Vectors and matrices

This page covers the operations most programs use every day: making vectors
and matrices, reading their entries, taking rows and columns, doing
arithmetic, and solving a square system. Everything here comes from one
import:

```scala mdoc:silent
import gale.linalg.*
```

## Make vectors and matrices

`Vec` makes a dense vector of `Double`s. `Matrix(rows, cols)` takes its values
in row-major order, so you write the matrix as you would on paper:

```scala mdoc
val x = Vec(1.0, 2.0, 3.0)

val a = Matrix(3, 3)(
  2.0, 1.0, 0.0,
  1.0, 3.0, 1.0,
  0.0, 1.0, 4.0
)

(x.length, a.rows, a.cols)
```

When each entry follows from its position, use `tabulate`. `zeros`, `fill`,
and `Matrix.eye` cover the other common cases:

```scala mdoc
val hilbert = Matrix.tabulate(3, 3)((i, j) => 1.0 / (i + j + 1))
val ones = Vec.fill(3)(1.0)
val identity = Matrix.eye(3)

hilbert.valuesRowMajor
```

The values are `DVec` and `DMat`. Neither has a readable `toString`, so the
examples on these pages print them with `toSeq` and `valuesRowMajor`.

## Read entries, rows, and columns

Indices start at zero. `a(i, j)` reads one entry and `x(i)` one element.
`row`, `col`, and `slice` select part of a matrix; slice bounds include the
first index and exclude the last.

```scala mdoc
a(1, 2)
x(0)
a.row(1).toSeq
a.col(2).toSeq
a.slice(0, 2, 1, 3).valuesRowMajor
```

`row`, `col`, `slice`, and the transpose `a.t` are views: they share the
original's storage instead of copying it. Sharing is safe because Gale's
vectors and matrices never change once made. To change an entry, ask for a new
value:

```scala mdoc
val b = a.updated(0, 0, 10.0)
(a(0, 0), b(0, 0))
```

`updated` copies the whole matrix, which is fine for a few changes. To fill
many entries, use `tabulate`, or a builder as described in
[Ownership](../advanced/ownership.md).

## Do arithmetic

`+` and `-` work entry by entry. `*` with a number scales; `*` with a vector
or matrix is the matrix product. `dot` and `norm2` give the inner product and
Euclidean length of vectors.

```scala mdoc
(a * x).toSeq
(a.t * a).valuesRowMajor
(a + identity).valuesRowMajor
(x * 2.0 - ones).toSeq
(x.dot(ones), x.norm2)
```

Each operation checks shapes first. Multiplying a 3 × 3 matrix by a vector of
length 2 throws a `LinAlgError` that names the mismatch, at the line that made
the mistake:

```scala mdoc
val mismatch =
  try
    a * Vec(1.0, 2.0)
    "no error"
  catch case e: LinAlgError => e.getClass.getSimpleName
```

## Work entry by entry

`*` between two matrices is the matrix product. For the entry-by-entry
product, or to apply a function to every entry, import `gale.syntax.all.*`
and use `pointwise`:

```scala mdoc
import gale.syntax.all.*

(a.pointwise * a).valuesRowMajor
a.pointwise.map(v => v * v + 1.0).valuesRowMajor
Numerics.exp(x).toSeq
```

`Numerics` also has `log`, `log1p`, `sigmoid`, and a stable `logSumExp`.

## Sum, average, and find extremes

A reduction without an axis covers every entry. `Axis.Rows` gives one value
per row and `Axis.Cols` one value per column.

```scala mdoc
val scores = Matrix(2, 3)(
  1.0, 4.0, 2.0,
  3.0, 0.0, 5.0
)

(scores.sum, scores.max, scores.argmax)
scores.sum(Axis.Rows).toSeq
scores.mean(Axis.Cols).toSeq
(x.sum, x.mean, x.argmax)
```

`argmax` on a matrix returns a `(row, column)` pair. To subtract each column's
mean from that column, build the matrix of means and subtract it:

```scala mdoc
val means = scores.mean(Axis.Cols)
val centered = scores - Matrix.tabulate(scores.rows, scores.cols)((_, j) => means(j))
centered.valuesRowMajor
```

[Reductions and norms](reductions-and-numerics.md) covers norms, ties and
NaN, exact sums, and softmax.

## Solve a square system

`a.solve(rhs)` solves *A x = b*. It returns an `Either`, because a singular
matrix has no unique solution and the caller should decide what that means:

```scala mdoc
val rhs = a * x
a.solve(rhs).map(_.toSeq)
```

```scala mdoc
val singular = Matrix(2, 2)(
  1.0, 2.0,
  2.0, 4.0
)
singular.solve(Vec(1.0, 1.0))
```

The second matrix has a zero determinant, so the result is a `Left` naming
the problem rather than a vector of meaningless numbers. Use
`.orThrow` when a failure should stop the program, as in a test or a script.

## Next steps

- [Dense systems and least squares](dense-systems.md) covers many right-hand
  sides, rectangular systems, and keeping a factorization to reuse.
- [Reductions and norms](reductions-and-numerics.md) covers sums, means,
  extrema, norms, and stable softmax.
- [Sparse matrices and operators](sparse-operators.md) covers matrices that
  are mostly zeros.
- [Core concepts](../core-concepts.md) explains views, failures, and backends
  in more depth.
