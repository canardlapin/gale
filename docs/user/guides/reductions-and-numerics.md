# Reductions, norms, and elementwise numerics

Gale vectors and matrices carry the everyday reductions as members: sums,
means, extrema and their positions, per-axis reductions, and norms. Named
elementwise functions and the log-domain operations (`logSumExp`, `softmax`,
`logSoftmax`) live in `gale.linalg.Numerics`. Everything here is portable
`gale-core` code and runs on the JVM and Scala.js.

## Reduce a vector

```scala mdoc
import gale.linalg.*

val x = Vec(3.0, -1.0, 4.0, -1.0, 5.0, -9.0)

assert(x.sum == 1.0)
assert(x.mean == 1.0 / 6.0)
assert(x.max == 5.0 && x.argmax == 4)
assert(x.min == -9.0 && x.argmin == 5)
assert(x.norm1 == 23.0)
assert(x.normInf == 9.0)
x.norm2
```

Views reduce in place: `x.slice(1, 4).sum`, a matrix row, or a column of a
row-major matrix (a strided view) needs no copy.

`argmax` and `argmin` return the **first** position of a tie. If any entry is
NaN, `max` and `min` are NaN and `argmax`/`argmin` point at the first NaN, so a
NaN cannot hide behind a larger value.

```scala mdoc
val tied = Vec(2.0, 7.0, 7.0, 1.0)
assert(tied.argmax == 1)

val withNaN = Vec(1.0, Double.NaN, 9.0)
assert(withNaN.argmax == 1 && withNaN.max.isNaN)
```

## Choose `sum` or `sumExact`

`sum` is the fast path. It keeps several partial sums, so its last bits may
differ from a left-to-right loop and between JVM and Scala.js. `sumExact`
computes the exact total and rounds it once, so it is the same on every
platform and for every order of the entries. Use it for cancellation-prone
data or when a result must be reproducible bit for bit.

```scala mdoc
val cancelling = Vec(1e100, 1.0, -1e100, 3.0)
cancelling.sum
assert(cancelling.sumExact == 4.0)
```

## Empty inputs

A sum of nothing is `0.0`, and so is every norm of an empty vector or matrix.
There is no meaningful mean, maximum, or position of an empty input, so those
throw the typed `LinAlgError.EmptyInput`:

```scala mdoc
val empty = Vec()
assert(empty.sum == 0.0 && empty.normInf == 0.0)

val failure =
  try
    empty.max
    None
  catch case error: LinAlgError.EmptyInput => Some(error.operation)
assert(failure.contains("DVec.max"))
```

## Reduce a matrix, whole or per axis

The whole-matrix forms reduce every entry; `argmax` returns `(row, col)` of the
first maximum in row-major order, whatever the storage layout.

Per-axis forms take an `Axis`, named after **what the result is indexed by**:
`Axis.Rows` gives one value per row, and `Axis.Cols` gives one value per column.
In Breeze terms, `a.sum(Axis.Rows)` is `sum(a(*, ::))` and `a.sum(Axis.Cols)`
is `sum(a(::, *))`. In NumPy terms they are `axis=1` and `axis=0`.

```scala mdoc
val a = Matrix(2, 3)(
  1.0, 2.0, 3.0,
  4.0, 5.0, 9.0
)

assert(a.sum == 24.0)
assert(a.argmax == (1, 2))
assert(a.sum(Axis.Rows).toSeq == Seq(6.0, 18.0))
assert(a.sum(Axis.Cols).toSeq == Seq(5.0, 7.0, 12.0))
assert(a.mean(Axis.Cols).toSeq == Seq(2.5, 3.5, 6.0))
assert(a.max(Axis.Rows).toSeq == Seq(3.0, 9.0))
```

## Matrix norms

`norm1` is the largest absolute column sum, `normInf` the largest absolute row
sum, and `normFrobenius` the square root of the sum of squares. The Frobenius
norm scales internally, so it does not overflow for entries near `1e300`.

```scala mdoc
val b = Matrix(2, 2)(
  1.0, -2.0,
  -3.0, 4.0
)
assert(b.norm1 == 6.0)
assert(b.normInf == 7.0)
assert(math.abs(b.normFrobenius - math.sqrt(30.0)) < 1e-15)

val huge = Matrix.tabulate(2, 2)((_, _) => 1e300)
assert(math.abs(huge.normFrobenius - 2e300) < 1e286)
```

## Elementwise functions

Elementwise operations are explicit in Gale. `a.pointwise.map(f)` (from
`gale.syntax.all.*`) applies any function. `Numerics.exp`, `log`, `log1p`,
`expm1`, and `sigmoid` are the fast paths for common ones: they apply the same
scalar function, so `Numerics.exp(a)` equals `a.pointwise.map(math.exp)` value
for value, but runs as a primitive loop. Results are new, owned row-major
matrices or contiguous vectors.

```scala mdoc
import gale.syntax.all.*

val z = Matrix(2, 2)(
  0.0, 1.0,
  -1.0, 2.0
)
assert(Numerics.exp(z).valuesRowMajor == z.pointwise.map(math.exp).valuesRowMajor)

val probabilities = Numerics.sigmoid(Vec(-1000.0, 0.0, 1000.0))
assert(probabilities.toSeq == Seq(0.0, 0.5, 1.0))
```

`sigmoid` is evaluated in an overflow-free form, so very large or very
negative inputs saturate cleanly to `1.0` and `0.0`.

## Log-sum-exp and softmax

`logSumExp`, `softmax`, and `logSoftmax` subtract the maximum first, so inputs
of magnitude 1000 neither overflow nor underflow to a useless result.

```scala mdoc
val logits = Vec(1000.0, 1000.0, 999.0)
val lse = Numerics.logSumExp(logits)
assert(math.abs(lse - (1000.0 + math.log(2.0 + math.exp(-1.0)))) < 1e-12)

val p = Numerics.softmax(logits)
assert(math.abs(p.sumExact - 1.0) < 1e-15)
val logP = Numerics.logSoftmax(logits)
assert(math.abs(logP(0) - (1000.0 - lse)) < 1e-12)
```

For matrices, the forms without an axis treat every entry as one input; the
`Axis` overloads give one distribution per row or per column:

```scala mdoc
val scores = Matrix(2, 3)(
  1.0, 2.0, 3.0,
  10.0, 10.0, 10.0
)
val perRow = Numerics.softmax(scores, Axis.Rows)
assert(math.abs(perRow.row(0).sum - 1.0) < 1e-15)
assert(perRow.row(1).toSeq.forall(v => math.abs(v - 1.0 / 3.0) < 1e-15))
Numerics.logSumExp(scores, Axis.Rows).toSeq
```

An empty `logSumExp` is `-Inf` (the log of an empty sum). If a line contains
NaN, contains `+Inf`, or is entirely `-Inf`, `softmax` and `logSoftmax` return
NaN for that whole line.

**Coming from Breeze:** Breeze's `softmax(v)` returns the scalar log-sum-exp.
Its Gale equivalent is `Numerics.logSumExp(v)`; Gale's `Numerics.softmax(v)` is
the normalized-exponential vector.

The complete NaN, empty-input, and determinism rules are in the
[numerical contract](../advanced/numerical-contract.md).
