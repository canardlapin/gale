# Axis margins (planned for 0.2)

**Unreleased API preview.** This guide requires the 0.2 development sources.
These methods are not available in Gale 0.1.0-M1 from
Maven Central; the published dependency in the getting-started guide cannot run
these examples.

A margin holds one value per row or one value per column. Its type retains that
choice, including when a square matrix makes both lengths identical. Gale never
infers a broadcast direction from vector length.

## Centering and row/column weights

```scala mdoc
import gale.linalg.*
import gale.syntax.all.*

val x = Matrix(2, 3)(1.0, 2.0, 3.0, 5.0, 8.0, 11.0)
val means: Margin.ByCol = x.byCol.mean
val centered = x - means
assert(centered.row(0).toSeq == Seq(-2.0, -3.0, -4.0))
assert(centered.row(1).toSeq == Seq(2.0, 3.0, 4.0))

val rowWeights = Vec(2.0, 3.0).along(Axis.Rows)
val weighted = x.pointwise * rowWeights
assert(weighted.row(1).toSeq == Seq(15.0, 24.0, 33.0))
```

`byRow` means one result per row; `byCol` means one per column. Both offer
`sum`, `mean`, `min`, and `max` with the existing reductions' empty/nonfinite
behavior. Existing calls such as `x.mean(Axis.Cols)` still return plain `DVec`.
Use `byAxis(axis)` if the axis is selected at runtime.
Use `Margin.ByRow` and `Margin.ByCol` for named row and column margin types.

Tagging `v.along(axis)` and unwrapping `margin.toDVec` do not copy elements.
Margin arithmetic and `map` retain the axis and create independent vector results.
There is one small wrapper per margin. A plain vector cannot broadcast, and a
margin cannot implicitly participate in matrix-vector multiplication:

```scala mdoc:fail
x * means
```

To request matrix-vector multiplication, write `x * means.toDVec`. For pointwise
multiplication/division, write `x.pointwise * means` or `x.pointwise / means`.
`x.pointwise + 2.0` and the corresponding scalar `-`, `*`, `/` apply the written
operation to every element.

## Applying supplied centers and scales

This preview has no `std` or variance reduction. Centering above is supported
directly; estimating z-scores from a dataset awaits an explicit estimator contract
(including the divisor and numerical behavior).

If centers and scales are already available, for example from instrument
calibration or a model fitted elsewhere, `standardize` applies
`(x - center) / scale` in one traversal. The numbers below are supplied calibration
parameters, not standard deviations estimated from `x`:

```scala mdoc
val calibrationCenters: Margin.ByCol = Vec(1.0, 2.0, 3.0).along(Axis.Cols)
val calibrationScales: Margin.ByCol = Vec(2.0, 3.0, 4.0).along(Axis.Cols)
val calibrated = x.standardize(calibrationCenters, calibrationScales)
assert(calibrated.row(0).toSeq == Seq(0.0, 0.0, 0.0))
assert(calibrated.row(1).toSeq == Seq(2.0, 2.0, 2.0))

val out = DMatBuilder.zeros(x.rows, x.cols)
x.standardizeInto(calibrationCenters, calibrationScales, out)
val result = out.result() // transfers storage and permanently closes out
assert(result.row(0).toSeq == calibrated.row(0).toSeq)
```

`affine(scale, shift)` and `affineInto` instead evaluate `(x * scale) + shift`.
Neither operation uses FMA or changes the arithmetic order. In particular,
standardization is not rewritten as a reciprocal-based affine transformation:
that rewrite can change rounding and overflow. Zero, negative, and nonfinite
scales follow ordinary IEEE arithmetic, including division by zero. Validate
scales yourself if your application requires positive finite values.

Both operands must have the same axis. This is rejected even on square matrices:

```scala mdoc:fail
Matrix.zeros(2, 2).standardize(
  Vec(1.0, 2.0).along(Axis.Rows),
  Vec(1.0, 2.0).along(Axis.Cols)
)
```

## Storing parameters without losing their axis

Fused operations deliberately require a shared axis type. Two independent fields
typed `Margin[?]` erase the evidence that their axes agree; a runtime axis check
alone does not make such a call compile:

```scala mdoc:fail
def erasedTransform(a: DMat, center: Margin[?], scale: Margin[?]): DMat =
  a.standardize(center, scale)
```

For fixed-axis fields, use `Margin.ByRow` or `Margin.ByCol`. For a reusable model,
keep both margins under one type parameter; the whole model can then be stored
with an unknown axis:

```scala mdoc
final case class Calibration[A <: Axis](center: Margin[A], scale: Margin[A]):
  def applyTo(a: DMat): DMat = a.standardize(center, scale)

val saved: Calibration[?] = Calibration(calibrationCenters, calibrationScales)
assert(saved.applyTo(x).row(1).toSeq == calibrated.row(1).toSeq)
```

When selecting an axis at runtime, use the same stable `axis` value for
`x.byAxis(axis)` and `v.along(axis)` so both retain the shared `axis.type`.
The same restriction applies to `affine` and both fused destination forms.
There is no wildcard overload: it would also let statically opposite axes compile.
Single-margin matrix arithmetic accepts `Margin[?]`, because it has no second
axis to compare.

## Destinations and limits

Use `x.addInto(margin, out)`, `x.subInto(margin, out)`,
`x.pointwise.mulInto(margin, out)`, or `x.pointwise.divInto(margin, out)` to reuse an
open builder before transferring its storage. All shape and alias checks happen
before writing. Destination contents are replaced, including NaN-filled storage.
No broadcast operation exposes a writable repeated-element view.

The output shape always equals the input matrix shape. Transposes, slices, and
strided vector views are accepted. Results are row-major, so traversal reads a
transposed input with a stride. With supplied parameters, a fused transform
allocates one result matrix buffer; subtracting and then dividing allocates two.
This is an allocation guarantee, not a general speedup claim.

A wrong length is an error, even when the
matrix has an empty dimension. Sparse matrices, tensors, singleton-matrix
expansion, and implicit outer products are outside this API. Axis types preserve
orientation; they cannot prove that you chose the correct statistic or that two
data sources use the same feature ordering.
