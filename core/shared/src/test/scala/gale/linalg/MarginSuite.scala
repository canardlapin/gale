package gale.linalg

import gale.platform.DoubleArray.*
import gale.syntax.all.*

class MarginSuite extends munit.FunSuite:
  private def same(actual: Double, expected: Double): Unit =
    if expected.isNaN then assert(actual.isNaN)
    else assertEquals(java.lang.Double.doubleToLongBits(actual), java.lang.Double.doubleToLongBits(expected))

  private def matrixSame(actual: DMat, expected: DMat): Unit =
    assertEquals(actual.shape, expected.shape)
    for i <- 0 until actual.rows; j <- 0 until actual.cols do same(actual(i, j), expected(i, j))

  private def layouts(a: DMat): Seq[DMat] =
    val padded = Matrix.tabulate(a.rows + 2, a.cols + 3)((i, j) =>
      if i > 0 && i <= a.rows && j > 0 && j <= a.cols then a(i - 1, j - 1) else -999.0
    )
    Seq(a, Matrix.tabulate(a.cols, a.rows)((i, j) => a(j, i)).t, padded.slice(1, a.rows + 1, 1, a.cols + 1))

  private def strided(v: DVec): DVec = Matrix.tabulate(v.length, 3)((i, j) => if j == 1 then v(i) else -777.0).col(1)

  test("square reductions retain distinct row and column meanings") {
    val a = Matrix(2, 2)(1.0, 3.0, 5.0, 11.0)
    val rows = a.byRow.mean
    val cols = a.byCol.mean
    assertEquals(rows.toDVec.toSeq, Seq(2.0, 8.0))
    assertEquals(cols.toDVec.toSeq, Seq(3.0, 7.0))
    matrixSame(a - rows, Matrix(2, 2)(-1.0, 1.0, -3.0, 3.0))
    matrixSame(a - cols, Matrix(2, 2)(-2.0, -4.0, 2.0, 4.0))
    matrixSame(a.pointwise * cols, Matrix(2, 2)(3.0, 21.0, 15.0, 77.0))
    assertEquals((a * cols.toDVec).toSeq, Seq(24.0, 92.0))
  }

  test("tagging is zero-copy; margin arithmetic retains orientation and owns results") {
    val v = strided(Vec(2.0, 4.0, 8.0)).slice(1, 3)
    val m = v.along(Axis.Cols)
    assert(m.toDVec eq v)
    val sum = (m + m).map(_ + 1.0) / 2.0
    assertEquals(sum.axis, Axis.Cols)
    assertEquals(sum.toDVec.toSeq, Seq(4.5, 8.5))
    assertEquals((m - 1.0).toDVec.toSeq, Seq(3.0, 7.0))
    assertEquals((m * 2.0).toDVec.toSeq, Seq(8.0, 16.0))
    assertEquals((m - m).toDVec.toSeq, Seq(0.0, 0.0))
    assert(!(sum.toDVec eq v))
    assertEquals(v.toSeq, Seq(4.0, 8.0))
  }

  test("tagged reductions preserve existing algorithms, including empty-line and nonfinite policy") {
    for
      a <- Seq(Matrix(2, 3)(1.0, 2.0, 3.0, 9.0, 5.0, 7.0), Matrix.zeros(0, 3), Matrix.zeros(3, 0));
      axis <- Seq(Axis.Rows, Axis.Cols)
    do
      val r = a.byAxis(axis)
      assertEquals(r.sum.toDVec.toSeq, a.sum(axis).toSeq)
      val lines = if axis == Axis.Rows then a.rows else a.cols
      val count = if axis == Axis.Rows then a.cols else a.rows
      if lines > 0 && count == 0 then
        intercept[LinAlgError.EmptyInput](r.mean)
        intercept[LinAlgError.EmptyInput](r.min)
        intercept[LinAlgError.EmptyInput](r.max)
      else
        assertEquals(r.mean.toDVec.toSeq, a.mean(axis).toSeq)
        assertEquals(r.min.toDVec.toSeq, a.min(axis).toSeq)
        assertEquals(r.max.toDVec.toSeq, a.max(axis).toSeq)
    val a = Matrix(2, 2)(1e308, 1e308, Double.NaN, Double.PositiveInfinity)
    same(a.byRow.mean(0), 1e308)
    assert(a.byRow.mean(1).isNaN)
    assert(a.byRow.min(1).isNaN)
  }

  test("all binary operations agree with explicit indexing across axes and storage layouts") {
    val rng = new scala.util.Random(72L)
    for rows <- Seq(0, 1, 2, 5); cols <- Seq(0, 1, 3, 5) do
      val base = Matrix.tabulate(rows, cols)((_, _) => rng.nextDouble() * 20.0 - 10.0)
      for a <- layouts(base); axis <- Seq(Axis.Rows, Axis.Cols) do
        val n = if axis == Axis.Rows then rows else cols
        val v = strided(Vec.tabulate(n)(i => i + 0.75))
        val margin = v.along(axis)
        val operations: Seq[(Double, Double) => Double] = Seq(_ + _, _ - _, _ * _, _ / _)
        val results = Seq(a + margin, a - margin, a.pointwise * margin, a.pointwise / margin)
        for (result, op) <- results.zip(operations) do
          matrixSame(result, Matrix.tabulate(rows, cols)((i, j) => op(a(i, j), v(if axis == Axis.Rows then i else j))))
        val out = DMatBuilder.zeros(rows, cols)
        a.addInto(margin, out)
        matrixSame(out.result(), results(0))
        val sub = DMatBuilder.zeros(rows, cols)
        a.subInto(margin, sub)
        matrixSame(sub.result(), results(1))
        val mul = DMatBuilder.zeros(rows, cols)
        a.pointwise.mulInto(margin, mul)
        matrixSame(mul.result(), results(2))
        val div = DMatBuilder.zeros(rows, cols)
        a.pointwise.divInto(margin, div)
        matrixSame(div.result(), results(3))
  }

  test("fused operations preserve explicit primitive order on strided inputs and both axes") {
    for
      base <- Seq(Matrix.tabulate(3, 4)((i, j) => i * 7.0 - j * 1.25), Matrix.zeros(0, 3), Matrix.zeros(3, 0));
      a <- layouts(base); axis <- Seq(Axis.Rows, Axis.Cols)
    do
      val n = if axis == Axis.Rows then a.rows else a.cols
      val center = strided(Vec.tabulate(n)(i => i + 0.5)).along(axis)
      val scale = strided(Vec.tabulate(n)(i => i + 2.0)).along(axis)
      val standard = Matrix.tabulate(a.rows, a.cols)((i, j) =>
        val k = if axis == Axis.Rows then i else j
        (a(i, j) - center(k)) / scale(k)
      )
      val affine = Matrix.tabulate(a.rows, a.cols)((i, j) =>
        val k = if axis == Axis.Rows then i else j
        (a(i, j) * scale(k)) + center(k)
      )
      matrixSame(a.standardize(center, scale), standard)
      matrixSame(a.affine(scale, center), affine)
      val out = DMatBuilder.zeros(a.rows, a.cols)
      out.fill(Double.NaN)
      a.standardizeInto(center, scale, out)
      matrixSame(out.result(), standard)
      val out2 = DMatBuilder.zeros(a.rows, a.cols)
      out2.fill(Double.NaN)
      a.affineInto(scale, center, out2)
      matrixSame(out2.result(), affine)
  }

  test("standardization is not reciprocal-affine arithmetic and never contracts an FMA") {
    val a = Matrix(1, 1)(10000000000000002.0)
    val center = Vec(10000000000000000.0).along(Axis.Cols)
    val scale = Vec(3.0).along(Axis.Cols)
    same(a.standardize(center, scale)(0, 0), 2.0 / 3.0)
    same(a.affine(scale.map(1.0 / _), center.zipMap(scale)((c, s) => -c / s))(0, 0), 0.5)
    val delta = math.pow(2.0, -27)
    val unfused = Matrix(1, 1)(1.0 + delta).affine(Vec(1.0 - delta).along(Axis.Rows), Vec(-1.0).along(Axis.Rows))
    same(unfused(0, 0), 0.0) // An FMA would retain -2^-54.
    same(Matrix(1, 1)(1e308).standardize(Vec(1e308).along(Axis.Rows), Vec(1e-308).along(Axis.Rows))(0, 0), 0.0)
  }

  test("IEEE values, signed zero, scalar operations and zero coefficients are literal arithmetic") {
    val values =
      Seq(0.0, -0.0, 1.0, -1.0, Double.MinPositiveValue, Double.PositiveInfinity, Double.NegativeInfinity, Double.NaN)
    val base = Matrix.tabulate(2, values.length)((_, j) => values(j))
    for a <- layouts(base) ++ Seq(Matrix.zeros(0, 3), Matrix.zeros(3, 0)); s <- values do
      val m = Vec.tabulate(a.rows)(_ => s).along(Axis.Rows)
      for (out, f) <- Seq[(DMat, Double => Double)](
          (a.pointwise + s, _ + s),
          (a.pointwise - s, _ - s),
          (a.pointwise * s, _ * s),
          (a.pointwise / s, _ / s),
          (a + m, _ + s),
          (a - m, _ - s),
          (a.pointwise * m, _ * s),
          (a.pointwise / m, _ / s)
        )
      do matrixSame(out, Matrix.tabulate(a.rows, a.cols)((i, j) => f(a(i, j))))
      matrixSame(a.pointwise * s, a * s)
    val zero = Vec(0.0, 0.0).along(Axis.Rows)
    assert(base.affine(zero, zero)(0, 5).isNaN)
    // Reciprocal scaling overflows even though literal division has a finite result.
    same((Matrix(1, 1)(Double.MinPositiveValue).pointwise / Double.MinPositiveValue)(0, 0), 1.0)
  }

  test("shape, axis, closed-builder and alias errors precede all destination writes") {
    val a = Matrix(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val good = Vec(1.0, 2.0, 3.0).along(Axis.Cols)
    val bad = Vec(1.0, 2.0).along(Axis.Cols)
    val out = DMatBuilder.zeros(2, 3)
    out.fill(123.0)
    intercept[LinAlgError.VectorLengthMismatch](a.subInto(bad, out))
    intercept[LinAlgError.VectorLengthMismatch](a.standardizeInto(good, bad, out))
    intercept[LinAlgError.VectorLengthMismatch](a.pointwise.mulInto(bad, out))
    for i <- 0 until 2; j <- 0 until 3 do same(out(i, j), 123.0)
    val wrongShape = DMatBuilder.zeros(3, 2)
    wrongShape.fill(456.0)
    intercept[LinAlgError.DimensionMismatch](a.addInto(good, wrongShape))
    same(wrongShape(0, 0), 456.0)
    val dynamicRows: Margin[Axis] = Margin.tag(Vec(1.0, 2.0), Axis.Rows: Axis)
    val dynamicCols: Margin[Axis] = Margin.tag(Vec(1.0, 2.0, 3.0), Axis.Cols: Axis)
    intercept[LinAlgError.InvalidArgument](dynamicRows.zipMap(dynamicCols)(_ + _))
    intercept[LinAlgError.InvalidArgument](a.affineInto(dynamicRows, dynamicCols, out))
    same(out(0, 0), 123.0)
    val aliased = new DMat(out.writableData, Offset(0), Rows(2), Cols(3), Stride(3), Stride(1))
    intercept[LinAlgError.UnsupportedOperation](aliased.subInto(good, out))
    val aliasVec = new DVec(out.writableData, Offset(0), Length(3), Stride(1)).along(Axis.Cols)
    intercept[LinAlgError.UnsupportedOperation](a.addInto(aliasVec, out))
    intercept[LinAlgError.UnsupportedOperation](a.affineInto(good, aliasVec, out))
    same(out(0, 0), 123.0)
    out.result()
    intercept[LinAlgError.UnsupportedOperation](a.standardizeInto(good, good, out))
  }

  test("outputs and repeated destination replacements retain independent immutable ownership") {
    val a = Matrix(2, 2)(1.0, 2.0, 3.0, 4.0)
    val margin = Vec(10.0, 20.0).along(Axis.Cols)
    val old = a + margin
    val out = a.toBuilder
    a.subInto(margin, out)
    a.addInto(margin, out)
    matrixSame(out.result(), old)
    matrixSame(a, Matrix(2, 2)(1.0, 2.0, 3.0, 4.0))
    matrixSame(old, Matrix(2, 2)(11.0, 22.0, 13.0, 24.0))
  }
