package userland

import gale.linalg.*
import gale.syntax.all.*

class MarginTypeSuite extends munit.FunSuite:
  test("reductions, tagging and margin transforms retain their singleton axis") {
    val a = Matrix(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    val rows: Margin.ByRow = a.byRow.mean
    val cols: Margin.ByCol = a.byCol.mean
    val tagged: Margin.ByCol = Vec(1.0, 2.0, 3.0).along(Axis.Cols)
    val transformed: Margin.ByCol = (cols + tagged).map(_ * 2.0) / 2.0
    assertEquals(rows.length, 2)
    assertEquals(transformed.length, 3)
    val oldReduction: DVec = a.mean(Axis.Cols)
    val oldProduct: DVec = a * oldReduction
    assertEquals(oldProduct.length, 2)
    val scaled: DMat = a.pointwise * transformed
    assertEquals(scaled.shape, a.shape)
    for axis <- Seq(Axis.Rows, Axis.Cols) do
      val mean = a.byAxis(axis).mean
      val ones = Vec.tabulate(mean.length)(_ => 1.0).along(axis)
      assertEquals(a.standardize(mean, ones).shape, a.shape)
  }

  test("plain vectors never broadcast and margins never silently select matvec") {
    assert(compileErrors("Matrix.zeros(2, 2) + Vec(1.0, 2.0)").nonEmpty)
    assert(compileErrors("Matrix.zeros(2, 2).pointwise * Vec(1.0, 2.0)").nonEmpty)
    assert(compileErrors("Matrix.zeros(2, 2) * Vec(1.0, 2.0).along(Axis.Cols)").nonEmpty)
    assert(compileErrors("val v: DVec = Vec(1.0).along(Axis.Rows)").nonEmpty)
    assert(compileErrors("val v: Margin[Axis] = Vec(1.0).along(Axis.Rows)").nonEmpty)
  }

  test("opposite axes cannot mix even when vector lengths happen to match") {
    assert(compileErrors("Vec(1.0).along(Axis.Rows) + Vec(1.0).along(Axis.Cols)").nonEmpty)
    assert(
      compileErrors("Matrix.zeros(1, 1).standardize(Vec(1.0).along(Axis.Rows), Vec(1.0).along(Axis.Cols))").nonEmpty
    )
    assert(compileErrors("Matrix.zeros(1, 1).affine(Vec(1.0).along(Axis.Rows), Vec(1.0).along(Axis.Cols))").nonEmpty)
    assert(compileErrors("new Margin(Axis.Rows, Vec(1.0))").nonEmpty)
  }

  test("stored parameters retain a shared axis; independently erased margins deliberately cannot mix") {
    final case class Parameters[A <: Axis](center: Margin[A], scale: Margin[A]):
      def applyTo(a: DMat): DMat = a.standardize(center, scale)

    val a = Matrix(2, 2)(1.0, 3.0, 5.0, 11.0)
    val stored: Parameters[?] = Parameters(a.byCol.mean, Vec(2.0, 4.0).along(Axis.Cols))
    assertEquals(stored.applyTo(a).row(0).toSeq, Seq(-1.0, -1.0))
    assertEquals(stored.applyTo(a).row(1).toSeq, Seq(1.0, 1.0))

    assert(compileErrors("def f(a: DMat, x: Margin[?], y: Margin[?]) = a.affine(x, y)").nonEmpty)
    assert(compileErrors("def f(a: DMat, x: Margin[?], y: Margin[?]) = a.standardize(x, y)").nonEmpty)
    assert(
      compileErrors("def f(a: DMat, x: Margin[?], y: Margin[?], out: DMatBuilder) = a.affineInto(x, y, out)").nonEmpty
    )
    assert(
      compileErrors(
        "def f(a: DMat, x: Margin[?], y: Margin[?], out: DMatBuilder) = a.standardizeInto(x, y, out)"
      ).nonEmpty
    )
  }
