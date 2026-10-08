package gale.linalg

import gale.syntax.all.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

class NumericsSuite extends ScalaCheckSuite:
  private val Nan = Double.NaN
  private val PInf = Double.PositiveInfinity
  private val NInf = Double.NegativeInfinity

  private def assertClose(actual: Double, expected: Double, tolerance: Double = 1e-14)(using munit.Location): Unit =
    assert(
      math.abs(actual - expected) <= tolerance * math.max(1.0, math.abs(expected)),
      s"actual=$actual expected=$expected"
    )

  private def assertBits(actual: Seq[Double], expected: Seq[Double], clue: String)(using munit.Location): Unit =
    assertEquals(actual.map(java.lang.Double.doubleToLongBits), expected.map(java.lang.Double.doubleToLongBits), clue)

  private def matrixLayouts(rows: Int, cols: Int, f: (Int, Int) => Double): Seq[(String, DMat)] =
    val rowMajor = Matrix.tabulate(rows, cols)(f)
    val colMajor = Matrix.tabulate(cols, rows)((j, i) => f(i, j)).t
    val slice = Matrix.tabulate(rows + 2, cols + 3)((i, j) =>
      if i >= 1 && i <= rows && j >= 2 && j < cols + 2 then f(i - 1, j - 2) else 1e9
    ).slice(1, rows + 1, 2, cols + 2)
    Seq("row-major" -> rowMajor, "col-major" -> colMajor, "slice" -> slice)

  private def vectorLayouts(xs: Seq[Double]): Seq[(String, DVec)] =
    val strided = Matrix.tabulate(xs.length, 2)((i, j) => if j == 0 then xs(i) else 42.0).col(0)
    Seq("contiguous" -> Vec(xs*), "strided" -> strided)

  private val specials = Seq(0.0, -0.0, 1e-300, -1e-300, 0.5, -0.5, 1.0, -1.0, 30.0, -30.0, 709.0, 710.0, -745.0,
    -746.0, 1000.0, -1000.0, PInf, NInf, Nan)

  private val elementwise: Seq[(String, DVec => DVec, DMat => DMat, Double => Double)] = Seq(
    ("exp", Numerics.exp, Numerics.exp, math.exp),
    ("log", Numerics.log, Numerics.log, math.log),
    ("log1p", Numerics.log1p, Numerics.log1p, math.log1p),
    ("expm1", Numerics.expm1, Numerics.expm1, math.expm1)
  )

  test("exp/log/log1p/expm1 equal pointwise.map of the same scalar function bit for bit") {
    for (name, onVec, onMat, f) <- elementwise do
      for (layout, x) <- vectorLayouts(specials) do
        assertBits(onVec(x).toSeq, specials.map(f), s"$name $layout")
      for rows <- 1 to 5; cols <- 1 to 5 do
        val cells = (i: Int, j: Int) => specials((i * 7 + j * 3) % specials.length)
        for (layout, a) <- matrixLayouts(rows, cols, cells) do
          val got = onMat(a)
          assertEquals((got.rows, got.cols), (rows, cols))
          assert(got.isContiguousRowMajor)
          assertBits(got.valuesRowMajor, a.pointwise.map(f).valuesRowMajor, s"$name $layout ${rows}x$cols")
  }

  test("sigmoid is overflow-free, saturates exactly, and matches the textbook formula in range") {
    val x = Vec(0.0, 1.0, -1.0, 20.0, -20.0, 709.0, -709.0, 800.0, -800.0, 1000.0, -1000.0, PInf, NInf, Nan)
    val s = Numerics.sigmoid(x).toSeq
    assertEquals(s.head, 0.5)
    assertClose(s(1), 1.0 / (1.0 + math.exp(-1.0)))
    assertClose(s(2), 1.0 / (1.0 + math.exp(1.0)))
    assertClose(s(3), 1.0 / (1.0 + math.exp(-20.0)))
    assertClose(s(4) / (1.0 / (1.0 + math.exp(20.0))), 1.0)
    assertEquals(s(5), 1.0)
    assertClose(s(6) / math.exp(-709.0), 1.0)
    assertEquals(s(7), 1.0)
    assert(s(8) == 0.0)
    assertEquals(s(9), 1.0)
    assertEquals(s(10), 0.0)
    assertEquals(s(11), 1.0)
    assertEquals(s(12), 0.0)
    assert(s(13).isNaN)
    for (layout, a) <- matrixLayouts(3, 4, (i, j) => (i - j) * 3.5) do
      assertBits(
        Numerics.sigmoid(a).valuesRowMajor,
        Numerics.sigmoid(Vec(a.valuesRowMajor*)).toSeq,
        layout
      )
  }

  test("logSumExp is stable for large magnitudes and follows the documented special-value rules") {
    assertClose(Numerics.logSumExp(Vec(1000.0, 1000.0)), 1000.0 + math.log(2.0))
    assertClose(Numerics.logSumExp(Vec(-1000.0, -1000.0)), -1000.0 + math.log(2.0))
    assertClose(Numerics.logSumExp(Vec(1000.0, -1000.0)), 1000.0)
    assertClose(Numerics.logSumExp(Vec(1.0, 2.0, 3.0)), math.log(math.exp(1.0) + math.exp(2.0) + math.exp(3.0)))
    assertEquals(Numerics.logSumExp(Vec()), NInf)
    assertEquals(Numerics.logSumExp(Vec(NInf, NInf, NInf)), NInf)
    assertEquals(Numerics.logSumExp(Vec(NInf, 0.0)), 0.0)
    assertEquals(Numerics.logSumExp(Vec(1.0, PInf, NInf)), PInf)
    assert(Numerics.logSumExp(Vec(1.0, Nan)).isNaN)
    assert(Numerics.logSumExp(Vec(PInf, Nan)).isNaN)
    assert(Numerics.logSumExp(Vec(Nan, PInf)).isNaN)
    for (layout, x) <- vectorLayouts(Seq(0.5, -2.0, 1000.0, 999.0, 3.0)) do
      assertClose(Numerics.logSumExp(x), 1000.0 + math.log1p(math.exp(-1.0)), 1e-15)
      assert(layout.nonEmpty)
  }

  test("softmax and logSoftmax are stable at +-1000 and consistent with each other") {
    assertEquals(Numerics.softmax(Vec(1000.0, 1000.0)).toSeq, Seq(0.5, 0.5))
    assertEquals(Numerics.softmax(Vec(-1000.0, -1000.0)).toSeq, Seq(0.5, 0.5))
    assertEquals(Numerics.softmax(Vec(1000.0, -1000.0)).toSeq, Seq(1.0, 0.0))
    assertEquals(Numerics.softmax(Vec(NInf, 3.0)).toSeq, Seq(0.0, 1.0))
    Numerics.logSoftmax(Vec(1000.0, 1000.0)).toSeq.foreach(v => assertClose(v, -math.log(2.0)))
    assertEquals(Numerics.logSoftmax(Vec(1000.0, -1000.0)).toSeq, Seq(0.0, -2000.0))
    assertEquals(Numerics.logSoftmax(Vec(NInf, 3.0)).toSeq, Seq(NInf, 0.0))
    for (layout, x) <- vectorLayouts(Seq(0.1, -3.0, 2.5, 7.0, 7.0, -0.25)) do
      val s = Numerics.softmax(x)
      val ls = Numerics.logSoftmax(x)
      assertClose(s.sumExact, 1.0, 1e-15)
      s.toSeq.zip(ls.toSeq).foreach((p, lp) => assertClose(p, math.exp(lp), 1e-15))
      assertEquals(s.argmax, x.argmax, layout)
  }

  test("softmax and logSoftmax: empty is empty, non-finite maxima give all NaN") {
    assertEquals(Numerics.softmax(Vec()).length, 0)
    assertEquals(Numerics.logSoftmax(Vec()).length, 0)
    for xs <- Seq(Seq(1.0, Nan, 2.0), Seq(1.0, PInf, 2.0), Seq(NInf, NInf), Seq(PInf, PInf)) do
      assert(Numerics.softmax(Vec(xs*)).toSeq.forall(_.isNaN), xs.toString)
      assert(Numerics.logSoftmax(Vec(xs*)).toSeq.forall(_.isNaN), xs.toString)
  }

  test("matrix logSumExp/softmax/logSoftmax: whole-matrix and per-axis forms on every layout") {
    val random = new scala.util.Random(11L)
    val cells = IndexedSeq.fill(4, 5)(random.nextGaussian() * 300.0)
    for (layout, a) <- matrixLayouts(4, 5, cells(_)(_)) do
      val flat = Vec(cells.flatten*)
      assertClose(Numerics.logSumExp(a), Numerics.logSumExp(flat), 1e-15)
      assertBits(Numerics.softmax(a).valuesRowMajor, Numerics.softmax(flat).toSeq, layout)
      assertBits(Numerics.logSoftmax(a).valuesRowMajor, Numerics.logSoftmax(flat).toSeq, layout)
      val rowsLse = Numerics.logSumExp(a, Axis.Rows).toSeq
      val colsLse = Numerics.logSumExp(a, Axis.Cols).toSeq
      assertEquals(rowsLse.length, 4)
      assertEquals(colsLse.length, 5)
      for i <- 0 until 4 do assertEquals(rowsLse(i), Numerics.logSumExp(Vec(cells(i)*)), layout)
      for j <- 0 until 5 do assertEquals(colsLse(j), Numerics.logSumExp(Vec(cells.map(_(j))*)), layout)
      val sRows = Numerics.softmax(a, Axis.Rows)
      val sCols = Numerics.softmax(a, Axis.Cols)
      val lsRows = Numerics.logSoftmax(a, Axis.Rows)
      val lsCols = Numerics.logSoftmax(a, Axis.Cols)
      for i <- 0 until 4 do
        assertBits(sRows.row(i).toSeq, Numerics.softmax(Vec(cells(i)*)).toSeq, s"$layout row $i")
        assertBits(lsRows.row(i).toSeq, Numerics.logSoftmax(Vec(cells(i)*)).toSeq, s"$layout row $i")
      for j <- 0 until 5 do
        val column = Vec(cells.map(_(j))*)
        assertBits(sCols.col(j).toSeq, Numerics.softmax(column).toSeq, s"$layout col $j")
        assertBits(lsCols.col(j).toSeq, Numerics.logSoftmax(column).toSeq, s"$layout col $j")
      assertClose(Numerics.logSoftmax(a).valuesRowMajor.head, cells(0)(0) - Numerics.logSumExp(flat), 1e-13)
  }

  test("matrix numerics on empty shapes") {
    val noRows = Matrix.zeros(0, 3)
    assertEquals(Numerics.exp(noRows).shape, noRows.shape)
    assertEquals(Numerics.softmax(noRows, Axis.Cols).shape, noRows.shape)
    assertEquals(Numerics.logSumExp(noRows), NInf)
    assertEquals(Numerics.logSumExp(noRows, Axis.Rows).length, 0)
    assertEquals(Numerics.logSumExp(noRows, Axis.Cols).toSeq, Seq(NInf, NInf, NInf))
    assertEquals(Numerics.softmax(Matrix.zeros(2, 0), Axis.Rows).shape, Matrix.zeros(2, 0).shape)
  }

  property("softmax is a probability vector and argmax-preserving") {
    forAll(Gen.nonEmptyListOf(Gen.choose(-800.0, 800.0))) { xs =>
      val x = Vec(xs*)
      val s = Numerics.softmax(x)
      assert(s.toSeq.forall(p => p >= 0.0 && p <= 1.0))
      assert(math.abs(s.sumExact - 1.0) <= 4.0 * xs.length * 2.3e-16)
      assertEquals(s.max, s(x.argmax))
      val lse = Numerics.logSumExp(x)
      assert(lse >= x.max && lse <= x.max + math.log(xs.length.toDouble) + 1e-12)
    }
  }
