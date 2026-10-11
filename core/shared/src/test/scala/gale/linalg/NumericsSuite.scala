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

  /** Bit-identical to the per-line vector kernel unless the axis is streamed
    * slice by slice, which reassociates the shifted sum; then within a few ulps.
    */
  private def assertLine(actual: Seq[Double], expected: Seq[Double], streamed: Boolean, clue: String)(using
      munit.Location
  ): Unit =
    if streamed then
      assertEquals(actual.length, expected.length, clue)
      actual.zip(expected).foreach { (x, y) =>
        assert((x.isNaN && y.isNaN) || x == y || math.abs(x - y) <= 1e-14 * math.max(1.0, math.abs(y)), s"$clue: $x vs $y")
      }
    else assertBits(actual, expected, clue)

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

  /** `sigmoid(x)` to 40 significant digits for finite `x`, from `java.math.BigDecimal`
    * series: `e^y = e^floor(y) * e^frac(y)`, Taylor for the fraction.
    */
  private object SigmoidReference:
    import java.math.{BigDecimal as BD, MathContext, RoundingMode}
    private val Mc = new MathContext(40, RoundingMode.HALF_EVEN)
    private val E =
      var e = BD.ONE
      var term = BD.ONE
      for k <- 1 until 60 do
        term = term.divide(BD.valueOf(k.toLong), Mc)
        e = e.add(term, Mc)
      e
    private val powers = Array.iterate(BD.ONE, 791)(_.multiply(E, Mc))

    private def expNonNegative(y: Double): BD =
      val n = math.floor(y).toInt
      val f = new BD(y).subtract(BD.valueOf(n.toLong))
      var s = BD.ONE
      var term = BD.ONE
      for k <- 1 until 45 do
        term = term.multiply(f, Mc).divide(BD.valueOf(k.toLong), Mc)
        s = s.add(term, Mc)
      s.multiply(powers(n), Mc)

    def apply(x: Double): BD =
      if x > 790.0 then BD.ONE // 1 - e^-790 is far within half an ulp of 1
      else if x < -790.0 then BD.ZERO // e^-790 is far below half of MinPositiveValue
      else if x >= 0.0 then BD.ONE.divide(BD.ONE.add(BD.ONE.divide(expNonNegative(x), Mc), Mc), Mc)
      else
        val t = BD.ONE.divide(expNonNegative(-x), Mc)
        t.divide(BD.ONE.add(t, Mc), Mc)

    /** `|got - sigmoid(x)|` in ulps of the correctly rounded result. */
    def ulpError(got: Double, x: Double): Double =
      val r = apply(x)
      val rounded = r.doubleValue
      val ulp = if rounded == 0.0 then Double.MinPositiveValue else math.ulp(rounded)
      new BD(got).subtract(r, Mc).abs.divide(new BD(ulp), Mc).doubleValue

  /** Measured worst case of the kernel against [[SigmoidReference]] over 10^6
    * points is 1.97 ulps (and 2.37 for the previous `exp(-|x|)` form); `Math.exp`
    * may differ by platform within its 1-ulp contract, hence the margin.
    */
  private val SigmoidUlps = 3.0

  test("sigmoid is within a few ulps of a 40-digit reference across both regimes and at edges") {
    val floor = gale.kernel.DoubleKernels.SigmoidDirectFloor
    val edges = Seq(0.0, -0.0, 1e-300, -1e-300, Double.MinPositiveValue, -Double.MinPositiveValue,
      java.lang.Double.MIN_NORMAL, -java.lang.Double.MIN_NORMAL, 1e-17, -1e-17, 0.5, -0.5, 1.0, -1.0, 36.7, -36.7,
      37.5, -37.5, floor, math.nextDown(floor), math.nextUp(floor), 709.0, -709.0, 709.79, -709.79, 744.0, -744.0,
      745.0, -745.0, 746.0, -746.0, 1000.0, -1000.0, Double.MaxValue, -Double.MaxValue)
    val rng = new scala.util.Random(20261011L)
    val random = Seq.fill(400)(rng.between(-40.0, 40.0)) ++ Seq.fill(200)(rng.between(-750.0, 750.0))
    val xs = edges ++ random
    for (layout, x) <- vectorLayouts(xs) do
      val s = Numerics.sigmoid(x).toSeq
      xs.zip(s).foreach { (v, g) =>
        val e = SigmoidReference.ulpError(g, v)
        assert(e <= SigmoidUlps, s"$layout sigmoid($v) = $g is $e ulps from the reference")
      }
  }

  test("sigmoid is exactly exp(x) once 1 + exp(x) rounds to 1, subnormal results included") {
    val xs = Seq(-37.5, -40.0, -100.0, -708.0, -709.0, -709.79, -720.0, -744.0, -745.0, -745.1, -746.0, -1000.0,
      -Double.MaxValue, Double.NegativeInfinity)
    assertBits(Numerics.sigmoid(Vec(xs*)).toSeq, xs.map(math.exp), "exp(x) tail")
    assert(Numerics.sigmoid(Vec(-744.0))(0) > 0.0 && Numerics.sigmoid(Vec(-744.0))(0) < java.lang.Double.MIN_NORMAL)
    assertBits(
      Numerics.sigmoid(Vec(0.0, -0.0, Double.MaxValue, Double.PositiveInfinity, Nan)).toSeq,
      Seq(0.5, 0.5, 1.0, 1.0, Nan),
      "specials"
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
      val streamRows = a.streamsAxis(Axis.Rows)
      val streamCols = a.streamsAxis(Axis.Cols)
      assertLine(rowsLse, (0 until 4).map(i => Numerics.logSumExp(Vec(cells(i)*))), streamRows, layout)
      assertLine(colsLse, (0 until 5).map(j => Numerics.logSumExp(Vec(cells.map(_(j))*))), streamCols, layout)
      val sRows = Numerics.softmax(a, Axis.Rows)
      val sCols = Numerics.softmax(a, Axis.Cols)
      val lsRows = Numerics.logSoftmax(a, Axis.Rows)
      val lsCols = Numerics.logSoftmax(a, Axis.Cols)
      for i <- 0 until 4 do
        assertLine(sRows.row(i).toSeq, Numerics.softmax(Vec(cells(i)*)).toSeq, streamRows, s"$layout row $i")
        assertLine(lsRows.row(i).toSeq, Numerics.logSoftmax(Vec(cells(i)*)).toSeq, streamRows, s"$layout row $i")
      for j <- 0 until 5 do
        val column = Vec(cells.map(_(j))*)
        assertLine(sCols.col(j).toSeq, Numerics.softmax(column).toSeq, streamCols, s"$layout col $j")
        assertLine(lsCols.col(j).toSeq, Numerics.logSoftmax(column).toSeq, streamCols, s"$layout col $j")
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

  test("per-axis log-domain forms keep the special-value rules on streamed and per-line layouts") {
    // Line 0: ordinary at +-1000; line 1: has NaN; line 2: has +Inf; line 3: all -Inf; line 4: -Inf and finite.
    val lines = IndexedSeq(
      IndexedSeq(1000.0, -1000.0, 999.0),
      IndexedSeq(1.0, Nan, 2.0),
      IndexedSeq(PInf, 1.0, PInf),
      IndexedSeq(NInf, NInf, NInf),
      IndexedSeq(NInf, 3.0, NInf)
    )
    // Lines as columns of a 3x5 matrix (Axis.Cols), and as rows of its transpose (Axis.Rows).
    val asColumns = matrixLayouts(3, 5, (i, j) => lines(j)(i))
    val asRows = matrixLayouts(5, 3, (i, j) => lines(i)(j))
    assert(asColumns.exists(_._2.streamsAxis(Axis.Cols)) && asRows.exists(_._2.streamsAxis(Axis.Rows)))
    for ((layout, a), axis) <- asColumns.map(_ -> Axis.Cols) ++ asRows.map(_ -> Axis.Rows) do
      val streamed = a.streamsAxis(axis)
      val expectedLse = lines.map(l => Numerics.logSumExp(Vec(l*)))
      assertLine(Numerics.logSumExp(a, axis).toSeq, expectedLse, streamed, s"$layout lse")
      val soft = Numerics.softmax(a, axis)
      val logSoft = Numerics.logSoftmax(a, axis)
      for (line, index) <- lines.zipWithIndex do
        def extract(m: DMat) = if axis == Axis.Cols then m.col(index).toSeq else m.row(index).toSeq
        assertLine(extract(soft), Numerics.softmax(Vec(line*)).toSeq, streamed, s"$layout softmax $index")
        assertLine(extract(logSoft), Numerics.logSoftmax(Vec(line*)).toSeq, streamed, s"$layout logSoftmax $index")
  }

  test("whole-matrix logSumExp on a strided view applies the special-value rules without copying") {
    for (cells, expected) <- Seq(
        ((i: Int, j: Int) => if (i, j) == (1, 2) then Nan else 1.0, Nan),
        ((i: Int, j: Int) => if (i, j) == (2, 0) then PInf else 1.0, PInf),
        ((_: Int, _: Int) => NInf, NInf),
        ((i: Int, j: Int) => if (i, j) == (0, 1) then 1000.0 else NInf, 1000.0)
      )
    do
      for (layout, a) <- matrixLayouts(3, 4, cells) do
        val got = Numerics.logSumExp(a)
        assert((got.isNaN && expected.isNaN) || got == expected, s"$layout: $got vs $expected")
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
