package gale.linalg

import gale.numeric.ExactSum
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

class ReductionsSuite extends ScalaCheckSuite:
  private val Nan = Double.NaN
  private val PInf = Double.PositiveInfinity
  private val NInf = Double.NegativeInfinity

  private def values(n: Int, seed: Long): Seq[Double] =
    val random = new scala.util.Random(seed)
    Seq.fill(n)(random.nextGaussian() * 10.0)

  /** Contiguous, strided (a column of a row-major matrix), and sliced views of `xs`. */
  private def layouts(xs: Seq[Double]): Seq[(String, DVec)] =
    val n = xs.length
    val contiguous = Vec(xs*)
    // Column 1 of an n×3 row-major matrix: stride 3.
    val strided = Matrix.tabulate(n, 3)((i, j) => if j == 1 then xs(i) else 99.0).col(1)
    val sliced = Vec((Seq(-7.0, 8.0) ++ xs ++ Seq(5.0))*).slice(2, 2 + n)
    val stridedSlice = Matrix.tabulate(n + 2, 2)((i, j) => if j == 0 && i >= 1 && i <= n then xs(i - 1) else 77.0)
      .col(0)
      .slice(1, n + 1)
    Seq("contiguous" -> contiguous, "strided" -> strided, "sliced" -> sliced, "strided slice" -> stridedSlice)

  private def assertClose(actual: Double, expected: Double, tolerance: Double = 1e-12)(using munit.Location): Unit =
    assert(
      math.abs(actual - expected) <= tolerance * math.max(1.0, math.abs(expected)),
      s"actual=$actual expected=$expected"
    )

  private def assertSameDouble(actual: Double, expected: Double)(using munit.Location): Unit =
    assert(
      (actual.isNaN && expected.isNaN) || actual == expected,
      s"actual=$actual expected=$expected"
    )

  private def firstIndex(xs: Seq[Double], better: (Double, Double) => Boolean): Int =
    val nan = xs.indexWhere(_.isNaN)
    if nan >= 0 then nan
    else
      var best = 0
      var i = 1
      while i < xs.length do
        if better(xs(i), xs(best)) then best = i
        i += 1
      best

  test("DVec reductions agree with naive loops on every layout and unroll tail") {
    for n <- 1 to 13; (name, x) <- layouts(values(n, n.toLong)) do
      val xs = x.toSeq
      assertEquals(x.length, n, name)
      assertClose(x.sum, xs.sum)
      assertClose(x.mean, xs.sum / n)
      assertEquals(x.max, xs.max, s"$name n=$n")
      assertEquals(x.min, xs.min, s"$name n=$n")
      assertEquals(x.argmax, xs.indexOf(xs.max), s"$name n=$n")
      assertEquals(x.argmin, xs.indexOf(xs.min), s"$name n=$n")
      assertClose(x.norm1, xs.map(math.abs).sum)
      assertEquals(x.normInf, xs.map(math.abs).max, s"$name n=$n")
      assertEquals(x.sumExact, xs.foldLeft(ExactSum.zero()) { (acc, v) => acc.add(v); acc }.value)
  }

  test("empty vectors: sums and norms are zero, order statistics and mean throw EmptyInput") {
    for (name, x) <- layouts(Seq.empty) do
      assertEquals(x.sum, 0.0, name)
      assertEquals(x.sumExact, 0.0, name)
      assertEquals(x.norm1, 0.0, name)
      assertEquals(x.normInf, 0.0, name)
      assertEquals(x.norm2, 0.0, name)
      for (op, f) <- Seq[(String, DVec => Any)](
          "mean" -> (_.mean),
          "max" -> (_.max),
          "min" -> (_.min),
          "argmax" -> (_.argmax),
          "argmin" -> (_.argmin)
        )
      do
        val error = intercept[LinAlgError.EmptyInput](f(x))
        assertEquals(error.operation, s"DVec.$op")
        assertEquals(error.getMessage, s"DVec.$op requires a non-empty input")
  }

  test("length-one vectors") {
    for (_, x) <- layouts(Seq(-2.5)) do
      assertEquals(x.sum, -2.5)
      assertEquals(x.mean, -2.5)
      assertEquals(x.max, -2.5)
      assertEquals(x.min, -2.5)
      assertEquals(x.argmax, 0)
      assertEquals(x.argmin, 0)
      assertEquals(x.norm1, 2.5)
      assertEquals(x.normInf, 2.5)
  }

  test("the first NaN wins argmax/argmin in every unroll lane, and max/min/sum/normInf propagate NaN") {
    for n <- 1 to 11; at <- 0 until n; (name, x) <- layouts(values(n, 31L).updated(at, Nan)) do
      assertEquals(x.argmax, at, s"$name n=$n at=$at")
      assertEquals(x.argmin, at, s"$name n=$n at=$at")
      assert(x.max.isNaN && x.min.isNaN && x.sum.isNaN && x.normInf.isNaN && x.norm1.isNaN, s"$name n=$n")
    // Two NaNs: the earlier one is reported even when a later lane sees its NaN first.
    for (_, x) <- layouts(Seq(1.0, 2.0, 3.0, Nan, 0.0, Nan, 4.0, 9.0)) do
      assertEquals(x.argmax, 3)
      assertEquals(x.argmin, 3)
  }

  test("argmax/argmin return the first occurrence of a tie across unroll lanes") {
    val tied = Seq(1.0, 5.0, 2.0, 5.0, 5.0, 0.0, 5.0, -3.0, -3.0)
    for (name, x) <- layouts(tied) do
      assertEquals(x.argmax, 1, name)
      assertEquals(x.argmin, 7, name)
    // A tie between lane 3 (index 3) and lane 0 of the next block (index 4).
    for (name, x) <- layouts(Seq(0.0, 0.0, 0.0, 9.0, 9.0, 0.0, 0.0, 0.0)) do assertEquals(x.argmax, 3, name)
    // Signed zeros compare equal, so the first one wins.
    for (_, x) <- layouts(Seq(-0.0, 0.0, -0.0, 0.0, 0.0)) do
      assertEquals(x.argmax, 0)
      assertEquals(x.argmin, 0)
  }

  test("infinities: extrema, all-infinite inputs, and IEEE sums") {
    for n <- 1 to 9 do
      for (name, x) <- layouts(Seq.fill(n)(NInf)) do
        assertEquals(x.argmax, 0, name)
        assertEquals(x.max, NInf)
        assertEquals(x.sum, NInf)
        assertEquals(x.normInf, PInf)
      for (name, x) <- layouts(Seq.fill(n)(PInf)) do
        assertEquals(x.argmin, 0, name)
        assertEquals(x.min, PInf)
    for (_, x) <- layouts(Seq(1.0, PInf, 3.0, NInf, 2.0)) do
      assertEquals(x.argmax, 1)
      assertEquals(x.argmin, 3)
      assert(x.sum.isNaN)
      assertEquals(x.normInf, PInf)
      assertEquals(x.norm1, PInf)
  }

  test("sumExact is correctly rounded and order-independent where sum is not") {
    val xs = Seq(1e100, 1.0, -1e100, 1e-30, 3.0)
    val expected = xs.foldLeft(ExactSum.zero()) { (acc, v) => acc.add(v); acc }.value
    assertEquals(expected, 4.0)
    for perm <- xs.permutations; (name, x) <- layouts(perm) do assertEquals(x.sumExact, 4.0, name)
    assert(Vec(1e308, 1e308, -1e308).sumExact == 1e308)
    assert(Vec(1.0, Nan).sumExact.isNaN)
    assertEquals(Vec(PInf, 1.0).sumExact, PInf)
  }

  property("sum stays within the fast-sum error bound of sumExact") {
    forAll(Gen.listOf(Gen.choose(-1e6, 1e6))) { xs =>
      val x = Vec(xs*)
      val bound = xs.length.toDouble * 2.3e-16 * xs.map(math.abs).sum + 1e-300
      assert(math.abs(x.sum - x.sumExact) <= bound)
      if xs.nonEmpty then
        assertEquals(x.max, xs.max)
        assertEquals(x.argmin, xs.indexOf(xs.min))
    }
  }

  test("vector norm2 reports +Inf for several infinities and NaN only for NaN entries") {
    assertEquals(Vec(PInf, PInf).norm2, PInf)
    assertEquals(Vec(1.0, NInf, 2.0, PInf).norm2, PInf)
    assert(Vec(PInf, Nan).norm2.isNaN)
    assert(Vec(Nan, PInf, PInf).norm2.isNaN)
    assertClose(Vec(3.0, 4.0).norm2, 5.0)
  }

  // ---------------------------------------------------------------------------
  // Matrices
  // ---------------------------------------------------------------------------

  /** Row-major, column-major (transpose of a row-major), strided slice, and the
    * transpose of a strided slice, each presenting the logical matrix `rows×cols`.
    */
  private def matrixLayouts(rows: Int, cols: Int, f: (Int, Int) => Double): Seq[(String, DMat)] =
    val rowMajor = Matrix.tabulate(rows, cols)(f)
    val colMajor = Matrix.tabulate(cols, rows)((j, i) => f(i, j)).t
    val slice = Matrix.tabulate(rows + 2, cols + 3)((i, j) =>
      if i >= 1 && i <= rows && j >= 2 && j < cols + 2 then f(i - 1, j - 2) else 1e9
    ).slice(1, rows + 1, 2, cols + 2)
    val sliceT = Matrix.tabulate(cols + 1, rows + 2)((j, i) =>
      if j >= 1 && i >= 2 then f(i - 2, j - 1) else 1e9
    ).slice(1, cols + 1, 2, rows + 2).t
    Seq("row-major" -> rowMajor, "col-major" -> colMajor, "slice" -> slice, "transposed slice" -> sliceT)

  private def naive(a: DMat): IndexedSeq[IndexedSeq[Double]] =
    IndexedSeq.tabulate(a.rows, a.cols)((i, j) => a(i, j))

  test("matrix layouts in this suite present the intended logical values") {
    val layouts = matrixLayouts(3, 4, (i, j) => 10.0 * i + j)
    assert(layouts(0)._2.isContiguousRowMajor)
    assert(layouts(1)._2.isContiguousColMajor)
    assert(!layouts(2)._2.isContiguousRowMajor && !layouts(2)._2.isContiguousColMajor)
    for (name, a) <- layouts do assertEquals(naive(a), naive(layouts(0)._2), name)
  }

  test("whole-matrix reductions agree with naive loops on every layout") {
    for rows <- 1 to 6; cols <- 1 to 6 do
      val random = new scala.util.Random(rows * 31L + cols)
      val cells = IndexedSeq.fill(rows, cols)(random.nextGaussian())
      for (name, a) <- matrixLayouts(rows, cols, cells(_)(_)) do
        val flat = cells.flatten
        val label = s"$name ${rows}x$cols"
        assertClose(a.sum, flat.sum)
        assertEquals(a.sumExact, Vec(flat*).sumExact, label)
        assertClose(a.mean, flat.sum / flat.length)
        assertEquals(a.max, flat.max, label)
        assertEquals(a.min, flat.min, label)
        val k = flat.indexOf(flat.max)
        assertEquals(a.argmax, (k / cols, k % cols), label)
        val m = flat.indexOf(flat.min)
        assertEquals(a.argmin, (m / cols, m % cols), label)
        assertClose(a.norm1, (0 until cols).map(j => (0 until rows).map(i => math.abs(cells(i)(j))).sum).max)
        assertClose(a.normInf, cells.map(_.map(math.abs).sum).max)
        assertClose(a.normFrobenius, math.sqrt(flat.map(v => v * v).sum))
  }

  test("matrix argmax/argmin use row-major first occurrence and first NaN, independent of layout") {
    // Ties at (0, 2) and (1, 0): row-major order puts (0, 2) first, column-major would not.
    val tied = (i: Int, j: Int) => if (i, j) == (0, 2) || (i, j) == (1, 0) then 7.0 else 0.0
    for (name, a) <- matrixLayouts(3, 3, tied) do
      assertEquals(a.argmax, (0, 2), name)
      assertEquals(a.argmin, (0, 0), name)
    val withNaN = (i: Int, j: Int) => if (i, j) == (0, 2) || (i, j) == (1, 0) then Nan else (i + j).toDouble
    for (name, a) <- matrixLayouts(3, 3, withNaN) do
      assertEquals(a.argmax, (0, 2), name)
      assertEquals(a.argmin, (0, 2), name)
      assert(a.max.isNaN && a.min.isNaN && a.sum.isNaN, name)
      assert(a.norm1.isNaN && a.normInf.isNaN && a.normFrobenius.isNaN, name)
  }

  test("per-axis reductions: Axis.Rows gives one value per row, Axis.Cols one per column") {
    val a = Matrix(2, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 9.0
    )
    assertEquals(a.sum(Axis.Rows).toSeq, Seq(6.0, 18.0))
    assertEquals(a.sum(Axis.Cols).toSeq, Seq(5.0, 7.0, 12.0))
    assertEquals(a.mean(Axis.Rows).toSeq, Seq(2.0, 6.0))
    assertEquals(a.mean(Axis.Cols).toSeq, Seq(2.5, 3.5, 6.0))
    assertEquals(a.max(Axis.Rows).toSeq, Seq(3.0, 9.0))
    assertEquals(a.min(Axis.Cols).toSeq, Seq(1.0, 2.0, 3.0))
    assertEquals(a.t.sum(Axis.Cols).toSeq, Seq(6.0, 18.0))
  }

  test("per-axis reductions agree with per-line vector reductions on every layout") {
    for rows <- 1 to 7; cols <- 1 to 7 do
      val random = new scala.util.Random(rows * 131L + cols)
      val cells = IndexedSeq.fill(rows, cols)(random.nextGaussian())
      for (name, a) <- matrixLayouts(rows, cols, cells(_)(_)) do
        val label = s"$name ${rows}x$cols"
        val rowLines = cells.map(r => Vec(r*))
        val colLines = (0 until cols).map(j => Vec(cells.map(_(j))*))
        def check(axis: Axis, lines: IndexedSeq[DVec]): Unit =
          val sums = a.sum(axis).toSeq
          assertEquals(sums.length, lines.length, label)
          sums.zip(lines).foreach((s, line) => assertClose(s, line.sumExact))
          a.mean(axis).toSeq.zip(lines).foreach((s, line) => assertClose(s, line.sumExact / line.length))
          assertEquals(a.max(axis).toSeq, lines.map(_.max), label)
          assertEquals(a.min(axis).toSeq, lines.map(_.min), label)
        check(Axis.Rows, rowLines)
        check(Axis.Cols, colLines)
  }

  test("per-axis NaN propagation is per line") {
    for (name, a) <- matrixLayouts(2, 3, (i, j) => if (i, j) == (1, 1) then Nan else (i * 3 + j).toDouble) do
      val maxRows = a.max(Axis.Rows).toSeq
      assertEquals(maxRows.head, 2.0, name)
      assert(maxRows(1).isNaN, name)
      val sumCols = a.sum(Axis.Cols).toSeq
      assertEquals(sumCols(0), 3.0, name)
      assert(sumCols(1).isNaN, name)
      assertEquals(sumCols(2), 7.0, name)
  }

  test("empty matrices: totals and norms are zero, order statistics throw, per-axis follows line emptiness") {
    for (rows, cols) <- Seq((0, 0), (0, 3), (3, 0)) do
      val a = Matrix.zeros(rows, cols)
      assertEquals(a.sum, 0.0)
      assertEquals(a.sumExact, 0.0)
      assertEquals(a.norm1, 0.0)
      assertEquals(a.normInf, 0.0)
      assertEquals(a.normFrobenius, 0.0)
      for (op, f) <- Seq[(String, DMat => Any)](
          "mean" -> (_.mean),
          "max" -> (_.max),
          "min" -> (_.min),
          "argmax" -> (_.argmax),
          "argmin" -> (_.argmin)
        )
      do assertEquals(intercept[LinAlgError.EmptyInput](f(a)).operation, s"DMat.$op")
    // 0×3: no rows to reduce; three empty columns.
    val noRows = Matrix.zeros(0, 3)
    assertEquals(noRows.sum(Axis.Rows).length, 0)
    assertEquals(noRows.max(Axis.Rows).length, 0)
    assertEquals(noRows.sum(Axis.Cols).toSeq, Seq(0.0, 0.0, 0.0))
    assertEquals(intercept[LinAlgError.EmptyInput](noRows.mean(Axis.Cols)).operation, "DMat.mean(Axis.Cols)")
    assertEquals(intercept[LinAlgError.EmptyInput](noRows.max(Axis.Cols)).operation, "DMat.max(Axis.Cols)")
    // 3×0: three empty rows.
    val noCols = Matrix.zeros(3, 0)
    assertEquals(noCols.sum(Axis.Rows).toSeq, Seq(0.0, 0.0, 0.0))
    assertEquals(intercept[LinAlgError.EmptyInput](noCols.min(Axis.Rows)).operation, "DMat.min(Axis.Rows)")
    assertEquals(noCols.mean(Axis.Cols).length, 0)
  }

  test("1x1 matrices") {
    val a = Matrix(1, 1)(-3.0)
    assertEquals(a.sum, -3.0)
    assertEquals(a.argmax, (0, 0))
    assertEquals(a.norm1, 3.0)
    assertEquals(a.normInf, 3.0)
    assertEquals(a.normFrobenius, 3.0)
    assertEquals(a.sum(Axis.Rows).toSeq, Seq(-3.0))
  }

  test("normFrobenius is overflow- and underflow-safe and IEEE for non-finite entries") {
    for (name, a) <- matrixLayouts(3, 4, (_, _) => 1e300) do
      assertClose(a.normFrobenius, math.sqrt(12.0) * 1e300, 1e-14)
      assert(name.nonEmpty)
    for (_, a) <- matrixLayouts(2, 2, (_, _) => 1e-300) do assertClose(a.normFrobenius / 1e-300, 2.0, 1e-14)
    for (_, a) <- matrixLayouts(2, 3, (i, j) => if i == 0 && j == 0 then 1e300 else 1e-300) do
      assertClose(a.normFrobenius / 1e300, 1.0, 1e-15)
    for (_, a) <- matrixLayouts(2, 2, (i, _) => if i == 0 then PInf else NInf) do assertEquals(a.normFrobenius, PInf)
    for (_, a) <- matrixLayouts(2, 2, (i, j) => if i == 0 then PInf else if j == 0 then Nan else 1.0) do
      assert(a.normFrobenius.isNaN)
    // Agreement with the vector 2-norm in the ordinary range.
    val random = new scala.util.Random(7L)
    val cells = IndexedSeq.fill(5, 6)(random.nextGaussian())
    for (_, a) <- matrixLayouts(5, 6, cells(_)(_)) do assertClose(a.normFrobenius, Vec(cells.flatten*).norm2, 1e-14)
  }

  test("matrix norm1/normInf are the max absolute column/row sums and are transposes of each other") {
    val a = Matrix(2, 3)(
      1.0, -2.0, 3.0,
      -4.0, 5.0, -16.0
    )
    assertEquals(a.norm1, 19.0)
    assertEquals(a.normInf, 25.0)
    assertEquals(a.t.norm1, 25.0)
    assertEquals(a.t.normInf, 19.0)
  }
