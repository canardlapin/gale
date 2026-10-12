package gale.backend.jvm.vector

import gale.backend.PureDenseDoubleKernel
import gale.platform.DoubleArray

import jdk.incubator.vector.DoubleVector

/** Parity for the row-major SIMD `gemv` tile body at fixed 2-, 4- and 8-lane species,
  * not only the host's preferred one. The eight-lane paths (the masked column step and
  * the rotate-and-add row-sum tree) are compiled only for an eight-lane species, so on
  * a two-lane host they run here through the Vector API's Java fallback.
  *
  * `A` is an interior view of a NaN-padded parent and `x`/`y` sit at offsets inside
  * NaN/sentinel padding, so any read past a row or past `x`, or any write past `y`,
  * shows up. The reference is the pure kernel on the same arrays; the SIMD tile
  * reassociates, so finite values agree within a tolerance; NaNs and infinities must
  * coincide.
  */
class VectorGemvRowMajorSuite extends munit.FunSuite:
  private val AbsTol = 1e-12
  private val RelTol = 1e-10
  private val Sentinel = 99.0

  // Species widths in doubles.
  private val Variants = Seq("2", "4", "8")

  private def runTiles(
      lanes: String,
      rows: Int,
      cols: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      x: Array[Double],
      xOffset: Int,
      beta: Double,
      y: Array[Double],
      yOffset: Int
  ): Unit =
    lanes match
      case "2" =>
        VectorDenseDoubleKernel.gemvRowMajorTiles(
          DoubleVector.SPECIES_128, rows, cols, alpha, a, aOffset, aRowStride, x, xOffset, beta, y, yOffset
        )
      case "4" =>
        VectorDenseDoubleKernel.gemvRowMajorTiles(
          DoubleVector.SPECIES_256, rows, cols, alpha, a, aOffset, aRowStride, x, xOffset, beta, y, yOffset
        )
      case "8" =>
        VectorDenseDoubleKernel.gemvRowMajorTiles(
          DoubleVector.SPECIES_512, rows, cols, alpha, a, aOffset, aRowStride, x, xOffset, beta, y, yOffset
        )
      case other => fail(s"unknown tile variant $other")

  private final case class Problem(
      rows: Int,
      cols: Int,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      x: Array[Double],
      xOffset: Int,
      y: Array[Double],
      yOffset: Int
  )

  /** `rows × cols` at row offset 2 and column offset 3 of a NaN-padded parent with
    * three spare columns; `x` at offset 3 and `y` at offset 2, both padded.
    */
  private def problem(
      rows: Int,
      cols: Int,
      fill: (Int, Int) => Double = (i, j) => ((i * 13 + j * 7) % 11 - 5).toDouble / 4.0 + i * 1e-3 - j * 1e-4,
      yInit: Int => Double = i => (i % 7).toDouble - 3.0
  ): Problem =
    val rowOff = 2
    val colOff = 3
    val parentCols = colOff + cols + 3
    val lead = 1
    val a = Array.fill(lead + (rowOff + rows + 1) * parentCols)(Double.NaN)
    for i <- 0 until rows; j <- 0 until cols do a(lead + (rowOff + i) * parentCols + colOff + j) = fill(i, j)
    val x = Array.fill(cols + 6)(Double.NaN)
    for j <- 0 until cols do x(3 + j) = ((j * 5) % 9 - 4).toDouble / 3.0 + j * 1e-3
    val y = Array.tabulate(rows + 4)(k => if k < 2 || k >= rows + 2 then Sentinel else yInit(k - 2))
    Problem(rows, cols, a, lead + rowOff * parentCols + colOff, parentCols, x, 3, y, 2)

  private def pureY(p: Problem, alpha: Double, beta: Double): Array[Double] =
    val y = p.y.clone()
    PureDenseDoubleKernel.gemv(
      p.rows, p.cols, alpha,
      DoubleArray.adopt(p.a), p.aOffset, p.aRowStride, 1,
      DoubleArray.adopt(p.x), p.xOffset, 1,
      beta, DoubleArray.adopt(y), p.yOffset, 1
    )
    y

  private def simdY(lanes: String, p: Problem, alpha: Double, beta: Double): Array[Double] =
    val y = p.y.clone()
    runTiles(lanes, p.rows, p.cols, alpha, p.a, p.aOffset, p.aRowStride, p.x, p.xOffset, beta, y, p.yOffset)
    y

  private def assertMatches(actual: Array[Double], expected: Array[Double], label: String)(using munit.Location): Unit =
    assertEquals(actual.length, expected.length)
    actual.indices.foreach { i =>
      val (a, e) = (actual(i), expected(i))
      if e.isNaN || a.isNaN then assert(a.isNaN && e.isNaN, s"$label: NaN mismatch at $i: simd=$a pure=$e")
      else if e.isInfinite || a.isInfinite then assertEquals(a, e, s"$label: infinity mismatch at $i")
      else
        val tolerance = AbsTol + RelTol * math.max(math.abs(a), math.abs(e))
        assert(math.abs(a - e) <= tolerance, s"$label: mismatch at $i: simd=$a pure=$e (|Δ|=${math.abs(a - e)})")
    }

  private def checkAllLanes(p: Problem, alpha: Double, beta: Double)(using munit.Location): Unit =
    val expected = pureY(p, alpha, beta)
    for lanes <- Variants do
      val label = s"${p.rows}x${p.cols} lanes=$lanes alpha=$alpha beta=$beta"
      val actual = simdY(lanes, p, alpha, beta)
      assertMatches(actual, expected, label)
      assertEquals(
        Seq(actual(0), actual(1), actual(p.rows + 2), actual(p.rows + 3)),
        Seq(Sentinel, Sentinel, Sentinel, Sentinel),
        s"$label: wrote outside y"
      )

  // Rows straddle the eight-row and four-row tiles (1..7 leftover rows); columns are
  // below one vector, exact multiples of 2/4/8 lanes, and every tail length.
  private val Shapes = Seq(
    (1, 1), (1, 7), (1, 16), (7, 1), (8, 1), (8, 5), (8, 8), (8, 13), (9, 50), (12, 3),
    (15, 17), (16, 64), (17, 67), (23, 131), (64, 50), (100, 7)
  )

  test("row-major tile parity at 2, 4 and 8 lanes over tails and leftover rows"):
    for (rows, cols) <- Shapes do checkAllLanes(problem(rows, cols), 1.0, 0.0)

  test("row-major tile parity with non-unit alpha and beta"):
    for (rows, cols) <- Shapes; (alpha, beta) <- Seq((-0.5, 2.0), (1.5, 1.0), (0.75, -1.25)) do
      checkAllLanes(problem(rows, cols), alpha, beta)

  test("row-major tile: n x 1, 1 x n and a 1000 x 50 design"):
    checkAllLanes(problem(37, 1), 1.0, 0.0)
    checkAllLanes(problem(1, 37), 1.0, 0.0)
    checkAllLanes(problem(1000, 50), 1.0, 0.0)

  test("row-major tile: beta = 0 overwrites a NaN y"):
    for (rows, cols) <- Shapes do
      val p = problem(rows, cols, yInit = _ => Double.NaN)
      checkAllLanes(p, 1.0, 0.0)
      for lanes <- Variants do
        val y = simdY(lanes, p, 1.0, 0.0)
        for i <- 0 until rows do assert(!y(p.yOffset + i).isNaN, s"${rows}x$cols lanes=$lanes row $i kept the prior NaN")

  test("row-major tile: alpha = 0 matches pure, with and without a NaN in A"):
    for (rows, cols) <- Shapes do
      checkAllLanes(problem(rows, cols), 0.0, 1.5)
      checkAllLanes(problem(rows, cols), 0.0, 0.0)
      checkAllLanes(problem(rows, cols, fill = (i, j) => if i == rows - 1 && j == cols - 1 then Double.NaN else 1.0), 0.0, 1.0)

  test("row-major tile: a NaN or infinity in A or x propagates to exactly its rows"):
    for (rows, cols) <- Shapes; poison <- Seq(Double.NaN, Double.PositiveInfinity) do
      // Last row and last column: inside the masked step or scalar tail on every species.
      val (pr, pc) = (rows - 1, cols - 1)
      val p = problem(rows, cols, fill = (i, j) => if i == pr && j == pc then poison else (i + j % 3).toDouble)
      checkAllLanes(p, 1.0, 0.0)
      for lanes <- Variants do
        val y = simdY(lanes, p, 1.0, 0.0)
        for i <- 0 until rows do
          val v = y(p.yOffset + i)
          if i == pr then assert(v.isNaN || v.isInfinite, s"${rows}x$cols lanes=$lanes: row $i lost $poison")
          else assert(v.isFinite, s"${rows}x$cols lanes=$lanes: row $i picked up $poison")
      // A NaN in x reaches every row.
      val q = problem(rows, cols)
      q.x(q.xOffset + pc) = Double.NaN
      checkAllLanes(q, 1.0, 0.0)
      for lanes <- Variants do
        val y = simdY(lanes, q, 1.0, 0.0)
        for i <- 0 until rows do assert(y(q.yOffset + i).isNaN, s"${rows}x$cols lanes=$lanes: row $i missed NaN x")

  test("row-major tile is deterministic run to run"):
    for (rows, cols) <- Shapes; lanes <- Variants do
      val p = problem(rows, cols)
      val first = simdY(lanes, p, 0.75, -1.25).map(java.lang.Double.doubleToRawLongBits)
      val second = simdY(lanes, p, 0.75, -1.25).map(java.lang.Double.doubleToRawLongBits)
      assertEquals(first.toSeq, second.toSeq, s"${rows}x$cols lanes=$lanes")
