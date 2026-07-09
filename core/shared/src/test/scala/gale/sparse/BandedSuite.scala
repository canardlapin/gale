package gale.sparse

import gale.linalg.*
import scala.util.Random

class BandedSuite extends munit.FunSuite:
  /** A dense matrix that is banded with lower/upper bandwidth `kl`/`ku`: band
    * entries are random nonzeros, everything outside the band is exactly zero.
    */
  private def bandedDense(rows: Int, cols: Int, kl: Int, ku: Int, rng: Random): DMat =
    Matrix.tabulate(rows, cols): (i, j) =>
      val d = i - j
      if d <= kl && d >= -ku then
        // A nonzero (offset keeps band edges detectable by fromDense).
        rng.between(1, 9).toDouble * (if (i + j) % 2 == 0 then 1.0 else -1.0)
      else 0.0

  private def denseDiagonal(a: DMat): DVec =
    val n = math.min(a.rows, a.cols)
    DVec.tabulate(n)(i => a(i, i))

  private def assertVecEquals(actual: DVec, expected: DVec): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assertEquals(actual(i), expected(i), s"at $i")
      i += 1

  private val shapes = Seq((5, 5), (4, 6), (6, 4))

  test("Banded.fromDense matches the dense matrix across shapes and bandwidths") {
    val rng = new Random(20260709L)
    for (rows, cols) <- shapes do
      val bandwidths = Seq((0, 0), (1, 0), (0, 1), (2, 1), (rows - 1, cols - 1))
      for (kl, ku) <- bandwidths do
        val dense = bandedDense(rows, cols, kl, ku, rng)
        val banded = Banded.fromDense(dense)
        val label = s"${rows}x$cols kl=$kl ku=$ku"

        // apply agrees at every position.
        var i = 0
        while i < rows do
          var j = 0
          while j < cols do
            assertEquals(banded(i, j), dense(i, j), s"$label apply($i,$j)")
            j += 1
          i += 1

        // toDense round-trips.
        assertEquals(banded.toDense().valuesRowMajor, dense.valuesRowMajor, s"$label toDense")

        // matvec and transpose-matvec agree.
        val x = DVec.tabulate(cols)(k => (k + 1).toDouble)
        assertVecEquals(banded * x, dense * x)
        val xt = DVec.tabulate(rows)(k => (k + 1).toDouble)
        val bandedT = MutableDVec.zeros(cols)
        banded.transposeApplyTo(xt, bandedT)
        assertVecEquals(bandedT.asVec, dense.t * xt)

        // diagonal, trace.
        assertVecEquals(banded.diagonal, denseDiagonal(dense))
        assertEquals(banded.trace, denseDiagonal(dense).toSeq.sum, s"$label trace")

        // transpose is a banded matrix equal to the dense transpose.
        assertEquals(banded.t.toDense().valuesRowMajor, dense.t.valuesRowMajor, s"$label t")
  }

  test("kl=ku=0 is a diagonal matrix") {
    val dense = Matrix.dense(3, 3)(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0)
    val banded = Banded.fromDense(dense)
    assertEquals(banded.kl, 0)
    assertEquals(banded.ku, 0)
    assertEquals(banded.nnz, 3)
    assertVecEquals(banded * Vec(1.0, 1.0, 1.0), Vec(2.0, 3.0, 4.0))
  }

  test("Banded.fromDiagonals places sub/main/super diagonals") {
    // 4x4 tridiagonal: main = [2,2,2,2], super(+1) = [-1,-1,-1], sub(-1) = [-1,-1,-1].
    val banded =
      Banded.fromDiagonals(
        4,
        4,
        Map(
          0 -> Seq(2.0, 2.0, 2.0, 2.0),
          1 -> Seq(-1.0, -1.0, -1.0),
          -1 -> Seq(-1.0, -1.0, -1.0)
        )
      )
    assertEquals(banded.kl, 1)
    assertEquals(banded.ku, 1)
    val dense = Matrix.dense(4, 4)(
      2.0, -1.0, 0.0, 0.0,
      -1.0, 2.0, -1.0, 0.0,
      0.0, -1.0, 2.0, -1.0,
      0.0, 0.0, -1.0, 2.0
    )
    assertEquals(banded.toDense().valuesRowMajor, dense.valuesRowMajor)
    assertVecEquals(banded * Vec(1.0, 2.0, 3.0, 4.0), dense * Vec(1.0, 2.0, 3.0, 4.0))
  }

  test("density reflects the stored band footprint") {
    val banded = Banded.fromDiagonals(3, 3, Map(0 -> Seq(1.0, 1.0, 1.0)))
    assertEquals(banded.nnz, 3)
    assertEqualsDouble(banded.density, 3.0 / 9.0, 1e-12)
  }
