package gale.parity

import breeze.linalg.DenseMatrix
import breeze.linalg.DenseVector
import breeze.linalg.norm as breezeNorm
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.syntax.all.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

/** Breeze differential tests for the operations that appear most often while
  * porting application code: construction, access, slicing, indexed selection,
  * persistent updates, norms, and elementwise transforms.
  *
  * The generators produce dimensions and deterministic data seeds. A failing
  * property therefore reports a small shape and seed that can be replayed with
  * [[ParitySupport.matrixData]] or [[ParitySupport.vectorData]].
  */
class EverydayOpsParitySuite extends ScalaCheckSuite:

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(40).withWorkers(1)

  private val dimensionGen = Gen.choose(1, 12)
  private val seedGen = Gen.choose(1L, 9_000_000L)
  private val scalarGen = Gen.chooseNum(-10.0, 10.0)

  private val matrixCaseGen: Gen[(Int, Int, Long)] =
    for
      rows <- dimensionGen
      cols <- dimensionGen
      seed <- seedGen
    yield (rows, cols, seed)

  private val sliceCaseGen: Gen[(Int, Int, Long, Int, Int, Int, Int)] =
    for
      rows <- dimensionGen
      cols <- dimensionGen
      seed <- seedGen
      rowFrom <- Gen.choose(0, rows - 1)
      rowUntil <- Gen.choose(rowFrom + 1, rows)
      colFrom <- Gen.choose(0, cols - 1)
      colUntil <- Gen.choose(colFrom + 1, cols)
    yield (rows, cols, seed, rowFrom, rowUntil, colFrom, colUntil)

  test("literal, zero, identity, and tabulated constructors match Breeze") {
    val literalGale = Matrix(2, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0
    )
    val literalBreeze = DenseMatrix(
      (1.0, 2.0, 3.0),
      (4.0, 5.0, 6.0)
    )
    assertMatClose(literalGale, literalBreeze, 0.0, "literal matrix")
    assertMatClose(Matrix.zeros(3, 4), DenseMatrix.zeros[Double](3, 4), 0.0, "zeros")
    assertMatClose(Matrix.eye(5), DenseMatrix.eye[Double](5), 0.0, "identity")

    val tabulatedGale = Matrix.tabulate(4, 3)((i, j) => i * 10.0 + j)
    val tabulatedBreeze = DenseMatrix.tabulate(4, 3)((i, j) => i * 10.0 + j)
    assertMatClose(tabulatedGale, tabulatedBreeze, 0.0, "tabulate")
  }

  property("row, column, and contiguous slice access match Breeze") {
    forAll(sliceCaseGen) {
      (sample: (Int, Int, Long, Int, Int, Int, Int)) =>
        val (rows, cols, seed, rowFrom, rowUntil, colFrom, colUntil) = sample
        val data = matrixData(rows, cols, seed)
        val gale = galeMatrix(data)
        val breeze = breezeMatrix(data)
        val expectedRow = DenseVector.tabulate(cols)(j => breeze(rowFrom, j))
        val expectedColumn = DenseVector.tabulate(rows)(i => breeze(i, colFrom))

        assertVecClose(
          gale.row(rowFrom),
          expectedRow,
          0.0,
          s"row rows=$rows cols=$cols seed=$seed row=$rowFrom"
        )
        assertVecClose(
          gale.col(colFrom),
          expectedColumn,
          0.0,
          s"column rows=$rows cols=$cols seed=$seed col=$colFrom"
        )
        assertMatClose(
          gale.slice(rowFrom, rowUntil, colFrom, colUntil),
          breeze(rowFrom until rowUntil, colFrom until colUntil),
          0.0,
          s"slice rows=$rows cols=$cols seed=$seed"
        )
    }
  }

  property("indexed gathers preserve order and repeated indices") {
    forAll(matrixCaseGen) { (sample: (Int, Int, Long)) =>
      val (rows, cols, seed) = sample
      val data = matrixData(rows, cols, seed)
      val gale = galeMatrix(data)
      val breeze = breezeMatrix(data)
      val rowIndices = IndexedSeq(rows - 1, 0, rows - 1)
      val colIndices = IndexedSeq(cols - 1, 0, cols - 1)

      val expectedRows =
        DenseMatrix.tabulate(rowIndices.length, cols)((i, j) => breeze(rowIndices(i), j))
      val expectedColumns =
        DenseMatrix.tabulate(rows, colIndices.length)((i, j) => breeze(i, colIndices(j)))

      assertMatClose(gale.gatherRows(rowIndices), expectedRows, 0.0, s"gather rows $sample")
      assertMatClose(gale.gatherColumns(colIndices), expectedColumns, 0.0, s"gather columns $sample")
    }
  }

  property("persistent matrix and vector updates match an updated Breeze copy") {
    forAll(matrixCaseGen, scalarGen) {
      (sample: (Int, Int, Long), replacement: Double) =>
        val (rows, cols, seed) = sample
        val data = matrixData(rows, cols, seed)
        val row = Math.floorMod(seed, rows.toLong).toInt
        val col = Math.floorMod(seed * 17L, cols.toLong).toInt

        val breeze = breezeMatrix(data)
        val updatedBreeze = breeze.copy
        updatedBreeze(row, col) = replacement
        assertMatClose(
          galeMatrix(data).updated(row, col, replacement),
          updatedBreeze,
          0.0,
          s"matrix updated $sample"
        )

        val vector = vectorData(cols, seed + 1L)
        val vectorIndex = Math.floorMod(seed, cols.toLong).toInt
        val updatedBreezeVector = breezeVector(vector).copy
        updatedBreezeVector(vectorIndex) = replacement
        assertVecClose(
          galeVector(vector).updated(vectorIndex, replacement),
          updatedBreezeVector,
          0.0,
          s"vector updated length=$cols seed=$seed"
        )
    }
  }

  property("vector slices, gathers, and Euclidean norms match Breeze") {
    forAll(dimensionGen, seedGen) { (length: Int, seed: Long) =>
      val data = vectorData(length, seed)
      val from = Math.floorMod(seed, length.toLong).toInt
      val until = length
      val indices = IndexedSeq(length - 1, 0, length - 1)
      val gale = galeVector(data)
      val breeze = breezeVector(data)

      val expectedSlice = DenseVector.tabulate(until - from)(i => breeze(from + i))
      val expectedGather = DenseVector.tabulate(indices.length)(i => breeze(indices(i)))
      assertVecClose(gale.slice(from, until), expectedSlice, 0.0, s"vector slice length=$length seed=$seed")
      assertVecClose(gale.gather(indices), expectedGather, 0.0, s"vector gather length=$length seed=$seed")
      assertScalarClose(gale.norm2, breezeNorm(breeze), 1e-14, s"vector norm length=$length seed=$seed")
    }
  }

  property("pointwise multiply, divide, and map match Breeze") {
    forAll(matrixCaseGen) { (sample: (Int, Int, Long)) =>
      val (rows, cols, seed) = sample
      val leftData = matrixData(rows, cols, seed)
      val rightData = matrixData(rows, cols, seed + 1L)
      val galeLeft = galeMatrix(leftData)
      val galeRight = galeMatrix(rightData)
      val breezeLeft = breezeMatrix(leftData)
      val breezeRight = breezeMatrix(rightData)

      assertMatClose(
        galeLeft.pointwise * galeRight,
        breezeLeft *:* breezeRight,
        0.0,
        s"pointwise multiply $sample"
      )
      assertMatClose(
        galeLeft.pointwise / galeRight,
        breezeLeft /:/ breezeRight,
        0.0,
        s"pointwise divide $sample"
      )
      assertMatClose(
        galeLeft.pointwise.map(math.tanh),
        breezeLeft.mapValues(math.tanh),
        0.0,
        s"pointwise map $sample"
      )
    }
  }

  property("diagonal addition and average symmetrization match Breeze expressions") {
    forAll(dimensionGen, seedGen, scalarGen) {
      (size: Int, seed: Long, diagonalShift: Double) =>
        val data = matrixData(size, size, seed)
        val gale = galeMatrix(data)
        val breeze = breezeMatrix(data)

        val shiftedBreeze = breeze.copy
        var i = 0
        while i < size do
          shiftedBreeze(i, i) = shiftedBreeze(i, i) + diagonalShift
          i += 1

        assertMatClose(
          gale.addToDiagonal(diagonalShift),
          shiftedBreeze,
          0.0,
          s"addToDiagonal size=$size seed=$seed"
        )
        assertMatClose(
          gale.symmetrizedAverage,
          (breeze + breeze.t) * 0.5,
          0.0,
          s"symmetrizedAverage size=$size seed=$seed"
        )
    }
  }
