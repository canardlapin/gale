package gale.laws

import gale.linalg.*
import gale.sparse.*
import munit.ScalaCheckSuite
import org.scalacheck.{Gen, Shrink}
import org.scalacheck.Prop.forAll

/** Replayable generators and structure-preserving shrinkers for cases that are
  * easy to lose when arbitrary tuples shrink outside their original bounds.
  */
class AdversarialGeneratorSuite extends ScalaCheckSuite:
  override def scalaCheckInitialSeed =
    "3IpKYLWqvse9f3GvOj9DgO4pqDgfpJ3mSfRFzU9i4yL="

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(40).withWorkers(1)

  private final case class ShapeCase(rows: Int, cols: Int, seed: Long)

  private given Shrink[ShapeCase] = Shrink: sample =>
    val rows =
      Shrink.shrink(sample.rows).filter(_ >= 0).map(value => sample.copy(rows = value))
    val cols =
      Shrink.shrink(sample.cols).filter(_ >= 0).map(value => sample.copy(cols = value))
    val seeds =
      Shrink.shrink(sample.seed).map(value => sample.copy(seed = value))
    rows.lazyAppendedAll(cols).lazyAppendedAll(seeds)

  private val shapeGen: Gen[ShapeCase] =
    for
      rows <- Gen.choose(0, 12)
      cols <- Gen.choose(0, 12)
      seed <- Gen.choose(0L, 9_000_000L)
    yield ShapeCase(rows, cols, seed)

  private final case class StrideCase(length: Int, gap: Int, seed: Long)

  private given Shrink[StrideCase] = Shrink: sample =>
    val lengths =
      Shrink.shrink(sample.length).filter(_ >= 0).map(value => sample.copy(length = value))
    val gaps =
      Shrink.shrink(sample.gap).filter(_ >= 1).map(value => sample.copy(gap = value))
    val seeds =
      Shrink.shrink(sample.seed).map(value => sample.copy(seed = value))
    lengths.lazyAppendedAll(gaps).lazyAppendedAll(seeds)

  private val strideGen: Gen[StrideCase] =
    for
      length <- Gen.choose(0, 12)
      gap <- Gen.choose(1, 5)
      seed <- Gen.choose(0L, 9_000_000L)
    yield StrideCase(length, gap, seed)

  private enum NonFiniteValue(val value: Double):
    case NaN              extends NonFiniteValue(Double.NaN)
    case PositiveInfinity extends NonFiniteValue(Double.PositiveInfinity)
    case NegativeInfinity extends NonFiniteValue(Double.NegativeInfinity)

  private final case class NonFiniteCase(
      rows: Int,
      cols: Int,
      row: Int,
      col: Int,
      nonFinite: NonFiniteValue
  ):
    def value: Double = nonFinite.value

  private given Shrink[NonFiniteCase] = Shrink: sample =>
    val rows =
      Shrink
        .shrink(sample.rows)
        .filter(_ >= 1)
        .map(value => sample.copy(rows = value, row = math.min(sample.row, value - 1)))
    val cols =
      Shrink
        .shrink(sample.cols)
        .filter(_ >= 1)
        .map(value => sample.copy(cols = value, col = math.min(sample.col, value - 1)))
    val row =
      Shrink.shrink(sample.row).filter(value => value >= 0 && value < sample.rows).map(value => sample.copy(row = value))
    val col =
      Shrink.shrink(sample.col).filter(value => value >= 0 && value < sample.cols).map(value => sample.copy(col = value))
    rows
      .lazyAppendedAll(cols)
      .lazyAppendedAll(row)
      .lazyAppendedAll(col)

  private val nonFiniteGen: Gen[NonFiniteCase] =
    for
      rows <- Gen.choose(1, 8)
      cols <- Gen.choose(1, 8)
      row <- Gen.choose(0, rows - 1)
      col <- Gen.choose(0, cols - 1)
      nonFinite <- Gen.oneOf(NonFiniteValue.values.toSeq)
    yield NonFiniteCase(rows, cols, row, col, nonFinite)

  property("degenerate and ordinary matrix shapes preserve transpose identity"):
    forAll(shapeGen): sample =>
      val random = new scala.util.Random(sample.seed)
      val values = IndexedSeq.fill(sample.rows * sample.cols)(random.nextDouble() * 2.0 - 1.0)
      val matrix = Matrix.dense(sample.rows, sample.cols, values)
      MatrixLaws.assertExact(matrix.t.t, matrix)

  property("strided views preserve logical values for zero and nonzero lengths"):
    forAll(strideGen): sample =>
      val random = new scala.util.Random(sample.seed)
      val cols = sample.gap + 1
      val matrix = Matrix.dense(
        sample.length,
        cols,
        IndexedSeq.fill(sample.length * cols)(random.nextDouble() * 2.0 - 1.0)
      )
      val strided = matrix.col(sample.gap)
      val contiguous = DVec.tabulate(sample.length)(row => matrix(row, sample.gap))
      VecLaws.assertClose(strided, contiguous, 0.0)

  property("strict sparse ingestion rejects every non-finite value without mutation"):
    forAll(nonFiniteGen): sample =>
      val strict =
        Sparse
          .cooChecked(sample.rows, sample.cols, SparseValuePolicy.RequireFinite)
          .toOption
          .get
      val result = strict.tryAdd(sample.row, sample.col, sample.value)
      assert(result.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
      assertEquals(strict.nnz, 0)

      val permissive = Sparse.cooChecked(sample.rows, sample.cols).toOption.get
      assert(permissive.tryAdd(sample.row, sample.col, sample.value).isRight)
      val observed = permissive.tryToCOO().toOption.get(sample.row, sample.col)
      if sample.value.isNaN then assert(observed.isNaN)
      else assertEquals(observed, sample.value)
