package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.LinearOperator
import gale.linalg.Matrix
import gale.linalg.MutableDVec
import scala.collection.mutable.ArrayBuffer

class GeneralizedBlockKernelsSuite extends munit.FunSuite:

  private def diagonalOperator(diagonal: IndexedSeq[Double]): DoubleLinearOperator =
    LinearOperator.fromFunction(diagonal.length, diagonal.length): (x, into) =>
      var i = 0
      while i < diagonal.length do
        into(i) = diagonal(i) * x(i)
        i += 1

  private def assertIdentity(matrix: DMat, tolerance: Double): Unit =
    assertEquals(matrix.rows, matrix.cols)
    var i = 0
    while i < matrix.rows do
      var j = 0
      while j < matrix.cols do
        val expected = if i == j then 1.0 else 0.0
        assertEqualsDouble(matrix(i, j), expected, tolerance, s"entry ($i,$j)")
        j += 1
      i += 1

  test("full-rank block is Cholesky-whitened in the B inner product") {
    val metric = diagonalOperator(IndexedSeq(1.0, 2.0, 3.0, 4.0))
    val candidates = Matrix.dense(4, 2)(
      1.0, 2.0,
      -1.0, 1.0,
      0.5, -2.0,
      3.0, 0.25
    )

    val block = GeneralizedBlockKernels.bOrthonormalize(candidates, metric).toOption.get
    assertEquals(block.cols, 2)
    assertIdentity(block.vectors.t * block.metricImages, 1e-11)
    assert(GeneralizedBlockKernels.metricOrthogonalityError(block) < 1e-11)
  }

  test("dependent initial columns are removed and deterministically replenished") {
    val metric = diagonalOperator(IndexedSeq(1.0, 2.0, 3.0, 4.0, 5.0))
    val repeated = DVec.fromSeq(Seq(1.0, -2.0, 0.5, 1.5, -1.0))
    val candidates = DMat.tabulate(5, 3)((i, _) => repeated(i))

    val first =
      GeneralizedBlockKernels.bOrthonormalizeAndReplenish(candidates, metric, targetColumns = 3).toOption.get
    val second =
      GeneralizedBlockKernels.bOrthonormalizeAndReplenish(candidates, metric, targetColumns = 3).toOption.get

    assertEquals(first.cols, 3)
    assertIdentity(first.vectors.t * first.metricImages, 1e-10)
    assertEquals(first.vectors.valuesRowMajor, second.vectors.valuesRowMajor)
  }

  test("orthogonal complement directions are rank-revealed and replenished deterministically") {
    val metric = diagonalOperator(IndexedSeq.fill(4)(1.0))
    val existingVectors = DMat.tabulate(4, 1)((row, _) => if row == 0 then 1.0 else 0.0)
    val existing = GeneralizedBlockKernels.MetricBlock(existingVectors, existingVectors)
    val repeated = DMat.tabulate(4, 2): (row, col) =>
      if row == 0 then 1.0
      else if row == 1 then (col + 1).toDouble
      else 0.0
    val candidates = GeneralizedBlockKernels.MetricBlock(repeated, repeated)

    val first = GeneralizedBlockKernels
      .bOrthonormalizeAgainstAndReplenish(candidates, metric, existing, targetColumns = 2)
      .toOption
      .get
    val second = GeneralizedBlockKernels
      .bOrthonormalizeAgainstAndReplenish(candidates, metric, existing, targetColumns = 2)
      .toOption
      .get
    val combined = GeneralizedBlockKernels.concatenate(existing, first).toOption.get

    assertEquals(first.cols, 2)
    assertIdentity(combined.vectors.t * combined.metricImages, 1e-11)
    assertEquals(first.vectors.valuesRowMajor, second.vectors.valuesRowMajor)
  }

  test("operator-retained mutable destinations cannot mutate a completed block image") {
    val retained = ArrayBuffer.empty[MutableDVec]
    val operator = new DoubleLinearOperator:
      def rows: Int = 3
      def cols: Int = 3
      def applyTo(x: DVec, into: MutableDVec): Unit =
        retained += into
        var i = 0
        while i < 3 do
          into(i) = (i + 1).toDouble * x(i)
          i += 1

    val image = GeneralizedBlockKernels
      .applyBlock(operator, Matrix.dense(3, 1)(2.0, 3.0, 4.0))
      .toOption
      .get
    val before = image.valuesRowMajor
    retained.foreach: destination =>
      var i = 0
      while i < destination.length do
        destination(i) = Double.NaN
        i += 1

    assertEquals(image.valuesRowMajor, before)
    assertEquals(image.valuesRowMajor.toSeq, Seq(2.0, 6.0, 12.0))
  }

  test("rank-revealing path reports encountered indefinite metric geometry") {
    val metric = diagonalOperator(IndexedSeq(1.0, -1.0, 2.0))
    val candidates = Matrix.eye(3)

    GeneralizedBlockKernels.bOrthonormalize(candidates, metric) match
      case Left(_: LinAlgError.NotPositiveDefinite) => ()
      case other                                    => fail(s"expected NotPositiveDefinite, got $other")
  }

  test("non-finite metric images fail before reaching projected arithmetic") {
    val metric = diagonalOperator(IndexedSeq(1.0, Double.NaN, 2.0))
    GeneralizedBlockKernels.bOrthonormalize(Matrix.eye(3), metric) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected InvalidArgument, got $other")
  }

  test("symmetricProjection averages roundoff asymmetry") {
    val basis = Matrix.eye(2)
    val image = Matrix.dense(2, 2)(
      2.0, 1.0 + 1e-12,
      1.0, 4.0
    )
    val projected = GeneralizedBlockKernels.symmetricProjection(basis, image).toOption.get
    assertEqualsDouble(projected(0, 1), 1.0 + 0.5e-12, 1e-15)
    assertEqualsDouble(projected(1, 0), projected(0, 1), 0.0)
  }

  test("small projected generalized solve returns the analytic pencil spectrum") {
    val projectedA = Matrix.dense(3, 3)(
      2.0, 0.0, 0.0,
      0.0, 6.0, 0.0,
      0.0, 0.0, 20.0
    )
    val projectedB = Matrix.dense(3, 3)(
      1.0, 0.0, 0.0,
      0.0, 2.0, 0.0,
      0.0, 0.0, 4.0
    )
    val result = GeneralizedBlockKernels
      .projectedGeneralizedEigen(projectedA, projectedB, 2, EigenOrder.SmallestAlgebraic)
      .toOption
      .get

    assertEqualsDouble(result.eigenvalues(0), 2.0, 1e-12)
    assertEqualsDouble(result.eigenvalues(1), 3.0, 1e-12)
    assert(result.diagnostics.orthogonalityError < 1e-12)
  }

  test("shape and tolerance violations are typed Left values") {
    val metric = diagonalOperator(IndexedSeq(1.0, 2.0, 3.0))
    assert(GeneralizedBlockKernels.bOrthonormalize(DMat.zeros(2, 1), metric).isLeft)
    assert(GeneralizedBlockKernels.bOrthonormalizeAndReplenish(DMat.zeros(3, 1), metric, 4).isLeft)
    assert(
      GeneralizedBlockKernels
        .bOrthonormalizeAndReplenish(DMat.zeros(3, 1), metric, 2, tolerance = Double.NaN)
        .isLeft
    )
  }
