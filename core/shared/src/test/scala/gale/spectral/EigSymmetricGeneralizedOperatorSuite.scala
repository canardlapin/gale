package gale.spectral

import gale.linalg.*
import gale.solvers.Preconditioner
import scala.collection.mutable.ArrayBuffer

class EigSymmetricGeneralizedOperatorSuite extends munit.FunSuite:

  private def diagonalOperator(diagonal: IndexedSeq[Double]): DoubleLinearOperator =
    LinearOperator.fromFunction(diagonal.length, diagonal.length): (x, into) =>
      var i = 0
      while i < diagonal.length do
        into(i) = diagonal(i) * x(i)
        i += 1

  private def diagonalMatrix(diagonal: IndexedSeq[Double]): DMat =
    DMat.tabulate(diagonal.length, diagonal.length): (row, col) =>
      if row == col then diagonal(row) else 0.0

  test("public typed operator facade agrees with the dense generalized solve") {
    val generalizedValues = IndexedSeq(-2.0, 0.5, 1.0, 3.0, 7.0, 12.0)
    val metricDiagonal = IndexedSeq(1.0, 2.0, 0.5, 3.0, 4.0, 1.5)
    val operatorDiagonal = generalizedValues
      .zip(metricDiagonal)
      .map:
        case (value, weight) => value * weight
    val operator = diagonalOperator(operatorDiagonal)
    val metric = diagonalOperator(metricDiagonal)

    val iterative = Eigen
      .eigSymmetricGeneralized(
        operator.assumeSymmetricOperator,
        metric.assumePositiveDefiniteOperator,
        generalizedValues.length,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(tolerance = 1e-10, maxIterations = 80)
      )
      .toOption
      .get
    val dense = Eigen
      .eigSymmetricGeneralized(
        diagonalMatrix(operatorDiagonal),
        diagonalMatrix(metricDiagonal),
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic)
      )
      .toOption
      .get

    assert(iterative.diagnostics.allConverged, iterative.diagnostics.toString)
    assertEquals(iterative.size, dense.size)
    var i = 0
    while i < dense.size do
      assertEqualsDouble(iterative.eigenvalues(i), dense.eigenvalues(i), 1e-9)
      assert(iterative.diagnostics.residuals(i) <= 1e-10)
      i += 1
    assert(iterative.diagnostics.orthogonalityError < 1e-10)
  }

  test("values-only mode preserves values and omits vectors") {
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0)
    val operator = diagonalOperator(values)
    val metric = diagonalOperator(IndexedSeq.fill(values.length)(1.0))
    val result = Eigen
      .eigSymmetricGeneralized(
        operator.assumeSymmetricOperator,
        metric.assumePositiveDefiniteOperator,
        values.length,
        EigenSelection.Count(2, EigenOrder.LargestAlgebraic),
        GeneralizedSpectralOptions(
          tolerance = 1e-10,
          maxIterations = 60,
          returnVectors = EigenVectors.ValuesOnly
        )
      )
      .toOption
      .get

    assert(result.diagnostics.allConverged)
    assertEquals(result.eigenvectors.rows, values.length)
    assertEquals(result.eigenvectors.cols, 0)
    assertEqualsDouble(result.eigenvalues(0), 8.0, 1e-9)
    assertEqualsDouble(result.eigenvalues(1), 16.0, 1e-9)
  }

  test("operator and preconditioner destination retention cannot mutate the result") {
    final class RetainingDiagonal(val diagonal: IndexedSeq[Double]) extends DoubleLinearOperator:
      val retained = ArrayBuffer.empty[MutableDVec]
      def rows: Int = diagonal.length
      def cols: Int = diagonal.length
      def applyTo(x: DVec, into: MutableDVec): Unit =
        retained += into
        var i = 0
        while i < diagonal.length do
          into(i) = diagonal(i) * x(i)
          i += 1
      def corrupt(): Unit =
        retained.foreach: destination =>
          var i = 0
          while i < destination.length do
            destination(i) = Double.NaN
            i += 1

    final class RetainingIdentity extends Preconditioner:
      val retained = ArrayBuffer.empty[MutableVec[Double]]
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        retained += into
        var i = 0
        while i < r.length do
          into(i) = r(i)
          i += 1
      def corrupt(): Unit =
        retained.foreach: destination =>
          var i = 0
          while i < destination.length do
            destination(i) = Double.NaN
            i += 1

    val a = new RetainingDiagonal(IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0))
    val b = new RetainingDiagonal(IndexedSeq.fill(6)(1.0))
    val preconditioner = new RetainingIdentity
    val result = Eigen
      .eigSymmetricGeneralized(
        a.assumeSymmetricOperator,
        b.assumePositiveDefiniteOperator,
        6,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(tolerance = 1e-10, maxIterations = 80),
        preconditioner
      )
      .toOption
      .get

    assert(result.diagnostics.allConverged)
    val valuesBefore = result.eigenvalues.toSeq
    val vectorsBefore = result.eigenvectors.valuesRowMajor
    a.corrupt()
    b.corrupt()
    preconditioner.corrupt()
    assertEquals(result.eigenvalues.toSeq, valuesBefore)
    assertEquals(result.eigenvectors.valuesRowMajor, vectorsBefore)
    assert(result.eigenvalues.toSeq.forall(_.isFinite))
    assert(result.eigenvectors.valuesRowMajor.forall(_.isFinite))
  }

  test("facade rejects malformed requests before applying either operator") {
    var applications = 0
    val counted = LinearOperator.fromFunction(4, 4): (x, into) =>
      applications += 1
      var i = 0
      while i < 4 do
        into(i) = x(i)
        i += 1
    val symmetric = counted.assumeSymmetricOperator
    val metric = counted.assumePositiveDefiniteOperator

    val denseOnly = Eigen.eigSymmetricGeneralized(
      symmetric,
      metric,
      4,
      EigenSelection.All
    )
    val badDimension = Eigen.eigSymmetricGeneralized(
      symmetric,
      metric,
      5,
      EigenSelection.Count(1, EigenOrder.SmallestAlgebraic)
    )
    val badOrder = Eigen.eigSymmetricGeneralized(
      symmetric,
      metric,
      4,
      EigenSelection.Count(1, EigenOrder.SmallestMagnitude)
    )
    val badOptions = Eigen.eigSymmetricGeneralized(
      symmetric,
      metric,
      4,
      EigenSelection.Count(1, EigenOrder.SmallestAlgebraic),
      GeneralizedSpectralOptions(tolerance = Double.NaN)
    )
    val badInitial = Eigen.eigSymmetricGeneralized(
      symmetric,
      metric,
      4,
      EigenSelection.Count(1, EigenOrder.SmallestAlgebraic),
      GeneralizedSpectralOptions(initialSubspace = Some(DMat.zeros(4, 2)))
    )
    val badVectors = Eigen.eigSymmetricGeneralized(
      symmetric,
      metric,
      4,
      EigenSelection.Count(1, EigenOrder.SmallestAlgebraic),
      GeneralizedSpectralOptions(returnVectors = EigenVectors.Left)
    )

    assert(denseOnly.isLeft)
    assert(badDimension.left.exists(_.isInstanceOf[LinAlgError.DimensionMismatch]))
    assert(badOrder.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(badOptions.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(badInitial.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(badVectors.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assertEquals(applications, 0)
  }

  test("residual convergence and extreme certification remain distinct") {
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
    val initial = DMat.tabulate(6, 2)((row, col) => if row == col then 1.0 else 0.0)
    val result = Eigen
      .eigSymmetricGeneralized(
        diagonalOperator(values).assumeSymmetricOperator,
        diagonalOperator(IndexedSeq.fill(6)(1.0)).assumePositiveDefiniteOperator,
        6,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(tolerance = 1e-12, initialSubspace = Some(initial))
      )
      .toOption
      .get

    assert(result.requireConverged.isRight)
    assert(
      result.requireExtremeCertified.left.exists(
        _.isInstanceOf[LinAlgError.SpectralExtremeNotCertified]
      )
    )
  }

  test("an exact invariant start explores its complement before accepting the requested end") {
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
    val unwantedLargest =
      DMat.tabulate(6, 2)((row, col) => if row == col + 4 then 1.0 else 0.0)
    val result = Eigen
      .eigSymmetricGeneralized(
        diagonalOperator(values).assumeSymmetricOperator,
        diagonalOperator(IndexedSeq.fill(6)(1.0)).assumePositiveDefiniteOperator,
        6,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(
          tolerance = 1e-10,
          maxIterations = 80,
          initialSubspace = Some(unwantedLargest)
        )
      )
      .toOption
      .get

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assertEqualsDouble(result.eigenvalues(0), 1.0, 1e-9)
    assertEqualsDouble(result.eigenvalues(1), 2.0, 1e-9)
    assert(result.diagnostics.iterations > 0)
  }
