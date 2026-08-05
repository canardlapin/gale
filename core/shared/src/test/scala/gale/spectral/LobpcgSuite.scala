package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.LinearOperator
import gale.linalg.Matrix
import gale.linalg.MutableVec
import gale.solvers.Preconditioner

class LobpcgSuite extends munit.FunSuite:

  private def diagonalOperator(diagonal: IndexedSeq[Double]): DoubleLinearOperator =
    LinearOperator.fromFunction(diagonal.length, diagonal.length): (x, into) =>
      var i = 0
      while i < diagonal.length do
        into(i) = diagonal(i) * x(i)
        i += 1

  private def solveDiagonal(
      generalizedValues: IndexedSeq[Double],
      metricDiagonal: IndexedSeq[Double],
      k: Int,
      order: EigenOrder,
      options: GeneralizedSpectralOptions,
      preconditioner: Preconditioner = Preconditioner.Identity
  ): Either[LinAlgError, EigenDecomposition] =
    val operatorDiagonal = generalizedValues
      .zip(metricDiagonal)
      .map:
        case (value, weight) => value * weight
    Lobpcg.solve(
      diagonalOperator(operatorDiagonal),
      diagonalOperator(metricDiagonal),
      generalizedValues.length,
      k,
      order,
      options,
      preconditioner
    )

  private def generalizedResidual(
      values: IndexedSeq[Double],
      metricDiagonal: IndexedSeq[Double],
      result: EigenDecomposition,
      column: Int
  ): Double =
    var sum = 0.0
    var row = 0
    while row < values.length do
      val residual =
        metricDiagonal(row) * (values(row) - result.eigenvalues(column)) * result.eigenvectors(row, column)
      sum += residual * residual
      row += 1
    math.sqrt(sum)

  private def metricOrthogonalityError(
      metricDiagonal: IndexedSeq[Double],
      vectors: DMat
  ): Double =
    var sum = 0.0
    var i = 0
    while i < vectors.cols do
      var j = 0
      while j < vectors.cols do
        var entry = 0.0
        var row = 0
        while row < vectors.rows do
          entry += vectors(row, i) * metricDiagonal(row) * vectors(row, j)
          row += 1
        val delta = if i == j then entry - 1.0 else entry
        sum += delta * delta
        j += 1
      i += 1
    math.sqrt(sum)

  test("smallest generalized eigenpairs converge with true residuals and B-orthonormal vectors") {
    val values = IndexedSeq(1.0, 2.0, 3.5, 5.0, 8.0, 13.0, 21.0, 34.0)
    val metric = IndexedSeq(1.0, 2.0, 0.5, 3.0, 4.0, 1.5, 2.5, 5.0)
    val result = solveDiagonal(
      values,
      metric,
      k = 2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(tolerance = 1e-9, maxIterations = 100)
    ).toOption.get

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assertEqualsDouble(result.eigenvalues(0), 1.0, 1e-8)
    assertEqualsDouble(result.eigenvalues(1), 2.0, 1e-8)
    assert(generalizedResidual(values, metric, result, 0) <= 1e-9)
    assert(generalizedResidual(values, metric, result, 1) <= 1e-9)
    assert(metricOrthogonalityError(metric, result.eigenvectors) < 1e-10)
  }

  test("largest algebraic selection is returned in ascending layout") {
    val values = IndexedSeq(-4.0, -1.0, 0.5, 2.0, 7.0, 11.0)
    val metric = IndexedSeq(1.0, 3.0, 2.0, 0.5, 4.0, 1.5)
    val result = solveDiagonal(
      values,
      metric,
      k = 2,
      EigenOrder.LargestAlgebraic,
      GeneralizedSpectralOptions(tolerance = 1e-9, maxIterations = 80)
    ).toOption.get

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assertEqualsDouble(result.eigenvalues(0), 7.0, 1e-8)
    assertEqualsDouble(result.eigenvalues(1), 11.0, 1e-8)
  }

  test("repeated eigenvalues are represented by a B-orthonormal invariant subspace") {
    val values = IndexedSeq(1.0, 1.0, 1.0, 4.0, 6.0, 9.0)
    val metric = IndexedSeq.fill(values.length)(1.0)
    val result = solveDiagonal(
      values,
      metric,
      k = 2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(tolerance = 1e-10, maxIterations = 60)
    ).toOption.get

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assertEqualsDouble(result.eigenvalues(0), 1.0, 1e-9)
    assertEqualsDouble(result.eigenvalues(1), 1.0, 1e-9)
    assert(metricOrthogonalityError(metric, result.eigenvectors) < 1e-10)
    var column = 0
    while column < result.size do
      var tailSquared = 0.0
      var row = 3
      while row < values.length do
        tailSquared += result.eigenvectors(row, column) * result.eigenvectors(row, column)
        row += 1
      assert(math.sqrt(tailSquared) < 1e-9, s"column $column escaped the repeated eigenspace")
      column += 1
  }

  test("dependent caller subspace is deterministically replenished") {
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0)
    val metric = IndexedSeq.fill(values.length)(1.0)
    val repeated = DVec.fromSeq(Seq(1.0, -1.0, 0.5, 2.0, -0.25))
    val initial = DMat.tabulate(values.length, 2)((row, _) => repeated(row))
    val options =
      GeneralizedSpectralOptions(tolerance = 1e-9, maxIterations = 60, initialSubspace = Some(initial))

    val first =
      solveDiagonal(values, metric, 2, EigenOrder.SmallestAlgebraic, options).toOption.get
    val second =
      solveDiagonal(values, metric, 2, EigenOrder.SmallestAlgebraic, options).toOption.get

    assert(first.diagnostics.allConverged, first.diagnostics.toString)
    assertEquals(first.eigenvalues.toSeq, second.eigenvalues.toSeq)
    assertEquals(first.eigenvectors.valuesRowMajor, second.eigenvectors.valuesRowMajor)
  }

  test("soft locking does not precondition an already converged Ritz pair") {
    var applications = 0
    val countingIdentity = new Preconditioner:
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        applications += 1
        var i = 0
        while i < r.length do
          into(i) = r(i)
          i += 1

    val initial = Matrix.dense(5, 2)(
      1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0
    )
    val result = solveDiagonal(
      IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0),
      IndexedSeq.fill(5)(1.0),
      2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(tolerance = 1e-12, maxIterations = 1, initialSubspace = Some(initial)),
      countingIdentity
    ).toOption.get

    assertEquals(applications, 1)
    assertEquals(result.diagnostics.iterations, 1)
  }

  test("zero preconditioner directions are replenished instead of causing breakdown") {
    val zero = new Preconditioner:
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        var i = 0
        while i < into.length do
          into(i) = 0.0
          i += 1

    val result = solveDiagonal(
      IndexedSeq(1.0, 2.0, 3.0, 4.0, 5.0),
      IndexedSeq.fill(5)(1.0),
      1,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(tolerance = 0.0, maxIterations = 1),
      zero
    ).toOption.get

    assertEquals(result.diagnostics.iterations, 1)
  }

  test("iteration exhaustion returns an honest zero-convergence Right") {
    val initial = Matrix.dense(4, 1)(1.0, 0.0, 1.0, 0.0)
    val result = solveDiagonal(
      IndexedSeq(1.0, 2.0, 4.0, 8.0),
      IndexedSeq.fill(4)(1.0),
      1,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(tolerance = 0.0, maxIterations = 0, initialSubspace = Some(initial))
    ).toOption.get

    assertEquals(result.size, 0)
    assertEquals(result.diagnostics.requested, 1)
    assertEquals(result.diagnostics.converged, 0)
    assertEquals(result.diagnostics.iterations, 0)
    assert(!result.diagnostics.allConverged)
  }

  test("indefinite metric geometry and illegal requests are typed Left values") {
    val indefinite = Lobpcg.solve(
      diagonalOperator(IndexedSeq(1.0, 2.0, 3.0)),
      diagonalOperator(IndexedSeq(1.0, -1.0, 2.0)),
      3,
      2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(
        initialSubspace = Some(DMat.tabulate(3, 2)((row, col) => if row == col then 1.0 else 0.0))
      ),
      Preconditioner.Identity
    )
    assert(indefinite.left.exists(_.isInstanceOf[LinAlgError.NotPositiveDefinite]))

    val illegalOrder = solveDiagonal(
      IndexedSeq(1.0, 2.0, 3.0),
      IndexedSeq.fill(3)(1.0),
      1,
      EigenOrder.SmallestMagnitude,
      GeneralizedSpectralOptions()
    )
    assert(illegalOrder.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
  }
