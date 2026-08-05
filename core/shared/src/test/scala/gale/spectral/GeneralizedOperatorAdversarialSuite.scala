package gale.spectral

import gale.linalg.*
import gale.solvers.Preconditioner
import scala.collection.mutable.ArrayBuffer

class GeneralizedOperatorAdversarialSuite extends munit.FunSuite:

  private final class CountingDiagonal(
      diagonal: IndexedSeq[Double],
      recordInputs: Boolean = false
  ) extends DoubleLinearOperator:
    var applications = 0
    val inputs = ArrayBuffer.empty[DVec]
    def rows: Int = diagonal.length
    def cols: Int = diagonal.length
    def applyTo(x: DVec, into: MutableDVec): Unit =
      applications += 1
      if recordInputs then inputs += x.copy
      var i = 0
      while i < diagonal.length do
        into(i) = diagonal(i) * x(i)
        i += 1

  private final class CountingIdentityPreconditioner extends Preconditioner:
    var applications = 0
    def solve(r: DVec, into: MutableVec[Double]): Unit =
      applications += 1
      var i = 0
      while i < r.length do
        into(i) = r(i)
        i += 1

  private def diagonalOperator(diagonal: IndexedSeq[Double]): DoubleLinearOperator =
    new CountingDiagonal(diagonal)

  private def solveDiagonal(
      values: IndexedSeq[Double],
      metricDiagonal: IndexedSeq[Double],
      k: Int,
      order: EigenOrder,
      options: GeneralizedSpectralOptions
  ): Either[LinAlgError, EigenDecomposition] =
    val operatorDiagonal = values.zip(metricDiagonal).map:
      case (value, weight) => value * weight
    Eigen.eigSymmetricGeneralized(
      diagonalOperator(operatorDiagonal).assumeSymmetricOperator,
      diagonalOperator(metricDiagonal).assumePositiveDefiniteOperator,
      values.length,
      EigenSelection.Count(k, order),
      options
    )

  private def randomMatrix(n: Int, seed: Long): DMat =
    val random = new scala.util.Random(seed)
    DMat.tabulate(n, n)((_, _) => random.nextDouble() * 2.0 - 1.0)

  private def randomSymmetric(n: Int, seed: Long): DMat =
    val raw = randomMatrix(n, seed)
    DMat.tabulate(n, n)((row, col) => 0.5 * (raw(row, col) + raw(col, row)))

  private def randomSpd(n: Int, seed: Long): DMat =
    val raw = randomMatrix(n, seed)
    val gram = raw * raw.t
    DMat.tabulate(n, n): (row, col) =>
      gram(row, col) + (if row == col then n.toDouble else 0.0)

  private def expectedCoordinateBasis(
      metricDiagonal: IndexedSeq[Double],
      columns: Int
  ): DMat =
    DMat.tabulate(metricDiagonal.length, columns): (row, col) =>
      if row == col then 1.0 / math.sqrt(metricDiagonal(row)) else 0.0

  private def frobenius(matrix: DMat): Double =
    var sum = 0.0
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        sum += matrix(row, col) * matrix(row, col)
        col += 1
      row += 1
    math.sqrt(sum)

  private def metricProjector(
      vectors: DMat,
      metricDiagonal: IndexedSeq[Double]
  ): DMat =
    val bx = DMat.tabulate(vectors.rows, vectors.cols): (row, col) =>
      metricDiagonal(row) * vectors(row, col)
    vectors * bx.t

  test("matrix-free LOBPCG differentially agrees with the dense generalized solver") {
    val n = 9
    val a = randomSymmetric(n, 42L)
    val b = randomSpd(n, 99L)
    val iterative = Eigen
      .eigSymmetricGeneralized(
        a.assumeSymmetricOperator,
        b.assumePositiveDefiniteOperator,
        n,
        EigenSelection.Count(3, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(tolerance = 1e-8, maxIterations = 200)
      )
      .toOption
      .get
    val dense = Eigen
      .eigSymmetricGeneralized(
        a,
        b,
        EigenSelection.Count(3, EigenOrder.SmallestAlgebraic)
      )
      .toOption
      .get

    assert(iterative.diagnostics.allConverged, iterative.diagnostics.toString)
    var i = 0
    while i < dense.size do
      assertEqualsDouble(iterative.eigenvalues(i), dense.eigenvalues(i), 1e-7)
      i += 1
  }

  test("repeated and tightly clustered roots are checked through metric projectors") {
    val metric = IndexedSeq(0.5, 2.0, 3.0, 1.0, 4.0, 1.5, 2.5)
    val cases = Seq(
      IndexedSeq(1.0, 1.0, 1.0, 4.0, 6.0, 8.0, 10.0),
      IndexedSeq(1.0, 1.0 + 1e-9, 1.0 + 2e-9, 4.0, 6.0, 8.0, 10.0)
    )
    cases.foreach: values =>
      val result = solveDiagonal(
        values,
        metric,
        3,
        EigenOrder.SmallestAlgebraic,
        GeneralizedSpectralOptions(tolerance = 1e-10, maxIterations = 120)
      ).toOption.get
      val expected = expectedCoordinateBasis(metric, columns = 3)
      val error =
        frobenius(
          metricProjector(result.eigenvectors, metric) -
            metricProjector(expected, metric)
        )

      assert(result.diagnostics.allConverged, result.diagnostics.toString)
      assert(error < 1e-8, s"projector error $error for $values")
  }

  test("strongly ill-conditioned diagonal SPD metric preserves the analytic pencil") {
    val values = IndexedSeq(1.0, 2.0, 3.0, 5.0, 8.0, 13.0)
    val metric =
      IndexedSeq(1e-12, 1e-8, 1e-4, 1.0, 1e4, 1e8)
    val initial =
      DMat.tabulate(6, 2)((row, col) => if row == col then 1.0 else 0.0)
    val result = solveDiagonal(
      values,
      metric,
      2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(
        tolerance = 1e-8,
        maxIterations = 20,
        initialSubspace = Some(initial)
      )
    ).toOption.get

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assertEqualsDouble(result.eigenvalues(0), 1.0, 1e-8)
    assertEqualsDouble(result.eigenvalues(1), 2.0, 1e-8)
    assert(result.diagnostics.orthogonalityError < 1e-8)
  }

  test("partial and zero convergence contain exactly the residual-passing pairs") {
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
    val metric = IndexedSeq.fill(values.length)(1.0)
    val partialInitial = Matrix.dense(6, 2)(
      1.0, 0.0,
      0.0, 0.0,
      0.0, 1.0,
      0.0, 1.0,
      0.0, 0.0,
      0.0, 0.0
    )
    val partial = solveDiagonal(
      values,
      metric,
      2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(
        tolerance = 1e-12,
        maxIterations = 0,
        initialSubspace = Some(partialInitial)
      )
    ).toOption.get
    val zero = solveDiagonal(
      values,
      metric,
      2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(tolerance = 0.0, maxIterations = 1)
    ).toOption.get

    assertEquals(partial.size, 1)
    assertEquals(partial.diagnostics.converged, partial.size)
    assertEqualsDouble(partial.eigenvalues(0), 1.0, 1e-12)
    assertEquals(zero.diagnostics.converged, zero.size)
    assert(zero.size < zero.diagnostics.requested)
    assert(zero.diagnostics.residuals.toSeq.forall(_ == 0.0))
  }

  test("indefinite and non-finite operator geometry are typed failures") {
    val initial = DMat.tabulate(4, 2)((row, col) => if row == col then 1.0 else 0.0)
    val indefinite = solveDiagonal(
      IndexedSeq(1.0, 2.0, 3.0, 4.0),
      IndexedSeq(1.0, -1.0, 2.0, 3.0),
      2,
      EigenOrder.SmallestAlgebraic,
      GeneralizedSpectralOptions(initialSubspace = Some(initial))
    )
    assert(indefinite.left.exists(_.isInstanceOf[LinAlgError.NotPositiveDefinite]))

    val nonFiniteA = LinearOperator.fromFunction(4, 4): (_, into) =>
      var i = 0
      while i < 4 do
        into(i) = Double.NaN
        i += 1
    val nonFiniteB = LinearOperator.fromFunction(4, 4): (_, into) =>
      var i = 0
      while i < 4 do
        into(i) = Double.PositiveInfinity
        i += 1
    val identity = diagonalOperator(IndexedSeq.fill(4)(1.0))
    val options =
      GeneralizedSpectralOptions(initialSubspace = Some(initial))

    val badA = Eigen.eigSymmetricGeneralized(
      nonFiniteA.assumeSymmetricOperator,
      identity.assumePositiveDefiniteOperator,
      4,
      EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
      options
    )
    val badB = Eigen.eigSymmetricGeneralized(
      identity.assumeSymmetricOperator,
      nonFiniteB.assumePositiveDefiniteOperator,
      4,
      EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
      options
    )
    assert(badA.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    assert(badB.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
  }

  test("pure zero-iteration work uses exactly cached AX and BX with no preconditioner") {
    val a = new CountingDiagonal(IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0))
    val b = new CountingDiagonal(IndexedSeq.fill(5)(1.0))
    val preconditioner = new CountingIdentityPreconditioner
    val initial = Matrix.dense(5, 2)(
      1.0, 0.0,
      1.0, 1.0,
      0.0, 1.0,
      1.0, 0.0,
      0.0, 1.0
    )

    val result = Eigen
      .eigSymmetricGeneralized(
        a.assumeSymmetricOperator,
        b.assumePositiveDefiniteOperator,
        5,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(
          tolerance = 0.0,
          maxIterations = 0,
          initialSubspace = Some(initial)
        ),
        preconditioner
      )
      .toOption
      .get

    assertEquals(a.applications, 2)
    assertEquals(b.applications, 2)
    assertEquals(preconditioner.applications, 0)
    assertEquals(result.diagnostics.iterations, 0)
  }

  test("one soft-locking iteration accounts independently for A, B, and preconditioning") {
    val a = new CountingDiagonal(IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0))
    val b = new CountingDiagonal(IndexedSeq.fill(5)(1.0))
    val preconditioner = new CountingIdentityPreconditioner
    val initial = Matrix.dense(5, 2)(
      1.0, 0.0,
      0.0, 0.0,
      0.0, 1.0,
      0.0, 1.0,
      0.0, 0.0
    )

    val _ = Eigen
      .eigSymmetricGeneralized(
        a.assumeSymmetricOperator,
        b.assumePositiveDefiniteOperator,
        5,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(
          tolerance = 1e-12,
          maxIterations = 1,
          initialSubspace = Some(initial)
        ),
        preconditioner
      )
      .toOption
      .get

    assertEquals(a.applications, 3)
    assertEquals(b.applications, 3)
    assertEquals(preconditioner.applications, 1)
  }

  test("matrix-free solve never materializes the operator through a canonical basis sweep") {
    val n = 12
    val a =
      new CountingDiagonal(
        IndexedSeq.tabulate(n)(i => (i + 1).toDouble),
        recordInputs = true
      )
    val b = new CountingDiagonal(IndexedSeq.fill(n)(1.0), recordInputs = true)
    val result = Eigen
      .eigSymmetricGeneralized(
        a.assumeSymmetricOperator,
        b.assumePositiveDefiniteOperator,
        n,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedSpectralOptions(tolerance = 1e-8, maxIterations = 100)
      )
      .toOption
      .get

    def canonicalIndex(vector: DVec): Option[Int] =
      var found = -1
      var i = 0
      while i < vector.length do
        if vector(i) == 1.0 then
          if found >= 0 then return None
          found = i
        else if vector(i) != 0.0 then return None
        i += 1
      Option.when(found >= 0)(found)

    val canonicalA = a.inputs.flatMap(canonicalIndex).toSet
    val canonicalB = b.inputs.flatMap(canonicalIndex).toSet
    assert(result.diagnostics.allConverged)
    assert(canonicalA.size < n, "A was materialized by a full coordinate sweep")
    assert(canonicalB.size < n, "B was materialized by a full coordinate sweep")
  }
