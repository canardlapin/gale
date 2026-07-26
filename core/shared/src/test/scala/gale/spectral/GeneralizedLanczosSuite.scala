package gale.spectral

import gale.linalg.*
import gale.solvers.SolverConfig
import gale.solvers.ToleranceMode

class GeneralizedLanczosSuite extends munit.FunSuite:

  private def diagonal(values: IndexedSeq[Double]): DMat =
    DMat.tabulate(values.length, values.length): (row, col) =>
      if row == col then values(row) else 0.0

  private def directSolve(
      a: DMat,
      b: DMat,
      k: Int,
      order: EigenOrder,
      options: GeneralizedLanczosOptions = GeneralizedLanczosOptions()
  ): Either[LinAlgError, EigenDecomposition] =
    val positiveMetric = b.assumePositiveDefinite
    val solve =
      MetricSolveOperator
        .bind(positiveMetric, LinearSolveOperator.direct(b.cholesky.orThrow))
        .orThrow
    Eigen.eigSymmetricGeneralizedLanczos(
      a.assumeSymmetric,
      solve,
      a.rows,
      EigenSelection.Count(k, order),
      options
    )

  private def randomMatrix(n: Int, seed: Long): DMat =
    val random = new scala.util.Random(seed)
    DMat.tabulate(n, n)((_, _) => random.nextDouble() * 2.0 - 1.0)

  private def randomSymmetric(n: Int, seed: Long): DMat =
    val raw = randomMatrix(n, seed)
    DMat.tabulate(n, n)((i, j) => 0.5 * (raw(i, j) + raw(j, i)))

  private def randomSpd(n: Int, seed: Long): DMat =
    val raw = randomMatrix(n, seed)
    val gram = raw * raw.t
    DMat.tabulate(n, n): (i, j) =>
      gram(i, j) + (if i == j then n.toDouble else 0.0)

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

  private def metricProjector(vectors: DMat, b: DMat): DMat =
    vectors * (b * vectors).t

  test("analytic diagonal pencil supports both algebraic ends") {
    val values = IndexedSeq(-3.0, -1.0, 0.5, 2.0, 4.0, 7.0, 11.0, 16.0)
    val metricDiagonal = IndexedSeq(0.5, 2.0, 1.5, 3.0, 4.0, 0.75, 2.5, 5.0)
    val b = diagonal(metricDiagonal)
    val a = diagonal(values.zip(metricDiagonal).map((value, weight) => value * weight))

    val smallest =
      directSolve(
        a,
        b,
        3,
        EigenOrder.SmallestAlgebraic,
        GeneralizedLanczosOptions(tolerance = 1e-10, maxIterations = 4)
      ).orThrow
    val largest =
      directSolve(
        a,
        b,
        3,
        EigenOrder.LargestAlgebraic,
        GeneralizedLanczosOptions(tolerance = 1e-10, maxIterations = 4)
      ).orThrow

    assert(smallest.diagnostics.allConverged, smallest.diagnostics.toString)
    assert(largest.diagnostics.allConverged, largest.diagnostics.toString)
    Seq(-3.0, -1.0, 0.5).zipWithIndex.foreach: (expected, i) =>
      assertEqualsDouble(smallest.eigenvalues(i), expected, 1e-9)
    Seq(7.0, 11.0, 16.0).zipWithIndex.foreach: (expected, i) =>
      assertEqualsDouble(largest.eigenvalues(i), expected, 1e-9)
    assert(smallest.diagnostics.orthogonalityError < 1e-9)
    assert(smallest.diagnostics.extremalityCertified)
    assert(smallest.diagnostics.innerSolve.exists(_.solves > 0))
  }

  test("generalized Lanczos differentially agrees with the dense generalized solver") {
    val n = 9
    val a = randomSymmetric(n, 42L)
    val b = randomSpd(n, 99L)
    val iterative =
      directSolve(
        a,
        b,
        3,
        EigenOrder.SmallestAlgebraic,
        GeneralizedLanczosOptions(tolerance = 1e-8, maxIterations = 4)
      ).orThrow
    val dense =
      Eigen
        .eigSymmetricGeneralized(
          a,
          b,
          EigenSelection.Count(3, EigenOrder.SmallestAlgebraic)
        )
        .orThrow

    assert(iterative.diagnostics.allConverged, iterative.diagnostics.toString)
    var i = 0
    while i < dense.size do
      assertEqualsDouble(iterative.eigenvalues(i), dense.eigenvalues(i), 1e-7)
      assert(iterative.diagnostics.residuals(i) <= 1e-8)
      i += 1
  }

  test("repeated roots are compared through the B-projector") {
    val values = IndexedSeq(1.0, 1.0, 1.0, 4.0, 6.0, 8.0, 10.0)
    val weights = IndexedSeq(0.5, 2.0, 3.0, 1.0, 4.0, 1.5, 2.5)
    val b = diagonal(weights)
    val a = diagonal(values.zip(weights).map((value, weight) => value * weight))
    val result =
      directSolve(
        a,
        b,
        3,
        EigenOrder.SmallestAlgebraic,
        GeneralizedLanczosOptions(tolerance = 1e-10, maxIterations = 4)
      ).orThrow
    val expected =
      DMat.tabulate(values.length, 3): (row, col) =>
        if row == col then 1.0 / math.sqrt(weights(row)) else 0.0
    val error =
      frobenius(
        metricProjector(result.eigenvectors, b) -
          metricProjector(expected, b)
      )

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assert(error < 1e-8, s"metric projector error $error")
  }

  test("rank-deficient initial blocks are replenished and thick restart converges") {
    val n = 30
    val values = IndexedSeq.tabulate(n)(i => (i + 1).toDouble)
    val b = Matrix.eye(n)
    val a = diagonal(values)
    val initial = DMat.tabulate(n, 2)((row, _) => if row == 0 then 1.0 else 0.0)
    val result =
      directSolve(
        a,
        b,
        2,
        EigenOrder.SmallestAlgebraic,
        GeneralizedLanczosOptions(
          tolerance = 1e-8,
          maxIterations = 20,
          subspaceDimension = Some(8),
          initialSubspace = Some(initial)
        )
      ).orThrow

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assertEqualsDouble(result.eigenvalues(0), 1.0, 1e-7)
    assertEqualsDouble(result.eigenvalues(1), 2.0, 1e-7)
    assert(result.diagnostics.iterations > 1)
    assert(!result.diagnostics.extremalityCertified)
  }

  test("block expansion respects a subspace cap that is not a block multiple") {
    val n = 25
    val values = IndexedSeq.tabulate(n)(i => (i + 1).toDouble)
    val result =
      directSolve(
        diagonal(values),
        Matrix.eye(n),
        4,
        EigenOrder.SmallestAlgebraic,
        GeneralizedLanczosOptions(
          tolerance = 1e-7,
          maxIterations = 12,
          subspaceDimension = Some(19)
        )
      ).orThrow

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    var i = 0
    while i < 4 do
      assertEqualsDouble(result.eigenvalues(i), (i + 1).toDouble, 1e-6)
      i += 1
  }

  test("common pencil scaling preserves the spectrum") {
    val values = IndexedSeq(1.0, 2.0, 3.5, 5.0, 8.0, 13.0, 21.0)
    val weights = IndexedSeq(0.5, 2.0, 3.0, 1.0, 4.0, 1.5, 2.5)
    val a = diagonal(values.zip(weights).map((value, weight) => value * weight))
    val b = diagonal(weights)
    val scale = 3.25
    val base =
      directSolve(a, b, 2, EigenOrder.SmallestAlgebraic).orThrow
    val scaledA =
      DMat.tabulate(a.rows, a.cols)((row, col) => scale * a(row, col))
    val scaledB =
      DMat.tabulate(b.rows, b.cols)((row, col) => scale * b(row, col))
    val scaled =
      directSolve(
        scaledA,
        scaledB,
        2,
        EigenOrder.SmallestAlgebraic
      ).orThrow

    assert(base.diagnostics.allConverged)
    assert(scaled.diagnostics.allConverged)
    var i = 0
    while i < base.size do
      assertEqualsDouble(scaled.eigenvalues(i), base.eigenvalues(i), 1e-8)
      i += 1
  }

  test("iterative metric solves expose exact aggregate inner work") {
    val values = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0)
    val weights = IndexedSeq(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
    val a = diagonal(values.zip(weights).map((value, weight) => value * weight))
    val b = diagonal(weights)
    val positiveMetric = b.assumePositiveDefinite
    val iterativeSolve =
      LinearSolveOperator
        .conjugateGradient(
          positiveMetric,
          SolverConfig(tolerance = 1e-12, maxIterations = 20),
          toleranceMode = ToleranceMode.RelativeToRhs
        )
        .orThrow
    val metricSolve =
      MetricSolveOperator.bind(positiveMetric, iterativeSolve).orThrow
    val result =
      Eigen
        .eigSymmetricGeneralizedLanczos(
          a.assumeSymmetric,
          metricSolve,
          values.length,
          EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
          GeneralizedLanczosOptions(tolerance = 1e-8, maxIterations = 4)
        )
        .orThrow
    val work = result.diagnostics.innerSolve.get

    assert(result.diagnostics.allConverged, result.diagnostics.toString)
    assert(work.solves > 0)
    assertEquals(work.converged, work.solves)
    assert(work.iterations > 0L)
    assertEquals(
      work.operatorApplications,
      work.iterations + work.solves.toLong
    )
    assert(work.worstResidualNorm.exists(_ <= 1e-10))
  }

  test("inner non-convergence is distinct from outer spectral exhaustion") {
    val n = 6
    val a = diagonal(IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0))
    val b = Matrix.eye(n)
    val positiveMetric = b.assumePositiveDefinite
    val failingSolve =
      LinearSolveOperator
        .conjugateGradient(
          positiveMetric,
          SolverConfig(tolerance = 0.0, maxIterations = 0)
        )
        .orThrow
    val metricSolve =
      MetricSolveOperator.bind(positiveMetric, failingSolve).orThrow
    val result =
      Eigen.eigSymmetricGeneralizedLanczos(
        a.assumeSymmetric,
        metricSolve,
        n,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedLanczosOptions(tolerance = 0.0, maxIterations = 3)
      )

    result match
      case Left(error: LinAlgError.InnerSolveDidNotConverge) =>
        assertEquals(error.outerIteration, 1)
        assertEquals(error.completedSolves, 1)
        assertEquals(error.innerIterations, 0L)
        assertEquals(error.operatorApplications, 1L)
        assert(error.residual > 0.0)
      case other => fail(s"expected InnerSolveDidNotConverge, got $other")
  }

  test("structural inner-solve failures retain their typed cause and outer location") {
    val n = 6
    val a = diagonal(IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0))
    val b = Matrix.eye(n)
    val positiveMetric = b.assumePositiveDefinite
    val failingSolve =
      LinearSolveOperator
        .backendProvided(n): _ =>
          Left(LinAlgError.SingularMatrix(3))
        .orThrow
    val metricSolve =
      MetricSolveOperator.bind(positiveMetric, failingSolve).orThrow
    val result =
      Eigen.eigSymmetricGeneralizedLanczos(
        a.assumeSymmetric,
        metricSolve,
        n,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedLanczosOptions(tolerance = 0.0, maxIterations = 3)
      )

    result match
      case Left(error: LinAlgError.InnerSolveFailed) =>
        assertEquals(error.outerIteration, 1)
        assertEquals(error.completedSolves, 0)
        assertEquals(error.innerIterations, 0L)
        assertEquals(error.operatorApplications, 0L)
        assert(error.failure.isInstanceOf[LinAlgError.SingularMatrix])
      case other => fail(s"expected InnerSolveFailed, got $other")
  }

  test("zero outer iterations returns residual-passing initial pairs and no inner work") {
    val n = 6
    val a = diagonal(IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0))
    val b = Matrix.eye(n)
    val initial =
      DMat.tabulate(n, 2)((row, col) => if row == col then 1.0 else 0.0)
    val result =
      directSolve(
        a,
        b,
        2,
        EigenOrder.SmallestAlgebraic,
        GeneralizedLanczosOptions(
          tolerance = 1e-12,
          maxIterations = 0,
          initialSubspace = Some(initial),
          returnVectors = EigenVectors.ValuesOnly
        )
      ).orThrow

    assert(result.diagnostics.allConverged)
    assertEquals(result.eigenvectors.cols, 0)
    assertEquals(result.diagnostics.iterations, 0)
    assertEquals(result.diagnostics.innerSolve, Some(LinearSolveSummary.Empty))
    assert(!result.diagnostics.extremalityCertified)
  }

  test("facade rejects invalid selection and options before applying the solve") {
    val a = Matrix.eye(4).assumeSymmetric
    val b = Matrix.eye(4).assumePositiveDefinite
    var solves = 0
    val rawSolve =
      LinearSolveOperator
        .backendProvided(4): rhs =>
          solves += 1
          Right(
            LinearSolveResult(
              rhs,
              LinearSolveDiagnostics(true, 0, None, 0L)
            )
          )
        .orThrow
    val metricSolve = MetricSolveOperator.bind(b, rawSolve).orThrow

    val badSelection =
      Eigen.eigSymmetricGeneralizedLanczos(
        a,
        metricSolve,
        4,
        EigenSelection.All
      )
    val badSubspace =
      Eigen.eigSymmetricGeneralizedLanczos(
        a,
        metricSolve,
        4,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedLanczosOptions(subspaceDimension = Some(2))
      )
    val badVectors =
      Eigen.eigSymmetricGeneralizedLanczos(
        a,
        metricSolve,
        4,
        EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
        GeneralizedLanczosOptions(returnVectors = EigenVectors.Left)
      )

    assert(badSelection.isLeft)
    assert(badSubspace.isLeft)
    assert(badVectors.isLeft)
    assertEquals(solves, 0)
  }
