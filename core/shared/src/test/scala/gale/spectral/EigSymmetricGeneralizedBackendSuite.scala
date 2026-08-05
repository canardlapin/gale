package gale.spectral

import gale.linalg.*
import gale.solvers.Preconditioner

class EigSymmetricGeneralizedBackendSuite extends munit.FunSuite:

  private val n = 6
  private val generalizedValues = IndexedSeq(1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
  private val metricDiagonal = IndexedSeq(1.0, 2.0, 0.5, 3.0, 4.0, 1.5)
  private val operatorDiagonal =
    generalizedValues
      .zip(metricDiagonal)
      .map:
        case (value, weight) => value * weight

  private def diagonalOperator(diagonal: IndexedSeq[Double]): DoubleLinearOperator =
    LinearOperator.fromFunction(diagonal.length, diagonal.length): (x, into) =>
      var i = 0
      while i < diagonal.length do
        into(i) = diagonal(i) * x(i)
        i += 1

  private val operator = diagonalOperator(operatorDiagonal)
  private val metric = diagonalOperator(metricDiagonal)
  private val selection =
    EigenSelection.Count(2, EigenOrder.SmallestAlgebraic)
  private val options =
    GeneralizedSpectralOptions(tolerance = 1e-10, maxIterations = 80)

  private def publicSolve(
      backend: SpectralBackend,
      requestOptions: GeneralizedSpectralOptions = options
  ): Either[LinAlgError, EigenDecomposition] =
    Eigen.eigSymmetricGeneralized(
      operator.assumeSymmetricOperator,
      metric.assumePositiveDefiniteOperator,
      n,
      selection,
      requestOptions
    )(using backend)

  private final class ReversingProvider extends SpectralBackend:
    var calls = 0
    def name: String = "reversing-generalized"
    def capabilities: Set[SpectralCapability] =
      Set(SpectralCapability.IterativeGeneralized)
    override def iterativeGeneralizedEigen(
        a: DoubleLinearOperator,
        b: DoubleLinearOperator,
        dimension: Int,
        k: Int,
        order: EigenOrder,
        requestOptions: GeneralizedSpectralOptions,
        preconditioner: Preconditioner
    ): Either[LinAlgError, RawIterativeGeneralizedEigen] =
      calls += 1
      Lobpcg
        .solve(
          a,
          b,
          dimension,
          k,
          order,
          requestOptions.copy(returnVectors = EigenVectors.Right),
          preconditioner
        )
        .map: result =>
          val reversed = (0 until result.size).reverse.toArray
          RawIterativeGeneralizedEigen(
            DVec.tabulate(result.size)(i => result.eigenvalues(reversed(i))),
            DMat.tabulate(dimension, result.size): (row, col) =>
              result.eigenvectors(row, reversed(col)) * (col + 2).toDouble,
            BackendConvergence(k, result.size, iterations = 17)
          )

  test("facade reorders, B-normalizes, and independently derives diagnostics") {
    val provider = new ReversingProvider
    val routed = publicSolve(provider).toOption.get
    val pure = publicSolve(SpectralBackend.none).toOption.get

    assertEquals(provider.calls, 1)
    assertEquals(routed.eigenvalues.toSeq, pure.eigenvalues.toSeq)
    assertEquals(routed.diagnostics.iterations, 17)
    assert(routed.diagnostics.allConverged)
    assert(routed.diagnostics.worstResidual <= options.tolerance)
    assert(routed.diagnostics.orthogonalityError < 1e-10)
    assertEqualsDouble(routed.eigenvalues(0), 1.0, 1e-9)
    assertEqualsDouble(routed.eigenvalues(1), 2.0, 1e-9)
  }

  test("values-only routing still validates hidden provider vectors then omits them") {
    val provider = new ReversingProvider
    val result = publicSolve(
      provider,
      options.copy(returnVectors = EigenVectors.ValuesOnly)
    ).toOption.get

    assertEquals(provider.calls, 1)
    assertEquals(result.eigenvectors.rows, n)
    assertEquals(result.eigenvectors.cols, 0)
    assert(result.diagnostics.allConverged)
    assert(result.diagnostics.worstResidual <= options.tolerance)
    assertEquals(result.diagnostics.orthogonalityError, 0.0)
  }

  test("provider Left is a decline and falls back to pure LOBPCG") {
    var calls = 0
    val provider = new SpectralBackend:
      def name: String = "declining-generalized"
      def capabilities: Set[SpectralCapability] =
        Set(SpectralCapability.IterativeGeneralized)
      override def iterativeGeneralizedEigen(
          a: DoubleLinearOperator,
          b: DoubleLinearOperator,
          dimension: Int,
          k: Int,
          order: EigenOrder,
          requestOptions: GeneralizedSpectralOptions,
          preconditioner: Preconditioner
      ): Either[LinAlgError, RawIterativeGeneralizedEigen] =
        calls += 1
        Left(LinAlgError.UnsupportedOperation("deliberate decline"))

    val routed = publicSolve(provider).toOption.get
    val pure = publicSolve(SpectralBackend.none).toOption.get
    assertEquals(calls, 1)
    assertEquals(routed.eigenvalues.toSeq, pure.eigenvalues.toSeq)
    assertEquals(routed.diagnostics.residuals.toSeq, pure.diagnostics.residuals.toSeq)
  }

  test("validation precedes provider invocation") {
    val provider = new ReversingProvider
    val badDimension = Eigen.eigSymmetricGeneralized(
      operator.assumeSymmetricOperator,
      metric.assumePositiveDefiniteOperator,
      n + 1,
      selection
    )(using provider)
    val badOrder = Eigen.eigSymmetricGeneralized(
      operator.assumeSymmetricOperator,
      metric.assumePositiveDefiniteOperator,
      n,
      EigenSelection.Count(2, EigenOrder.SmallestMagnitude)
    )(using provider)
    val badInitial = Eigen.eigSymmetricGeneralized(
      operator.assumeSymmetricOperator,
      metric.assumePositiveDefiniteOperator,
      n,
      selection,
      options.copy(initialSubspace = Some(DMat.zeros(n, 1)))
    )(using provider)

    assert(badDimension.isLeft)
    assert(badOrder.isLeft)
    assert(badInitial.isLeft)
    assertEquals(provider.calls, 0)
  }

  test("IterativeGeneralized is not inferred from dense capabilities") {
    var iterativeCalls = 0
    val provider = new SpectralBackend:
      def name: String = "dense-only"
      def capabilities: Set[SpectralCapability] =
        Set(SpectralCapability.DenseSymmetricEigen)
      override def iterativeGeneralizedEigen(
          a: DoubleLinearOperator,
          b: DoubleLinearOperator,
          dimension: Int,
          k: Int,
          order: EigenOrder,
          requestOptions: GeneralizedSpectralOptions,
          preconditioner: Preconditioner
      ): Either[LinAlgError, RawIterativeGeneralizedEigen] =
        iterativeCalls += 1
        fail("dense capability must not route the iterative generalized seam")

    val result = publicSolve(provider).toOption.get
    assert(result.diagnostics.allConverged)
    assertEquals(iterativeCalls, 0)
  }

  test("zero-convergence raw result remains an honest Right") {
    val provider = rawProvider(
      RawIterativeGeneralizedEigen(
        DVec.zeros(0),
        DMat.zeros(n, 0),
        BackendConvergence(requested = 2, converged = 0, iterations = 9)
      )
    )
    val result = publicSolve(provider).toOption.get

    assertEquals(result.size, 0)
    assertEquals(result.diagnostics.requested, 2)
    assertEquals(result.diagnostics.converged, 0)
    assertEquals(result.diagnostics.iterations, 9)
    assert(!result.diagnostics.allConverged)
  }

  test("malformed raw shapes and non-finite entries fail loudly") {
    val wrongCount = rawProvider(
      RawIterativeGeneralizedEigen(
        DVec.fromSeq(Seq(1.0)),
        DMat.zeros(n, 1),
        BackendConvergence(requested = 2, converged = 2, iterations = 1)
      )
    )
    val _ = intercept[LinAlgError.InvalidArgument](publicSolve(wrongCount))

    val wrongRows = rawProvider(
      RawIterativeGeneralizedEigen(
        DVec.fromSeq(Seq(1.0, 2.0)),
        DMat.zeros(n - 1, 2),
        BackendConvergence(requested = 2, converged = 2, iterations = 1)
      )
    )
    val _ = intercept[LinAlgError.InvalidArgument](publicSolve(wrongRows))

    val nonFinite = rawProvider(
      RawIterativeGeneralizedEigen(
        DVec.fromSeq(Seq(Double.NaN, 2.0)),
        DMat.tabulate(n, 2)((row, col) => if row == col then 1.0 else 0.0),
        BackendConvergence(requested = 2, converged = 2, iterations = 1)
      )
    )
    val _ = intercept[LinAlgError.InvalidArgument](publicSolve(nonFinite))
  }

  test("false convergence and non-B-orthogonal Ritz blocks fail loudly") {
    val falseConvergence = rawProvider(
      RawIterativeGeneralizedEigen(
        DVec.fromSeq(Seq(999.0)),
        DMat.tabulate(n, 1)((row, _) => if row == 0 then 1.0 else 0.0),
        BackendConvergence(requested = 2, converged = 1, iterations = 1)
      )
    )
    val _ = intercept[LinAlgError.InvalidArgument](publicSolve(falseConvergence))

    val duplicate = DMat.tabulate(n, 2)((row, _) => if row == 0 then 1.0 else 0.0)
    val nonOrthogonal = rawProvider(
      RawIterativeGeneralizedEigen(
        DVec.fromSeq(Seq(1.0, 1.0)),
        duplicate,
        BackendConvergence(requested = 2, converged = 2, iterations = 1)
      )
    )
    val _ = intercept[LinAlgError.InvalidArgument](publicSolve(nonOrthogonal))
  }

  test("backend final diagnostics account for exactly one A and B application per returned pair") {
    var operatorApplications = 0
    var metricApplications = 0
    val countedOperator =
      LinearOperator.fromFunction(n, n): (x, into) =>
        operatorApplications += 1
        var i = 0
        while i < n do
          into(i) = operatorDiagonal(i) * x(i)
          i += 1
    val countedMetric =
      LinearOperator.fromFunction(n, n): (x, into) =>
        metricApplications += 1
        var i = 0
        while i < n do
          into(i) = metricDiagonal(i) * x(i)
          i += 1
    var preconditionerApplications = 0
    val countedPreconditioner = new Preconditioner:
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        preconditionerApplications += 1
        var i = 0
        while i < r.length do
          into(i) = r(i)
          i += 1
    val provider = rawProvider(
      RawIterativeGeneralizedEigen(
        DVec.fromSeq(Seq(2.0, 1.0)),
        DMat.tabulate(n, 2): (row, col) =>
          if col == 0 && row == 1 then 3.0
          else if col == 1 && row == 0 then 2.0
          else 0.0,
        BackendConvergence(requested = 2, converged = 2, iterations = 7)
      )
    )

    val result = Eigen
      .eigSymmetricGeneralized(
        countedOperator.assumeSymmetricOperator,
        countedMetric.assumePositiveDefiniteOperator,
        n,
        selection,
        options,
        countedPreconditioner
      )(using provider)
      .toOption
      .get

    assert(result.diagnostics.allConverged)
    assertEquals(operatorApplications, 2)
    assertEquals(metricApplications, 2)
    assertEquals(preconditionerApplications, 0)
  }

  private def rawProvider(raw: RawIterativeGeneralizedEigen): SpectralBackend =
    new SpectralBackend:
      def name: String = "raw-fixture"
      def capabilities: Set[SpectralCapability] =
        Set(SpectralCapability.IterativeGeneralized)
      override def iterativeGeneralizedEigen(
          a: DoubleLinearOperator,
          b: DoubleLinearOperator,
          dimension: Int,
          k: Int,
          order: EigenOrder,
          requestOptions: GeneralizedSpectralOptions,
          preconditioner: Preconditioner
      ): Either[LinAlgError, RawIterativeGeneralizedEigen] =
        Right(raw)
