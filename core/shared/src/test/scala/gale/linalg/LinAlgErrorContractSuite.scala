package gale.linalg

class LinAlgErrorContractSuite extends munit.FunSuite:
  test("shape and rank errors retain typed payloads and stable context"):
    val expected = Matrix.zeros(2, 3).shape
    val actual = Matrix.zeros(4, 1).shape
    val dimension = LinAlgError.DimensionMismatch(expected, actual)
    assertEquals(dimension.expected, expected)
    assertEquals(dimension.actual, actual)
    assertEquals(dimension.getMessage, "dimension mismatch: expected 2x3, got 4x1")

    val rank = LinAlgError.RankDeficient(2, 5)
    assertEquals(rank.rank, 2)
    assertEquals(rank.cols, 5)
    assertEquals(rank.getMessage, "matrix is rank deficient: rank 2 for 5 columns")

  test("convergence errors distinguish outer, inner, and certification context"):
    val outer = LinAlgError.DidNotConverge(12, 3.5e-6)
    assertEquals(outer.iterations, 12)
    assertEquals(outer.residual, 3.5e-6)
    assert(outer.getMessage.startsWith("solver did not converge after 12 iterations; residual="))

    val inner = LinAlgError.InnerSolveDidNotConverge(
      outerIteration = 7,
      completedSolves = 3,
      innerIterations = 42L,
      residual = 1.25e-4,
      operatorApplications = 51L
    )
    assertEquals(inner.outerIteration, 7)
    assertEquals(inner.completedSolves, 3)
    assertEquals(inner.innerIterations, 42L)
    assertEquals(inner.operatorApplications, 51L)
    assert(inner.getMessage.contains("outer iteration 7"))
    assert(inner.getMessage.contains("42 inner iterations"))
    assert(inner.getMessage.contains("operatorApplications=51"))

    val uncertified = LinAlgError.SpectralExtremeNotCertified(19, 8.0e-9)
    assertEquals(uncertified.iterations, 19)
    assertEquals(uncertified.residual, 8.0e-9)
    assert(uncertified.getMessage.contains("membership in the requested extreme is not certified"))

  test("nested failures retain the original typed cause"):
    val cause = LinAlgError.NotPositiveDefinite(4)
    val failure = LinAlgError.InnerSolveFailed(
      outerIteration = 2,
      completedSolves = 1,
      innerIterations = 9L,
      operatorApplications = 11L,
      failure = cause
    )
    assert(failure.failure eq cause)
    assertEquals(failure.outerIteration, 2)
    assertEquals(failure.getMessage, "inner solve failed at outer iteration 2 after 1 completed solves and 9 inner iterations; operatorApplications=11: matrix is not positive definite at leading minor 4")

  test("free-form boundary errors preserve their supplied diagnostic"):
    val invalid = LinAlgError.InvalidArgument("tolerance must be finite and non-negative")
    val unsupported = LinAlgError.UnsupportedRepresentation("complex sparse values are not supported")
    assertEquals(invalid.getMessage, "tolerance must be finite and non-negative")
    assertEquals(unsupported.getMessage, "complex sparse values are not supported")
