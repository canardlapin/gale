package gale.spectral

import gale.TestAccess
import gale.linalg.*

class DenseSymmetricWorkspaceSuite extends munit.FunSuite:
  private def randomSymmetric(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    val raw = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        val value = rng.nextDouble() * 2.0 - 1.0
        raw(i)(j) = value
        raw(j)(i) = value
        j += 1
      i += 1
    Matrix.tabulate(n, n)((r, c) => raw(r)(c))

  private def projector(vectors: DMat): DMat =
    vectors * vectors.t

  private def assertMatrixApprox(actual: DMat, expected: DMat, tolerance: Double = 1e-11): Unit =
    assertEquals(actual.shape, expected.shape)
    var r = 0
    while r < actual.rows do
      var c = 0
      while c < actual.cols do
        assert(
          math.abs(actual(r, c) - expected(r, c)) <= tolerance,
          s"($r,$c): ${actual(r, c)} != ${expected(r, c)}"
        )
        c += 1
      r += 1

  test("symmetric eigen reports checked exact values-only and vector requirements") {
    assertEquals(
      Eigen.symmetricScratchRequirement(7, EigenVectors.ValuesOnly).map(_.doubleElements),
      Right(56)
    )
    assertEquals(
      Eigen.symmetricScratchRequirement(7, EigenVectors.Right).map(_.doubleElements),
      Right(7)
    )
    assertEquals(Eigen.symmetricScratchRequirement(0, EigenVectors.ValuesOnly).map(_.doubleElements), Right(0))
    assert(Eigen.symmetricScratchRequirement(-1, EigenVectors.ValuesOnly).isLeft)
    assert(Eigen.symmetricScratchRequirement(50000, EigenVectors.ValuesOnly).isLeft)
    assert(Eigen.symmetricScratchRequirement(7, EigenVectors.Left).isLeft)
  }

  test("values-only workspace route is deterministic, exact-sized, and reuses scratch") {
    val a = randomSymmetric(9, 2026072001L)
    val ordinary = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    val requirement = Eigen.symmetricScratchRequirement(9, EigenVectors.ValuesOnly).toOption.get
    val workspace = DenseWorkspace.empty

    val first = Eigen.eigSymmetricWith(a, EigenSelection.All, EigenVectors.ValuesOnly, workspace).toOption.get
    assertEquals(workspace.doubleCapacity, requirement.doubleElements)
    assertEquals(workspace.indexCapacity, 0)
    assertEquals(first.eigenvalues.toSeq, ordinary.eigenvalues.toSeq)
    val firstValues = first.eigenvalues.toSeq
    val backing = TestAccess.workBacking(workspace)

    val second = Eigen.eigSymmetricWith(a, EigenSelection.All, EigenVectors.ValuesOnly, workspace).toOption.get
    assert(TestAccess.sameStorage(backing, TestAccess.workBacking(workspace)))
    assertEquals(second.eigenvalues.toSeq, firstValues)
    assertEquals(first.eigenvalues.toSeq, firstValues)
  }

  test("vector workspace route preserves the repeated-eigenspace projector and owns results") {
    val diagonal = Array(5.0, 5.0, 3.0, 2.0, 1.0)
    val a = Matrix.tabulate(5, 5)((r, c) => if r == c then diagonal(r) else 0.0)
    val workspace = DenseWorkspace.empty
    val first = Eigen.eigSymmetricWith(a, EigenSelection.All, workspace).toOption.get
    assertEquals(workspace.doubleCapacity, 5)
    assert(first.diagnostics.worstResidual < 1e-12)
    assert(first.diagnostics.orthogonalityError < 1e-12)
    val values = first.eigenvalues.toSeq
    val firstProjector = projector(first.eigenvectors)
    val backing = TestAccess.workBacking(workspace)

    val second = Eigen.eigSymmetricWith(a, EigenSelection.All, workspace).toOption.get
    assert(TestAccess.sameStorage(backing, TestAccess.workBacking(workspace)))
    assertEquals(second.eigenvalues.toSeq, values)
    assertEquals(first.eigenvalues.toSeq, values)
    assertMatrixApprox(projector(second.eigenvectors), firstProjector)
    assertMatrixApprox(projector(first.eigenvectors), firstProjector)
  }

  test("workspace facade preserves dense validation as typed errors") {
    val workspace = DenseWorkspace.empty
    val rectangular = Matrix.zeros(2, 3)
    assert(Eigen.eigSymmetricWith(rectangular, EigenSelection.All, workspace).isLeft)
    val square = Matrix.eye(3)
    assert(
      Eigen
        .eigSymmetricWith(square, EigenSelection.Count(0, EigenOrder.LargestAlgebraic), workspace)
        .isLeft
    )
    assert(Eigen.eigSymmetricWith(square, EigenSelection.All, EigenVectors.LeftAndRight, workspace).isLeft)
    assertEquals(workspace.doubleCapacity, 0)
  }
