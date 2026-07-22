package gale.solvers

import gale.linalg.*

class CgWorkspaceSuite extends munit.FunSuite:
  private val matrix = Matrix.dense(4, 4)(
    5.0, -1.0, 0.2, 0.0,
    -1.0, 4.0, -0.5, 0.1,
    0.2, -0.5, 3.0, -0.4,
    0.0, 0.1, -0.4, 2.5
  )

  test("cgWith matches cg with relative tolerance and preconditioning"):
    val truth = Vec(0.5, -1.2, 2.0, 0.75)
    val rhs = matrix * truth
    val config = SolverConfig(tolerance = 1e-12, maxIterations = 50)
    val preconditioner = Preconditioner.Jacobi(matrix)

    val owning = cg(
      matrix,
      rhs,
      config,
      preconditioner,
      toleranceMode = ToleranceMode.RelativeToRhs
    )
    val workspace = CgWorkspace(4)
    val reused = cgWith(
      matrix,
      rhs,
      workspace,
      config,
      preconditioner,
      toleranceMode = ToleranceMode.RelativeToRhs
    )

    assert(owning.converged)
    assert(reused.converged)
    assertEquals(reused.iterations, owning.iterations)
    assertEqualsDouble(reused.residual, owning.residual, 1e-15)
    assert(residualNorm(matrix, reused.solution, rhs) <= 1e-11)
    assert(vectorDifference(reused.solution, truth) <= 1e-11)

  test("workspace reuse overwrites the view but not an explicit copy"):
    val workspace = CgWorkspace(4)
    val firstTruth = Vec(1.0, 2.0, 3.0, 4.0)
    val secondTruth = Vec(-2.0, 0.5, 1.5, -1.0)
    val config = SolverConfig(tolerance = 1e-12, maxIterations = 50)

    cgWith(matrix, matrix * firstTruth, workspace, config)
    val alias = workspace.solution
    val copy = workspace.solutionCopy
    assert(vectorDifference(alias, firstTruth) <= 1e-11)

    cgWith(matrix, matrix * secondTruth, workspace, config)
    assert(vectorDifference(alias, secondTruth) <= 1e-11)
    assert(vectorDifference(copy, firstTruth) <= 1e-11)
    assert(vectorDifference(copy, alias) > 1.0)

  test("exact initial solution converges without an iteration"):
    val truth = Vec(0.25, -0.5, 0.75, 1.25)
    val rhs = matrix * truth
    val workspace = CgWorkspace(4)

    cgWith(
      matrix,
      rhs,
      workspace,
      SolverConfig(tolerance = 1e-12, maxIterations = 20),
      initial = Some(truth)
    )

    assert(workspace.converged)
    assertEquals(workspace.iterations, 0)
    assert(workspace.residual <= 1e-14)
    assert(vectorDifference(workspace.solution, truth) == 0.0)

  test("workspace and initial dimensions are checked"):
    intercept[LinAlgError.VectorLengthMismatch]:
      cgWith(matrix, Vec(1.0, 2.0, 3.0, 4.0), CgWorkspace(3))

    intercept[LinAlgError.VectorLengthMismatch]:
      cgWith(
        matrix,
        Vec(1.0, 2.0, 3.0, 4.0),
        CgWorkspace(4),
        initial = Some(Vec(1.0, 2.0, 3.0))
      )

  test("zero iteration budget reports the initial residual without stale state"):
    val workspace = CgWorkspace(4)
    val rhs = Vec(1.0, -2.0, 0.5, 3.0)
    cgWith(matrix, matrix * Vec(1.0, 1.0, 1.0, 1.0), workspace)

    cgWith(matrix, rhs, workspace, SolverConfig(tolerance = 1e-30, maxIterations = 0))

    assert(!workspace.converged)
    assertEquals(workspace.iterations, 0)
    assertEqualsDouble(workspace.residual, rhs.norm2, 1e-15)
    assertEqualsDouble(workspace.solution.norm2, 0.0, 0.0)

  test("matrix-free destination application is sufficient for cgWith"):
    val diagonal = Array(1.0, 2.0, 4.0, 8.0)
    val operator = LinearOperator.fromFunction(4, 4): (input, output) =>
      var i = 0
      while i < 4 do
        output(i) = diagonal(i) * input(i)
        i += 1
    val truth = Vec(2.0, -1.0, 0.5, 3.0)
    val rhs = Vec.tabulate(4)(i => diagonal(i) * truth(i))
    val workspace = CgWorkspace(4)

    cgWith(operator, rhs, workspace, SolverConfig(tolerance = 1e-13, maxIterations = 10))

    assert(workspace.converged)
    assert(vectorDifference(workspace.solution, truth) <= 1e-12)

  private def residualNorm(operator: DoubleLinearOperator, solution: DVec, rhs: DVec): Double =
    val actual = operator * solution
    vectorDifference(actual, rhs)

  private def vectorDifference(left: DVec, right: DVec): Double =
    var sum = 0.0
    var i = 0
    while i < left.length do
      val difference = left(i) - right(i)
      sum += difference * difference
      i += 1
    math.sqrt(sum)
