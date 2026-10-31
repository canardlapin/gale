package gale.laws

import gale.linalg.*
import gale.solvers.*

/** Conformance suite exercising the [[SolverLaws]] bundle against concrete systems solved by CG and preconditioned
  * GMRES.
  */
class SolverLawSuite extends munit.FunSuite:
  private val a = Matrix.dense(2, 2)(4.0, 1.0, 1.0, 3.0)
  private val b = Vec(1.0, 2.0)

  test("CG converges, solves the system, and agrees with the direct solve") {
    val tol = 1e-10
    val result = cg(a, b, SolverConfig(tolerance = tol))
    SolverLaws.convergedWithinTolerance(result, tol)
    SolverLaws.solvesSystem(a, result.x, b, 1e-8)
    SolverLaws.solutionsAgree(result.x, a.solve(b).orThrow)
  }

  for scale <- Seq(1e-12, 1.0, 1e12) do
    test(s"Jacobi GMRES solves the original system at matrix scale $scale") {
      val matrix = Matrix.dense(2, 2)(4.0 * scale, scale, 0.0, 3.0 * scale)
      val tol = 1e-10
      val result = gmres(matrix, b, SolverConfig(tolerance = tol), Preconditioner.Jacobi(matrix))
      SolverLaws.convergedWithinTolerance(result, tol)
      SolverLaws.solvesSystem(matrix, result.x, b, tol)
      SolverLaws.solutionsAgree(result.x, matrix.solve(b).orThrow)
    }
