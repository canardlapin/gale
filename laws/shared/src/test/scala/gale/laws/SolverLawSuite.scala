package gale.laws

import gale.linalg.*
import gale.solvers.*

/** Conformance suite exercising the [[SolverLaws]] bundle against a concrete
  * symmetric-positive-definite system solved by CG.
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
