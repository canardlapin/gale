package gale.laws

import gale.linalg.*
import gale.solvers.*
import munit.Assertions

/** Reusable laws for the iterative and least-squares solvers, expressed against
  * the public API.
  */
object SolverLaws extends Assertions:
  /** A converged result reports a residual within its own tolerance (with a
    * small relative slack for the final-step estimate).
    */
  def convergedWithinTolerance(result: SolverResult, tolerance: Double): Unit =
    assert(result.converged, s"solver did not converge (residual ${result.residual})")
    assert(
      result.residual <= tolerance * (1.0 + 1e-6) + 1e-30,
      s"residual ${result.residual} exceeds tolerance $tolerance"
    )

  /** The solution actually solves the square system: `||A x - b|| <= tolerance`. */
  def solvesSystem(a: DoubleLinearOperator, x: DVec, b: DVec, tolerance: Double): Unit =
    val ax = a * x
    val residual = (ax - b).norm2
    assert(residual <= tolerance, s"||A x - b|| = $residual exceeds tolerance $tolerance")

  /** Two solutions agree elementwise within a relative tolerance — e.g. an
    * iterative least-squares solution against a direct QR least-squares solution.
    */
  def solutionsAgree(actual: DVec, expected: DVec, rel: Double = 1e-6): Unit =
    VecLaws.assertCloseRel(actual, expected, rel)
