package gale.solvers

import gale.linalg.*
import gale.sparse.*

class IterativeSolversSuite extends munit.FunSuite:
  test("CG solves a sparse symmetric positive-definite system") {
    val A =
      Sparse
        .coo(4, 4)
        .add(0, 0, 4.0)
        .add(0, 1, -1.0)
        .add(1, 0, -1.0)
        .add(1, 1, 4.0)
        .add(1, 2, -1.0)
        .add(2, 1, -1.0)
        .add(2, 2, 4.0)
        .add(2, 3, -1.0)
        .add(3, 2, -1.0)
        .add(3, 3, 3.0)
        .toCSR()
    val truth = Vec(1.0, 2.0, 3.0, 4.0)
    val b = A * truth
    val result = cg(A, b, SolverConfig(tolerance = 1e-12, maxIterations = 100), Preconditioner.Jacobi(A))

    assert(result.converged)
    assert(norm(result.x - truth) < 1e-9)
  }

  test("BiCGSTAB and GMRES solve a nonsymmetric dense system") {
    val A = Matrix.dense(2, 2)(
      4.0, 1.0,
      2.0, 3.0
    )
    val truth = Vec(1.0, -2.0)
    val b = A * truth
    val config = SolverConfig(tolerance = 1e-12, maxIterations = 50, restart = 4)

    val bi = bicgstab(A, b, config)
    val gm = gmres(A, b, config)

    assert(bi.converged)
    assert(gm.converged)
    assert(norm(bi.x - truth) < 1e-9)
    assert(norm(gm.x - truth) < 1e-9)
  }

  test("CGNR normal-equation solve recovers full-rank least-squares coefficients") {
    val A = Matrix.dense(3, 2)(
      1.0, 0.0,
      1.0, 1.0,
      1.0, 2.0
    )
    val truth = Vec(1.0, 2.0)
    val b = A * truth
    val result = cgnr(A, b, SolverConfig(tolerance = 1e-12, maxIterations = 100))

    assert(result.converged)
    assert(norm(result.x - truth) < 1e-9)
  }

  test("CG returns nonconvergence diagnostics when iteration budget is too small") {
    val A = Matrix.dense(2, 2)(
      2.0, 0.0,
      0.0, 3.0
    )
    val b = Vec(1.0, 1.0)
    val result = cg(A, b, SolverConfig(tolerance = 1e-30, maxIterations = 0))

    assert(!result.converged)
    intercept[LinAlgError.DidNotConverge] {
      result.orThrow
    }
  }
