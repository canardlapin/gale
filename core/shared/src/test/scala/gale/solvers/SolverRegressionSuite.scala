package gale.solvers

import gale.linalg.*

class SolverRegressionSuite extends munit.FunSuite:
  private def spd2: DMat =
    Matrix.dense(2, 2)(
      2.0, 0.0,
      0.0, 3.0
    )

  // Item 8: a non-empty initial guess whose length does not match the system
  // must be rejected, not silently replaced by a zero vector.
  test("cg rejects an initial guess of the wrong length") {
    intercept[LinAlgError.VectorLengthMismatch] {
      cg(spd2, Vec(1.0, 1.0), SolverConfig(), Preconditioner.Identity, Some(Vec(1.0, 2.0, 3.0)))
    }
  }

  test("bicgstab rejects an initial guess of the wrong length") {
    intercept[LinAlgError.VectorLengthMismatch] {
      bicgstab(spd2, Vec(1.0, 1.0), SolverConfig(), initial = Some(Vec(1.0, 2.0, 3.0)))
    }
  }

  test("gmres rejects an initial guess of the wrong length") {
    intercept[LinAlgError.VectorLengthMismatch] {
      gmres(spd2, Vec(1.0, 1.0), SolverConfig(), initial = Some(Vec(1.0, 2.0, 3.0)))
    }
  }

  // M8a: a length-0 initial for a nonzero-dimension system is a dimension error
  // like any other mismatch — it is no longer a "no guess" sentinel.
  test("a length-0 initial for a nonzero-dimension system is a dimension error") {
    intercept[LinAlgError.VectorLengthMismatch] {
      cg(spd2, Vec(1.0, 1.0), SolverConfig(), Preconditioner.Identity, Some(DVec.zeros(0)))
    }
  }

  test("a correctly sized initial guess is still accepted") {
    val truth = Vec(0.5, 1.0 / 3.0)
    val result =
      cg(spd2, Vec(1.0, 1.0), SolverConfig(tolerance = 1e-12, maxIterations = 50), Preconditioner.Identity, Some(Vec(0.4, 0.3)))
    assert(result.converged)
    assert(norm(result.x - truth) < 1e-9)
  }

  // Item 10: a zero diagonal makes Symmetric Gauss-Seidel divide by zero. Reject
  // it at construction rather than producing silent NaNs during solve.
  test("SymmetricGaussSeidel rejects a zero diagonal at construction") {
    val singular = Matrix.dense(2, 2)(
      0.0, 1.0,
      1.0, 3.0
    )
    intercept[LinAlgError.SingularMatrix] {
      Preconditioner.SymmetricGaussSeidel(singular)
    }
  }

  // Item 10: bicgstab and gmres must reject a mismatched right-hand side up front
  // like cg does, rather than failing deep inside a vector subtraction.
  test("bicgstab rejects a right-hand side of the wrong length") {
    intercept[LinAlgError.DimensionMismatch] {
      bicgstab(spd2, Vec(1.0, 2.0, 3.0), SolverConfig())
    }
  }

  test("gmres rejects a right-hand side of the wrong length") {
    intercept[LinAlgError.DimensionMismatch] {
      gmres(spd2, Vec(1.0, 2.0, 3.0), SolverConfig())
    }
  }

  // Item 9: an exact eigenvector triggers happy breakdown at the first Arnoldi
  // step; GMRES must converge there and report ~1 iteration, not inflate to the
  // full inner limit.
  test("gmres converges in one step on an exact eigenvector, without inflating iterations") {
    val id = Matrix.eye(4)
    val b = Vec(1.0, 0.0, 0.0, 0.0)
    val result = gmres(id, b, SolverConfig(tolerance = 1e-10, maxIterations = 100, restart = 30))
    assert(result.converged)
    assert(result.iterations <= 2, s"iterations=${result.iterations}")
  }

  // Item 9: a defective operator whose Krylov space cannot represent the solution
  // produces zero basis vectors and rank-deficient least-squares systems. The old
  // code swallowed those failures and spun to maxIterations; the solver must
  // instead detect the breakdown and stop promptly.
  test("gmres stops promptly on a Krylov breakdown instead of spinning to maxIterations") {
    // A = [[0, 1], [0, 0]]; b = e1 lies in an invariant subspace that cannot
    // represent the true solution, so GMRES cannot converge.
    val defective =
      LinearOperator.fromFunction(2, 2): (x, into) =>
        into(0) = x(1)
        into(1) = 0.0
    val b = Vec(1.0, 0.0)
    val result = gmres(defective, b, SolverConfig(tolerance = 1e-10, maxIterations = 100, restart = 30))
    assert(!result.converged)
    assert(result.iterations < 10, s"iterations inflated to ${result.iterations}")
  }
