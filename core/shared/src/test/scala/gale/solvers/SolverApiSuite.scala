package gale.solvers

import gale.linalg.*
import gale.sparse.*

class SolverApiSuite extends munit.FunSuite:
  test("SolverResult exposes convergence and orThrow") {
    val x = Vec(1.0, 2.0)
    val ok = SolverResult.Converged(x, iterations = 3, residual = 1e-12)
    val bad = SolverResult.NotConverged(x, iterations = 10, residual = 1.0)

    assert(ok.converged)
    assertEquals(ok.orThrow.toSeq, Seq(1.0, 2.0))
    intercept[LinAlgError.DidNotConverge] {
      bad.orThrow
    }
  }

  test("Identity and Jacobi preconditioners solve diagonal approximations") {
    val r = Vec(2.0, 6.0)
    val identityOut = Preconditioner.Identity(r)
    val A = Sparse.diagonal(2.0, 3.0)
    val jacobiOut = Preconditioner.Jacobi(A)(r)

    assertEquals(identityOut.toSeq, Seq(2.0, 6.0))
    assertEquals(jacobiOut.toSeq, Seq(1.0, 2.0))
  }

  test("BlockJacobi and SymmetricGaussSeidel constructors are usable") {
    val A = Matrix.dense(2, 2)(
      4.0, 1.0,
      1.0, 3.0
    )
    val r = Vec(1.0, 2.0)

    assertEquals(Preconditioner.BlockJacobi(A, blockSize = 1)(r).length, 2)
    assertEquals(Preconditioner.SymmetricGaussSeidel(A)(r).length, 2)
  }

  test("BlockJacobi solves each diagonal block exactly and ignores off-block entries") {
    // Blocks [0,1], [2,3], [4]. Off-block entries are set to 9 to prove the
    // preconditioner reads only the diagonal blocks.
    val A = Matrix.dense(5, 5)(
      2.0, 1.0, 9.0, 9.0, 9.0,
      1.0, 2.0, 9.0, 9.0, 9.0,
      9.0, 9.0, 3.0, 1.0, 9.0,
      9.0, 9.0, 1.0, 3.0, 9.0,
      9.0, 9.0, 9.0, 9.0, 5.0
    )
    val r = Vec(1.0, 1.0, 1.0, 1.0, 1.0)
    val out = Preconditioner.BlockJacobi(A, blockSize = 2)(r).toSeq

    // [[2,1],[1,2]]^{-1} (1,1) = (1/3, 1/3); [[3,1],[1,3]]^{-1} (1,1) = (1/4, 1/4);
    // [5]^{-1} 1 = 1/5.
    val expected = Array(1.0 / 3, 1.0 / 3, 0.25, 0.25, 0.2)
    var i = 0
    while i < expected.length do
      assert(math.abs(out(i) - expected(i)) < 1e-12, s"index $i: ${out(i)} != ${expected(i)}")
      i += 1
  }

  test("BlockJacobi with blockSize >= n is an exact solve of the whole matrix") {
    val A = Matrix.dense(3, 3)(
      4.0, 1.0, 0.0,
      1.0, 3.0, 1.0,
      0.0, 1.0, 2.0
    )
    val r = Vec(1.0, 2.0, 3.0)
    val out = Preconditioner.BlockJacobi(A, blockSize = 5)(r).toSeq
    val exact = A.solve(r).orThrow.toSeq

    var i = 0
    while i < exact.length do
      assert(math.abs(out(i) - exact(i)) < 1e-12, s"index $i: ${out(i)} != ${exact(i)}")
      i += 1
  }

  test("BlockJacobi preconditioning beats point Jacobi on a block-tridiagonal SPD system") {
    // Strong within-block coupling (9) versus a diagonal of 10 makes each 2x2
    // block far from diagonal; block Jacobi inverts those blocks exactly while
    // point Jacobi cannot, so CG must converge in strictly fewer iterations.
    val n = 20
    val A = Matrix.tabulate(n, n) { (i, j) =>
      if i == j then 10.0
      else if i / 2 == j / 2 then 9.0 // same 2x2 diagonal block
      else if math.abs(i - j) == 2 then -0.4 // adjacent-block coupling
      else 0.0
    }
    val truth = Vec.tabulate(n)(i => (i + 1).toDouble)
    val b = A * truth
    val config = SolverConfig(tolerance = 1e-10, maxIterations = 500)

    val jacobi = cg(A, b, config, Preconditioner.Jacobi(A))
    val block = cg(A, b, config, Preconditioner.BlockJacobi(A, blockSize = 2))

    assert(jacobi.converged, s"jacobi residual=${jacobi.residual}")
    assert(block.converged, s"block residual=${block.residual}")
    assert(
      block.iterations < jacobi.iterations,
      s"block=${block.iterations} jacobi=${jacobi.iterations}"
    )
  }
