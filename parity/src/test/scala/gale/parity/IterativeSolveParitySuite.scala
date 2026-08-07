package gale.parity

import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.solvers.*
import gale.sparse.*

/** Solution-equivalence for Gale iterative solvers versus Breeze dense `\`.
  *
  * Breeze 2.1 has no first-class `cg` / `bicgstab` / `gmres` / `lsqr` in
  * `breeze.linalg`, so this suite does '''not''' claim algorithm-vs-algorithm
  * parity. It checks that common migrate-and-solve workloads — sparse or dense
  * operators solved iteratively in Gale — recover the same solution (within
  * tolerance) as Breeze's dense direct solve on the same coefficient matrix and
  * right-hand side. That is a replaceability contract for migrants, not a claim
  * that the Krylov methods match netlib LAPACK path-for-path.
  *
  * Fixtures stay well-conditioned so both sides converge cleanly; ill-conditioned
  * and nonconvergent cases remain covered by Gale's core solver suites.
  */
class IterativeSolveParitySuite extends munit.FunSuite:

  private val solveTol = 1e-8
  private val config =
    SolverConfig(tolerance = 1e-12, maxIterations = 500, restart = 40)

  private def galeCsr(data: Array[Array[Double]]): CSR =
    val builder = Sparse.coo(data.length, data(0).length)
    var i = 0
    while i < data.length do
      var j = 0
      while j < data(0).length do
        if data(i)(j) != 0.0 then builder.add(i, j, data(i)(j))
        j += 1
      i += 1
    builder.toCSR()

  private def assertConvergedClose(
      result: SolverResult,
      breezeX: breeze.linalg.DenseVector[Double],
      clue: String
  ): Unit =
    assert(result.converged, s"$clue: Gale solver did not converge (iters=${result.iterations}, residual=${result.residual})")
    assertVecClose(result.x, breezeX, solveTol, clue)

  // ---------------------------------------------------------------------------
  // CG on SPD sparse systems
  // ---------------------------------------------------------------------------

  test("CG on SPD CSR matches Breeze dense backslash") {
    for n <- List(8, 16, 32); seed <- List(1L, 2L) do
      val data = spd(n, seed)
      val gA = galeCsr(data)
      val bA = breezeMatrix(data)
      val xData = vectorData(n, seed * 11 + 3)
      val truth = galeVector(xData)
      val gb = gA * truth
      val bb = bA * breezeVector(xData)

      val result = cg(gA, gb, config, Preconditioner.Jacobi(gA))
      val breezeX = bA \ bb
      assertConvergedClose(result, breezeX, s"CG SPD n=$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // BiCGSTAB / GMRES on nonsymmetric systems
  // ---------------------------------------------------------------------------

  test("BiCGSTAB on diagonally dominant CSR matches Breeze dense backslash") {
    for n <- List(10, 20); seed <- List(1L, 2L, 3L) do
      val data = diagonallyDominant(n, seed)
      val gA = galeCsr(data)
      val bA = breezeMatrix(data)
      val xData = vectorData(n, seed * 13 + 5)
      val gb = gA * galeVector(xData)
      val bb = bA * breezeVector(xData)

      val result = bicgstab(gA, gb, config)
      val breezeX = bA \ bb
      assertConvergedClose(result, breezeX, s"BiCGSTAB n=$n seed=$seed")
  }

  test("GMRES on diagonally dominant dense matches Breeze dense backslash") {
    for n <- List(8, 15); seed <- List(4L, 5L) do
      val data = diagonallyDominant(n, seed)
      val gA = galeMatrix(data)
      val bA = breezeMatrix(data)
      val xData = vectorData(n, seed * 17 + 7)
      val gb = gA * galeVector(xData)
      val bb = bA * breezeVector(xData)

      val result = gmres(gA, gb, config)
      val breezeX = bA \ bb
      assertConvergedClose(result, breezeX, s"GMRES n=$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Least-squares iterative paths
  // ---------------------------------------------------------------------------

  test("LSQR on tall full-rank systems matches Breeze dense backslash") {
    for (m, n, seed) <- Seq((20, 6, 21L), (30, 8, 22L), (16, 5, 23L)) do
      val aData = matrixData(m, n, seed)
      // Make columns better conditioned: add a diagonal boost on the leading n rows.
      var i = 0
      while i < n do
        aData(i)(i) += 3.0
        i += 1
      val gA = galeMatrix(aData)
      val bA = breezeMatrix(aData)
      val xData = vectorData(n, seed * 19 + 9)
      val gb = gA * galeVector(xData)
      val bb = bA * breezeVector(xData)

      val result = lsqr(gA, gb, config)
      val breezeX = bA \ bb
      assertConvergedClose(result, breezeX, s"LSQR ${m}x$n seed=$seed")
  }

  test("CGNR on tall full-rank systems matches Breeze dense backslash") {
    for (m, n, seed) <- Seq((18, 5, 31L), (24, 7, 32L)) do
      val aData = matrixData(m, n, seed)
      var i = 0
      while i < n do
        aData(i)(i) += 4.0
        i += 1
      val gA = galeMatrix(aData)
      val bA = breezeMatrix(aData)
      val xData = vectorData(n, seed * 23 + 11)
      val gb = gA * galeVector(xData)
      val bb = bA * breezeVector(xData)

      val result = cgnr(gA, gb, config)
      val breezeX = bA \ bb
      assertConvergedClose(result, breezeX, s"CGNR ${m}x$n seed=$seed")
  }
