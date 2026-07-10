package gale.parity

import breeze.linalg.CSCMatrix
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.sparse.*

/** Sparse / banded parity: gale's `Banded`, `CSR`, `CSC`, and `Diagonal` matvec and
  * transpose-matvec versus the equivalent `breeze.linalg.CSCMatrix` operations, on
  * the same sparsity pattern. These are exact same-arithmetic sums of products, so
  * they must agree to `1e-12`.
  *
  * gale exposes no public sparse '''solve''' (its solvers take dense `DMat` or a
  * matrix-free `DoubleLinearOperator`, not a stored sparse factorization), so there
  * is no sparse-solve surface to compare here; dense solve parity is covered by
  * `FactorizationParitySuite`.
  */
class BandedSparseParitySuite extends munit.FunSuite:

  private val tol = 1e-12

  /** A `rows × cols` array whose nonzeros lie within bandwidth `[-kl, ku]`. */
  private def bandedData(rows: Int, cols: Int, kl: Int, ku: Int, seed: Long): Array[Array[Double]] =
    val rng = new scala.util.Random(seed)
    Array.tabulate(rows, cols): (i, j) =>
      val d = i - j
      if d <= kl && d >= -ku then rng.nextDouble() * 2.0 - 1.0 else 0.0

  /** A `rows × cols` array with roughly `density` fraction of nonzeros. */
  private def sparseData(rows: Int, cols: Int, density: Double, seed: Long): Array[Array[Double]] =
    val rng = new scala.util.Random(seed)
    Array.tabulate(rows, cols): (_, _) =>
      if rng.nextDouble() < density then rng.nextDouble() * 2.0 - 1.0 else 0.0

  private def breezeCsc(data: Array[Array[Double]]): CSCMatrix[Double] =
    CSCMatrix.tabulate(data.length, data(0).length)((i, j) => data(i)(j))

  /** A breeze CSC of the transpose (built directly, avoiding `.t` view ambiguity). */
  private def breezeCscTransposed(data: Array[Array[Double]]): CSCMatrix[Double] =
    CSCMatrix.tabulate(data(0).length, data.length)((i, j) => data(j)(i))

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

  // ---------------------------------------------------------------------------
  // Banded
  // ---------------------------------------------------------------------------

  test("Banded matvec / transpose-matvec vs breeze CSCMatrix") {
    for n <- List(6, 16, 40); (kl, ku) <- List((1, 1), (2, 0), (0, 2), (3, 2)); seed <- List(1L, 2L) do
      val data = bandedData(n, n, kl, ku, seed)
      val gA   = Sparse.banded(galeMatrix(data))
      val bA   = breezeCsc(data)
      val bAt  = breezeCscTransposed(data)

      val xData = vectorData(n, seed * 43 + 1)
      val gx    = galeVector(xData)
      val bx    = breezeVector(xData)

      assertVecClose(gA * gx, bA * bx, tol, s"Banded A·x n=$n band=($kl,$ku) seed=$seed")
      assertVecClose(gA.t * gx, bAt * bx, tol, s"Banded Aᵀ·x n=$n band=($kl,$ku) seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // CSR / CSC
  // ---------------------------------------------------------------------------

  test("CSR matvec / transpose-matvec vs breeze CSCMatrix (rectangular)") {
    for (m, n) <- List((8, 8), (12, 7), (7, 15)); seed <- List(1L, 2L, 3L) do
      val data = sparseData(m, n, 0.3, seed)
      val gA   = galeCsr(data)
      val bA   = breezeCsc(data)
      val bAt  = breezeCscTransposed(data)

      val xF = galeVector(vectorData(n, seed * 47 + 1))
      val xT = galeVector(vectorData(m, seed * 53 + 2))
      val bxF = breezeVector(vectorData(n, seed * 47 + 1))
      val bxT = breezeVector(vectorData(m, seed * 53 + 2))

      assertVecClose(gA * xF, bA * bxF, tol, s"CSR A·x ${m}x$n seed=$seed")
      assertVecClose(gA.t * xT, bAt * bxT, tol, s"CSR Aᵀ·x ${m}x$n seed=$seed")
  }

  test("CSC matvec / transpose-matvec vs breeze CSCMatrix (rectangular)") {
    for (m, n) <- List((8, 8), (12, 7), (7, 15)); seed <- List(4L, 5L) do
      val data = sparseData(m, n, 0.3, seed)
      val gA   = galeCsr(data).toCSC
      val bA   = breezeCsc(data)
      val bAt  = breezeCscTransposed(data)

      val xF = galeVector(vectorData(n, seed * 59 + 1))
      val xT = galeVector(vectorData(m, seed * 61 + 2))
      val bxF = breezeVector(vectorData(n, seed * 59 + 1))
      val bxT = breezeVector(vectorData(m, seed * 61 + 2))

      assertVecClose(gA * xF, bA * bxF, tol, s"CSC A·x ${m}x$n seed=$seed")
      assertVecClose(gA.t * xT, bAt * bxT, tol, s"CSC Aᵀ·x ${m}x$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Diagonal
  // ---------------------------------------------------------------------------

  test("Diagonal matvec vs breeze CSCMatrix") {
    for n <- List(4, 10, 25); seed <- List(1L, 2L) do
      val diag = vectorData(n, seed * 67 + 1)
      val data = Array.tabulate(n, n)((i, j) => if i == j then diag(i) else 0.0)
      val gA   = Sparse.diagonal(diag.toIndexedSeq*)
      val bA   = breezeCsc(data)

      val xData = vectorData(n, seed * 71 + 3)
      assertVecClose(gA * galeVector(xData), bA * breezeVector(xData), tol, s"Diagonal A·x n=$n seed=$seed")
  }
