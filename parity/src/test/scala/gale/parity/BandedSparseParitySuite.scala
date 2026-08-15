package gale.parity

import breeze.linalg.CSCMatrix
import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.sparse.*

/** Sparse / banded parity: gale's `Banded`, `CSR`, `CSC`, `COO`, `Diagonal`, and
  * structured constructors (`identity` / `zero` / `permutation`) versus the
  * equivalent `breeze.linalg.CSCMatrix` operations — matvec, transpose-matvec,
  * `+`/`−`/scale, sparse×dense products, and inspect helpers (`apply` / `row` /
  * `col` / `t` / `trace` / `toDense`). These are exact same-arithmetic sums of
  * products, so they must agree to `1e-12`.
  *
  * gale exposes no public sparse '''solve''' (its solvers take dense `DMat` or a
  * matrix-free `DoubleLinearOperator`, not a stored sparse factorization), so there
  * is no sparse-solve surface to compare here; dense solve parity is covered by
  * `FactorizationParitySuite`.
  */
class BandedSparseParitySuite extends munit.FunSuite:

  private val tol = 1e-12

  /** Compare a gale sparse matrix to a Breeze CSC by probing every stored and a
    * sample of structural zero locations (elementwise, not just nnz pattern).
    */
  private def assertSparseClose(g: SparseMatrix[Double], b: CSCMatrix[Double], tol: Double, clue: String): Unit =
    assertEquals((g.rows, g.cols), (b.rows, b.cols), clue)
    var i = 0
    while i < g.rows do
      var j = 0
      while j < g.cols do
        assertScalarClose(g(i, j), b(i, j), tol, s"$clue ($i,$j)")
        j += 1
      i += 1

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

  /** Row `i` of a dense array as a Breeze vector (CSC has no `::` row slice). */
  private def breezeRow(data: Array[Array[Double]], i: Int): BDV[Double] =
    BDV(data(i).clone())

  /** Column `j` of a dense array as a Breeze vector. */
  private def breezeCol(data: Array[Array[Double]], j: Int): BDV[Double] =
    BDV.tabulate(data.length)(i => data(i)(j))

  private def breezeTrace(data: Array[Array[Double]]): Double =
    var s = 0.0
    var i = 0
    while i < data.length do
      s += data(i)(i)
      i += 1
    s

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

  test("Diagonal matvec / transpose-matvec vs breeze CSCMatrix") {
    for n <- List(4, 10, 25); seed <- List(1L, 2L) do
      val diag = vectorData(n, seed * 67 + 1)
      val data = Array.tabulate(n, n)((i, j) => if i == j then diag(i) else 0.0)
      val gA   = Sparse.diagonal(diag.toIndexedSeq*)
      val bA   = breezeCsc(data)
      val bAt  = breezeCscTransposed(data)

      val xData = vectorData(n, seed * 71 + 3)
      val gx = galeVector(xData)
      val bx = breezeVector(xData)
      assertVecClose(gA * gx, bA * bx, tol, s"Diagonal A·x n=$n seed=$seed")
      assertVecClose(gA.t * gx, bAt * bx, tol, s"Diagonal Aᵀ·x n=$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // CSR / CSC arithmetic (+, −, scalar *)
  // ---------------------------------------------------------------------------

  test("CSR add/subtract/scale vs breeze CSCMatrix") {
    for (m, n) <- List((8, 8), (10, 6), (5, 12)); seed <- List(1L, 2L) do
      val aData = sparseData(m, n, 0.35, seed)
      val bData = sparseData(m, n, 0.35, seed * 17 + 3)
      val gA = galeCsr(aData)
      val gB = galeCsr(bData)
      val bA = breezeCsc(aData)
      val bB = breezeCsc(bData)
      val alpha = 1.7 - seed.toDouble * 0.1

      assertSparseClose(gA + gB, bA + bB, tol, s"CSR A+B ${m}x$n seed=$seed")
      assertSparseClose(gA - gB, bA - bB, tol, s"CSR A−B ${m}x$n seed=$seed")
      assertSparseClose(gA * alpha, bA * alpha, tol, s"CSR αA ${m}x$n seed=$seed")
  }

  test("CSC add/subtract/scale vs breeze CSCMatrix") {
    for (m, n) <- List((8, 8), (9, 5)); seed <- List(3L, 4L) do
      val aData = sparseData(m, n, 0.35, seed)
      val bData = sparseData(m, n, 0.35, seed * 19 + 5)
      val gA = galeCsr(aData).toCSC
      val gB = galeCsr(bData).toCSC
      val bA = breezeCsc(aData)
      val bB = breezeCsc(bData)
      val alpha = -0.5 * seed.toDouble

      assertSparseClose(gA + gB, bA + bB, tol, s"CSC A+B ${m}x$n seed=$seed")
      assertSparseClose(gA - gB, bA - bB, tol, s"CSC A−B ${m}x$n seed=$seed")
      assertSparseClose(gA * alpha, bA * alpha, tol, s"CSC αA ${m}x$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // COO and structured sparse types
  // ---------------------------------------------------------------------------

  private def galeCoo(data: Array[Array[Double]]): COO =
    val builder = Sparse.coo(data.length, data(0).length)
    var i = 0
    while i < data.length do
      var j = 0
      while j < data(0).length do
        if data(i)(j) != 0.0 then builder.add(i, j, data(i)(j))
        j += 1
      i += 1
    builder.toCOO()

  test("COO matvec / transpose-matvec vs breeze CSCMatrix (rectangular)") {
    for (m, n) <- List((8, 8), (12, 7), (7, 15)); seed <- List(1L, 2L) do
      val data = sparseData(m, n, 0.3, seed)
      val gA = galeCoo(data)
      val bA = breezeCsc(data)
      val bAt = breezeCscTransposed(data)

      val xF = galeVector(vectorData(n, seed * 73 + 1))
      val xT = galeVector(vectorData(m, seed * 79 + 2))
      val bxF = breezeVector(vectorData(n, seed * 73 + 1))
      val bxT = breezeVector(vectorData(m, seed * 79 + 2))

      assertVecClose(gA * xF, bA * bxF, tol, s"COO A·x ${m}x$n seed=$seed")
      assertVecClose(gA.t * xT, bAt * bxT, tol, s"COO Aᵀ·x ${m}x$n seed=$seed")
  }

  test("Identity, Zero, and Permutation matvec match breeze CSCMatrix") {
    for n <- List(4, 9, 16); seed <- List(1L, 2L) do
      val xData = vectorData(n, seed * 101 + 1)
      val gx = galeVector(xData)
      val bx = breezeVector(xData)

      val identity = Array.tabulate(n, n)((i, j) => if i == j then 1.0 else 0.0)
      assertVecClose(Sparse.identity(n) * gx, breezeCsc(identity) * bx, tol, s"Identity A·x n=$n seed=$seed")
      assertVecClose(Sparse.identity(n).t * gx, breezeCsc(identity) * bx, tol, s"Identity Aᵀ·x n=$n seed=$seed")

      val zero = Array.fill(n, n)(0.0)
      assertVecClose(Sparse.zero(n, n) * gx, breezeCsc(zero) * bx, tol, s"Zero A·x n=$n seed=$seed")

      val rng = new scala.util.Random(seed)
      val columnsByRow = rng.shuffle((0 until n).toList).toArray
      val permData = Array.tabulate(n, n)((i, j) => if columnsByRow(i) == j then 1.0 else 0.0)
      val gP = Sparse.permutation(columnsByRow.toIndexedSeq*)
      val bP = breezeCsc(permData)
      val bPt = breezeCscTransposed(permData)
      assertVecClose(gP * gx, bP * bx, tol, s"Permutation A·x n=$n seed=$seed")
      assertVecClose(gP.t * gx, bPt * bx, tol, s"Permutation Aᵀ·x n=$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Structured constructors
  // ---------------------------------------------------------------------------

  test("Sparse.identity / zero / permutation match Breeze CSC images") {
    for n <- List(3, 8) do
      val gId = Sparse.identity(n)
      val bId = breezeCsc(Array.tabulate(n, n)((i, j) => if i == j then 1.0 else 0.0))
      assertSparseClose(gId, bId, tol, s"identity n=$n")
      val xData = vectorData(n, n * 101L + 1)
      assertVecClose(gId * galeVector(xData), bId * breezeVector(xData), tol, s"identity A·x n=$n")

      val gZ = Sparse.zero(n, n + 2)
      val bZ = breezeCsc(Array.fill(n, n + 2)(0.0))
      assertSparseClose(gZ, bZ, tol, s"zero ${n}x${n + 2}")
      assertVecClose(
        gZ * galeVector(vectorData(n + 2, n * 103L + 2)),
        bZ * breezeVector(vectorData(n + 2, n * 103L + 2)),
        tol,
        s"zero A·x ${n}x${n + 2}"
      )

    // columnsByRow = [2,0,1] ⇒ into(row) = x(columnsByRow(row))
    val colsByRow = IndexedSeq(2, 0, 1)
    val gP = Sparse.permutation(colsByRow*)
    val pData = Array.tabulate(3, 3)((i, j) => if colsByRow(i) == j then 1.0 else 0.0)
    val bP = breezeCsc(pData)
    assertSparseClose(gP, bP, tol, "permutation structure")
    val px = vectorData(3, 77L)
    assertVecClose(gP * galeVector(px), bP * breezeVector(px), tol, "permutation A·x")
    assertVecClose(gP.t * galeVector(px), breezeCscTransposed(pData) * breezeVector(px), tol, "permutation Aᵀ·x")
  }

  // ---------------------------------------------------------------------------
  // Sparse × dense products
  // ---------------------------------------------------------------------------

  test("CSR * dense matrix matches Breeze CSC * DenseMatrix") {
    for (m, k, n, seed) <- Seq((8, 5, 3, 41L), (6, 6, 4, 42L), (10, 4, 2, 43L)) do
      val aData = sparseData(m, k, 0.4, seed)
      val bData = matrixData(k, n, seed * 9 + 1)
      val gA = galeCsr(aData)
      val gB = galeMatrix(bData)
      val bA = breezeCsc(aData)
      val bB = breezeMatrix(bData)
      assertMatClose(gA * gB, bA * bB, tol, s"CSR*DMat ${m}x$k · ${k}x$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Inspect surface
  // ---------------------------------------------------------------------------

  test("CSR/CSC apply, row, col, transpose, trace, toDense vs Breeze CSC") {
    for (m, n, seed) <- Seq((7, 7, 51L), (6, 9, 52L)) do
      val data = sparseData(m, n, 0.35, seed)
      // Ensure a few diagonals exist so trace is exercised on square cases.
      if m == n then
        var d = 0
        while d < m do
          data(d)(d) = 1.5 + d * 0.1
          d += 1
      val gCsr = galeCsr(data)
      val gCsc = gCsr.toCSC
      val bA = breezeCsc(data)
      val bAt = breezeCscTransposed(data)
      val denseRef = BDM.tabulate(m, n)((r, c) => data(r)(c))

      assertSparseClose(gCsr, bA, tol, s"CSR apply ${m}x$n seed=$seed")
      assertSparseClose(gCsc, bA, tol, s"CSC apply ${m}x$n seed=$seed")
      assertSparseClose(gCsr.t, bAt, tol, s"CSR.t ${m}x$n seed=$seed")
      assertSparseClose(gCsc.t, bAt, tol, s"CSC.t ${m}x$n seed=$seed")

      var i = 0
      while i < m do
        assertVecClose(gCsr.row(i), breezeRow(data, i), tol, s"CSR.row($i) ${m}x$n seed=$seed")
        assertVecClose(gCsc.row(i), breezeRow(data, i), tol, s"CSC.row($i) ${m}x$n seed=$seed")
        i += 1
      var j = 0
      while j < n do
        assertVecClose(gCsr.col(j), breezeCol(data, j), tol, s"CSR.col($j) ${m}x$n seed=$seed")
        assertVecClose(gCsc.col(j), breezeCol(data, j), tol, s"CSC.col($j) ${m}x$n seed=$seed")
        j += 1

      if m == n then
        assertScalarClose(gCsr.trace, breezeTrace(data), tol, s"CSR.trace n=$n seed=$seed")
        assertScalarClose(gCsc.trace, breezeTrace(data), tol, s"CSC.trace n=$n seed=$seed")

      assertMatClose(gCsr.toDense(), denseRef, tol, s"CSR.toDense ${m}x$n seed=$seed")
      assertMatClose(gCsc.toDense(), denseRef, tol, s"CSC.toDense ${m}x$n seed=$seed")
  }
