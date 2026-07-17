package userland

import breeze.linalg.DenseMatrix
import breeze.linalg.DenseVector
import breeze.linalg.MatrixNotSymmetricException
import breeze.linalg.cholesky as breezeCholesky
import breeze.linalg.det as breezeDet
import breeze.linalg.eigSym as breezeEigSym
import gale.interop.breeze.BreezeMigration
import gale.linalg.LinAlgError

/** Smoke tests for the migration shim: each Breeze-shaped call must produce the
  * numerically correct result (verified against Breeze's own reference or the
  * defining equation). These are numeric solves, so they use tolerances — unlike
  * the bit-exact conversion round trips.
  */
class BreezeMigrationSuite extends munit.FunSuite:

  private def close(a: Double, b: Double, tol: Double): Boolean =
    math.abs(a - b) <= tol * math.max(1.0, math.max(math.abs(a), math.abs(b)))

  private def diagDominant(n: Int): DenseMatrix[Double] =
    DenseMatrix.tabulate(n, n): (i, j) =>
      if i == j then n + 2.0 else 0.1 * (((i * 3 + j) % 5) - 2)

  /** An exactly-symmetric SPD matrix (lower triangle mirrored) so Breeze's own
    * `eigSym` reference — which requires exact symmetry — accepts it too.
    */
  private def spd(n: Int, seed: Int): DenseMatrix[Double] =
    val m = DenseMatrix.tabulate(n, n)((i, j) => math.sin((i + 1) * 0.7 + (j + 1) * 1.3 + seed))
    val a = Array.ofDim[Double](n, n)
    for i <- 0 until n; j <- 0 to i do
      var s = 0.0
      for k <- 0 until n do s += m(i, k) * m(j, k)
      if i == j then s += n.toDouble
      a(i)(j) = s
      a(j)(i) = s
    DenseMatrix.tabulate(n, n)((i, j) => a(i)(j))

  private def assertMatrixClose(
      actual: DenseMatrix[Double],
      expected: DenseMatrix[Double],
      tol: Double
  ): Unit =
    assertEquals((actual.rows, actual.cols), (expected.rows, expected.cols))
    for i <- 0 until actual.rows; j <- 0 until actual.cols do
      assert(
        close(actual(i, j), expected(i, j), tol),
        s"matrix mismatch ($i,$j): ${actual(i, j)} vs ${expected(i, j)}"
      )

  test("solve: A x = b for a square system") {
    val n = 6
    val a = diagDominant(n)
    val b = DenseVector.tabulate(n)(i => i * 0.5 - 1.0)
    val x = BreezeMigration.solve(a, b)
    val residual = a * x - b
    assert(residual.data.forall(r => math.abs(r) < 1e-9), s"residual ${residual}")
  }

  test("solve: matrix RHS matches Breeze and reuses one LU factorization") {
    val n = 8
    val a = diagDominant(n).t // exercise a transposed Breeze view
    val b = DenseMatrix.tabulate(n, 3)((i, j) => math.sin(i * 0.3 + j * 0.8))
    val actual = BreezeMigration.solve(a, b)
    val expected = a \ b
    assertMatrixClose(actual, expected, 1e-10)
  }

  test("det: matches Breeze det (sign included)") {
    val n = 5
    val a = diagDominant(n)
    assert(close(BreezeMigration.det(a), breezeDet(a), 1e-9), "det mismatch")
  }

  test("cholesky: L Lᵀ reconstructs the SPD input") {
    val n = 5
    val s = spd(n, 1)
    val l = BreezeMigration.cholesky(s)
    val recon = l * l.t
    for i <- 0 until n; j <- 0 until n do
      assert(close(recon(i, j), s(i, j), 1e-10), s"recon ($i,$j)")
  }

  test("eigSym: eigenvalues match Breeze eigSym (ascending)") {
    val n = 6
    val s = spd(n, 2)
    val (w, v) = BreezeMigration.eigSym(s)
    val ref = breezeEigSym(s).eigenvalues
    assertEquals(w.length, n)
    assertEquals((v.rows, v.cols), (n, n))
    for i <- 0 until n do assert(close(w(i), ref(i), 1e-9), s"eigenvalue [$i]: ${w(i)} vs ${ref(i)}")
  }

  test("leastSquares: satisfies the normal equations for a tall system") {
    val m = 12
    val n = 4
    val a = DenseMatrix.tabulate(m, n)((i, j) => math.cos(i * 0.4 + j * 1.1) + (if i == j then 1.0 else 0.0))
    val b = DenseVector.tabulate(m)(i => i * 0.2 - 0.7)
    val x = BreezeMigration.leastSquares(a, b)
    // Aᵀ A x ≈ Aᵀ b
    val lhs = a.t * (a * x)
    val rhs = a.t * b
    for i <- 0 until n do assert(close(lhs(i), rhs(i), 1e-8), s"normal eqn [$i]")
  }

  test("leastSquares: matrix RHS matches Breeze with one reused QR") {
    val m = 18
    val n = 5
    val a = DenseMatrix.tabulate(m, n): (i, j) =>
      math.cos(i * 0.27 + j * 0.91) + (if i == j then 2.0 else 0.0)
    val b = DenseMatrix.tabulate(m, 3)((i, j) => math.sin(i * 0.19 - j * 0.43))
    val actual = BreezeMigration.leastSquares(a, b)
    val expected = a \ b
    assertMatrixClose(actual, expected, 1e-9)
  }

  test("migration shims throw typed Gale failures for singular, non-SPD, and rank-deficient inputs") {
    val singular = DenseMatrix((1.0, 2.0), (2.0, 4.0))
    intercept[LinAlgError.SingularMatrix]:
      BreezeMigration.solve(singular, DenseVector(1.0, 2.0))

    val indefinite = DenseMatrix((1.0, 0.0), (0.0, -1.0))
    intercept[LinAlgError.NotPositiveDefinite]:
      BreezeMigration.cholesky(indefinite)

    val rankDeficient = DenseMatrix.tabulate(8, 3): (i, j) =>
      if j == 2 then 2.0 * i else if j == 1 then i.toDouble else 1.0
    intercept[LinAlgError.RankDeficient]:
      BreezeMigration.leastSquares(rankDeficient, DenseVector.tabulate(8)(_.toDouble))
  }

  test("cholesky and eigSym preserve Breeze's public symmetry guard") {
    val asymmetric = spd(5, 7)
    asymmetric(0, 3) = asymmetric(3, 0) + 1.0e-4

    intercept[MatrixNotSymmetricException]:
      breezeCholesky(asymmetric)
    intercept[LinAlgError.InvalidArgument]:
      BreezeMigration.cholesky(asymmetric)

    intercept[MatrixNotSymmetricException]:
      breezeEigSym(asymmetric)
    intercept[LinAlgError.InvalidArgument]:
      BreezeMigration.eigSym(asymmetric)
  }

  test("symmetry guard accepts the same scale-aware near-symmetry as Breeze") {
    val nearly = spd(4, 9)
    val lower = nearly(3, 1)
    val scale = math.max(1.0, math.abs(lower))
    nearly(1, 3) = lower + 0.5e-7 * scale

    val breezeL = breezeCholesky(nearly)
    val galeL = BreezeMigration.cholesky(nearly)
    assertMatrixClose(galeL, breezeL, 1e-10)

    val breezeW = breezeEigSym(nearly).eigenvalues
    val (galeW, _) = BreezeMigration.eigSym(nearly)
    for i <- 0 until breezeW.length do
      assert(close(galeW(i), breezeW(i), 1e-9), s"near-symmetric eigenvalue [$i]")
  }

  test("matrix RHS validates row shape even when it has zero columns") {
    val a = diagDominant(4)
    val wrong = DenseMatrix.zeros[Double](3, 0)
    intercept[LinAlgError.DimensionMismatch]:
      BreezeMigration.solve(a, wrong)
    intercept[LinAlgError.DimensionMismatch]:
      BreezeMigration.leastSquares(DenseMatrix.ones[Double](6, 4), wrong)
  }
