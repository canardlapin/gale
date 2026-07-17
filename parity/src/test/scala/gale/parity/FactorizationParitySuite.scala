package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.LU as BreezeLU
import breeze.linalg.cholesky
import breeze.linalg.cond
import breeze.linalg.det
import breeze.linalg.qr
import breeze.linalg.rank
import gale.linalg.*
import gale.parity.ParitySupport.*

/** Factorization and linear-solve parity versus Breeze.
  *
  * Factorizations are not unique — pivoting order and reflector/sign conventions
  * differ between gale's pure kernels and Breeze's LAPACK backend — so these tests
  * compare '''invariants''' rather than raw factors: determinant values (sign
  * included), solution vectors, `A = L Lᵀ` / `A = Q R` reconstructions, and
  * sign-free products such as `Rᵀ R = Aᵀ A`.
  */
class FactorizationParitySuite extends munit.FunSuite:

  private val solveTol = 1e-10
  private val detTol   = 1e-11
  private val reconTol = 1e-11
  private val cholTol  = 1e-10

  // ---------------------------------------------------------------------------
  // Determinant
  // ---------------------------------------------------------------------------

  test("det: gale vs breeze on well-conditioned matrices") {
    for n <- List(4, 8, 15); seed <- List(1L, 2L, 3L) do
      val data = diagonallyDominant(n, seed)
      val gd   = galeMatrix(data).det.orThrow
      val bd   = det(breezeMatrix(data))
      assertScalarClose(gd, bd, detTol, s"det n=$n seed=$seed")
  }

  test("det: sign parity under a row swap (negative determinant)") {
    for n <- List(4, 9); seed <- List(5L, 6L) do
      val data = diagonallyDominant(n, seed)
      val tmp  = data(0); data(0) = data(1); data(1) = tmp // swap two rows -> flips sign
      val gd   = galeMatrix(data).det.orThrow
      val bd   = det(breezeMatrix(data))
      assert(gd < 0.0 == bd < 0.0, s"det sign disagreement n=$n seed=$seed: gale=$gd breeze=$bd")
      assertScalarClose(gd, bd, detTol, s"det(swapped) n=$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Solve A x = b
  // ---------------------------------------------------------------------------

  test("solve: gale A.solve(b) vs breeze A \\ b") {
    for n <- List(5, 12, 25); seed <- List(1L, 2L, 3L) do
      val aData = diagonallyDominant(n, seed)
      val bData = vectorData(n, seed * 41 + 1)
      val gx    = galeMatrix(aData).solve(galeVector(bData)).orThrow
      val bx    = breezeMatrix(aData) \ breezeVector(bData)
      assertVecClose(gx, bx, solveTol, s"solve n=$n seed=$seed")
  }

  test("solve: LU-reused solve matches a fresh breeze solve for several RHS") {
    val n  = 20
    val a  = diagonallyDominant(n, 99L)
    val lu = galeMatrix(a).lu.orThrow
    val ba = breezeMatrix(a)
    for r <- 0 until 4 do
      val bData = vectorData(n, 500L + r)
      val gx    = lu.solve(galeVector(bData)).orThrow
      val bx    = ba \ breezeVector(bData)
      assertVecClose(gx, bx, solveTol, s"LU solve rhs=$r")
  }

  // ---------------------------------------------------------------------------
  // LU reconstruction (A = P⁻¹ L U)
  // ---------------------------------------------------------------------------

  test("LU: both reconstruct the original A") {
    for n <- List(5, 11, 20); seed <- List(1L, 2L) do
      val data = diagonallyDominant(n, seed)

      // gale: packed L\U + pivots. pivots(i) is the original row now in row i, so
      // (L·U)(i, j) must equal A(pivots(i), j).
      val lu   = galeMatrix(data).lu.orThrow
      val piv  = lu.pivots
      val lMat = Matrix.tabulate(n, n)((i, j) => if i == j then 1.0 else if j < i then lu.packed(i, j) else 0.0)
      val uMat = Matrix.tabulate(n, n)((i, j) => if j >= i then lu.packed(i, j) else 0.0)
      val recon = lMat * uMat
      var i = 0
      while i < n do
        var j = 0
        while j < n do
          assertScalarClose(recon(i, j), data(piv(i))(j), reconTol, s"gale LU recon n=$n seed=$seed ($i,$j)")
          j += 1
        i += 1

      // breeze: P·L·U == A.
      val luB = BreezeLU(breezeMatrix(data))
      assertBreezeMatClose(luB.P * luB.L * luB.U, breezeMatrix(data), reconTol, s"breeze LU recon n=$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Cholesky (A = L Lᵀ, lower factor is unique)
  // ---------------------------------------------------------------------------

  test("cholesky: lower factor matches elementwise and reconstructs A") {
    for n <- List(4, 10, 20); seed <- List(1L, 2L, 3L) do
      val data = spd(n, seed)
      val gL   = galeMatrix(data).cholesky.orThrow.lower
      val bL   = cholesky(breezeMatrix(data))
      // Cholesky is unique for SPD with positive diagonal: compare the lower
      // triangles directly (breeze may leave nonzero garbage above the diagonal,
      // so only the lower triangle is compared).
      var i = 0
      while i < n do
        var j = 0
        while j <= i do
          assertScalarClose(gL(i, j), bL(i, j), cholTol, s"chol L n=$n seed=$seed ($i,$j)")
          j += 1
        i += 1
      // Reconstruction invariant for both.
      assertMatClose(gL * gL.t, breezeMatrix(data), reconTol, s"gale chol recon n=$n seed=$seed")
      // breeze's cholesky zeroes the strict upper triangle, so bL·bLᵀ reconstructs A.
      assertBreezeMatClose(bL * bL.t, breezeMatrix(data), reconTol, s"breeze chol recon n=$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // QR (Q R = A; compare sign-free invariants)
  // ---------------------------------------------------------------------------

  test("QR: reconstruction, orthonormality, and Rᵀ R = Aᵀ A parity") {
    for (m, n) <- List((5, 5), (8, 5), (20, 12)); seed <- List(1L, 2L) do
      val data = matrixData(m, n, seed)
      val ga   = galeMatrix(data)
      val ba   = breezeMatrix(data)

      val gqr = ga.qr
      val bqr = qr(ba) // complete QR: Q is m×m, R is m×n

      // Reconstruction Q·R == A for both.
      assertMatClose(gqr.q * gqr.r, ba, reconTol, s"gale QR recon ${m}x$n seed=$seed")
      assertBreezeMatClose(bqr.q * bqr.r, ba, reconTol, s"breeze QR recon ${m}x$n seed=$seed")

      // Orthonormality Qᵀ Q == I for both.
      assertMatClose(gqr.q.t * gqr.q, BDM.eye[Double](m), reconTol, s"gale Q orthonormal ${m}x$n seed=$seed")
      assertBreezeMatClose(bqr.q.t * bqr.q, BDM.eye[Double](m), reconTol, s"breeze Q orthonormal ${m}x$n seed=$seed")

      // Sign-free cross-library parity: Rᵀ R == Aᵀ A regardless of reflector signs.
      val ata = ba.t * ba
      assertMatClose(gqr.r.t * gqr.r, ata, reconTol, s"gale RᵀR=AᵀA ${m}x$n seed=$seed")
      assertBreezeMatClose(bqr.r.t * bqr.r, ata, reconTol, s"breeze RᵀR=AᵀA ${m}x$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Least squares (tall systems; the QR path as a SOLVER vs an independent ref)
  // ---------------------------------------------------------------------------

  test("least squares: gale leastSquares vs breeze A \\ b on tall systems") {
    for (m, n) <- List((10, 4), (25, 12), (40, 7)); seed <- List(1L, 2L) do
      val aData = matrixData(m, n, seed)
      val bData = vectorData(m, seed * 17 + 3)
      val gx    = galeMatrix(aData).leastSquares(galeVector(bData)).orThrow
      val bx    = breezeMatrix(aData) \ breezeVector(bData)
      assertVecClose(gx, bx, solveTol, s"least squares ${m}x$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Rank / conditioning: shared semantic subset
  // ---------------------------------------------------------------------------

  test("rankEstimate agrees with Breeze rank on clear full-rank and deficient cases") {
    val full = matrixData(12, 7, 71L)
    val deficient = Array.tabulate(12, 7): (i, j) =>
      // Use an exactly zero column, not a floating linear combination. Breeze's
      // SVD rank cutoff and Gale's QR cutoff are deliberately different, so a
      // numerically tiny singular value is outside the claimed overlap.
      if j == 6 then 0.0 else full(i)(j)

    for (label, data) <- Seq("full" -> full, "deficient" -> deficient) do
      val galeRank = galeMatrix(data).rankEstimate
      val breezeRank = rank(breezeMatrix(data))
      assertEquals(galeRank, breezeRank, s"rank mismatch for $label matrix")
  }

  test("conditionEstimate agrees with Breeze cond where 1-norm and 2-norm coincide") {
    // Breeze `cond` is exact SVD 2-norm conditioning. Gale deliberately returns
    // a Hager/Higham 1-norm estimate. Diagonal matrices are the honest overlap:
    // both condition numbers are max(|d|)/min(|d|).
    for diagonal <- Seq(
      Array(0.5, 1.0, 4.0, 25.0),
      Array(-1000.0, 7.0, -2.0, 0.25)
    ) do
      val data = Array.tabulate(diagonal.length, diagonal.length): (i, j) =>
        if i == j then diagonal(i) else 0.0
      val galeCond = galeMatrix(data).conditionEstimate.orThrow
      val breezeCond = cond(breezeMatrix(data))
      assertScalarClose(galeCond, breezeCond, 1e-11, s"diagonal condition ${diagonal.mkString(",")}")
  }
