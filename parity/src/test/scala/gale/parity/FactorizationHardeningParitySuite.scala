package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.LU as BreezeLU
import breeze.linalg.MatrixEmptyException
import breeze.linalg.MatrixSingularException
import breeze.linalg.NotConvergedException
import breeze.linalg.cholesky
import breeze.linalg.cond
import breeze.linalg.det
import breeze.linalg.inv
import breeze.linalg.pinv as breezePinv
import breeze.linalg.qr
import breeze.linalg.svd as breezeSvd
import gale.linalg.*
import gale.parity.ParitySupport.*

/** Factorization hardening versus Breeze (W5.3 of the Breeze competitiveness plan).
  *
  * The original [[FactorizationParitySuite]] uses well-conditioned inputs with
  * `n ≤ 25`. This suite adds what a future blocked kernel (blocked LU at
  * `n ≳ 96`, blocked Cholesky, trsm-routed multi-RHS solves) must survive: sizes
  * above every planned blocking threshold, `1 × 1` and empty inputs, exactly
  * singular inputs, ill-conditioned inputs, and non-contiguous views.
  *
  * '''Tolerances.''' Every forward comparison is `relDiff ≤ c · n · ε · κ` with
  * `c = 32`, `ε = 2⁻⁵²` and `κ` Breeze's SVD 2-norm condition number of the
  * actual input. That is the standard first-order forward-error bound for a
  * backward-stable solver (`‖δx‖/‖x‖ ≲ κ · γₙ`) applied to '''both''' libraries,
  * so their difference is at most twice it; the factor `c` absorbs the
  * 2-norm/∞-norm mismatch and modest pivot growth. Least squares uses `κ²`
  * (the LS sensitivity term). Where `κ ε ≳ 1` the forward error is meaningless,
  * so only the norm-wise backward error `η ≤ c · n · ε` is asserted, for both.
  *
  * '''Failure parity.''' Singular inputs are compared by ''failure'', not by
  * pivot index: gale returns a typed `Left`, Breeze throws.
  *
  * '''Behavioural divergences (asserted, not hidden).'''
  *  - `det` of an exactly singular matrix: gale `Left(SingularMatrix)`, Breeze
  *    returns `±0.0`.
  *  - Empty `0 × 0`: gale factors and solves it (`det = 1`, empty solutions);
  *    Breeze `det`/`inv`/`LU` throw `IndexOutOfBoundsException`, `cholesky`
  *    throws `MatrixEmptyException` and `qr` rejects it in `DGEQRF`.
  *  - Rank-deficient tall least squares: gale `Left(RankDeficient)`; Breeze `\`
  *    returns an unchecked non-solution, and Breeze `pinv` does not truncate the
  *    rounding-level singular value of the zero column (entries ~1e15), so gale's `pinv` is compared with `V Σ⁺ Uᵀ` built
  *    from Breeze's SVD instead.
  */
class FactorizationHardeningParitySuite extends munit.FunSuite:

  private val c = 32.0
  private val largeSizes = List(96, 128, 200, 300)

  private def forwardTol(n: Int, kappa: Double): Double = c * n * Eps * kappa

  private def kappa(data: Array[Array[Double]]): Double = cond(breezeMatrix(data))

  // ---------------------------------------------------------------------------
  // Above the blocking thresholds
  // ---------------------------------------------------------------------------

  test("large n LU: reconstruction, det, vector and multi-RHS solve vs breeze") {
    for (n, seed) <- largeSizes.zip(List(11L, 12L, 13L, 14L)) do
      val data = unitDiagonallyDominant(n, seed)
      val k = kappa(data)
      val tol = forwardTol(n, k)
      val lu = galeMatrix(data).lu.orThrow
      val packed = lu.packed
      val lower = Matrix.tabulate(n, n)((i, j) => if i == j then 1.0 else if j < i then packed(i, j) else 0.0)
      val upper = Matrix.tabulate(n, n)((i, j) => if j >= i then packed(i, j) else 0.0)
      val recon = lower * upper
      val reconTol = c * n * Eps * maxAbs(data)
      var worst = 0.0
      for i <- 0 until n; j <- 0 until n do
        worst = math.max(worst, math.abs(recon(i, j) - data(lu.pivots(i))(j)))
      assertBelow(worst, reconTol, s"gale LU recon n=$n")
      val luB = BreezeLU(breezeMatrix(data))
      val bRecon = luB.P * luB.L * luB.U
      assertBelow(relDiff(galeMatrix(data), bRecon) * maxAbs(data), reconTol, s"breeze LU recon n=$n")

      val gd = lu.det.orThrow
      val bd = det(breezeMatrix(data))
      assertBelow(math.abs(gd - bd) / math.abs(bd), tol, s"det n=$n (gale=$gd breeze=$bd)")

      val bData = vectorData(n, seed + 100)
      assertBelow(relDiff(lu.solve(galeVector(bData)).orThrow, breezeMatrix(data) \ breezeVector(bData)), tol, s"solve n=$n")
      val rhs = matrixData(n, 8, seed + 200)
      assertBelow(
        relDiff(lu.solve(galeMatrix(rhs)).orThrow, breezeMatrix(data) \ breezeMatrix(rhs)),
        tol,
        s"multi-RHS solve n=$n"
      )
  }

  test("large n Cholesky: lower factor and vector/multi-RHS solve vs breeze") {
    for (n, seed) <- largeSizes.zip(List(21L, 22L, 23L, 24L)) do
      val data = spd(n, seed)
      val tol = forwardTol(n, kappa(data))
      val factor = galeMatrix(data).cholesky.orThrow
      val bL = cholesky(breezeMatrix(data))
      val gLower = Matrix.tabulate(n, n)((i, j) => if j <= i then factor.lower(i, j) else 0.0)
      val bLower = BDM.tabulate(n, n)((i, j) => if j <= i then bL(i, j) else 0.0)
      assertBelow(relDiff(gLower, bLower), tol, s"chol L n=$n")
      val bData = vectorData(n, seed + 100)
      assertBelow(relDiff(factor.solve(galeVector(bData)).orThrow, breezeMatrix(data) \ breezeVector(bData)), tol, s"chol solve n=$n")
      val rhs = matrixData(n, 8, seed + 200)
      assertBelow(
        relDiff(factor.solve(galeMatrix(rhs)).orThrow, breezeMatrix(data) \ breezeMatrix(rhs)),
        tol,
        s"chol multi-RHS n=$n"
      )
  }

  test("large n inverse: gale solve(I) vs breeze inv") {
    for (n, seed) <- largeSizes.zip(List(31L, 32L, 33L, 34L)) do
      val data = unitDiagonallyDominant(n, seed)
      val tol = forwardTol(n, kappa(data))
      assertBelow(relDiff(galeMatrix(data).solve(Matrix.eye(n)).orThrow, inv(breezeMatrix(data))), tol, s"inv n=$n")
  }

  test("large n QR: reconstruction, orthonormality, and RᵀR = AᵀA") {
    for (n, seed) <- largeSizes.zip(List(41L, 42L, 43L, 44L)); m <- List(n, math.min(4 * n, 600)) do
      val data = matrixData(m, n, seed)
      val ba = breezeMatrix(data)
      val gqr = galeMatrix(data).qr
      val unitTol = c * m * Eps
      assertBelow(relDiff(gqr.q * gqr.r, ba), unitTol, s"QR recon ${m}x$n")
      assertBelow(relDiff(gqr.q.t * gqr.q, BDM.eye[Double](m)), unitTol, s"Q orthonormal ${m}x$n")
      val ata = ba.t * ba
      assertBelow(relDiff(gqr.r.t * gqr.r, ata), unitTol, s"RᵀR=AᵀA ${m}x$n")
      val bR = qr(ba).r
      val bRtR = bR.t * bR
      assertBelow(relDiff(galeMatrix(Array.tabulate(n, n)((i, j) => bRtR(i, j))), ata), unitTol, s"breeze RᵀR ${m}x$n")
  }

  test("large n least squares (tall, capped at 600 rows) vs breeze A \\ b") {
    for (n, seed) <- largeSizes.zip(List(51L, 52L, 53L, 54L)) do
      val m = math.min(4 * n, 600)
      val data = matrixData(m, n, seed)
      val k = kappa(data)
      val tol = c * m * Eps * k * k
      val bData = vectorData(m, seed + 100)
      assertBelow(
        relDiff(galeMatrix(data).leastSquares(galeVector(bData)).orThrow, breezeMatrix(data) \ breezeVector(bData)),
        tol,
        s"lstsq ${m}x$n κ=$k"
      )
      val rhs = matrixData(m, 4, seed + 200)
      assertBelow(
        relDiff(galeMatrix(data).leastSquares(galeMatrix(rhs)).orThrow, breezeMatrix(data) \ breezeMatrix(rhs)),
        tol,
        s"lstsq multi-RHS ${m}x$n"
      )
  }

  // ---------------------------------------------------------------------------
  // 1 × 1 and empty
  // ---------------------------------------------------------------------------

  test("1x1: LU, det, solve, inverse, Cholesky, QR and least squares match breeze") {
    for value <- List(3.5, -0.25, 1e-300, 7e250) do
      val data = Array(Array(value))
      val ga = galeMatrix(data)
      val ba = breezeMatrix(data)
      assertScalarClose(ga.det.orThrow, det(ba), 0.0, s"det [$value]")
      assertVecClose(ga.solve(Vec(2.0)).orThrow, ba \ BDV(2.0), Eps, s"solve [$value]")
      assertMatClose(ga.solve(Matrix(1, 3)(1.0, -2.0, 4.0)).orThrow, ba \ BDM((1.0, -2.0, 4.0)), Eps, s"solve RHS [$value]")
      assertMatClose(ga.solve(Matrix.eye(1)).orThrow, inv(ba), Eps, s"inv [$value]")
      val r = ga.qr.r(0, 0)
      assertScalarClose(math.abs(r), math.abs(qr(ba).r(0, 0)), 0.0, s"|R| [$value]")
      assertVecClose(ga.leastSquares(Vec(2.0)).orThrow, ba \ BDV(2.0), Eps, s"lstsq [$value]")
      if value > 0.0 then
        assertScalarClose(ga.cholesky.orThrow.lower(0, 0), cholesky(ba).apply(0, 0), Eps, s"chol [$value]")
      else
        assert(ga.cholesky.isLeft, s"gale chol must fail on [$value]")
        intercept[NotConvergedException](cholesky(ba))
  }

  test("empty 0x0: gale's documented results; breeze agrees only on backslash") {
    val ge = Matrix.zeros(0, 0)
    val be = BDM.zeros[Double](0, 0)
    // gale: factorizations of the empty matrix succeed; det is the empty product.
    assertEquals(ge.lu.orThrow.size, 0)
    assertEquals(ge.det.orThrow, 1.0)
    assertEquals(ge.solve(Vec.zeros(0)).orThrow.length, 0)
    val gx = ge.solve(Matrix.zeros(0, 3)).orThrow
    assertEquals((gx.rows, gx.cols), (0, 3))
    assertEquals(ge.cholesky.orThrow.solve(Vec.zeros(0)).orThrow.length, 0)
    assertEquals(ge.qr.r.rows, 0)
    assertEquals(ge.leastSquares(Vec.zeros(0)).orThrow.length, 0)
    assertEquals(Matrix.zeros(3, 0).leastSquares(Vec.zeros(3)).orThrow.length, 0)
    // breeze: backslash agrees; the factorizations throw.
    assertVecClose(ge.solve(Vec.zeros(0)).orThrow, be \ BDV.zeros[Double](0), 0.0, "empty solve")
    assertVecClose(
      Matrix.zeros(3, 0).leastSquares(Vec.zeros(3)).orThrow,
      BDM.zeros[Double](3, 0) \ BDV.zeros[Double](3),
      0.0,
      "3x0 lstsq"
    )
    intercept[IndexOutOfBoundsException](det(be))
    intercept[IndexOutOfBoundsException](inv(be))
    intercept[IndexOutOfBoundsException](BreezeLU(be))
    intercept[MatrixEmptyException](cholesky(be))
    intercept[IllegalArgumentException](qr(be))
  }

  // ---------------------------------------------------------------------------
  // Singular inputs: both fail
  // ---------------------------------------------------------------------------

  /** Zero column `j`: every elimination update leaves it exactly zero, so the
    * pivot search at step `j` (or earlier) finds an exact zero in both libraries.
    */
  private def zeroColumn(n: Int, seed: Long): Array[Array[Double]] =
    val a = matrixData(n, n, seed)
    val j = (seed % n).toInt
    a.foreach(row => row(j) = 0.0)
    a

  /** Two identical rows of small integers: elimination applies the same exact
    * operations to both, so their difference is exactly zero.
    */
  private def duplicatedRows(n: Int, seed: Long): Array[Array[Double]] =
    val rng = new scala.util.Random(seed)
    val a = Array.fill(n, n)((rng.nextInt(9) - 4).toDouble)
    a(n - 1) = a(0).clone()
    a

  test("singular LU: gale Left, breeze solve/inv throw MatrixSingularException, breeze det is 0") {
    val cases =
      List(3, 8, 96, 200).map(n => s"zero column n=$n" -> zeroColumn(n, n.toLong + 5)) ++
        List(3, 6, 12).map(n => s"duplicated rows n=$n" -> duplicatedRows(n, n.toLong))
    for (label, data) <- cases do
      val n = data.length
      val ga = galeMatrix(data)
      val ba = breezeMatrix(data)
      ga.lu match
        case Left(_: LinAlgError.SingularMatrix) => ()
        case other                               => fail(s"$label: gale lu = $other")
      assert(ga.det.isLeft, s"$label: gale det must be Left")
      assert(ga.solve(galeVector(vectorData(n, 1L))).isLeft, s"$label: gale solve must be Left")
      assert(ga.solve(Matrix.eye(n)).isLeft, s"$label: gale matrix solve must be Left")
      intercept[MatrixSingularException](ba \ breezeVector(vectorData(n, 1L)))
      intercept[MatrixSingularException](ba \ BDM.eye[Double](n))
      intercept[MatrixSingularException](inv(ba))
      assertEquals(math.abs(det(ba)), 0.0, s"$label: breeze det")
  }

  test("not positive definite: gale Cholesky Left, breeze throws NotConvergedException") {
    for n <- List(4, 40, 128) do
      val pd = spd(n, n.toLong)
      // Indefinite: shift past the largest eigenvalue bound (Gershgorin).
      val shift = pd.map(_.map(math.abs).sum).max + 1.0
      val indefinite = Array.tabulate(n, n)((i, j) => if i == j then pd(i)(j) - shift else pd(i)(j))
      // Semidefinite with an exactly zero diagonal pivot: zero row and column.
      val z = n / 2
      val semidefinite = Array.tabulate(n, n)((i, j) => if i == z || j == z then 0.0 else pd(i)(j))
      for (label, data) <- List("indefinite" -> indefinite, "zero row/col" -> semidefinite) do
        galeMatrix(data).cholesky match
          case Left(_: LinAlgError.NotPositiveDefinite) => ()
          case other                                    => fail(s"$label n=$n: gale cholesky = $other")
        intercept[NotConvergedException](cholesky(breezeMatrix(data)))
  }

  private def svdPseudoInverse(a: BDM[Double]): BDM[Double] =
    val decomposition = breezeSvd.reduced(a)
    val sigma = decomposition.S
    val cutoff = math.max(a.rows, a.cols) * Eps * sigma(0)
    val inverted = sigma.map(s => if s > cutoff then 1.0 / s else 0.0)
    val v = decomposition.Vt.t
    v * breeze.linalg.diag(inverted) * decomposition.U.t

  test("rank-deficient tall least squares: gale Left; pinv matches breeze pinv") {
    for (m, n, seed) <- List((12, 5, 61L), (40, 12, 62L), (150, 100, 63L)) do
      val data = matrixData(m, n, seed)
      data.foreach(row => row(n / 2) = 0.0)
      val ga = galeMatrix(data)
      ga.leastSquares(galeVector(vectorData(m, seed))) match
        case Left(_: LinAlgError.RankDeficient) => ()
        case other                              => fail(s"${m}x$n: gale lstsq = $other")
      // Minimum-norm parity: gale's SVD pseudo-inverse vs V Σ⁺ Uᵀ assembled from
      // breeze's SVD with the same max(m, n)·ε·σ₁ cutoff. (breeze's own `pinv` keeps
      // the rounding-level singular value, so its entries blow up to ~1e15.)
      val reduced = data.map(row => row.patch(n / 2, Nil, 1))
      val tol = c * m * Eps * kappa(reduced)
      assertBelow(relDiff(ga.pinv.orThrow, svdPseudoInverse(breezeMatrix(data))), tol, s"pinv ${m}x$n")
      assert(breeze.linalg.max(breeze.numerics.abs(breezePinv(breezeMatrix(data)))) > 1e8, s"breeze pinv ${m}x$n")
  }

  // ---------------------------------------------------------------------------
  // Ill-conditioned inputs
  // ---------------------------------------------------------------------------

  test("Hilbert 6..12: LU solve/det/inverse forward error where κε < 1e-3, backward error always") {
    for n <- 6 to 12 do
      val data = hilbert(n)
      val k = kappa(data)
      val bData = Array.tabulate(n)(i => 1.0 + i.toDouble)
      val gx = galeMatrix(data).solve(galeVector(bData)).orThrow
      val bx = breezeMatrix(data) \ breezeVector(bData)
      val backTol = c * n * Eps
      assertBelow(backwardError(data, gx(_), bData), backTol, s"gale η H$n")
      assertBelow(backwardError(data, bx(_), bData), backTol, s"breeze η H$n")
      if k * Eps < 1e-3 then
        val tol = forwardTol(n, k)
        assertBelow(relDiff(gx, bx), tol, s"solve H$n κ=$k")
        val gd = galeMatrix(data).det.orThrow
        val bd = det(breezeMatrix(data))
        assertBelow(math.abs(gd - bd) / math.abs(bd), tol, s"det H$n κ=$k")
        assertBelow(relDiff(galeMatrix(data).solve(Matrix.eye(n)).orThrow, inv(breezeMatrix(data))), tol, s"inv H$n κ=$k")
  }

  test("Hilbert 6..10: Cholesky factor, reconstruction and solve backward error") {
    for n <- 6 to 10 do
      val data = hilbert(n)
      val k = kappa(data)
      val factor = galeMatrix(data).cholesky.orThrow
      val bL = cholesky(breezeMatrix(data))
      val gLower = Matrix.tabulate(n, n)((i, j) => if j <= i then factor.lower(i, j) else 0.0)
      val bLower = BDM.tabulate(n, n)((i, j) => if j <= i then bL(i, j) else 0.0)
      assertBelow(relDiff(gLower * gLower.t, breezeMatrix(data)), c * n * Eps, s"gale LLᵀ H$n")
      if k * Eps < 1e-3 then assertBelow(relDiff(gLower, bLower), forwardTol(n, k), s"chol L H$n κ=$k")
      val bData = Array.tabulate(n)(i => if i % 2 == 0 then 1.0 else -1.0)
      val gx = factor.solve(galeVector(bData)).orThrow
      assertBelow(backwardError(data, gx(_), bData), c * n * Eps, s"gale chol η H$n")
  }

  test("prescribed spectrum κ = 1e4, 1e8, 1e12: LU/Cholesky solve, det and inverse scale with κε") {
    for kappaTarget <- List(1e4, 1e8, 1e12); n <- List(16, 64) do
      val data = withSpectrum(geometricSpectrum(n, kappaTarget), (n + kappaTarget.toLong % 97).toLong)
      val k = kappa(data)
      val tol = forwardTol(n, k)
      val clue = s"n=$n κ=$k"
      val bData = vectorData(n, 7L)
      val bx = breezeMatrix(data) \ breezeVector(bData)
      val gLu = galeMatrix(data).solve(galeVector(bData)).orThrow
      val gChol = galeMatrix(data).cholesky.orThrow.solve(galeVector(bData)).orThrow
      assertBelow(relDiff(gLu, bx), tol, s"LU solve $clue")
      assertBelow(relDiff(gChol, bx), tol, s"Cholesky solve $clue")
      assertBelow(backwardError(data, gLu(_), bData), c * n * Eps, s"LU η $clue")
      assertBelow(backwardError(data, gChol(_), bData), c * n * Eps, s"Cholesky η $clue")
      val bd = det(breezeMatrix(data))
      if math.abs(bd) > 1e-280 then // det = Πλ underflows for n = 64 at large κ
        val gd = galeMatrix(data).det.orThrow
        assertBelow(math.abs(gd - bd) / math.abs(bd), tol, s"det $clue")
      assertBelow(relDiff(galeMatrix(data).solve(Matrix.eye(n)).orThrow, inv(breezeMatrix(data))), tol, s"inv $clue")
  }

  // ---------------------------------------------------------------------------
  // Non-contiguous views (compared against breeze on materialized copies)
  // ---------------------------------------------------------------------------

  test("views: LU solve/det, Cholesky and QR on transposed and strided views") {
    for n <- List(7, 40, 130) do
      val general = unitDiagonallyDominant(n, n.toLong)
      val tol = forwardTol(n, kappa(general))
      val bData = vectorData(n, 3L)
      val bx = breezeMatrix(general) \ breezeVector(bData)
      val bd = det(breezeMatrix(general))
      val spdData = spd(n, n.toLong + 1)
      val cholTol = forwardTol(n, kappa(spdData))
      val bL = cholesky(breezeMatrix(spdData))
      val bLower = BDM.tabulate(n, n)((i, j) => if j <= i then bL(i, j) else 0.0)
      for ((label, view), (_, spdView)) <- galeViews(general).zip(galeViews(spdData)) do
        val clue = s"$label n=$n"
        val lu = view.lu.orThrow
        assertBelow(relDiff(lu.solve(galeVector(bData)).orThrow, bx), tol, s"LU solve $clue")
        assertBelow(relDiff(lu.solve(stridedVector(bData)).orThrow, bx), tol, s"LU solve strided RHS $clue")
        assertBelow(math.abs(lu.det.orThrow - bd) / math.abs(bd), tol, s"det $clue")
        val rhs = matrixData(n, 3, 4L)
        for (rhsLabel, rhsView) <- galeViews(rhs) do
          assertBelow(
            relDiff(view.solve(rhsView).orThrow, breezeMatrix(general) \ breezeMatrix(rhs)),
            tol,
            s"multi-RHS ($rhsLabel) $clue"
          )
        val factor = spdView.cholesky.orThrow
        val gLower = Matrix.tabulate(n, n)((i, j) => if j <= i then factor.lower(i, j) else 0.0)
        assertBelow(relDiff(gLower, bLower), cholTol, s"chol L $clue")
        assertBelow(relDiff(factor.solve(stridedVector(bData)).orThrow, breezeMatrix(spdData) \ breezeVector(bData)), cholTol, s"chol solve $clue")
        val gqr = view.qr
        val ba = breezeMatrix(general)
        assertBelow(relDiff(gqr.q * gqr.r, ba), c * n * Eps, s"QR recon $clue")
        assertBelow(relDiff(gqr.r.t * gqr.r, ba.t * ba), c * n * Eps, s"RᵀR $clue")
  }

  test("views: tall least squares on transposed and strided views") {
    for (m, n) <- List((9, 4), (60, 25), (260, 120)) do
      val data = matrixData(m, n, m.toLong)
      val k = kappa(data)
      val bData = vectorData(m, 5L)
      val bx = breezeMatrix(data) \ breezeVector(bData)
      for (label, view) <- galeViews(data) do
        assertBelow(relDiff(view.leastSquares(stridedVector(bData)).orThrow, bx), c * m * Eps * k * k, s"lstsq $label ${m}x$n")
  }
