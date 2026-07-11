package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.qr
import breeze.linalg.svd as breezeSvd
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.spectral.*

/** Parity for the two additions that postdate the original harness: partial SVD
  * (`Svds.svd`, Golub–Kahan–Lanczos) versus Breeze's dense `svd` extremes, and
  * the blocked / compact-WY QR path (engaged when `min(m, n) > 96`) through the
  * same invariants the unblocked QR parity uses.
  *
  * Both libraries return singular values '''descending''', so values compare
  * elementwise. Singular vectors are defined up to sign for simple singular
  * values — compared via `|⟨u_gale, u_breeze⟩| ≈ 1` per column — and the random
  * rectangular fixtures here have well-separated spectra.
  */
class SvdQrParitySuite extends munit.FunSuite:

  private val svdTol   = 1e-7
  private val vecTol   = 1e-6
  private val qrTol    = 1e-9
  private val solveTol = 1e-9

  private def galeColumn(m: DMat, j: Int): IndexedSeq[Double] =
    (0 until m.rows).map(m(_, j))

  private def breezeColumn(b: BDM[Double], j: Int): IndexedSeq[Double] =
    (0 until b.rows).map(b(_, j))

  private def absDot(x: IndexedSeq[Double], y: IndexedSeq[Double]): Double =
    math.abs(x.zip(y).map(_ * _).sum)

  // ---------------------------------------------------------------------------
  // Partial SVD vs Breeze dense SVD extremes
  // ---------------------------------------------------------------------------

  test("partial SVD (Largest): top-k matches breeze svd on tall and wide inputs") {
    for (m, n, seed) <- Seq((30, 12, 21L), (12, 30, 22L), (40, 15, 23L)) do
      val data = matrixData(m, n, seed)
      val k = 4
      val g = Svds.svd(galeMatrix(data), SingularSelection.Count(k, SingularOrder.Largest)).orThrow
      assert(g.diagnostics.allConverged, s"not converged ${m}x$n seed=$seed: ${g.diagnostics}")
      val b = breezeSvd(breezeMatrix(data))
      // Both descending; gale returns k, breeze the full spectrum.
      var i = 0
      while i < k do
        assertScalarClose(g.singularValues(i), b.singularValues(i), svdTol, s"sigma($i) ${m}x$n seed=$seed")
        i += 1
      // Simple singular values: vectors match up to sign, per column.
      i = 0
      while i < k do
        val uDot = absDot(galeColumn(g.u, i), breezeColumn(b.leftVectors, i))
        assert(math.abs(uDot - 1.0) < vecTol, s"u($i) misaligned ($uDot) ${m}x$n seed=$seed")
        // gale stores Vᵀ rows; breeze rightVectors is Vᵀ as well.
        val gv = (0 until g.vt.cols).map(g.vt(i, _))
        val bv = (0 until b.rightVectors.cols).map(b.rightVectors(i, _))
        val vDot = absDot(gv, bv)
        assert(math.abs(vDot - 1.0) < vecTol, s"v($i) misaligned ($vDot) ${m}x$n seed=$seed")
        i += 1
  }

  test("partial SVD (Smallest): bottom-k values match breeze svd tail") {
    for (m, n, seed) <- Seq((20, 8, 31L), (8, 20, 32L)) do
      val data = matrixData(m, n, seed)
      val p = math.min(m, n)
      val k = 2
      val g = Svds.svd(galeMatrix(data), SingularSelection.Count(k, SingularOrder.Smallest)).orThrow
      assert(g.diagnostics.allConverged, s"not converged ${m}x$n seed=$seed: ${g.diagnostics}")
      val b = breezeSvd(breezeMatrix(data))
      // gale returns the k smallest, descending; breeze's tail is those values.
      var i = 0
      while i < k do
        assertScalarClose(g.singularValues(i), b.singularValues(p - k + i), svdTol, s"sigma tail($i) ${m}x$n seed=$seed")
        i += 1
  }

  // ---------------------------------------------------------------------------
  // Blocked / compact-WY QR (min(m, n) > 96 engages the panel path)
  // ---------------------------------------------------------------------------

  test("blocked QR: reconstruction, orthonormality, RtR=AtA on blocked shapes") {
    for (m, n, seed) <- Seq((120, 110, 41L), (200, 120, 42L), (110, 128, 43L)) do
      val data = matrixData(m, n, seed)
      val ga = galeMatrix(data)
      val ba = breezeMatrix(data)
      val gqr = ga.qr

      assertMatClose(gqr.q * gqr.r, ba, qrTol, s"blocked QR recon ${m}x$n seed=$seed")
      assertMatClose(gqr.q.t * gqr.q, BDM.eye[Double](m), qrTol, s"blocked Q orthonormal ${m}x$n seed=$seed")
      val ata = ba.t * ba
      assertMatClose(gqr.r.t * gqr.r, ata, qrTol * 100, s"blocked RtR=AtA ${m}x$n seed=$seed")
      // Cross-check R against breeze's (sign-free): both RᵀR equal AᵀA.
      val bqr = qr(ba)
      assertBreezeMatClose(bqr.r.t * bqr.r, ata, qrTol * 100, s"breeze RtR=AtA ${m}x$n seed=$seed")
  }

  test("blocked least squares: gale vs breeze backslash above the panel threshold") {
    for (m, n, seed) <- Seq((240, 120, 51L), (400, 100, 52L)) do
      val aData = matrixData(m, n, seed)
      val bData = vectorData(m, seed * 13 + 7)
      val gx = galeMatrix(aData).leastSquares(galeVector(bData)).orThrow
      val bx = breezeMatrix(aData) \ breezeVector(bData)
      assertVecClose(gx, bx, solveTol, s"blocked lstsq ${m}x$n seed=$seed")
  }
