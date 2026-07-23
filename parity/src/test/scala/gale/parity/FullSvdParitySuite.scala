package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.kron as breezeKron
import breeze.linalg.pinv as breezePinv
import breeze.linalg.svd as breezeSvd
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.spectral.*

/** Parity for the full/economy dense SVD (`Svds.svd` with
  * `SingularSelection.All`, the Golub–Kahan–Reinsch bidiagonal kernel), the
  * Moore–Penrose `pinv`, and the Kronecker `kron` against Breeze.
  *
  * Both libraries return singular values '''descending''', so values compare
  * elementwise across the whole spectrum. Breeze's `svd` returns the '''full'''
  * square `U` (m×m) / `Vᵀ` (n×n) while gale returns the economy factors, so
  * vectors compare over the first `min(m, n)` columns/rows — up to sign, via
  * `|⟨·,·⟩| ≈ 1` per column, exactly as `SvdQrParitySuite` does, and only for
  * singular values whose relative gap to their neighbours is resolvable (a
  * clustered pair's individual vectors are not comparable across libraries).
  * `pinv` compares elementwise on full-rank fixtures (both libraries' default
  * cutoffs agree that nothing is truncated); rank-deficient `pinv` behaviour is
  * pinned by the core Moore–Penrose suite instead, keeping cutoff-policy
  * differences out of the parity claim. `kron` is plain arithmetic and
  * compares essentially exactly.
  */
class FullSvdParitySuite extends munit.FunSuite:

  private val svdTol  = 1e-9
  private val vecTol  = 1e-7
  private val pinvTol = 1e-8
  private val kronTol = 1e-12

  private def galeColumn(m: DMat, j: Int): IndexedSeq[Double] =
    (0 until m.rows).map(m(_, j))

  private def breezeColumn(b: BDM[Double], j: Int): IndexedSeq[Double] =
    (0 until b.rows).map(b(_, j))

  private def absDot(x: IndexedSeq[Double], y: IndexedSeq[Double]): Double =
    math.abs(x.zip(y).map(_ * _).sum)

  /** Relative gap of σ_i to its nearest neighbour, for vector comparability. */
  private def relativeGap(s: IndexedSeq[Double], i: Int): Double =
    val sigmaMax = math.max(s.head, 1e-300)
    val up = if i > 0 then (s(i - 1) - s(i)) / sigmaMax else 1.0
    val down = if i < s.length - 1 then (s(i) - s(i + 1)) / sigmaMax else 1.0
    math.min(up, down)

  // ---------------------------------------------------------------------------
  // Full SVD vs Breeze dense SVD
  // ---------------------------------------------------------------------------

  test("full SVD: all singular values match breeze on tall, wide, and square inputs") {
    for (m, n, seed) <- Seq((30, 12, 121L), (12, 30, 122L), (25, 25, 123L), (40, 15, 124L)) do
      val data = matrixData(m, n, seed)
      val p = math.min(m, n)
      val g = galeMatrix(data).svd.orThrow
      assert(g.diagnostics.allConverged, s"not converged ${m}x$n seed=$seed: ${g.diagnostics}")
      val b = breezeSvd(breezeMatrix(data))
      assertEquals(g.size, p, s"economy count ${m}x$n seed=$seed")
      var i = 0
      while i < p do
        assertScalarClose(g.singularValues(i), b.singularValues(i), svdTol, s"sigma($i) ${m}x$n seed=$seed")
        i += 1
  }

  test("full SVD: singular vectors match breeze up to sign on resolvable values") {
    for (m, n, seed) <- Seq((30, 12, 121L), (12, 30, 122L), (25, 25, 123L)) do
      val data = matrixData(m, n, seed)
      val p = math.min(m, n)
      val g = galeMatrix(data).svd.orThrow
      val b = breezeSvd(breezeMatrix(data))
      val sigmas = (0 until p).map(g.singularValues(_))
      var i = 0
      while i < p do
        // Vectors are only individually comparable when σ_i is separated from
        // its neighbours; the random fixtures resolve essentially all values.
        if relativeGap(sigmas, i) > 1e-3 then
          val uDot = absDot(galeColumn(g.u, i), breezeColumn(b.leftVectors, i))
          assert(math.abs(uDot - 1.0) < vecTol, s"u($i) misaligned ($uDot) ${m}x$n seed=$seed")
          // gale stores Vᵀ rows; breeze rightVectors is Vᵀ as well.
          val gv = (0 until g.vt.cols).map(g.vt(i, _))
          val bv = (0 until b.rightVectors.cols).map(b.rightVectors(i, _))
          val vDot = absDot(gv, bv)
          assert(math.abs(vDot - 1.0) < vecTol, s"v($i) misaligned ($vDot) ${m}x$n seed=$seed")
        i += 1
  }

  // ---------------------------------------------------------------------------
  // pinv vs breeze pinv
  // ---------------------------------------------------------------------------

  test("pinv: full-rank tall, wide, and square match breeze elementwise") {
    for (m, n, seed) <- Seq((10, 6, 131L), (6, 10, 132L), (8, 8, 133L)) do
      val data = matrixData(m, n, seed)
      val g = galeMatrix(data).pinv.orThrow
      val b = breezePinv(breezeMatrix(data))
      assertMatClose(g, b, pinvTol, s"pinv ${m}x$n seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // kron vs breeze kron
  // ---------------------------------------------------------------------------

  test("kron: matches breeze on rectangular operands") {
    for (m1, n1, m2, n2, seed) <- Seq((3, 2, 2, 4, 141L), (2, 5, 4, 3, 142L), (1, 4, 3, 1, 143L)) do
      val aData = matrixData(m1, n1, seed)
      val bData = matrixData(m2, n2, seed * 7 + 1)
      val g = galeMatrix(aData).kron(galeMatrix(bData))
      val b = breezeKron(breezeMatrix(aData), breezeMatrix(bData))
      assertMatClose(g, b, kronTol, s"kron ${m1}x$n1 (x) ${m2}x$n2 seed=$seed")
  }
