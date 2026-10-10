package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.kron as breezeKron
import breeze.linalg.pinv as breezePinv
import breeze.linalg.svd as breezeSvd
import gale.linalg.*
import gale.golden.GoldenTolerance.{svdFactorRatios, svdProjectorRatio}
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

  /** A deterministic Householder basis, independent of either library's QR. */
  private def householder(n: Int, offset: Int): Array[Array[Double]] =
    val v = Array.tabulate(n)(i => 1.0 + ((i + offset) % 7).toDouble)
    val norm = math.sqrt(v.map(x => x * x).sum)
    val q = v.map(_ / norm)
    Array.tabulate(n, n)((i, j) => (if i == j then 1.0 else 0.0) - 2.0 * q(i) * q(j))

  private def planted(m: Int, n: Int, sigma: Array[Double]): Array[Array[Double]] =
    val u = householder(m, 1)
    val v = householder(n, 3)
    Array.tabulate(m, n)((i, j) => sigma.indices.map(k => u(i)(k) * sigma(k) * v(j)(k)).sum)

  private def checkFullFactors(data: Array[Array[Double]], clue: String): Unit =
    val m = data.length
    val n = data(0).length
    val p = math.min(m, n)
    val b = breezeSvd(breezeMatrix(data))
    for (name, units) <- svdFactorRatios(data(_)(_), m, n, b.singularValues.toArray, (i, j) => b.leftVectors(i, j), (i, j) => b.rightVectors(i, j)) do
      assert(units <= 32.0, s"$clue Breeze $name: $units dimension*eps units")
    val g = galeMatrix(data).svd.fold(e => fail(s"$clue Gale SVD: $e"), identity)
    assert(g.diagnostics.allConverged, clue)
    assertEquals((g.size, g.u.rows, g.u.cols, g.vt.rows, g.vt.cols), (p, m, p, p, n))
    val unit = math.max(m, n) * Eps * b.singularValues(0)
    for k <- 0 until p do assertBelow(math.abs(g.singularValues(k) - b.singularValues(k)), 32.0 * unit, s"$clue sigma[$k]")
    for (label, s, u, vt) <- Seq(
      ("gale", g.singularValues.toArray, g.u.apply, g.vt.apply),
      ("breeze", b.singularValues.toArray, (i: Int, j: Int) => b.leftVectors(i, j), (i: Int, j: Int) => b.rightVectors(i, j))
    ) do
      for (name, units) <- svdFactorRatios(data(_)(_), m, n, s, u, vt) do
        assert(units <= 32.0, s"$clue $label $name: $units dimension*eps units")
    val s = b.singularValues.toArray
    assert(svdProjectorRatio(g.u.apply, (i, j) => b.leftVectors(i, j), s, m, math.max(m, n)) <= 32.0, s"$clue U subspaces")
    assert(svdProjectorRatio((i, j) => g.vt(j, i), (i, j) => b.rightVectors(j, i), s, n, math.max(m, n)) <= 32.0, s"$clue V subspaces")

  test("benchmark-size full SVD: both libraries reconstruct and have orthonormal factors") {
    for (m, n) <- Seq((128, 128), (256, 256), (192, 128), (128, 192)) do
      checkFullFactors(matrixData(m, n, 801L + m + n), s"${m}x$n")
  }

  for scale <- Seq(1e-200, 1.0, 1e200); (m, n) <- Seq((12, 8), (8, 12)) do
    test(s"clustered graded full SVD ${m}x$n scale=$scale retains factors") {
      val sigma = Array(8.0, 3.0, 3.0, 1.0, 1e-3, 1e-6, 1e-9, 1e-12)
      checkFullFactors(planted(m, n, sigma.map(_ * scale)), s"graded ${m}x$n scale=$scale")
    }

  test("rank-deficient pinv matches Breeze SVD with the same explicit cutoff") {
    // Breeze pinv's exact-zero cutoff differs. Build its reference from its SVD,
    // retaining exactly sigma > 1e-10*sigmaMax, the same policy passed to Gale.
    for (m, n) <- Seq((128, 64), (64, 128)) do
      val data = planted(m, n, Array(8.0, 3.0, 3.0, 1.0, 0.25, 1e-14))
      checkFullFactors(data, s"rank5 ${m}x$n")
      val g = galeMatrix(data).pinv(SvdCutoff.Relative(1e-10)).orThrow
      val b = breezeSvd(breezeMatrix(data))
      val retained = (0 until math.min(m, n)).filter(k => b.singularValues(k) > 1e-10 * b.singularValues(0))
      assertEquals(retained.size, 5)
      val expected = BDM.tabulate(n, m)((i, j) => retained.map(k => b.rightVectors(k, i) * b.leftVectors(j, k) / b.singularValues(k)).sum)
      assertMatClose(g, expected, 32.0 * math.max(m, n) * Eps * 32.0, s"explicit-cutoff ${m}x$n")
      val a = galeMatrix(data)
      val ax = a * g
      val xa = g * a
      for (left, right, label) <- Seq((ax * a, a, "AXA"), (xa * g, g, "XAX"), (ax, ax.t, "AX symmetric"), (xa, xa.t, "XA symmetric")) do
        var error = 0.0
        var norm = 0.0
        for i <- 0 until left.rows; j <- 0 until left.cols do
          error = math.hypot(error, left(i, j) - right(i, j))
          norm = math.hypot(norm, right(i, j))
        assert(error / norm <= 32.0 * math.max(m, n) * Eps * 32.0, s"${m}x$n $label: ${error / norm}")
  }

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

  /** `‖V_c V_cᵀ − W_c W_cᵀ‖_∞` over cluster columns. */
  private def projectorDiff(gv: DMat, bv: BDM[Double], cols: Range): Double =
    val n = gv.rows
    var worst = 0.0
    var i = 0
    while i < n do
      var j = 0
      while j < n do
        var pg = 0.0
        var pb = 0.0
        for c <- cols do
          pg += gv(i, c) * gv(j, c)
          pb += bv(i, c) * bv(j, c)
        worst = math.max(worst, math.abs(pg - pb))
        j += 1
      i += 1
    worst

  test("full SVD: clustered singular vectors match breeze as subspace projectors") {
    // Plant A = U Σ Vᵀ with a repeated σ = 3 so individual vectors are not
    // comparable across libraries; the 2-D projectors must still agree.
    val rngU = new scala.util.Random(301L)
    val rngV = new scala.util.Random(302L)
    val m = 10
    val n = 6
    val uRaw = Array.tabulate(m, m)((_, _) => rngU.nextDouble() * 2.0 - 1.0)
    val vRaw = Array.tabulate(n, n)((_, _) => rngV.nextDouble() * 2.0 - 1.0)
    val uQ = galeMatrix(uRaw).qr.q
    val vQ = galeMatrix(vRaw).qr.q
    val sigmas = IndexedSeq(7.0, 3.0, 3.0, 1.5, 0.4)
    val planted = Array.tabulate(m, n): (i, j) =>
      var sum = 0.0
      var k = 0
      while k < sigmas.length do
        sum += uQ(i, k) * sigmas(k) * vQ(j, k)
        k += 1
      sum
    val g = galeMatrix(planted).svd.orThrow
    val b = breezeSvd(breezeMatrix(planted))
    val cluster = 1 until 3
    val uDiff = projectorDiff(g.u, b.leftVectors, cluster)
    assert(uDiff < 1e-8, s"U cluster projector mismatch $uDiff")
    // gale stores Vᵀ; breeze rightVectors is also Vᵀ — compare V = (Vᵀ)ᵀ.
    val gV = Matrix.tabulate(g.vt.cols, g.vt.rows)((i, j) => g.vt(j, i))
    val bV = BDM.tabulate(b.rightVectors.cols, b.rightVectors.rows)((i, j) => b.rightVectors(j, i))
    val vDiff = projectorDiff(gV, bV, cluster)
    assert(vDiff < 1e-8, s"V cluster projector mismatch $vDiff")
    var i = 0
    while i < sigmas.length do
      assertScalarClose(g.singularValues(i), b.singularValues(i), svdTol, s"clustered sigma($i)")
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

  test("pinv: gale pinv(A)*b matches breeze A \\ b on full-rank square and tall systems") {
    for (m, n, seed) <- Seq((8, 8, 134L), (12, 5, 135L), (20, 7, 136L)) do
      val aData = matrixData(m, n, seed)
      val bData = vectorData(m, seed * 11 + 3)
      val gx = (galeMatrix(aData).pinv.orThrow * galeVector(bData))
      val bx = breezeMatrix(aData) \ breezeVector(bData)
      assertVecClose(gx, bx, pinvTol, s"pinv*b ${m}x$n seed=$seed")
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
