package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.MatrixEmptyException
import breeze.linalg.eigSym
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.spectral.*

/** Dense symmetric eigensolver hardening versus `breeze.linalg.eigSym` (W5.3).
  *
  * Covers sizes above the planned blocked-tridiagonalization threshold
  * (`n ≳ 64`), `1 × 1`, ill-conditioned spectra, and non-contiguous views.
  *
  * '''Tolerances.''' Symmetric eigenvalues are perfectly conditioned in the
  * absolute sense (Weyl): a backward-stable solver satisfies
  * `|λ̂ᵢ − λᵢ| ≲ γₙ ‖A‖₂`, so both libraries are compared with the absolute bound
  * `c · n · ε · ‖A‖₂` (`c = 32`) — '''independent of κ'''. Eigenvectors are only
  * determined up to the spectral gap (Davis–Kahan: `sin θ ≲ γₙ ‖A‖₂ / sep`), so
  * the spectrum is split into clusters separated by more than `√ε ‖A‖₂`, and each
  * cluster's invariant subspace is compared with
  * `‖G − B Bᵀ G‖_F ≤ √k · c · n · ε · ‖A‖₂ / sep` (`k` the cluster size) — sign- and rotation-invariant.
  *
  * '''Divergence.''' On `0 × 0` gale returns an empty decomposition (the
  * `lu`/`qr`/`cholesky` precedent) while Breeze's `eigSym` throws
  * `MatrixEmptyException`.
  */
class SpectralHardeningParitySuite extends munit.FunSuite:

  private val c = 32.0

  private def galeEig(a: DMat): EigenDecomposition =
    Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.Right).orThrow

  /** Compare gale's full decomposition of `view` with Breeze's of `data`. */
  private def assertEigParity(view: DMat, data: Array[Array[Double]], clue: String): Unit =
    val n = data.length
    val gd = galeEig(view)
    val es = eigSym(breezeMatrix(data))
    val bValues = (0 until n).map(es.eigenvalues(_))
    val norm2 = bValues.map(math.abs).max
    val valueTol = c * n * Eps * norm2
    assertEquals(gd.size, n, clue)
    var i = 0
    while i < n do
      assertBelow(math.abs(gd.eigenvalues(i) - bValues(i)), valueTol, s"$clue λ[$i]")
      i += 1
    for (cluster, sep) <- spectralClusters(bValues, math.sqrt(Eps) * norm2) do
      // Frobenius Davis–Kahan: ‖sin Θ‖_F ≤ √k · ‖E‖₂ / sep; a cluster spanning the
      // whole spectrum (sep = ∞) is the full space, where only roundoff remains.
      val bound = math.sqrt(cluster.length.toDouble) * c * n * Eps * (if sep.isInfinite then 1.0 else norm2 / sep)
      assertBelow(subspaceSine(gd.eigenvectors, es.eigenvectors, cluster), bound, s"$clue subspace $cluster sep=$sep")
    // Independent of Breeze: gale's own residual ‖A V − V Λ‖ is backward-stable small.
    val residual = galeMatrix(data) * gd.eigenvectors - gd.eigenvectors * Matrix.tabulate(n, n)((r, k) =>
      if r == k then gd.eigenvalues(r) else 0.0
    )
    var worst = 0.0
    for r <- 0 until n; k <- 0 until n do worst = math.max(worst, math.abs(residual(r, k)))
    assertBelow(worst, valueTol, s"$clue residual")

  test("large n eigSym: values (Weyl bound) and cluster subspaces (Davis–Kahan bound)") {
    for (n, seed) <- List(96, 128, 200, 300).zip(List(71L, 72L, 73L, 74L)) do
      val data = symmetric(n, seed)
      assertEigParity(galeMatrix(data), data, s"symmetric n=$n")
  }

  test("large n eigSym: repeated eigenvalues above the blocking threshold") {
    for (n, seed) <- List(96 -> 1L, 128 -> 2L, 160 -> 3L) do
      // Three nonzero eigenvalues of multiplicity n/4 plus a distinct tail.
      val spectrum = Array.tabulate(n)(i => if i < 3 * n / 4 then 1.0 + (i % 3).toDouble else 4.0 + i.toDouble / n)
      val data = withSpectrum(spectrum, seed)
      assertEigParity(galeMatrix(data), data, s"repeated n=$n seed=$seed")
  }

  test("repeated eigenvalues AT ZERO above the threshold (QL deflation regression)") {
    // Each case exhausted the 30-sweep cap (Left(DidNotConverge)) while the QL
    // deflation test was purely local (bd-01M4EKP8CWDHZ2W1WFBRQ64Z9X).
    def twoLevel(n: Int) = Array.tabulate(n)(i => if i < n / 2 then 0.0 else 1.0)
    def threeLevel(n: Int) = Array.tabulate(n)(i => (i % 3).toDouble)
    def withTail(n: Int) = Array.tabulate(n)(i => if i < 3 * n / 4 then (i % 3).toDouble else 3.0 + i.toDouble / n)
    val cases = List(
      ("two-level", twoLevel(96), 5L),
      ("three-level", threeLevel(96), 5L),
      ("three-level+tail", withTail(128), 128L),
      ("three-level+tail", withTail(128), 2L),
      ("three-level+tail", withTail(160), 160L),
      ("three-level+tail", withTail(160), 3L)
    )
    for (label, spectrum, seed) <- cases do
      val data = withSpectrum(spectrum, seed)
      assertEigParity(galeMatrix(data), data, s"$label n=${spectrum.length} seed=$seed")
  }

  test("n = 128 values-only and workspace routes: exact vs the ordinary route, Weyl bound vs breeze") {
    val n = 128
    val data = symmetric(n, 75L)
    val a = galeMatrix(data)
    val bValues = (0 until n).map(eigSym(breezeMatrix(data)).eigenvalues(_))
    val valueTol = c * n * Eps * bValues.map(math.abs).max
    def values(d: EigenDecomposition) = (0 until d.size).map(d.eigenvalues(_))
    val ordinaryValues = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).orThrow
    val ordinaryVectors = galeEig(a)
    val workspace = DenseWorkspace.empty
    val workspaceValues = Eigen.eigSymmetricWith(a, EigenSelection.All, EigenVectors.ValuesOnly, workspace).orThrow
    val workspaceVectors = Eigen.eigSymmetricWith(a, EigenSelection.All, workspace).orThrow
    // Ordinary and workspace routes share one kernel per mode (bit-pinning audit,
    // W1.5), so they agree exactly; values-only vs vectors may differ by ulps.
    assertEquals(values(workspaceValues), values(ordinaryValues), "values-only: workspace vs ordinary")
    assertEquals(values(workspaceVectors), values(ordinaryVectors), "vectors: workspace vs ordinary values")
    for r <- 0 until n; k <- 0 until n do
      assertEquals(workspaceVectors.eigenvectors(r, k), ordinaryVectors.eigenvectors(r, k), s"vectors: V($r,$k)")
    for (label, d) <- List("values-only" -> ordinaryValues, "vectors" -> ordinaryVectors) do
      for i <- 0 until n do assertBelow(math.abs(values(d)(i) - bValues(i)), valueTol, s"$label λ[$i]")
  }

  test("1x1 eigSym matches breeze") {
    for value <- List(-2.5, 0.0, 1e-300, 4e200) do
      val data = Array(Array(value))
      val gd = galeEig(galeMatrix(data))
      val es = eigSym(breezeMatrix(data))
      assertEquals(gd.eigenvalues(0), es.eigenvalues(0), s"λ [$value]")
      assertEquals(math.abs(gd.eigenvectors(0, 0)), 1.0, s"v [$value]")
  }

  test("empty 0x0: gale returns an empty decomposition, breeze throws MatrixEmptyException") {
    val gd = galeEig(Matrix.zeros(0, 0))
    assertEquals((gd.size, gd.eigenvectors.rows, gd.eigenvectors.cols), (0, 0, 0))
    intercept[MatrixEmptyException](eigSym(BDM.zeros[Double](0, 0)))
  }

  test("ill-conditioned spectra: κ = 1e12 SPD, indefinite, and Hilbert 6..12") {
    for n <- List(16, 64) do
      val spd12 = withSpectrum(geometricSpectrum(n, 1e12), n.toLong + 3)
      assertEigParity(galeMatrix(spd12), spd12, s"geometric κ=1e12 n=$n")
      val signed = geometricSpectrum(n, 1e12).zipWithIndex.map((v, i) => if i % 2 == 0 then v else -v)
      val indefinite = withSpectrum(signed, n.toLong + 4)
      assertEigParity(galeMatrix(indefinite), indefinite, s"indefinite κ=1e12 n=$n")
    for n <- 6 to 12 do
      val h = hilbert(n)
      assertEigParity(galeMatrix(h), h, s"Hilbert n=$n")
  }

  test("views: eigSym on transposed and strided views matches breeze on the copy") {
    for n <- List(6, 33, 100) do
      val data = symmetric(n, n.toLong + 9)
      for (label, view) <- galeViews(data) do assertEigParity(view, data, s"$label n=$n")
  }
