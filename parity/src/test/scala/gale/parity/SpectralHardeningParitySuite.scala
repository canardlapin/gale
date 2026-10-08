package gale.parity

import breeze.linalg.DenseMatrix as BDM
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
  * '''Divergence.''' Breeze's `eigSym` throws `MatrixEmptyException` on `0 × 0`.
  * gale's behaviour on `0 × 0` is not asserted here: `Eigen.eigSymmetric` currently
  * throws `ArrayIndexOutOfBoundsException` instead of a typed result (reported to
  * the coordinator; not covered until its intended semantics are decided).
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
    for n <- List(96, 128, 160); seed <- List(1L, 2L, 3L) do
      // Three nonzero eigenvalues of multiplicity n/4 plus a distinct tail. (A
      // multiple eigenvalue AT ZERO currently makes gale's QL kernel exhaust its
      // sweep cap for n ≳ 96 — a reported kernel bug, not covered here yet.)
      val spectrum = Array.tabulate(n)(i => if i < 3 * n / 4 then 1.0 + (i % 3).toDouble else 4.0 + i.toDouble / n)
      val data = withSpectrum(spectrum, seed)
      assertEigParity(galeMatrix(data), data, s"repeated n=$n seed=$seed")
  }

  test("1x1 eigSym matches breeze") {
    for value <- List(-2.5, 0.0, 1e-300, 4e200) do
      val data = Array(Array(value))
      val gd = galeEig(galeMatrix(data))
      val es = eigSym(breezeMatrix(data))
      assertEquals(gd.eigenvalues(0), es.eigenvalues(0), s"λ [$value]")
      assertEquals(math.abs(gd.eigenvectors(0, 0)), 1.0, s"v [$value]")
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
