package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.eigSym
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.spectral.*

/** Symmetric spectral parity versus `breeze.linalg.eigSym`.
  *
  * Both gale (tridiagonal QL/QR) and Breeze (LAPACK `dsyev`, lower triangle) return
  * eigenvalues '''ascending''', so eigenvalues are compared elementwise. Eigenvectors
  * are only defined up to sign (simple eigenvalues) or up to an orthogonal rotation
  * (repeated/clustered eigenvalues), so vectors are compared sign-aware for well-
  * separated eigenvalues and via the '''subspace projector''' `‖V Vᵀ − W Wᵀ‖` for
  * clusters. The iterative Lanczos path is checked against Breeze's dense extremes.
  */
class SpectralParitySuite extends munit.FunSuite:

  private val eigTol   = 1e-9
  private val vecTol   = 1e-8
  private val projTol  = 1e-8
  private val lanczTol = 1e-7

  // ---------------------------------------------------------------------------
  // Helpers
  // ---------------------------------------------------------------------------

  private def galeEig(data: Array[Array[Double]]): EigenDecomposition =
    Eigen.eigSymmetric(galeMatrix(data), EigenSelection.All, EigenVectors.Right).orThrow

  private def galeValues(d: EigenDecomposition): IndexedSeq[Double] =
    (0 until d.size).map(d.eigenvalues(_))

  private def breezeValues(es: eigSym.DenseEigSym): IndexedSeq[Double] =
    (0 until es.eigenvalues.length).map(es.eigenvalues(_))

  /** Partition ascending eigenvalue indices into clusters where consecutive gaps
    * fall at or below `gap` (a repeated/degenerate group is one cluster).
    */
  private def clusters(values: IndexedSeq[Double], gap: Double): List[Range] =
    val out = List.newBuilder[Range]
    var start = 0
    var i = 1
    while i < values.length do
      if values(i) - values(i - 1) > gap then
        out += (start until i)
        start = i
      i += 1
    out += (start until values.length)
    out.result()

  /** `‖V_c V_cᵀ − W_c W_cᵀ‖_∞` over the cluster columns `cols`: an orthogonal-
    * rotation-invariant (and sign-invariant) measure that gale's and Breeze's
    * eigenvectors span the same subspace.
    */
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

  /** Flip a column's sign so its largest-magnitude entry is positive. */
  private def signNormalized(col: IndexedSeq[Double]): IndexedSeq[Double] =
    var idx = 0
    var i = 1
    while i < col.length do
      if math.abs(col(i)) > math.abs(col(idx)) then idx = i
      i += 1
    if col(idx) < 0.0 then col.map(-_) else col

  // ---------------------------------------------------------------------------
  // Dense: All selection
  // ---------------------------------------------------------------------------

  test("dense eigSym All: eigenvalues match, eigenvectors span the same subspaces") {
    for n <- List(5, 12, 24); seed <- List(1L, 2L, 3L) do
      val data = symmetric(n, seed)
      val gd   = galeEig(data)
      val es   = eigSym(breezeMatrix(data))

      // Eigenvalues, ascending, elementwise.
      val gv = galeValues(gd)
      val bv = breezeValues(es)
      assertEquals(gv.length, bv.length)
      gv.zip(bv).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, eigTol, s"eigval n=$n seed=$seed [$i]")
      }

      // Eigenvectors: subspace parity per detected cluster (robust to sign and to
      // rotation within near-degenerate groups).
      for cluster <- clusters(gv, 1e-6) do
        val d = projectorDiff(gd.eigenvectors, es.eigenvectors, cluster)
        assert(d < projTol, s"subspace mismatch n=$n seed=$seed cluster=$cluster diff=$d")
  }

  // ---------------------------------------------------------------------------
  // Dense: Count selection
  // ---------------------------------------------------------------------------

  test("dense eigSym Count: top-k / bottom-k are breeze's extreme k") {
    for n <- List(10, 20); seed <- List(1L, 2L) do
      val data = symmetric(n, seed)
      val k    = 4
      val full = breezeValues(eigSym(breezeMatrix(data))) // ascending

      val top = Eigen.eigSymmetric(galeMatrix(data), EigenSelection.Count(k, EigenOrder.LargestAlgebraic), EigenVectors.ValuesOnly).orThrow
      val bot = Eigen.eigSymmetric(galeMatrix(data), EigenSelection.Count(k, EigenOrder.SmallestAlgebraic), EigenVectors.ValuesOnly).orThrow

      galeValues(top).zip(full.takeRight(k)).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, eigTol, s"top-k n=$n seed=$seed [$i]")
      }
      galeValues(bot).zip(full.take(k)).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, eigTol, s"bottom-k n=$n seed=$seed [$i]")
      }
  }

  // ---------------------------------------------------------------------------
  // Dense: fixed spectra
  // ---------------------------------------------------------------------------

  test("dense eigSym: distinct well-separated spectrum, sign-aware eigenvectors") {
    val n        = 8
    val spectrum = Array.tabulate(n)(i => math.pow(2.0, i.toDouble)) // 1,2,4,...,128
    for seed <- List(1L, 2L, 3L) do
      val data = withSpectrum(spectrum, seed)
      val gd   = galeEig(data)
      val es   = eigSym(breezeMatrix(data))

      // Recovered eigenvalues equal the prescribed (ascending) spectrum.
      galeValues(gd).zip(spectrum.toIndexedSeq).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, eigTol, s"fixed spectrum gale [$i] seed=$seed")
      }
      breezeValues(es).zip(spectrum.toIndexedSeq).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, eigTol, s"fixed spectrum breeze [$i] seed=$seed")
      }

      // Each eigenvector is simple and well-separated: sign-normalize and compare.
      var c = 0
      while c < n do
        val g = signNormalized((0 until n).map(gd.eigenvectors(_, c)))
        val b = signNormalized((0 until n).map(es.eigenvectors(_, c)))
        var i = 0
        while i < n do
          assertScalarClose(g(i), b(i), vecTol, s"eigvec seed=$seed col=$c row=$i")
          i += 1
        c += 1
  }

  test("dense eigSym: repeated eigenvalue, subspace parity") {
    val spectrum = Array(1.0, 1.0, 1.0, 5.0, 9.0) // triple at 1.0
    val n        = spectrum.length
    for seed <- List(1L, 2L, 3L) do
      val data = withSpectrum(spectrum, seed)
      val gd   = galeEig(data)
      val es   = eigSym(breezeMatrix(data))

      galeValues(gd).zip(spectrum.toIndexedSeq).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, eigTol, s"repeated spectrum [$i] seed=$seed")
      }
      // Individual columns in the triple are not comparable (arbitrary rotation),
      // but the degenerate subspace must agree.
      val d = projectorDiff(gd.eigenvectors, es.eigenvectors, 0 until 3)
      assert(d < projTol, s"degenerate subspace mismatch seed=$seed diff=$d")
  }

  // ---------------------------------------------------------------------------
  // Iterative: Lanczos vs breeze extremes
  // ---------------------------------------------------------------------------

  test("Lanczos top-k / bottom-k vs breeze eigSym extremes") {
    val n        = 30
    val spectrum = Array.tabulate(n)(i => (i + 1).toDouble) // 1..30, well separated
    val k        = 5
    for seed <- List(1L, 2L) do
      val data = withSpectrum(spectrum, seed)
      val op   = galeMatrix(data).asLinearOperator
      val full = breezeValues(eigSym(breezeMatrix(data))) // ascending

      val top = Eigen.eigSymmetric(op, n, EigenSelection.Count(k, EigenOrder.LargestAlgebraic)).orThrow
      val bot = Eigen.eigSymmetric(op, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic)).orThrow

      assert(top.diagnostics.allConverged, s"top-k Lanczos did not converge seed=$seed")
      assert(bot.diagnostics.allConverged, s"bottom-k Lanczos did not converge seed=$seed")

      galeValues(top).zip(full.takeRight(k)).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, lanczTol, s"Lanczos top-k seed=$seed [$i]")
      }
      galeValues(bot).zip(full.take(k)).zipWithIndex.foreach { case ((x, y), i) =>
        assertScalarClose(x, y, lanczTol, s"Lanczos bottom-k seed=$seed [$i]")
      }
  }
