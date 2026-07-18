package gale.spectral

import gale.linalg.DMat
import gale.linalg.LinAlgError
import gale.linalg.Matrix

/** Tests for the full/economy dense SVD — `Svds.svd` with
  * `SingularSelection.All` (and the now-served `Count(k = min(m, n))` edge),
  * backed by the Householder-bidiagonalization + Golub–Kahan–Reinsch kernel
  * (`DenseSvdKernel`). Checks known spectra, reconstruction, orthonormality,
  * the descending non-negative layout, rank on deficient inputs, degenerate
  * shapes, and agreement with the partial Golub–Kahan–Lanczos path.
  */
class DenseSvdSuite extends munit.FunSuite:

  // --- helpers ---------------------------------------------------------------

  private def rectDiag(m: Int, n: Int, d: Seq[Double]): DMat =
    Matrix.tabulate(m, n)((i, j) => if i == j && i < d.length then d(i) else 0.0)

  private def randomMat(m: Int, n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(m, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  private def values(s: SVD): IndexedSeq[Double] =
    (0 until s.size).map(s.singularValues(_))

  private def frob(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  private def orthoError(cols: DMat): Double =
    frob(cols.t * cols - Matrix.eye(cols.cols))

  /** `‖A − U·diag(σ)·Vᵀ‖_F` — the full reconstruction error. */
  private def reconstructionError(a: DMat, s: SVD): Double =
    val p = s.size
    val d = Matrix.tabulate(p, p)((i, j) => if i == j then s.singularValues(i) else 0.0)
    frob(a - s.u * d * s.vt)

  private def gramSpectrum(a: DMat): IndexedSeq[Double] =
    val gram = a.t * a
    val eig = Eigen.eigSymmetric(gram, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    (0 until eig.size)
      .map(i => math.sqrt(math.max(eig.eigenvalues(i), 0.0)))
      .reverse
      .take(math.min(a.rows, a.cols))

  private def assertDescendingNonNegative(v: IndexedSeq[Double]): Unit =
    var i = 0
    while i < v.length do
      assert(v(i) >= 0.0, s"negative sigma at $i: $v")
      if i > 0 then assert(v(i) <= v(i - 1) + 1e-12, s"not descending at $i: $v")
      i += 1

  private def assertClose(a: IndexedSeq[Double], b: IndexedSeq[Double], tol: Double): Unit =
    assertEquals(a.length, b.length, s"length $a vs $b")
    a.zip(b).foreach { case (x, y) => assert(math.abs(x - y) < tol, s"$x != $y in $a vs $b") }

  private def assertEconomy(s: SVD, m: Int, n: Int): Unit =
    val p = math.min(m, n)
    assertEquals(s.size, p)
    assertEquals(s.u.rows, m)
    assertEquals(s.u.cols, p)
    assertEquals(s.vt.rows, p)
    assertEquals(s.vt.cols, n)

  private def fullChecks(a: DMat, s: SVD, tol: Double): Unit =
    assertEconomy(s, a.rows, a.cols)
    assert(s.diagnostics.allConverged, s"${s.diagnostics}")
    assertDescendingNonNegative(values(s))
    assert(reconstructionError(a, s) < tol, s"reconstruction ${reconstructionError(a, s)}")
    assert(orthoError(s.u) < tol, s"U ortho ${orthoError(s.u)}")
    assert(orthoError(s.vt.t) < tol, s"V ortho ${orthoError(s.vt.t)}")

  // --- known spectra ---------------------------------------------------------

  test("full svd tall diagonal: known values, economy shapes, reconstruction") {
    val a = rectDiag(8, 5, Seq(9.0, 7.0, 5.0, 3.0, 1.0))
    val s = Svds.svd(a, SingularSelection.All).toOption.get
    fullChecks(a, s, 1e-10)
    assertClose(values(s), IndexedSeq(9.0, 7.0, 5.0, 3.0, 1.0), 1e-10)
    assertEquals(s.rank, 5)
    assertEquals(s.diagnostics.rank, Some(5))
    assertEquals(s.diagnostics.iterations, 0)
  }

  test("full svd wide diagonal: known values via the transpose orientation") {
    val a = rectDiag(4, 9, Seq(6.0, 4.0, 2.0, 1.0))
    val s = Svds.svd(a, SingularSelection.All).toOption.get
    fullChecks(a, s, 1e-10)
    assertClose(values(s), IndexedSeq(6.0, 4.0, 2.0, 1.0), 1e-10)
  }

  test("full svd of a negated diagonal keeps sigma non-negative") {
    val a = rectDiag(5, 3, Seq(-4.0, -2.0, -1.0))
    val s = Svds.svd(a, SingularSelection.All).toOption.get
    fullChecks(a, s, 1e-10)
    assertClose(values(s), IndexedSeq(4.0, 2.0, 1.0), 1e-10)
  }

  // --- random matrices -------------------------------------------------------

  test("full svd random tall: reconstruction, orthonormality, Gram spectrum") {
    val a = randomMat(20, 8, 31L)
    val s = Svds.svd(a, SingularSelection.All).toOption.get
    fullChecks(a, s, 1e-9)
    assertClose(values(s), gramSpectrum(a), 1e-7)
    assert(s.diagnostics.worstResidual < 1e-9, s"${s.diagnostics}")
  }

  test("full svd random wide: reconstruction, orthonormality, Gram spectrum") {
    val a = randomMat(7, 18, 42L)
    val s = Svds.svd(a, SingularSelection.All).toOption.get
    fullChecks(a, s, 1e-9)
    assertClose(values(s), gramSpectrum(a), 1e-7)
  }

  test("full svd random square") {
    val a = randomMat(12, 12, 57L)
    val s = Svds.svd(a, SingularSelection.All).toOption.get
    fullChecks(a, s, 1e-9)
    assertClose(values(s), gramSpectrum(a), 1e-7)
  }

  // --- agreement with the partial (Golub–Kahan–Lanczos) path -----------------

  test("full svd agrees with the partial path's top-k on deterministic fixtures") {
    for (m, n, seed) <- Seq((16, 9, 101L), (9, 16, 102L)) do
      val a = randomMat(m, n, seed)
      val k = 3
      val full = Svds.svd(a, SingularSelection.All).toOption.get
      val partial = Svds.svd(a, SingularSelection.Count(k, SingularOrder.Largest)).toOption.get
      assert(partial.diagnostics.allConverged, s"partial not converged ${m}x$n: ${partial.diagnostics}")
      var i = 0
      while i < k do
        assert(
          math.abs(full.singularValues(i) - partial.singularValues(i)) < 1e-7,
          s"sigma($i) full=${full.singularValues(i)} partial=${partial.singularValues(i)} ${m}x$n seed=$seed"
        )
        i += 1
  }

  // --- rank deficiency -------------------------------------------------------

  test("full svd rank-deficient: trailing zeros and honest rank") {
    // rank 2 by construction: two outer products.
    val x1 = (0 until 6).map(i => 1.0 + 0.3 * i)
    val y1 = (0 until 5).map(j => 2.0 - 0.4 * j)
    val x2 = (0 until 6).map(i => math.sin(i + 1.0))
    val y2 = (0 until 5).map(j => math.cos(j + 0.5))
    val a = Matrix.tabulate(6, 5)((i, j) => 3.0 * x1(i) * y1(j) + x2(i) * y2(j))
    val s = Svds.svd(a, SingularSelection.All).toOption.get
    assertEconomy(s, 6, 5)
    assertDescendingNonNegative(values(s))
    assert(reconstructionError(a, s) < 1e-9)
    assert(values(s)(2) < 1e-10 * values(s)(0), s"expected numerical zero, got ${values(s)}")
    assertEquals(s.rank, 2)
    assertEquals(s.diagnostics.rank, Some(2))
  }

  test("full svd of the zero matrix: all-zero values, rank 0, orthonormal factors") {
    val s = Svds.svd(DMat.zeros(4, 3), SingularSelection.All).toOption.get
    assertEconomy(s, 4, 3)
    assertClose(values(s), IndexedSeq(0.0, 0.0, 0.0), 1e-15)
    assertEquals(s.rank, 0)
    assert(orthoError(s.u) < 1e-12)
    assert(orthoError(s.vt.t) < 1e-12)
  }

  // --- degenerate shapes -----------------------------------------------------

  test("full svd 1x1, mx1, 1xn edges") {
    val one = Svds.svd(Matrix.dense(1, 1)(-3.0), SingularSelection.All).toOption.get
    assertEconomy(one, 1, 1)
    assertClose(values(one), IndexedSeq(3.0), 1e-15)
    assert(reconstructionError(Matrix.dense(1, 1)(-3.0), one) < 1e-15)

    val col = Matrix.dense(4, 1)(1.0, 2.0, 2.0, 4.0) // ‖·‖₂ = 5
    val sc = Svds.svd(col, SingularSelection.All).toOption.get
    assertEconomy(sc, 4, 1)
    assertClose(values(sc), IndexedSeq(5.0), 1e-12)
    fullChecks(col, sc, 1e-12)

    val row = Matrix.dense(1, 4)(2.0, 1.0, 2.0, 4.0) // ‖·‖₂ = 5
    val sr = Svds.svd(row, SingularSelection.All).toOption.get
    assertEconomy(sr, 1, 4)
    assertClose(values(sr), IndexedSeq(5.0), 1e-12)
    fullChecks(row, sr, 1e-12)
  }

  // --- selection routing -----------------------------------------------------

  test("Count(k = min(m,n)) routes to the dense kernel, either order") {
    val a = randomMat(10, 6, 77L)
    val all = Svds.svd(a, SingularSelection.All).toOption.get
    val largest = Svds.svd(a, SingularSelection.Count(6, SingularOrder.Largest)).toOption.get
    val smallest = Svds.svd(a, SingularSelection.Count(6, SingularOrder.Smallest)).toOption.get
    assertClose(values(largest), values(all), 1e-14)
    assertClose(values(smallest), values(all), 1e-14)
    assertEconomy(largest, 10, 6)
  }

  test("full svd values-only: empty factors, values correct") {
    val a = rectDiag(6, 4, Seq(8.0, 4.0, 2.0, 1.0))
    val s = Svds.svd(a, SingularSelection.All, EigenVectors.ValuesOnly).toOption.get
    assertEquals(s.u.cols, 0)
    assertEquals(s.vt.rows, 0)
    assertClose(values(s), IndexedSeq(8.0, 4.0, 2.0, 1.0), 1e-10)
    assertEquals(s.diagnostics.orthogonalityError, 0.0)
  }

  // --- structural violations -------------------------------------------------

  test("full svd structural violations return Left") {
    val a = randomMat(6, 4, 5L)
    Svds.svd(a, SingularSelection.All, EigenVectors.Left) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected InvalidArgument for Left vectors, got $other")
    assert(Svds.svd(a, SingularSelection.Count(5, SingularOrder.Largest)).isLeft) // k > min(6,4)
    assert(Svds.svd(DMat.zeros(0, 3), SingularSelection.All).isLeft)
  }

  // --- DMat facade -----------------------------------------------------------

  test("DMat.svd facade is the full economy SVD") {
    val a = randomMat(9, 5, 88L)
    val s = a.svd.toOption.get
    fullChecks(a, s, 1e-9)
    assertClose(values(s), values(Svds.svd(a, SingularSelection.All).toOption.get), 1e-14)
  }

  test("a DenseSvd-capable provider Left falls back to the pure kernel (S8 policy, no new failure mode)") {
    // Advertises the capability but declines every call (the trait default is
    // Left(UnsupportedOperation)) — the S7 seam must then compute pure, exactly
    // like Eigen.symmetricSpectrum, never surface the provider's Left.
    object DecliningSvdBackend extends SpectralBackend:
      val name = "declining-svd"
      val capabilities = Set(SpectralCapability.DenseSvd)
    val a = randomMat(7, 4, 21L)
    val viaDeclining = Svds.svd(a, SingularSelection.All)(using DecliningSvdBackend).toOption.get
    val pure = Svds.svd(a, SingularSelection.All)(using SpectralBackend.none).toOption.get
    assertClose(values(viaDeclining), values(pure), 1e-15) // same pure kernel → bit-identical
    fullChecks(a, viaDeclining, 1e-9)
  }
