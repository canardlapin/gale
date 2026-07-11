package gale.spectral

import gale.linalg.*

/** Tests for the QZ boundary seam: the projective `(α, β)` sorter
  * (`Eigen.generalizedIndices`), the `GeneralizedEigenDecomposition` packing
  * invariants, and — the money test — a FAKE in-package backend driving the
  * `eigGeneralizedNonsymmetric` facade end-to-end (hand-built raw QZ output →
  * facade canonicalizes → sealed result with correct ordering, packing, and
  * diagnostics), proving the boundary works before any real backend exists.
  */
class GeneralizedEigenSuite extends munit.FunSuite:

  private def assertVec(actual: DVec, expected: Seq[Double])(using munit.Location): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < expected.length do
      assertEqualsDouble(actual(i), expected(i), 1e-12)
      i += 1

  // --- projective (α, β) sorter ----------------------------------------------

  test("generalizedIndices: finite ordering by descending |α/β|") {
    // ratios 3, 1, 5 → descending 5, 3, 1 → indices [2, 0, 1].
    val idx = Eigen.generalizedIndices(Vec(3.0, 1.0, 5.0), Vec(0.0, 0.0, 0.0), Vec(1.0, 1.0, 1.0), EigenOrder.LargestMagnitude)
    assertEquals(idx.toSeq, Seq(2, 0, 1))
  }

  test("generalizedIndices: multiple infinites tie-break by descending |α|, then index") {
    // β=0 at 0 (|α|=2) and 1 (|α|=5); finite ratio 3 at 2. LargestMagnitude:
    // infinites first (|α| desc → 1 then 0), then the finite.
    val idx = Eigen.generalizedIndices(Vec(2.0, 5.0, 3.0), Vec(0.0, 0.0, 0.0), Vec(0.0, 0.0, 1.0), EigenOrder.LargestMagnitude)
    assertEquals(idx.toSeq, Seq(1, 0, 2))
  }

  test("generalizedIndices: conjugate pair kept adjacent while reordering") {
    // real ratio 1 at 0; pair |3±i|/1 ≈ 3.16 at (1,2). LargestMagnitude puts the
    // pair first (both columns, positive-imag member leading), then the real.
    val idx = Eigen.generalizedIndices(Vec(1.0, 3.0, 3.0), Vec(0.0, 1.0, -1.0), Vec(1.0, 1.0, 1.0), EigenOrder.LargestMagnitude)
    assertEquals(idx.toSeq, Seq(1, 2, 0))
  }

  test("generalizedIndices: SmallestMagnitude puts finites ascending, infinites last") {
    // β=0 at 0; finite ratios 5 at 1, 3 at 2. Smallest → 3, 5, then infinite.
    val idx = Eigen.generalizedIndices(Vec(2.0, 5.0, 3.0), Vec(0.0, 0.0, 0.0), Vec(0.0, 1.0, 1.0), EigenOrder.SmallestMagnitude)
    assertEquals(idx.toSeq, Seq(2, 1, 0))
  }

  // --- packing invariants ----------------------------------------------------

  test("GeneralizedEigenDecomposition: valid packing (real, infinite, pair) constructs") {
    val d = new GeneralizedEigenDecomposition(
      Vec(3.0, 2.0, 0.0, 0.0),
      Vec(0.0, 0.0, 1.0, -1.0),
      Vec(1.0, 0.0, 1.0, 1.0),
      DMat.zeros(4, 0),
      None,
      SpectralDiagnostics(4, 4, DVec.zeros(4), 0.0, 0)
    )
    assertEquals(d.size, 4)
    assert(!d.isInfinite(0) && d.isInfinite(1))
    assert(d.isRealPair(1)) // infinite is real
    assert(!d.isRealPair(2))
    assertEquals(d.eigenvalue(2), GeneralizedEigenvalue(Complex(0.0, 1.0), 1.0))
  }

  test("GeneralizedEigenDecomposition: complex-infinite (β=0, alphaIm≠0) is rejected") {
    intercept[IllegalArgumentException] {
      new GeneralizedEigenDecomposition(
        Vec(1.0, 1.0),
        Vec(3.0, -3.0),
        Vec(0.0, 0.0), // β=0 with a nonzero alphaIm — not representable
        DMat.zeros(2, 0),
        None,
        SpectralDiagnostics(2, 2, DVec.zeros(2), 0.0, 0)
      )
    }
  }

  test("GeneralizedEigenDecomposition: a split pair (unequal β) is rejected") {
    intercept[IllegalArgumentException] {
      new GeneralizedEigenDecomposition(
        Vec(2.0, 2.0),
        Vec(1.0, -1.0),
        Vec(1.0, 2.0), // conjugate α but different β — mispacked
        DMat.zeros(2, 0),
        None,
        SpectralDiagnostics(2, 2, DVec.zeros(2), 0.0, 0)
      )
    }
  }

  // --- the money test: fake backend → facade end-to-end ----------------------

  /** A stateless fake QZ backend returning a hand-built raw spectrum for the test
    * pencil below: real λ=3, complex pair λ=±i, and one infinite (β=0) eigenvalue.
    */
  private def fakeQz: SpectralBackend =
    new SpectralBackend:
      def name: String = "fake-qz"
      def capabilities: Set[SpectralCapability] = Set(SpectralCapability.GeneralizedNonsymmetricEigen)
      override def generalizedNonsymmetricEigen(
          x: DMat,
          y: DMat,
          vectors: EigenVectors
      ): Either[LinAlgError, RawGeneralizedEigen] =
        Right(
          RawGeneralizedEigen(
            alphaRe = Vec(3.0, 0.0, 0.0, 2.0),
            alphaIm = Vec(0.0, 1.0, -1.0, 0.0),
            beta = Vec(1.0, 1.0, 1.0, 0.0),
            // Columns: e0 (λ=3); (0,1,0,0)/(0,0,-1,0) for the ±i pair; e3 (infinite).
            rightPacked = Matrix.dense(4, 4)(
              1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, -1.0, 0.0,
              0.0, 0.0, 0.0, 1.0
            ),
            leftPacked = None,
            schur = None
          )
        )

  test("eigGeneralizedNonsymmetric: fake QZ backend drives the facade end-to-end") {
    // Pencil A x = λ B x with A, B block diagonal so the fake's raw output is an
    // exact eigendecomposition (residuals vanish).
    val a = Matrix.dense(4, 4)(
      3.0, 0.0, 0.0, 0.0,
      0.0, 0.0, -1.0, 0.0,
      0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 2.0
    )
    val b = Matrix.dense(4, 4)(
      1.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0,
      0.0, 0.0, 0.0, 0.0
    )
    val d = Eigen.eigGeneralizedNonsymmetric(a, b)(using fakeQz).toOption.get
    assertEquals(d.size, 4)

    // Canonical order (LargestMagnitude): infinite (∞) first, then λ=3, then ±i.
    assert(d.isInfinite(0))
    assert(d.isRealPair(0)) // infinite owns a single real column
    assertEquals(d.eigenvalue(0), GeneralizedEigenvalue(Complex(2.0, 0.0), 0.0))
    assert(!d.isInfinite(1))
    assertEquals(d.eigenvalue(1), GeneralizedEigenvalue(Complex(3.0, 0.0), 1.0))
    assertEquals(d.eigenvalue(1).value, Complex(3.0, 0.0))
    // The conjugate pair, adjacent, positive-imaginary first, equal β.
    assert(!d.isRealPair(2))
    assertEquals(d.eigenvalue(2), GeneralizedEigenvalue(Complex(0.0, 1.0), 1.0))
    assertEquals(d.eigenvalue(3), GeneralizedEigenvalue(Complex(0.0, -1.0), 1.0))

    // Eigenvectors decode from the reordered packing.
    val (vr0, vi0) = d.eigenvector(0)
    assertVec(vr0, Seq(0.0, 0.0, 0.0, 1.0))
    assertVec(vi0, Seq(0.0, 0.0, 0.0, 0.0))
    val (vr2, vi2) = d.eigenvector(2)
    assertVec(vr2, Seq(0.0, 1.0, 0.0, 0.0))
    assertVec(vi2, Seq(0.0, 0.0, -1.0, 0.0))

    // Diagnostics: dense one-shot, rank = number of finite eigenvalues, residuals
    // (facade-derived homogeneous ‖βAx − αBx‖) vanish.
    assertEquals(d.diagnostics.requested, 4)
    assertEquals(d.diagnostics.converged, 4)
    assert(d.diagnostics.allConverged)
    assertEquals(d.diagnostics.orthogonalityError, 0.0)
    assertEquals(d.diagnostics.rank, Some(3))
    assert(d.diagnostics.worstResidual < 1e-12, s"residual ${d.diagnostics.worstResidual}")
  }

  // --- residual math guards (nonzero, hand-computed) -------------------------

  /** A fake QZ backend returning a fixed raw payload (ignores its inputs). */
  private def fakeQzRaw(raw: RawGeneralizedEigen): SpectralBackend =
    new SpectralBackend:
      def name: String = "fake-qz-raw"
      def capabilities: Set[SpectralCapability] = Set(SpectralCapability.GeneralizedNonsymmetricEigen)
      override def generalizedNonsymmetricEigen(x: DMat, y: DMat, v: EigenVectors): Either[LinAlgError, RawGeneralizedEigen] =
        Right(raw)

  test("generalizedQzResiduals: complex pair with a deliberately wrong vector gives ‖√50‖") {
    // A = [[0,-1],[1,0]], B = I, eigenvalues ±i. The backend returns the WRONG
    // eigenvector (v_re=(2,0), v_im=(0,3)) so the homogeneous residual is nonzero:
    //   realPart = A·v_re + B·v_im = (0,2)+(0,3) = (0,5)
    //   imagPart = A·v_im − B·v_re = (-3,0)-(2,0) = (-5,0)  ⇒ ‖·‖ = √50.
    val a = Matrix.dense(2, 2)(0.0, -1.0, 1.0, 0.0)
    val b = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val raw = RawGeneralizedEigen(
      Vec(0.0, 0.0),
      Vec(1.0, -1.0),
      Vec(1.0, 1.0),
      Matrix.dense(2, 2)(2.0, 0.0, 0.0, 3.0), // col0 = (2,0), col1 = (0,3)
      None,
      None
    )
    val d = Eigen.eigGeneralizedNonsymmetric(a, b)(using fakeQzRaw(raw)).toOption.get
    assertEqualsDouble(d.diagnostics.residuals(0), math.sqrt(50.0), 1e-9)
    assertEqualsDouble(d.diagnostics.residuals(1), math.sqrt(50.0), 1e-9)
  }

  test("generalizedQzResiduals: infinite eigenvalue with x ∉ null(B) gives |α|·‖Bx‖") {
    // n=1: A=[[0]], B=[[1]]. Backend claims β=0 (infinite), α=2, x=(1). x is not in
    // null(B), so residual = ‖−2·B x‖ = 2·‖(1)‖ = 2.
    val a = Matrix.dense(1, 1)(0.0)
    val b = Matrix.dense(1, 1)(1.0)
    val raw = RawGeneralizedEigen(Vec(2.0), Vec(0.0), Vec(0.0), Matrix.dense(1, 1)(1.0), None, None)
    val d = Eigen.eigGeneralizedNonsymmetric(a, b)(using fakeQzRaw(raw)).toOption.get
    assert(d.isInfinite(0))
    assertEqualsDouble(d.diagnostics.residuals(0), 2.0, 1e-9)
  }

  test("assembleQz: negative β is rejected as a backend contract violation") {
    val a = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val raw = RawGeneralizedEigen(Vec(1.0, 1.0), Vec(0.0, 0.0), Vec(1.0, -1.0), DMat.zeros(2, 0), None, None)
    Eigen.eigGeneralizedNonsymmetric(a, a)(using fakeQzRaw(raw)) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected InvalidArgument for β<0, got $other")
  }

  test("left residual honesty: a wrong LEFT vector surfaces in worstResidual, not a deceptive 0") {
    // A = diag(3,5), B = I. Left-only result: the backend returns the correct left
    // vector for λ=5 (residual 0) but a WRONG one, (1,1), for λ=3:
    //   ‖1·Aᵀ(1,1) − 3·Bᵀ(1,1)‖ = ‖(3,5)-(3,3)‖ = ‖(0,2)‖ = 2.
    val a = Matrix.dense(2, 2)(3.0, 0.0, 0.0, 5.0)
    val b = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val raw = RawGeneralizedEigen(
      Vec(3.0, 5.0),
      Vec(0.0, 0.0),
      Vec(1.0, 1.0),
      DMat.zeros(2, 0), // no right vectors — a left-only result
      Some(Matrix.dense(2, 2)(1.0, 0.0, 1.0, 1.0)), // col0 = (1,1) for λ=3, col1 = (0,1) for λ=5
      None
    )
    val d = Eigen.eigGeneralizedNonsymmetric(a, b, EigenVectors.Left)(using fakeQzRaw(raw)).toOption.get
    // Canonical order puts λ=5 (residual 0) first, λ=3 (residual 2) second.
    assertEqualsDouble(d.eigenvalue(0).value.re, 5.0, 1e-12)
    assertEqualsDouble(d.diagnostics.residuals(1), 2.0, 1e-9)
    assertEqualsDouble(d.diagnostics.worstResidual, 2.0, 1e-9)
    // Left vectors present, right vectors absent.
    d.leftEigenvector(1)
    intercept[LinAlgError.UnsupportedOperation](d.eigenvector(0))
  }

  // --- S6: rank-deficient GSVD routing ---------------------------------------

  /** A fake rank-deficient GSVD backend returning raw `c`/`s` that classify to
    * Finite / Infinite / Zero (values-only; empty factor matrices).
    */
  private def fakeRankDeficientGsvd: SpectralBackend =
    new SpectralBackend:
      def name: String = "fake-gsvd"
      def capabilities: Set[SpectralCapability] = Set(SpectralCapability.RankDeficientGsvd)
      override def rankDeficientGsvd(x: DMat, y: DMat, wantVectors: Boolean): Either[LinAlgError, RawGsvd] =
        // ratios: 0.6/0.8 = 0.75 (Finite), 1/0 = ∞ (Infinite), 0/1 = 0 (Zero).
        Right(RawGsvd(DMat.zeros(x.rows, 0), DMat.zeros(y.rows, 0), DMat.zeros(x.cols, 0), Vec(0.6, 1.0, 0.0), Vec(0.8, 0.0, 1.0)))

  test("Svds.gsvd: rank-deficient pencil routes to a capable backend, else Left(RankDeficient)") {
    // m+p = 2 < n = 3 ⇒ rank-deficient stacked pencil.
    val a = Matrix.tabulate(1, 3)((_, _) => 1.0)
    val b = Matrix.tabulate(1, 3)((_, _) => 2.0)
    // With `none` (default), the shipped Left(RankDeficient) stands — unchanged.
    Svds.gsvd(a, b) match
      case Left(_: LinAlgError.RankDeficient) => ()
      case other                              => fail(s"expected RankDeficient with none, got $other")
    // With a capable backend, the facade canonicalizes the raw factors into a
    // GeneralizedSVD (Infinite-first / Zero-last descending order).
    val g = Svds.gsvd(a, b)(using fakeRankDeficientGsvd).toOption.get
    assertEquals(g.size, 3)
    assertEquals(g.values(0), GeneralizedSingularValue.Infinite)
    assertEquals(g.values(2), GeneralizedSingularValue.Zero)
    g.values(1) match
      case GeneralizedSingularValue.Finite(r) => assert(math.abs(r - 0.75) < 1e-9, s"ratio $r")
      case other                              => fail(s"expected Finite, got $other")
    // Descending ratio: +∞ ≥ 0.75 ≥ 0.
    assert(g.ratio(0) >= g.ratio(1) && g.ratio(1) >= g.ratio(2))
  }

  test("eigGeneralizedNonsymmetric: no backend → Left(UnsupportedOperation); structural Lefts validate first") {
    val a = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 2.0)
    val b = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    // Default given is `none`: unsupported.
    Eigen.eigGeneralizedNonsymmetric(a, b) match
      case Left(_: LinAlgError.UnsupportedOperation) => ()
      case other                                     => fail(s"expected UnsupportedOperation, got $other")
    // Shape validation happens before the backend is consulted, even with a capable one.
    Eigen.eigGeneralizedNonsymmetric(Matrix.dense(2, 3)(1.0, 0, 0, 0, 1.0, 0), b)(using fakeQz) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected NonSquareMatrix, got $other")
    Eigen.eigGeneralizedNonsymmetric(a, Matrix.dense(3, 3)(Seq.fill(9)(0.0)*))(using fakeQz) match
      case Left(_: LinAlgError.DimensionMismatch) => ()
      case other                                  => fail(s"expected DimensionMismatch, got $other")
  }
