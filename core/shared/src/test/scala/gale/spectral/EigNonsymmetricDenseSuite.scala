package gale.spectral

import gale.linalg.DMat
import gale.linalg.LinAlgError
import gale.linalg.Matrix

/** Tests for the dense nonsymmetric facade `Eigen.eigNonsymmetric(a, selection,
  * vectors)` — full Francis-QR solve with the selection realized as a permutation
  * of the full spectrum. Selection membership, the canonical ordering guarantee
  * (§ 2), the never-split boundary-pair rule, and the structural `Left`s are all
  * exercised here; the packing invariant is additionally enforced by the
  * `NonsymmetricEigenDecomposition` constructor on every returned result.
  */
class EigNonsymmetricDenseSuite extends munit.FunSuite:

  // --- fixtures --------------------------------------------------------------

  /** Block-diagonal matrix with a fully-known spectrum:
    *   6 (real),  1 ± 3i  (|λ| = √10),  −4 (real),  0.5 ± 0.5i (|λ| = √0.5).
    * The two rotation blocks make the complex eigenvectors non-trivial, so the
    * real-Schur packing is genuinely exercised, while the values stay exact.
    */
  private val blockDiag: DMat = Matrix.dense(6, 6)(
    6.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, -3.0, 0.0, 0.0, 0.0,
    0.0, 3.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -4.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.5, -0.5,
    0.0, 0.0, 0.0, 0.0, 0.5, 0.5
  )

  private def randomReal(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  private def round(x: Double): Double =
    math.round(x * 1e6) / 1e6

  private def valueSet(d: NonsymmetricEigenDecomposition): Set[(Double, Double)] =
    (0 until d.size).map(i => (round(d.eigenvalue(i).re), round(d.eigenvalue(i).im))).toSet

  private def frobenius(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  /** ‖A v − λ v‖ for the i-th eigenpair in complex arithmetic (A real). */
  private def eigenpairResidual(a: DMat, d: NonsymmetricEigenDecomposition, i: Int): Double =
    val lambda = d.eigenvalue(i)
    val (vr, vi) = d.eigenvector(i)
    val realPart = (a * vr) - (vr * lambda.re - vi * lambda.im)
    val imagPart = (a * vi) - (vi * lambda.re + vr * lambda.im)
    math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))

  /** Assert conjugate pairs are adjacent, positive-imaginary member first, with
    * exact conjugate symmetry — a positive check on top of the constructor's.
    */
  private def assertPacking(d: NonsymmetricEigenDecomposition): Unit =
    var i = 0
    while i < d.size do
      val lam = d.eigenvalue(i)
      if lam.im > 0.0 then
        assert(i + 1 < d.size, s"positive-imag eigenvalue $i has no successor")
        val nxt = d.eigenvalue(i + 1)
        assertEquals(nxt.im, -lam.im, s"pair $i imaginary parts not exact negatives")
        assertEquals(nxt.re, lam.re, s"pair $i real parts differ")
      i += 1

  // --- All / ordering --------------------------------------------------------

  test("dense All: full spectrum in descending-magnitude canonical order") {
    val d = Eigen.eigNonsymmetric(blockDiag, EigenSelection.All).toOption.get
    assertEquals(d.size, 6)
    // Layout: 6, −4, then the √10 pair, then the √0.5 pair (magnitudes 6,4,√10,√0.5).
    assertEquals(d.eigenvalue(0), Complex(6.0, 0.0))
    assertEquals(d.eigenvalue(1), Complex(-4.0, 0.0))
    assertEqualsDouble(d.eigenvalue(2).re, 1.0, 1e-9)
    assertEqualsDouble(d.eigenvalue(2).im, 3.0, 1e-9)
    assertEqualsDouble(d.eigenvalue(4).re, 0.5, 1e-9)
    assertEqualsDouble(math.abs(d.eigenvalue(4).im), 0.5, 1e-9)
    assertPacking(d)
    assert(d.diagnostics.allConverged)
    var i = 0
    while i < d.size do
      assert(eigenpairResidual(blockDiag, d, i) < 1e-9, s"eigenpair $i residual")
      i += 1
  }

  // --- selection membership for the four legal orders ------------------------

  test("dense Count LargestMagnitude: membership by |λ|") {
    val k1 = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(1, EigenOrder.LargestMagnitude)).toOption.get
    assertEquals(valueSet(k1), Set((6.0, 0.0)))
    assertEquals(k1.eigenvalue(0), Complex(6.0, 0.0))

    val k2 = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(2, EigenOrder.LargestMagnitude)).toOption.get
    assertEquals(valueSet(k2), Set((6.0, 0.0), (-4.0, 0.0)))
  }

  test("dense Count SmallestMagnitude: membership by |λ|, layout ascending") {
    // k=2 selects the √0.5 pair (smallest magnitude); it is a single unit.
    val d = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(2, EigenOrder.SmallestMagnitude)).toOption.get
    assertEquals(valueSet(d), Set((0.5, 0.5), (0.5, -0.5)))
    // Ascending-magnitude layout ⇒ the pair is first, positive-imag member first.
    assertEquals(d.eigenvalue(0), Complex(0.5, 0.5))
    assertEquals(d.eigenvalue(1), Complex(0.5, -0.5))
  }

  test("dense Count LargestRealPart: membership by re(λ)") {
    val k1 = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(1, EigenOrder.LargestRealPart)).toOption.get
    assertEquals(valueSet(k1), Set((6.0, 0.0)))
    // k=2 boundary falls in the re=1 pair ⇒ whole pair kept ⇒ 3 values.
    val k2 = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(2, EigenOrder.LargestRealPart)).toOption.get
    assertEquals(valueSet(k2), Set((6.0, 0.0), (1.0, 3.0), (1.0, -3.0)))
    assertEquals(k2.size, 3)
  }

  test("dense Count SmallestRealPart: membership by re(λ)") {
    val k1 = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(1, EigenOrder.SmallestRealPart)).toOption.get
    assertEquals(valueSet(k1), Set((-4.0, 0.0)))
    val k2 = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(2, EigenOrder.SmallestRealPart)).toOption.get
    assertEquals(valueSet(k2), Set((-4.0, 0.0), (0.5, 0.5), (0.5, -0.5)))
  }

  // --- boundary-pair rule lock -----------------------------------------------

  test("dense Count boundary pair: k=3 keeps the whole straddling pair (returns k+1)") {
    // LargestMagnitude accumulates 6 (1), −4 (2); the 3rd/4th are the √10 pair,
    // whose two members share |λ| = √10. The never-split rule returns all four.
    val d = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(3, EigenOrder.LargestMagnitude)).toOption.get
    assertEquals(d.size, 4)
    assertEquals(valueSet(d), Set((6.0, 0.0), (-4.0, 0.0), (1.0, 3.0), (1.0, -3.0)))
    assertPacking(d)
    // A dense one-shot solve: requested == converged == returned, allConverged true.
    assertEquals(d.diagnostics.requested, 4)
    assertEquals(d.diagnostics.converged, 4)
    assert(d.diagnostics.allConverged)
  }

  test("dense Count boundary pair: SmallestMagnitude k=1 returns the whole smallest pair") {
    val d = Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(1, EigenOrder.SmallestMagnitude)).toOption.get
    assertEquals(d.size, 2)
    assertEquals(valueSet(d), Set((0.5, 0.5), (0.5, -0.5)))
  }

  // --- ordering tie-break ----------------------------------------------------

  test("dense tie-break: equal-|λ| pair precedes the real by descending real part") {
    // Spectrum: −√5 (real, placed first in the Schur order) and 1 ± 2i, all with
    // |λ| = √5. The real eigenvalue is set to −hypot(1,2) so the magnitude tie is
    // EXACT (the 2×2 block yields re=1, im=2 exactly; the 1×1 keeps its entry), so
    // the § 2 tie-break — descending real part — is what orders them: the pair
    // (re=1) must come before the real (re≈−2.236). orderUnits actively reorders
    // the kernel's real-first Schur output.
    val c = math.hypot(1.0, 2.0)
    val a = Matrix.dense(3, 3)(
      -c, 0.0, 0.0,
      0.0, 1.0, -2.0,
      0.0, 2.0, 1.0
    )
    val d = Eigen.eigNonsymmetric(a, EigenSelection.All).toOption.get
    assertEquals(d.eigenvalue(0), Complex(1.0, 2.0))
    assertEquals(d.eigenvalue(1), Complex(1.0, -2.0))
    assertEqualsDouble(d.eigenvalue(2).re, -c, 1e-12)
    assertEqualsDouble(d.eigenvalue(2).im, 0.0, 1e-12)
    assertPacking(d)

    // The boundary keeps the pair whole: Count(1) picks the top unit (the pair, by
    // the tie-break over the equal-magnitude real) and returns both members.
    val one = Eigen.eigNonsymmetric(a, EigenSelection.Count(1, EigenOrder.LargestMagnitude)).toOption.get
    assertEquals(one.size, 2)
    assertEquals(valueSet(one), Set((1.0, 2.0), (1.0, -2.0)))
  }

  // --- values-only -----------------------------------------------------------

  test("dense values-only: eigenvalues present, eigenvectors absent") {
    val d = Eigen.eigNonsymmetric(blockDiag, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    assertEquals(d.size, 6)
    assertEquals(valueSet(d), Set((6.0, 0.0), (-4.0, 0.0), (1.0, 3.0), (1.0, -3.0), (0.5, 0.5), (0.5, -0.5)))
    // isRealPair still works from the SoA parts.
    assert(d.isRealPair(0))
    assert(!d.isRealPair(2))
    // Vectors were not computed.
    intercept[LinAlgError.UnsupportedOperation](d.eigenvector(0))
  }

  test("dense values-only agrees with the vector solve on the spectrum") {
    val a = randomReal(9, 2027L)
    val full = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Right).toOption.get
    val vo = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    assertEquals(vo.size, full.size)
    var i = 0
    while i < full.size do
      assertEquals(vo.eigenvalue(i), full.eigenvalue(i), s"eigenvalue $i")
      i += 1
  }

  // --- residuals through the typed accessors ---------------------------------

  test("dense random matrix: eigenpairs satisfy A v = λ v and diagnostics residuals align") {
    val n = 12
    val a = randomReal(n, 20260711L)
    val d = Eigen.eigNonsymmetric(a, EigenSelection.All).toOption.get
    assertEquals(d.size, n)
    assertPacking(d)
    val scale = math.max(frobenius(a), 1.0)
    var i = 0
    while i < n do
      val res = eigenpairResidual(a, d, i)
      assert(res < 1e-8 * scale, s"eigenpair $i residual $res")
      // The recorded per-eigenvalue residual matches the independently computed one.
      assert(math.abs(d.diagnostics.residuals(i) - res) < 1e-9 * scale, s"diagnostics residual $i")
      i += 1
  }

  // --- structural Lefts ------------------------------------------------------

  test("dense structural violations return the mapped Left") {
    // non-square.
    Eigen.eigNonsymmetric(Matrix.tabulate(3, 4)((_, _) => 1.0), EigenSelection.All) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected NonSquareMatrix, got $other")

    // IndexRange / ValueInterval are symmetric-only.
    assert(Eigen.eigNonsymmetric(blockDiag, EigenSelection.IndexRange(0, 2)).isLeft)
    assert(Eigen.eigNonsymmetric(blockDiag, EigenSelection.ValueInterval(-1.0, 5.0)).isLeft)
    Eigen.eigNonsymmetric(blockDiag, EigenSelection.IndexRange(0, 2)) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected InvalidArgument, got $other")

    // Illegal (symmetric-only) orders.
    for order <- Seq(EigenOrder.LargestAlgebraic, EigenOrder.SmallestAlgebraic, EigenOrder.BothEnds) do
      Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(2, order)) match
        case Left(_: LinAlgError.InvalidArgument) => ()
        case other                                => fail(s"expected InvalidArgument for $order, got $other")

    // Bad k.
    assert(Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(0, EigenOrder.LargestMagnitude)).isLeft)
    assert(Eigen.eigNonsymmetric(blockDiag, EigenSelection.Count(7, EigenOrder.LargestMagnitude)).isLeft)

    // Left / left-and-right vectors are deferred (kernel computes right vectors only).
    for flag <- Seq(EigenVectors.Left, EigenVectors.LeftAndRight) do
      Eigen.eigNonsymmetric(blockDiag, EigenSelection.All, flag) match
        case Left(_: LinAlgError.UnsupportedOperation) => ()
        case other                                     => fail(s"expected UnsupportedOperation for $flag, got $other")
  }
