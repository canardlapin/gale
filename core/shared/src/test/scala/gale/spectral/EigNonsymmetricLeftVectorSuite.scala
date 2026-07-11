package gale.spectral

import gale.linalg.DMat
import gale.linalg.LinAlgError
import gale.linalg.Matrix

/** Tests for '''left''' eigenvectors of the dense nonsymmetric path
  * (`Eigen.eigNonsymmetric(a, selection, EigenVectors.Left | LeftAndRight)`),
  * recovered from the right-eigenvector matrix via the biorthogonal `V⁻¹` route.
  * Convention: `wᴴ A = λ wᴴ`, unit 2-norm, real-Schur SoA packing (positive-imag
  * member first), read through `leftEigenvector`.
  */
class EigNonsymmetricLeftVectorSuite extends munit.FunSuite:

  // --- fixtures --------------------------------------------------------------

  /** Companion of (x−1)(x−2)(x−3): eigenvalues 1, 2, 3 (real). */
  private val companion: DMat = Matrix.dense(3, 3)(
    0.0, 0.0, 6.0,
    1.0, 0.0, -11.0,
    0.0, 1.0, 6.0
  )

  /** Scaled rotation with eigenvalues 1 ± √3·i. */
  private val rotation: DMat =
    val s3 = math.sqrt(3.0)
    Matrix.dense(2, 2)(1.0, -s3, s3, 1.0)

  /** Block diagonal: real 5 and a 2×2 with eigenvalues 1 ± 2i. */
  private val mixed: DMat = Matrix.dense(3, 3)(
    5.0, 0.0, 0.0,
    0.0, 1.0, -2.0,
    0.0, 2.0, 1.0
  )

  private def randomReal(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

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

  /** ‖wᴴ A − λ wᴴ‖ for the i-th left eigenpair, in real arithmetic via Aᵀ. */
  private def leftResidual(a: DMat, d: NonsymmetricEigenDecomposition, i: Int): Double =
    val lam = d.eigenvalue(i)
    val (wr, wi) = d.leftEigenvector(i)
    val at = a.t
    val realPart = (at * wr) - (wr * lam.re) - (wi * lam.im)
    val imagPart = (wi * lam.re) - (at * wi) - (wr * lam.im)
    math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))

  /** The complex inner product `wᵢᴴ vⱼ = (wr − i·wi)·(vr + i·vj)` as `(re, im)`. */
  private def leftRightInner(d: NonsymmetricEigenDecomposition, i: Int, j: Int): Double =
    val (wr, wi) = d.leftEigenvector(i)
    val (vr, vi) = d.eigenvector(j)
    val re = wr.dot(vr) + wi.dot(vi)
    val im = wr.dot(vi) - wi.dot(vr)
    math.sqrt(re * re + im * im)

  private def assertUnitNorm(d: NonsymmetricEigenDecomposition): Unit =
    var i = 0
    while i < d.size do
      val (wr, wi) = d.leftEigenvector(i)
      val nrm = math.sqrt(wr.dot(wr) + wi.dot(wi))
      assert(math.abs(nrm - 1.0) < 1e-9, s"left vector $i norm $nrm")
      i += 1

  // --- wᴴA = λwᴴ on fixed spectra --------------------------------------------

  test("left eigenvectors satisfy wᴴA = λwᴴ (companion, rotation, mixed)") {
    for a <- Seq(companion, rotation, mixed) do
      val d = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.LeftAndRight).toOption.get
      val scale = math.max(frobenius(a), 1.0)
      assertUnitNorm(d)
      var i = 0
      while i < d.size do
        assert(leftResidual(a, d, i) < 1e-10 * scale, s"left residual $i for $a")
        i += 1
  }

  test("left eigenvectors of a random matrix satisfy wᴴA = λwᴴ") {
    for seed <- Seq(1L, 42L, 2027L) do
      val a = randomReal(9, seed)
      val d = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Left).toOption.get
      val scale = math.max(frobenius(a), 1.0)
      assertUnitNorm(d)
      var i = 0
      while i < d.size do
        assert(leftResidual(a, d, i) < 1e-7 * scale, s"left residual $i (seed $seed) = ${leftResidual(a, d, i)}")
        i += 1
  }

  // --- biorthogonality (distinct eigenvalues) --------------------------------

  test("biorthogonality: wᵢᴴvⱼ ≈ 0 for i≠j, nonzero for i=j (unit norm, not biorthonormal)") {
    // Mixed spectrum {5, 1±2i} — all distinct — so left/right vectors biorthogonalize.
    val d = Eigen.eigNonsymmetric(mixed, EigenSelection.All, EigenVectors.LeftAndRight).toOption.get
    val n = d.size
    var i = 0
    while i < n do
      var j = 0
      while j < n do
        val inner = leftRightInner(d, i, j)
        if i == j then assert(inner > 1e-6, s"wᵢᴴvᵢ should be nonzero at $i (got $inner)")
        else assert(inner < 1e-9, s"wᵢᴴvⱼ should vanish at ($i,$j) (got $inner)")
        j += 1
      i += 1
    // Normalization is unit 2-norm, NOT biorthonormal: some wᵢᴴvᵢ ≠ 1.
    assertUnitNorm(d)
  }

  // --- packing ---------------------------------------------------------------

  test("left vector packing: conjugate pairs decode with negated imaginary part") {
    val d = Eigen.eigNonsymmetric(mixed, EigenSelection.All, EigenVectors.Left).toOption.get
    var i = 0
    while i < d.size do
      val lam = d.eigenvalue(i)
      if lam.im > 0.0 then
        val (wrA, wiA) = d.leftEigenvector(i)
        val (wrB, wiB) = d.leftEigenvector(i + 1)
        // Conjugate: same real part, negated imaginary part.
        var r = 0
        while r < wrA.length do
          assert(math.abs(wrB(r) - wrA(r)) < 1e-12, s"pair $i left real parts differ at $r")
          assert(math.abs(wiB(r) + wiA(r)) < 1e-12, s"pair $i left imag parts not negated at $r")
          r += 1
      i += 1
  }

  // --- flag consistency ------------------------------------------------------

  test("Left, Right, LeftAndRight agree on eigenvalues; each exposes only its vectors") {
    val a = randomReal(6, 7L)
    val right = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Right).toOption.get
    val left = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Left).toOption.get
    val both = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.LeftAndRight).toOption.get

    // Eigenvalues identical (same kernel solve underneath).
    var i = 0
    while i < right.size do
      assertEquals(left.eigenvalue(i), right.eigenvalue(i), s"eigenvalue $i")
      assertEquals(both.eigenvalue(i), right.eigenvalue(i), s"eigenvalue $i")
      i += 1

    // Right-only: right present, left absent.
    right.eigenvector(0)
    intercept[LinAlgError.UnsupportedOperation](right.leftEigenvector(0))
    // Left-only: left present, right absent.
    left.leftEigenvector(0)
    intercept[LinAlgError.UnsupportedOperation](left.eigenvector(0))
    // LeftAndRight: both present, and the left/right vectors match the one-sided
    // solves elementwise (same kernel solve + same V⁻¹ recovery underneath).
    both.eigenvector(0)
    both.leftEigenvector(0)
    assertVecEquals(both.leftEigenvector(0)._1, left.leftEigenvector(0)._1)
    assertVecEquals(both.eigenvector(0)._1, right.eigenvector(0)._1)
  }

  private def assertVecEquals(a: gale.linalg.DVec, b: gale.linalg.DVec): Unit =
    assertEquals(a.length, b.length)
    var i = 0
    while i < a.length do
      assertEquals(a(i), b(i), s"element $i")
      i += 1

  // --- selection lockstep ----------------------------------------------------

  test("selection permutes left columns in lockstep with right and eigenvalues") {
    val a = randomReal(10, 909L)
    val k = 4
    val d = Eigen.eigNonsymmetric(a, EigenSelection.Count(k, EigenOrder.LargestMagnitude), EigenVectors.LeftAndRight).toOption.get
    val scale = math.max(frobenius(a), 1.0)
    // If a left column had drifted out of lockstep with its eigenvalue, wᴴA = λwᴴ
    // would fail. Both residual sides tiny ⇒ the pair-aware permutation moved left
    // and right columns together.
    var i = 0
    while i < d.size do
      assert(leftResidual(a, d, i) < 1e-7 * scale, s"left residual $i under selection")
      i += 1
    // The diagnostics residual is the worse of the two sides — also tiny.
    assert(d.diagnostics.worstResidual < 1e-7 * scale)
  }

  // --- defective inputs ------------------------------------------------------

  test("defective matrices reject left vectors with Left(SingularMatrix)") {
    val j2 = Matrix.dense(2, 2)(2.0, 1.0, 0.0, 2.0) // Jordan block, one eigenvector
    val j3 = Matrix.dense(3, 3)(
      2.0, 1.0, 0.0,
      0.0, 2.0, 1.0,
      0.0, 0.0, 2.0
    )
    val nilpotent = Matrix.dense(2, 2)(0.0, 1.0, 0.0, 0.0)
    for (name, a) <- Seq(("J2", j2), ("J3", j3), ("nilpotent", nilpotent)) do
      for flag <- Seq(EigenVectors.Left, EigenVectors.LeftAndRight) do
        Eigen.eigNonsymmetric(a, EigenSelection.All, flag) match
          case Left(_: LinAlgError.SingularMatrix) => ()
          case other                               => fail(s"$name/$flag: expected SingularMatrix, got $other")
      // Right-only requests still succeed — right eigenvectors exist for a defective A;
      // the guard must fire only when left vectors were requested.
      val right = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Right)
      assert(right.isRight, s"$name: Right-only should succeed, got $right")
      val vo = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly)
      assert(vo.isRight, s"$name: ValuesOnly should succeed, got $vo")
  }

  test("near-defective (distinct eigenvalues 2, 2+1e-8) is accepted with honest left vectors") {
    // Distinct eigenvalues ⇒ diagonalizable ⇒ passes both guards (empirically the
    // embedded-LU pivot ratio ~1e-8, well above the ~ε defective floor; left residual
    // ~0). This pins the accepted side of the threshold.
    val a = Matrix.dense(2, 2)(2.0, 1.0, 0.0, 2.0 + 1e-8)
    Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.Left) match
      case Right(d) =>
        assertEquals(d.size, 2)
        var i = 0
        while i < d.size do
          assert(leftResidual(a, d, i) < 1e-6 * math.max(frobenius(a), 1.0), s"near-defective left residual $i")
          i += 1
      case other => fail(s"expected Right (near-defective is diagonalizable), got $other")
  }

  // --- unchanged paths -------------------------------------------------------

  test("values-only and right-only are unchanged; Arnoldi still rejects left") {
    val a = randomReal(5, 3L)
    val vo = Eigen.eigNonsymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    intercept[LinAlgError.UnsupportedOperation](vo.leftEigenvector(0))
    intercept[LinAlgError.UnsupportedOperation](vo.eigenvector(0))

    // The iterative path has no left-vector notion.
    val opts = SpectralOptions(returnVectors = EigenVectors.Left)
    Eigen.eigNonsymmetric(a, 5, EigenSelection.Count(1, EigenOrder.LargestMagnitude), opts) match
      case Left(_: LinAlgError.UnsupportedOperation) => ()
      case other                                     => fail(s"expected UnsupportedOperation, got $other")
    val opts2 = SpectralOptions(returnVectors = EigenVectors.LeftAndRight)
    assert(Eigen.eigNonsymmetric(a, 5, EigenSelection.Count(1, EigenOrder.LargestMagnitude), opts2).isLeft)
  }
