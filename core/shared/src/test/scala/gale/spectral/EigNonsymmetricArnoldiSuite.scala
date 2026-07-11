package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.LinearOperator
import gale.linalg.Matrix

/** Tests for the iterative partial nonsymmetric eigensolver
  * `Eigen.eigNonsymmetric(op, n, Count, options)` — Arnoldi with full
  * reorthogonalization over a growing subspace. Convergence is checked against the
  * dense solve's extremes and against a matrix-free operator with a known complex
  * spectrum; pair adjacency, exact conjugate symmetry, the partial/zero-convergence
  * boundaries, and the structural `Left`s are exercised too.
  */
class EigNonsymmetricArnoldiSuite extends munit.FunSuite:

  // --- fixtures --------------------------------------------------------------

  /** Block-upper-triangular matrix (`n = 2·mags.length`): 2×2 scaled-rotation
    * blocks `mags(i)·[[c,-s],[s,c]]` on the diagonal (eigenvalues `mags(i)·e^{±iθ}`)
    * plus a strict block-upper-triangular random perturbation. The spectrum is the
    * union of the diagonal blocks' pairs (upper-triangular ⇒ perturbation-invariant
    * eigenvalues), and the perturbation makes it genuinely non-normal.
    */
  private def blockUpperTriangular(mags: Array[Double], theta: Double, fill: Double, seed: Long): DMat =
    val n = 2 * mags.length
    val c = math.cos(theta)
    val s = math.sin(theta)
    val rng = new scala.util.Random(seed)
    val noise = Array.ofDim[Double](n, n)
    var r = 0
    while r < n do
      var cc = 0
      while cc < n do
        noise(r)(cc) = rng.nextDouble() * 2.0 - 1.0
        cc += 1
      r += 1
    Matrix.tabulate(n, n): (row, col) =>
      val br = row / 2
      val bc = col / 2
      if br == bc then
        val m = mags(br)
        (row % 2, col % 2) match
          case (0, 0) => m * c
          case (0, 1) => -m * s
          case (1, 0) => m * s
          case _      => m * c
      else if bc > br then fill * noise(row)(col)
      else 0.0

  /** Matrix-free direct sum of scaled planar rotations: block `b` maps
    * `(x_{2b}, x_{2b+1})` by `mags(b)·[[c,-s],[s,c]]`. Known complex spectrum
    * `mags(b)·e^{±iθ}`; never materialized.
    */
  private def blockRotationOp(mags: Array[Double], theta: Double): DoubleLinearOperator =
    val n = 2 * mags.length
    val c = math.cos(theta)
    val s = math.sin(theta)
    LinearOperator.fromFunction(n, n): (x, into) =>
      var b = 0
      while b < mags.length do
        val x0 = x(2 * b)
        val x1 = x(2 * b + 1)
        into(2 * b) = mags(b) * (c * x0 - s * x1)
        into(2 * b + 1) = mags(b) * (s * x0 + c * x1)
        b += 1

  private def randomReal(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  /** ‖A v − λ v‖ for the i-th eigenpair, via the operator and the decoded vectors. */
  private def opResidual(op: DoubleLinearOperator, d: NonsymmetricEigenDecomposition, i: Int): Double =
    val lam = d.eigenvalue(i)
    val (vr, vi) = d.eigenvector(i)
    val realPart = (op * vr) - (vr * lam.re - vi * lam.im)
    val imagPart = (op * vi) - (vi * lam.re + vr * lam.im)
    math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))

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

  /** Every eigenvalue of `arn` matches some eigenvalue of `dense` within `tol`. */
  private def subsetOf(arn: NonsymmetricEigenDecomposition, dense: NonsymmetricEigenDecomposition, tol: Double): Boolean =
    (0 until arn.size).forall: i =>
      val la = arn.eigenvalue(i)
      (0 until dense.size).exists: j =>
        val ld = dense.eigenvalue(j)
        math.abs(la.re - ld.re) < tol && math.abs(la.im - ld.im) < tol

  // --- vs the dense solve ----------------------------------------------------

  test("arnoldi: largest-magnitude extremes match the dense solve") {
    val mags = Array.tabulate(15)(i => 20.0 / (i + 1)) // 20, 10, 6.67, 5, ... n = 30
    val a = blockUpperTriangular(mags, theta = 0.4, fill = 0.3, seed = 5150L)
    val n = a.rows
    val dense = Eigen.eigNonsymmetric(a, EigenSelection.All).toOption.get

    val two = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude)).toOption.get
    assert(two.diagnostics.allConverged, s"k=2 not converged: ${two.diagnostics}")
    assertEquals(two.size, 2)
    assert(subsetOf(two, dense, 1e-6), "k=2 values do not match the dense extremes")
    assertPacking(two)
    // The dominant pair really is the largest-magnitude pair of the spectrum.
    assertEqualsDouble(two.eigenvalue(0).magnitude, 20.0, 1e-6)

    val four = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(4, EigenOrder.LargestMagnitude)).toOption.get
    assert(four.diagnostics.allConverged, s"k=4 not converged: ${four.diagnostics}")
    assertEquals(four.size, 4)
    assert(subsetOf(four, dense, 1e-6), "k=4 values do not match the dense extremes")
    // Top two magnitudes are 20 and 10.
    assertEqualsDouble(four.eigenvalue(0).magnitude, 20.0, 1e-6)
    assertEqualsDouble(four.eigenvalue(2).magnitude, 10.0, 1e-6)
  }

  test("arnoldi: largest-real-part matches the dense solve") {
    val mags = Array.tabulate(15)(i => 20.0 / (i + 1))
    val a = blockUpperTriangular(mags, theta = 0.4, fill = 0.3, seed = 6161L)
    val n = a.rows
    val dense = Eigen.eigNonsymmetric(a, EigenSelection.All).toOption.get
    // Real part = |λ|·cos θ with θ fixed, so largest real part = largest magnitude.
    val d = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(4, EigenOrder.LargestRealPart)).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assert(subsetOf(d, dense, 1e-6), "largest-real-part values do not match dense")
  }

  // --- matrix-free operator with a known spectrum ----------------------------

  test("arnoldi: matrix-free operator resolves the dominant complex pair") {
    val mags = Array(30.0, 8.0, 5.0, 3.0, 2.0, 1.0) // n = 12
    val theta = 0.5
    val op = blockRotationOp(mags, theta)
    val n = 2 * mags.length

    val d = Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude)).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assertEquals(d.size, 2)
    assertPacking(d)
    // Analytic dominant pair 30·e^{±iθ}, positive-imag member first.
    val lam = d.eigenvalue(0)
    assertEqualsDouble(lam.re, 30.0 * math.cos(theta), 1e-7)
    assertEqualsDouble(lam.im, 30.0 * math.sin(theta), 1e-7)
    assertEquals(d.eigenvalue(1), Complex(lam.re, -lam.im))
    // Eigenvector residual through the typed accessors.
    assert(opResidual(op, d, 0) < 1e-7, "dominant eigenvector residual too large")
    assert(d.diagnostics.worstResidual < 1e-7)
  }

  test("arnoldi: values-only omits eigenvectors but still converges") {
    val mags = Array(30.0, 8.0, 5.0, 3.0, 2.0, 1.0)
    val op = blockRotationOp(mags, 0.5)
    val n = 2 * mags.length
    val options = SpectralOptions(returnVectors = EigenVectors.ValuesOnly)
    val d = Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude), options).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assertEqualsDouble(d.eigenvalue(0).magnitude, 30.0, 1e-7)
    intercept[LinAlgError.UnsupportedOperation](d.eigenvector(0))
  }

  test("arnoldi: explicit start vector converges to the dominant pair") {
    val mags = Array(30.0, 8.0, 5.0, 3.0, 2.0, 1.0)
    val op = blockRotationOp(mags, 0.5)
    val n = 2 * mags.length
    val start = DVec.tabulate(n)(i => if i % 3 == 0 then 1.0 else 0.4)
    val options = SpectralOptions(startVector = Some(start))
    val d = Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude), options).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assertEqualsDouble(d.eigenvalue(0).magnitude, 30.0, 1e-7)
  }

  // --- mixed real + complex spectrum (real branch + iterative boundary rule) --

  test("arnoldi: mixed real+complex spectrum exercises the real branch and the k+1 boundary") {
    // Block-diagonal operator with spectrum {10 (real), 6 ± 2i, 3 (real)}, driven
    // through the iterative path (a DMat is a DoubleLinearOperator). Magnitudes:
    // 10, √40 ≈ 6.32, 3 — so LargestMagnitude order is [10, 6±2i, 3].
    val a: DMat = Matrix.dense(4, 4)(
      10.0, 0.0, 0.0, 0.0,
      0.0, 6.0, -2.0, 0.0,
      0.0, 2.0, 6.0, 0.0,
      0.0, 0.0, 0.0, 3.0
    )
    val n = 4

    // Count(2, LM): 10 (1) then the 6±2i pair straddles the boundary ⇒ whole pair
    // kept ⇒ size 3, and requested tracks the never-split target so allConverged
    // stays a clean true.
    val two = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude)).toOption.get
    assertEquals(two.size, 3)
    assertEquals(two.diagnostics.requested, 3)
    assertEquals(two.diagnostics.converged, 3)
    assert(two.diagnostics.allConverged, s"not converged: ${two.diagnostics}")
    assertPacking(two)
    assertEqualsDouble(two.eigenvalue(0).re, 10.0, 1e-9)
    assertEqualsDouble(two.eigenvalue(0).im, 0.0, 1e-12)
    assert(two.isRealPair(0))
    assertEqualsDouble(two.eigenvalue(1).re, 6.0, 1e-9)
    assertEqualsDouble(two.eigenvalue(1).im, 2.0, 1e-9)
    // Conjugate symmetry is exact by construction (constructor-enforced).
    assertEquals(two.eigenvalue(2), Complex(two.eigenvalue(1).re, -two.eigenvalue(1).im))

    // Count(1, LM): a lone REAL eigenvalue through the ritzReal branch — exactly
    // [10], real, unit residual tiny.
    val one = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(1, EigenOrder.LargestMagnitude)).toOption.get
    assertEquals(one.size, 1)
    assertEquals(one.diagnostics.requested, 1)
    assert(one.diagnostics.allConverged, s"not converged: ${one.diagnostics}")
    assert(one.isRealPair(0))
    assertEqualsDouble(one.eigenvalue(0).re, 10.0, 1e-9)
    assertEqualsDouble(one.eigenvalue(0).im, 0.0, 1e-12)
    assert(opResidual(a, one, 0) < 1e-9, "real eigenvector residual too large")
  }

  // --- convergence boundaries ------------------------------------------------

  test("arnoldi: partial convergence returns Right with converged < requested") {
    val n = 40
    val a = randomReal(n, 909L)
    val k = 6
    // Minimal subspace + a single build + a punishing tolerance: one build cannot
    // resolve all k dominant pairs, so the result is partial (never a Left).
    val options = SpectralOptions(tolerance = 1e-13, maxIterations = 1, subspaceDimension = Some(k + 2))
    val d = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestMagnitude), options).toOption.get
    assert(d.diagnostics.requested >= k, s"requested ${d.diagnostics.requested}")
    assert(d.diagnostics.converged < d.diagnostics.requested, s"expected partial, got ${d.diagnostics}")
    assert(!d.diagnostics.allConverged)
    assertEquals(d.size, d.diagnostics.converged)
    assertPacking(d)
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  test("arnoldi: zero convergence returns Right(empty), never a Left") {
    val n = 30
    val a = randomReal(n, 8888L)
    val options = SpectralOptions(tolerance = 0.0, maxIterations = 1, subspaceDimension = Some(6))
    val d = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestMagnitude), options).toOption.get
    assertEquals(d.size, 0)
    assertEquals(d.diagnostics.converged, 0)
    assert(!d.diagnostics.allConverged)
    assert(d.diagnostics.requested >= 3)
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  test("arnoldi: shift-invert / target is deferred") {
    val op = blockRotationOp(Array(30.0, 8.0, 5.0, 3.0), 0.5)
    val n = 8
    val target = Some(SpectralTarget.ShiftInvert(0.5, LinearSolvePlan.Direct))
    Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude), SpectralOptions(), target) match
      case Left(_: LinAlgError.UnsupportedOperation) => ()
      case other                                     => fail(s"expected UnsupportedOperation, got $other")
  }

  // --- structural Lefts ------------------------------------------------------

  test("arnoldi: structural violations return the mapped Left") {
    val mags = Array(30.0, 8.0, 5.0, 3.0, 2.0, 1.0)
    val op = blockRotationOp(mags, 0.5)
    val n = 2 * mags.length
    val opts = SpectralOptions()

    // k out of range (k <= 0 and k >= n-1).
    assert(Eigen.eigNonsymmetric(op, n, EigenSelection.Count(0, EigenOrder.LargestMagnitude), opts).isLeft)
    assert(Eigen.eigNonsymmetric(op, n, EigenSelection.Count(n - 1, EigenOrder.LargestMagnitude), opts).isLeft)

    // Illegal (symmetric-only) orders.
    for order <- Seq(EigenOrder.LargestAlgebraic, EigenOrder.SmallestAlgebraic, EigenOrder.BothEnds) do
      Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, order), opts) match
        case Left(_: LinAlgError.InvalidArgument) => ()
        case other                                => fail(s"expected InvalidArgument for $order, got $other")

    // Non-Count selection is dense-only.
    assert(Eigen.eigNonsymmetric(op, n, EigenSelection.All, opts).isLeft)

    // Square operator whose shape disagrees with n -> DimensionMismatch.
    Eigen.eigNonsymmetric(op, n + 1, EigenSelection.Count(2, EigenOrder.LargestMagnitude), opts) match
      case Left(_: LinAlgError.DimensionMismatch) => ()
      case other                                  => fail(s"expected DimensionMismatch, got $other")

    // Genuinely non-square operator -> NonSquareMatrix.
    val rect = LinearOperator.fromFunction(n, n + 1)((_, _) => ())
    Eigen.eigNonsymmetric(rect, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude), opts) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected NonSquareMatrix, got $other")

    // Left / left-and-right vectors are deferred.
    for flag <- Seq(EigenVectors.Left, EigenVectors.LeftAndRight) do
      val leftOpts = SpectralOptions(returnVectors = flag)
      Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude), leftOpts) match
        case Left(_: LinAlgError.UnsupportedOperation) => ()
        case other                                     => fail(s"expected UnsupportedOperation for $flag, got $other")

    // Bad start vector: wrong length, then zero norm.
    val badLen = SpectralOptions(startVector = Some(DVec.zeros(n - 1)))
    assert(Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude), badLen).isLeft)
    val zero = SpectralOptions(startVector = Some(DVec.zeros(n)))
    assert(Eigen.eigNonsymmetric(op, n, EigenSelection.Count(2, EigenOrder.LargestMagnitude), zero).isLeft)
  }
