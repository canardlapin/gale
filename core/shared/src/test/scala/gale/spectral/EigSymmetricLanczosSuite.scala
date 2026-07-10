package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.LinearOperator
import gale.linalg.Matrix

/** Tests for the iterative partial symmetric eigensolver
  * `Eigen.eigSymmetric(op, n, Count, options)` — Lanczos with full
  * reorthogonalization. Convergence is checked against the dense solve's extremes.
  */
class EigSymmetricLanczosSuite extends munit.FunSuite:

  private def randomSymmetric(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        val v = rng.nextDouble() * 2.0 - 1.0
        a(i)(j) = v
        a(j)(i) = v
        j += 1
      i += 1
    Matrix.tabulate(n, n)((r, c) => a(r)(c))

  private def laplacian(n: Int): DMat =
    Matrix.tabulate(n, n): (i, j) =>
      if i == j then 2.0
      else if math.abs(i - j) == 1 then -1.0
      else 0.0

  /** Matrix-free diagonal operator `diag(d)` — a well-separated spectrum to
    * exercise the operator-only path (plain Lanczos resolves separated extremes
    * without the implicit restarting that clustered spectra would need).
    */
  private def diagonalOperator(d: Array[Double]): DoubleLinearOperator =
    val n = d.length
    LinearOperator.fromFunction(n, n): (x, into) =>
      var i = 0
      while i < n do
        into(i) = d(i) * x(i)
        i += 1

  private def values(d: EigenDecomposition): IndexedSeq[Double] =
    (0 until d.eigenvalues.length).map(d.eigenvalues(_))

  private def denseSpectrum(a: DMat): IndexedSeq[Double] =
    val d = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    (0 until d.size).map(d.eigenvalues(_))

  private def assertAscending(v: IndexedSeq[Double]): Unit =
    var i = 1
    while i < v.length do
      assert(v(i) >= v(i - 1) - 1e-10, s"not ascending at $i: $v")
      i += 1

  private def assertClose(a: IndexedSeq[Double], b: IndexedSeq[Double], tol: Double): Unit =
    assertEquals(a.length, b.length, s"length $a vs $b")
    a.zip(b).foreach { case (x, y) => assert(math.abs(x - y) < tol, s"$x != $y in $a vs $b") }

  test("lanczos: largest/smallest algebraic match the dense extremes") {
    val n = 60
    val a = randomSymmetric(n, 3131L)
    val full = denseSpectrum(a)
    val k = 6

    val largest = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestAlgebraic)).toOption.get
    assert(largest.diagnostics.allConverged, s"largest not converged: ${largest.diagnostics}")
    assertAscending(values(largest))
    assertClose(values(largest), full.takeRight(k), 1e-6)
    assert(largest.diagnostics.worstResidual < 1e-6)

    val smallest = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic)).toOption.get
    assert(smallest.diagnostics.allConverged, s"smallest not converged: ${smallest.diagnostics}")
    assertClose(values(smallest), full.take(k), 1e-6)
  }

  test("lanczos: BothEnds returns floor/ceil split of both extremes") {
    val n = 50
    val a = randomSymmetric(n, 5150L)
    val full = denseSpectrum(a)
    val k = 6
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.BothEnds)).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    val expected = (full.take(k / 2) ++ full.takeRight((k + 1) / 2)).sorted
    assertClose(values(d), expected, 1e-6)
  }

  test("lanczos: largest magnitude on a sign-straddling spectrum") {
    val n = 50
    val a = randomSymmetric(n, 2026L)
    val full = denseSpectrum(a)
    val k = 5
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestMagnitude)).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    val expected = full.sortBy(v => -math.abs(v)).take(k).sorted
    assertClose(values(d), expected, 1e-6)
  }

  test("lanczos: matrix-free diagonal operator resolves the smallest eigenvalues") {
    val n = 40
    val d = Array.tabulate(n)(i => (i + 1).toDouble) // eigenvalues 1..n
    val op = diagonalOperator(d)
    val k = 4
    val result = Eigen.eigSymmetric(op, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic)).toOption.get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(1.0, 2.0, 3.0, 4.0), 1e-7)
    // Eigenvectors satisfy op y = θ y.
    var c = 0
    while c < k do
      val y = result.eigenvectors.col(c)
      val r = op * y - y * result.eigenvalues(c)
      assert(r.norm2 < 1e-6)
      c += 1
  }

  test("lanczos: partial convergence returns Right with converged < requested") {
    val n = 40
    val a = randomSymmetric(n, 909L)
    val k = 6
    // Minimal subspace + a single restart + a punishing tolerance: one build
    // cannot resolve all k pairs, so the result is partial (never a Left).
    val options = SpectralOptions(tolerance = 1e-13, maxIterations = 1, subspaceDimension = Some(k + 1))
    val result = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestAlgebraic), options)
    val d = result.toOption.get
    assertEquals(d.diagnostics.requested, k)
    assert(d.diagnostics.converged < k, s"expected partial, got ${d.diagnostics.converged}")
    assert(!d.diagnostics.allConverged)
    assertEquals(d.size, d.diagnostics.converged)
    // requireConverged is the caller's opt-in to fail-fast.
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  test("lanczos: values-only omits eigenvectors but still converges") {
    val n = 50
    val a = randomSymmetric(n, 6161L)
    val k = 4
    val options = SpectralOptions(returnVectors = EigenVectors.ValuesOnly)
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestAlgebraic), options).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assertEquals(d.eigenvectors.cols, 0)
    assertClose(values(d), denseSpectrum(a).takeRight(k), 1e-6)
  }

  test("lanczos: explicit start vector") {
    val n = 40
    val a = randomSymmetric(n, 7171L)
    val k = 4
    val start = DVec.tabulate(n)(i => if i % 2 == 0 then 1.0 else 0.5)
    val options = SpectralOptions(startVector = Some(start))
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic), options).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assertClose(values(d), denseSpectrum(a).take(k), 1e-6)
  }

  test("lanczos: structural violations return Left") {
    val n = 20
    val a = laplacian(n)
    val opts = SpectralOptions()

    // k out of range.
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(0, EigenOrder.LargestAlgebraic), opts).isLeft)
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(n, EigenOrder.LargestAlgebraic), opts).isLeft)
    // illegal order for symmetric.
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestRealPart), opts).isLeft)
    // dimension disagreement.
    assert(Eigen.eigSymmetric(a, n + 1, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), opts).isLeft)
    // illegal vector flag.
    val leftVecs = SpectralOptions(returnVectors = EigenVectors.Left)
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), leftVecs).isLeft)
    // wrong-length start vector.
    val badStart = SpectralOptions(startVector = Some(DVec.zeros(n - 1)))
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), badStart).isLeft)
    // zero start vector.
    val zeroStart = SpectralOptions(startVector = Some(DVec.zeros(n)))
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), zeroStart).isLeft)
  }

  test("lanczos: repeated eigenvalue collapses to one copy (documented contract)") {
    // Characterizing test for the scaladoc'd multiplicity caveat: single-vector
    // Lanczos resolves at most one eigenvector per DISTINCT eigenvalue. For
    // diag(5,5,3,2,1) and k=3 largest-algebraic, the returned set is the top-3 of
    // the DISTINCT spectrum {2,3,5} — not the multiset [3,5,5] — and the result
    // still reports allConverged (each returned pair has a genuinely tiny
    // residual). If this ever starts returning the true multiset, the caveat can
    // be deleted; update this test rather than "fixing" the solver silently.
    val d = Array(5.0, 5.0, 3.0, 2.0, 1.0)
    val op = diagonalOperator(d)
    val result = Eigen.eigSymmetric(op, 5, EigenSelection.Count(3, EigenOrder.LargestAlgebraic)).toOption.get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(2.0, 3.0, 5.0), 1e-9)
  }

  test("lanczos: clustered spectrum converges via subspace growth (iterations > 1)") {
    // A 1e-6-wide cluster at the bottom is unresolvable in the initial minimal
    // subspace; the growth strategy (not a lossy restart) must separate it.
    val n = 40
    val d = Array.tabulate(n)(i => if i == 0 then 1.0 else if i == 1 then 1.0 + 1e-6 else (i + 1).toDouble)
    val op = diagonalOperator(d)
    val k = 2
    val options = SpectralOptions(tolerance = 1e-9, maxIterations = 50, subspaceDimension = Some(k + 2))
    val result = Eigen.eigSymmetric(op, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic), options).toOption.get
    assert(result.diagnostics.allConverged, s"cluster not resolved: ${result.diagnostics}")
    assert(result.diagnostics.iterations > 1, s"expected growth, converged in ${result.diagnostics.iterations} build(s)")
    assertClose(values(result), IndexedSeq(1.0, 1.0 + 1e-6), 1e-8)
  }

  test("lanczos: zero convergence returns Right(empty), never a Left") {
    // An unsatisfiable tolerance: no pair can converge, so the result is an empty
    // Right with converged == 0 — the parity doc's zero-convergence boundary.
    val n = 30
    val a = randomSymmetric(n, 8888L)
    val options = SpectralOptions(tolerance = 0.0, maxIterations = 1, subspaceDimension = Some(5))
    val result = Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), options)
    val d = result.toOption.get
    assertEquals(d.size, 0)
    assertEquals(d.diagnostics.converged, 0)
    assert(!d.diagnostics.allConverged)
    assertEquals(d.diagnostics.requested, 3)
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  test("lanczos: shift-invert / target is deferred") {
    val n = 20
    val a = laplacian(n)
    val target = Some(SpectralTarget.ShiftInvert(0.5, LinearSolvePlan.Direct))
    val result = Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.SmallestAlgebraic), SpectralOptions(), target)
    result match
      case Left(_: LinAlgError.UnsupportedOperation) => ()
      case other                                     => fail(s"expected UnsupportedOperation, got $other")
  }
