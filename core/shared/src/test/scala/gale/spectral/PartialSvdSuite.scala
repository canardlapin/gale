package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.linalg.MutableDVec

/** Tests for the partial SVD facade `Svd.svd` — Golub–Kahan–Lanczos
  * bidiagonalization solved via the Jordan–Wielandt augmented tridiagonal.
  * Convergence tests use well-separated spectra (the honest scope for plain GKL,
  * as with the Lanczos eigensolver); accuracy is checked by two-sided residuals,
  * U/V orthonormality, and cross-checking against the Gram matrix's eigenvalues.
  */
class PartialSvdSuite extends munit.FunSuite:

  // --- helpers ---------------------------------------------------------------

  private def rectDiag(m: Int, n: Int, d: Seq[Double]): DMat =
    Matrix.tabulate(m, n)((i, j) => if i == j && i < d.length then d(i) else 0.0)

  private def randomMat(m: Int, n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(m, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  /** A rectangular diagonal operator `diag(d)` implementing both `A·x` and `Aᵀ·y`,
    * to exercise the matrix-free path with a genuine (non-`DMat`) operator.
    */
  private def diagOperator(m: Int, n: Int, d: Array[Double]): DoubleLinearOperator =
    new DoubleLinearOperator:
      def rows: Int = m
      def cols: Int = n
      def applyTo(x: DVec, into: MutableDVec): Unit =
        var i = 0
        while i < m do
          into(i) = if i < d.length && i < n then d(i) * x(i) else 0.0
          i += 1
      override def transposeApplyTo(y: DVec, into: MutableDVec): Unit =
        var j = 0
        while j < n do
          into(j) = if j < d.length && j < m then d(j) * y(j) else 0.0
          j += 1

  private def values(s: SVD): IndexedSeq[Double] =
    (0 until s.size).map(s.singularValues(_))

  private def gramSpectrum(a: DMat): IndexedSeq[Double] =
    // Singular values of A are sqrt of the eigenvalues of AᵀA, descending.
    val gram = a.t * a
    val eig = Eigen.eigSymmetric(gram, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    (0 until eig.size).map(i => math.sqrt(math.max(eig.eigenvalues(i), 0.0))).reverse

  private def assertDescending(v: IndexedSeq[Double]): Unit =
    var i = 1
    while i < v.length do
      assert(v(i) <= v(i - 1) + 1e-12, s"not descending at $i: $v")
      i += 1

  private def assertClose(a: IndexedSeq[Double], b: IndexedSeq[Double], tol: Double): Unit =
    assertEquals(a.length, b.length, s"length $a vs $b")
    a.zip(b).foreach { case (x, y) => assert(math.abs(x - y) < tol, s"$x != $y in $a vs $b") }

  /** Worst two-sided residual max_i(max(‖A v_i − σ u_i‖, ‖Aᵀ u_i − σ v_i‖)). */
  private def worstTripletResidual(a: DMat, s: SVD): Double =
    var worst = 0.0
    var i = 0
    while i < s.size do
      val u = s.u.col(i)
      val v = s.vt.row(i)
      val sigma = s.singularValues(i)
      val rV = (a * v - u * sigma).norm2
      val rU = (a.t * u - v * sigma).norm2
      worst = math.max(worst, math.max(rV, rU))
      i += 1
    worst

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

  // --- fixed spectra ---------------------------------------------------------

  test("svd tall diagonal: known singular values, both orders, descending") {
    val a = rectDiag(8, 5, Seq(9.0, 7.0, 5.0, 3.0, 1.0))
    val largest = Svds.svd(a, SingularSelection.Count(3, SingularOrder.Largest)).toOption.get
    assert(largest.diagnostics.allConverged, s"largest: ${largest.diagnostics}")
    assertDescending(values(largest))
    assertClose(values(largest), IndexedSeq(9.0, 7.0, 5.0), 1e-8)
    assert(worstTripletResidual(a, largest) < 1e-8)
    assert(orthoError(largest.u) < 1e-9)
    assert(orthoError(largest.vt.t) < 1e-9)

    val smallest = Svds.svd(a, SingularSelection.Count(3, SingularOrder.Smallest)).toOption.get
    assert(smallest.diagnostics.allConverged, s"smallest: ${smallest.diagnostics}")
    assertDescending(values(smallest)) // still descending, even for Smallest membership
    assertClose(values(smallest), IndexedSeq(5.0, 3.0, 1.0), 1e-8)
    assert(worstTripletResidual(a, smallest) < 1e-8)
  }

  test("svd wide diagonal: min(m,n) bounds and known singular values") {
    val a = rectDiag(4, 9, Seq(6.0, 4.0, 2.0, 1.0))
    val d = Svds.svd(a, SingularSelection.Count(2, SingularOrder.Largest)).toOption.get
    assert(d.diagnostics.allConverged, s"${d.diagnostics}")
    assertClose(values(d), IndexedSeq(6.0, 4.0), 1e-8)
    assertEquals(d.u.rows, 4)
    assertEquals(d.vt.cols, 9)
    assert(worstTripletResidual(a, d) < 1e-8)
  }

  // --- random matrices: Gram cross-check -------------------------------------

  test("svd random tall: largest values match the Gram spectrum, residuals small") {
    val a = randomMat(20, 8, 31L)
    val k = 4
    val d = Svds.svd(a, SingularSelection.Count(k, SingularOrder.Largest)).toOption.get
    assert(d.diagnostics.allConverged, s"${d.diagnostics}")
    assertDescending(values(d))
    assertClose(values(d), gramSpectrum(a).take(k), 1e-6)
    assert(worstTripletResidual(a, d) < 1e-6)
    assert(orthoError(d.u) < 1e-8)
    assert(orthoError(d.vt.t) < 1e-8)
  }

  test("svd random wide: largest values match the Gram spectrum") {
    val a = randomMat(7, 18, 42L)
    val k = 3
    val d = Svds.svd(a, SingularSelection.Count(k, SingularOrder.Largest)).toOption.get
    assert(d.diagnostics.allConverged, s"${d.diagnostics}")
    assertClose(values(d), gramSpectrum(a).take(k), 1e-6)
    assert(worstTripletResidual(a, d) < 1e-6)
    // Orthonormality holds on the wide (transposed) path too.
    assert(orthoError(d.u) < 1e-8)
    assert(orthoError(d.vt.t) < 1e-8)
  }

  // --- rank deficiency -------------------------------------------------------

  test("svd rank-deficient: near-zero singular value and rank diagnostic") {
    // Smallest singular value is a genuine near-zero; the augmented-tridiagonal
    // form (not normal equations) keeps it usable. Rank counts σ > tol·σmax among
    // the returned values, so the tiny σ is excluded.
    val a = rectDiag(6, 4, Seq(5.0, 3.0, 1.0, 1e-8))
    val options = SpectralOptions(tolerance = 1e-6)
    val d = Svds.svd(a, 6, 4, SingularSelection.Count(2, SingularOrder.Smallest), options).toOption.get
    assert(d.diagnostics.allConverged, s"${d.diagnostics}")
    assertDescending(values(d))
    assert(math.abs(values(d)(0) - 1.0) < 1e-6, s"${values(d)}")
    assert(math.abs(values(d)(1) - 1e-8) < 1e-9, s"${values(d)}")
    assertEquals(d.rank, 1) // only σ = 1 exceeds tol·σmax; the 1e-8 is deficient
    assertEquals(d.diagnostics.rank, Some(1))
  }

  // --- convergence semantics -------------------------------------------------

  test("svd partial convergence: Right with converged < requested, requireConverged Left") {
    val a = randomMat(20, 12, 71L)
    val k = 4
    // Minimal subspace + one build + a punishing tolerance -> not all converge.
    val options = SpectralOptions(tolerance = 1e-13, maxIterations = 1, subspaceDimension = Some(k + 1))
    val d = Svds.svd(a, 20, 12, SingularSelection.Count(k, SingularOrder.Largest), options).toOption.get
    assertEquals(d.diagnostics.requested, k)
    assert(d.diagnostics.converged < k, s"expected partial, got ${d.diagnostics.converged}")
    assert(!d.diagnostics.allConverged)
    assertEquals(d.size, d.diagnostics.converged)
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  test("svd zero convergence: Right(empty), never a Left") {
    val a = randomMat(15, 10, 88L)
    val options = SpectralOptions(tolerance = 0.0, maxIterations = 1, subspaceDimension = Some(4))
    val result = Svds.svd(a, 15, 10, SingularSelection.Count(3, SingularOrder.Largest), options)
    val d = result.toOption.get
    assertEquals(d.size, 0)
    assertEquals(d.diagnostics.converged, 0)
    assert(!d.diagnostics.allConverged)
    assertEquals(d.diagnostics.requested, 3)
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  // --- matrix-free operator --------------------------------------------------

  test("svd matrix-free rectangular operator resolves its diagonal") {
    val d = Array(8.0, 6.0, 4.0, 2.0)
    val op = diagOperator(9, 4, d)
    val result = Svds.svd(op, 9, 4, SingularSelection.Count(3, SingularOrder.Largest)).toOption.get
    assert(result.diagnostics.allConverged, s"${result.diagnostics}")
    assertClose(values(result), IndexedSeq(8.0, 6.0, 4.0), 1e-7)
    // op y = σ x residual via the operator directly.
    var i = 0
    while i < result.size do
      val u = result.u.col(i)
      val v = result.vt.row(i)
      val av = MutableDVec.zeros(9)
      op.applyTo(v, av)
      av.axpyInPlace(-result.singularValues(i), u)
      assert(av.asVec.norm2 < 1e-6)
      i += 1
  }

  test("svd matrix-free WIDE operator exercises the transpose-oriented path") {
    // m < n forces GKL onto Aᵀ via transposeOp; a genuine operator (asymmetric
    // applyTo/transposeApplyTo, non-DMat) is exactly where a transpose-wiring or
    // U/V-swap bug would hide.
    val d = Array(9.0, 5.0, 3.0)
    val op = diagOperator(4, 9, d)
    val result = Svds.svd(op, 4, 9, SingularSelection.Count(3, SingularOrder.Largest)).toOption.get
    assert(result.diagnostics.allConverged, s"${result.diagnostics}")
    assertClose(values(result), IndexedSeq(9.0, 5.0, 3.0), 1e-7)
    // Two-sided residuals against the ORIGINAL operator: A v = σ u and Aᵀ u = σ v.
    var i = 0
    while i < result.size do
      val u = result.u.col(i)
      val v = result.vt.row(i)
      val av = MutableDVec.zeros(4)
      op.applyTo(v, av)
      av.axpyInPlace(-result.singularValues(i), u)
      assert(av.asVec.norm2 < 1e-6, s"A v != sigma u at $i")
      val atu = MutableDVec.zeros(9)
      op.transposeApplyTo(u, atu)
      atu.axpyInPlace(-result.singularValues(i), v)
      assert(atu.asVec.norm2 < 1e-6, s"At u != sigma v at $i")
      i += 1
  }

  test("svd exact rank deficiency: pending-right breakdown keeps the trailing coupling") {
    // A true null space (rank 3 < p = 4) plus a start vector with a null-space
    // component forces the α-breakdown-with-pending-right path: the projected
    // problem must be the rectangular mEff×(mEff+1) bidiagonal (odd-dimensional
    // augmented form). Truncating it square silently solves the wrong matrix and
    // converges nothing — the original latent bug this test pins.
    val d = Array(8.0, 6.0, 4.0)
    val tall = Svds.svd(diagOperator(9, 4, d), 9, 4, SingularSelection.Count(3, SingularOrder.Largest)).toOption.get
    assert(tall.diagnostics.allConverged, s"tall rank-deficient: ${tall.diagnostics}")
    assertClose(values(tall), IndexedSeq(8.0, 6.0, 4.0), 1e-9)
    assert(tall.diagnostics.worstResidual < 1e-9)
  }

  // --- values-only + explicit start ------------------------------------------

  test("svd values-only: no vectors, values still correct") {
    val a = rectDiag(10, 6, Seq(7.0, 5.0, 3.0, 2.0, 1.0, 0.5))
    val d = Svds.svd(a, SingularSelection.Count(3, SingularOrder.Largest), EigenVectors.ValuesOnly).toOption.get
    assert(d.diagnostics.allConverged)
    assertEquals(d.u.cols, 0)
    assertEquals(d.vt.rows, 0)
    assertClose(values(d), IndexedSeq(7.0, 5.0, 3.0), 1e-8)
  }

  test("svd explicit start vector") {
    val a = randomMat(16, 8, 99L)
    val start = DVec.tabulate(8)(i => if i % 2 == 0 then 1.0 else 0.5)
    val options = SpectralOptions(startVector = Some(start))
    val d = Svds.svd(a, 16, 8, SingularSelection.Count(3, SingularOrder.Largest), options).toOption.get
    assert(d.diagnostics.allConverged, s"${d.diagnostics}")
    assertClose(values(d), gramSpectrum(a).take(3), 1e-6)
  }

  // --- structural violations -------------------------------------------------

  test("svd structural violations return Left") {
    val a = randomMat(10, 6, 5L)

    // k out of range: k <= 0 and k > min(m,n).
    assert(Svds.svd(a, SingularSelection.Count(0, SingularOrder.Largest)).isLeft)
    assert(Svds.svd(a, SingularSelection.Count(7, SingularOrder.Largest)).isLeft) // k > min(10,6)=6
    // k = min(m,n) and All are served by the dense bidiagonal kernel now.
    assert(Svds.svd(a, SingularSelection.Count(6, SingularOrder.Largest)).isRight)
    assert(Svds.svd(a, SingularSelection.All).isRight)
    // The operator (matrix-free) path stays Count-only: All is still rejected.
    Svds.svd(a.asLinearOperator, 10, 6, SingularSelection.All, SpectralOptions()) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected InvalidArgument for operator All, got $other")
    // illegal vector flags.
    assert(Svds.svd(a, SingularSelection.Count(2, SingularOrder.Largest), EigenVectors.Left).isLeft)
    assert(Svds.svd(a, SingularSelection.Count(2, SingularOrder.Largest), EigenVectors.LeftAndRight).isLeft)
    // operator dimension disagreement -> DimensionMismatch (rectangular is legal).
    Svds.svd(a, 11, 6, SingularSelection.Count(2, SingularOrder.Largest), SpectralOptions()) match
      case Left(_: LinAlgError.DimensionMismatch) => ()
      case other                                  => fail(s"expected DimensionMismatch, got $other")
    // bad start vectors.
    val badLen = SpectralOptions(startVector = Some(DVec.zeros(5)))
    assert(Svds.svd(a, 10, 6, SingularSelection.Count(2, SingularOrder.Largest), badLen).isLeft)
    val zeroStart = SpectralOptions(startVector = Some(DVec.zeros(6)))
    assert(Svds.svd(a, 10, 6, SingularSelection.Count(2, SingularOrder.Largest), zeroStart).isLeft)
  }
