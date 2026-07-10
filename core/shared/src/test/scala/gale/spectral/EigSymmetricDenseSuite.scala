package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix

/** Tests for the dense symmetric eigendecomposition facade `Eigen.eigSymmetric`:
  * selection slicing consistency, fixed spectra, and boundary validation.
  */
class EigSymmetricDenseSuite extends munit.FunSuite:

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

  private def values(d: EigenDecomposition): IndexedSeq[Double] =
    (0 until d.eigenvalues.length).map(d.eigenvalues(_))

  private def fullSpectrum(a: DMat): IndexedSeq[Double] =
    values(Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get)

  private def assertAscending(v: IndexedSeq[Double]): Unit =
    var i = 1
    while i < v.length do
      assert(v(i) >= v(i - 1) - 1e-12, s"not ascending at $i")
      i += 1

  private def assertClose(a: IndexedSeq[Double], b: IndexedSeq[Double], tol: Double): Unit =
    assertEquals(a.length, b.length, s"length $a vs $b")
    a.zip(b).zipWithIndex.foreach { case ((x, y), i) =>
      assert(math.abs(x - y) < tol, s"[$i] $x != $y")
    }

  test("dense All: full ascending spectrum, small residual and orthogonality") {
    val n = 10
    val a = laplacian(n)
    val d = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.Right).toOption.get
    assertEquals(d.size, n)
    assertAscending(values(d))
    // analytic Laplacian spectrum
    val analytic = (0 until n).map(k => 2.0 - 2.0 * math.cos((k + 1) * math.Pi / (n + 1)))
    assertClose(values(d), analytic, 1e-10)
    assert(d.diagnostics.worstResidual < 1e-10)
    assert(d.diagnostics.orthogonalityError < 1e-10)
    assert(d.diagnostics.allConverged)
  }

  test("dense Count: each order is a slice of the full ascending spectrum") {
    val n = 14
    val a = randomSymmetric(n, 11L)
    val full = fullSpectrum(a)
    val k = 5

    def selected(order: EigenOrder): IndexedSeq[Double] =
      values(Eigen.eigSymmetric(a, EigenSelection.Count(k, order), EigenVectors.ValuesOnly).toOption.get)

    assertClose(selected(EigenOrder.SmallestAlgebraic), full.take(k), 1e-12)
    assertClose(selected(EigenOrder.LargestAlgebraic), full.takeRight(k), 1e-12)

    // BothEnds: floor(k/2) low + ceil(k/2) high.
    val both = full.take(k / 2) ++ full.takeRight((k + 1) / 2)
    assertClose(selected(EigenOrder.BothEnds), both.sorted, 1e-12)

    // Magnitude orders picked independently, output ascending.
    val byMagDesc = full.sortBy(v => -math.abs(v)).take(k).sorted
    val byMagAsc = full.sortBy(v => math.abs(v)).take(k).sorted
    assertClose(selected(EigenOrder.LargestMagnitude), byMagDesc, 1e-12)
    assertClose(selected(EigenOrder.SmallestMagnitude), byMagAsc, 1e-12)
  }

  test("dense Count: eigenvectors are the matching columns of the full solve") {
    val n = 8
    val a = randomSymmetric(n, 22L)
    val d = Eigen.eigSymmetric(a, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), EigenVectors.Right).toOption.get
    assertEquals(d.size, 3)
    assert(d.diagnostics.worstResidual < 1e-10)
    assert(d.diagnostics.orthogonalityError < 1e-10)
    // A v = λ v for each returned pair.
    var c = 0
    while c < 3 do
      val v = d.eigenvectors.col(c)
      val r = a * v - v * d.eigenvalues(c)
      assert(r.norm2 < 1e-9)
      c += 1
  }

  test("dense IndexRange: ascending-rank window") {
    val n = 12
    val a = randomSymmetric(n, 33L)
    val full = fullSpectrum(a)
    val d = Eigen.eigSymmetric(a, EigenSelection.IndexRange(3, 7), EigenVectors.ValuesOnly).toOption.get
    assertClose(values(d), full.slice(3, 8), 1e-12)
  }

  test("dense ValueInterval: half-open (lower, upper]") {
    val n = 12
    val a = randomSymmetric(n, 44L)
    val full = fullSpectrum(a)
    val lower = -0.5
    val upper = 1.5
    val d = Eigen.eigSymmetric(a, EigenSelection.ValueInterval(lower, upper), EigenVectors.ValuesOnly).toOption.get
    val expected = full.filter(v => v > lower && v <= upper)
    assertClose(values(d), expected, 1e-12)
  }

  test("dense ValuesOnly: empty eigenvectors, zero residuals") {
    val n = 6
    val a = randomSymmetric(n, 55L)
    val d = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    assertEquals(d.eigenvectors.cols, 0)
    assertEquals(d.diagnostics.worstResidual, 0.0)
    assertEquals(d.diagnostics.orthogonalityError, 0.0)
  }

  test("dense: structural and selection violations return Left") {
    val a = randomSymmetric(6, 66L)
    val rect = Matrix.dense(2, 3)(1.0, 0.0, 0.0, 0.0, 1.0, 0.0)

    assert(Eigen.eigSymmetric(rect, EigenSelection.All, EigenVectors.Right).isLeft)
    assertEquals(
      Eigen.eigSymmetric(a, EigenSelection.Count(0, EigenOrder.LargestAlgebraic), EigenVectors.Right).left.toOption.get.getClass,
      classOf[LinAlgError.InvalidArgument]
    )
    assert(Eigen.eigSymmetric(a, EigenSelection.Count(7, EigenOrder.LargestAlgebraic), EigenVectors.Right).isLeft)
    assert(Eigen.eigSymmetric(a, EigenSelection.Count(3, EigenOrder.LargestRealPart), EigenVectors.Right).isLeft)
    assert(Eigen.eigSymmetric(a, EigenSelection.Count(3, EigenOrder.SmallestRealPart), EigenVectors.Right).isLeft)
    assert(Eigen.eigSymmetric(a, EigenSelection.IndexRange(-1, 2), EigenVectors.Right).isLeft)
    assert(Eigen.eigSymmetric(a, EigenSelection.IndexRange(2, 6), EigenVectors.Right).isLeft)
    assert(Eigen.eigSymmetric(a, EigenSelection.ValueInterval(2.0, 1.0), EigenVectors.Right).isLeft)
    assert(Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.Left).isLeft)
    assert(Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.LeftAndRight).isLeft)
  }

  test("dense: default overload computes vectors") {
    val a = randomSymmetric(5, 77L)
    val d = Eigen.eigSymmetric(a, EigenSelection.All).toOption.get
    assertEquals(d.eigenvectors.cols, 5)
  }

  test("dense: residual diagnostics measure the symmetrized (decomposed) matrix") {
    // Only the lower triangle is read; the strict upper triangle is garbage here.
    // The decomposition is of the symmetrized matrix, so the reported residuals
    // must be tiny even though `a * v` with the garbage triangle would not be.
    val n = 6
    val sym = randomSymmetric(n, 4242L)
    val garbled = Matrix.tabulate(n, n)((i, j) => if i >= j then sym(i, j) else 999.0 + i + j)
    val d = Eigen.eigSymmetric(garbled, EigenSelection.All).toOption.get
    val reference = Eigen.eigSymmetric(sym, EigenSelection.All).toOption.get
    assert(d.diagnostics.worstResidual < 1e-10, s"residual ${d.diagnostics.worstResidual}")
    var i = 0
    while i < n do
      assert(math.abs(d.eigenvalues(i) - reference.eigenvalues(i)) < 1e-12)
      i += 1
  }
