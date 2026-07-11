package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.linalg.TriangularSolve

/** Tests for the dense generalized symmetric-definite solver
  * `Eigen.eigSymmetricGeneralized(A, B, selection, vectors)` — `A x = λ B x` with
  * `B` SPD, reduced via Cholesky to the standard symmetric problem. Checks the
  * analytic diagonal spectrum, true generalized residuals, `B`-orthonormality
  * (`XᵀBX = I`), ascending order, a cross-check against the naively-formed reduced
  * matrix, selection slicing, the lower-triangle read convention, and the
  * structural `Left`s.
  */
class EigSymmetricGeneralizedSuite extends munit.FunSuite:

  // --- helpers ---------------------------------------------------------------

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

  /** A well-conditioned SPD matrix `M Mᵀ + n·I`. */
  private def randomSPD(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    val m = Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)
    val mmt = m * m.t
    Matrix.tabulate(n, n)((i, j) => if i == j then mmt(i, j) + n.toDouble else mmt(i, j))

  private def mirror(m: DMat): DMat =
    val n = m.rows
    Matrix.tabulate(n, n)((i, j) => if i >= j then m(i, j) else m(j, i))

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

  private def values(d: EigenDecomposition): IndexedSeq[Double] =
    (0 until d.size).map(d.eigenvalues(_))

  private def assertAscending(v: IndexedSeq[Double]): Unit =
    var i = 1
    while i < v.length do
      assert(v(i) >= v(i - 1) - 1e-10, s"not ascending at $i: $v")
      i += 1

  /** ‖A x − λ B x‖ for the c-th returned generalized eigenpair. */
  private def generalizedResidual(aSym: DMat, bSym: DMat, d: EigenDecomposition, c: Int): Double =
    val x = d.eigenvectors.col(c)
    ((aSym * x) - (bSym * x) * d.eigenvalues(c)).norm2

  /** ‖XᵀBX − I‖_F of the returned eigenvectors. */
  private def bOrthoError(d: EigenDecomposition, bSym: DMat): Double =
    val g = d.eigenvectors.t * (bSym * d.eigenvectors)
    frobenius(g - Matrix.eye(g.rows))

  // --- analytic diagonal spectrum --------------------------------------------

  test("generalized: diagonal A, B give λ = a_i / b_i") {
    // A = diag(2,6,3), B = diag(1,2,3) ⇒ λ = {2, 3, 1}, ascending {1, 2, 3}.
    val a = Matrix.tabulate(3, 3)((i, j) => if i == j then Array(2.0, 6.0, 3.0)(i) else 0.0)
    val b = Matrix.tabulate(3, 3)((i, j) => if i == j then Array(1.0, 2.0, 3.0)(i) else 0.0)
    val d = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).toOption.get
    assertEquals(d.size, 3)
    val vs = values(d)
    assertAscending(vs)
    assert(math.abs(vs(0) - 1.0) < 1e-12, s"$vs")
    assert(math.abs(vs(1) - 2.0) < 1e-12, s"$vs")
    assert(math.abs(vs(2) - 3.0) < 1e-12, s"$vs")
    // B-orthonormal eigenvectors, tiny generalized residuals.
    assert(bOrthoError(d, b) < 1e-12, s"XᵀBX≠I: ${bOrthoError(d, b)}")
    var c = 0
    while c < 3 do
      assert(generalizedResidual(a, b, d, c) < 1e-12, s"residual $c")
      c += 1
  }

  // --- random SPD B + random symmetric A -------------------------------------

  test("generalized: random pencil — residuals tiny, B-orthonormal, ascending") {
    val n = 10
    val a = randomSymmetric(n, 4242L)
    val b = randomSPD(n, 9001L)
    val d = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).toOption.get
    assertEquals(d.size, n)
    assertAscending(values(d))
    val scale = math.max(frobenius(a), frobenius(b))
    var c = 0
    while c < n do
      assert(generalizedResidual(a, b, d, c) < 1e-8 * scale, s"residual $c = ${generalizedResidual(a, b, d, c)}")
      c += 1
    assert(bOrthoError(d, b) < 1e-10, s"XᵀBX≠I: ${bOrthoError(d, b)}")
    // The diagnostics residuals agree with the independently computed ones.
    c = 0
    while c < n do
      assert(math.abs(d.diagnostics.residuals(c) - generalizedResidual(a, b, d, c)) < 1e-12)
      c += 1
    assert(d.diagnostics.allConverged)
  }

  test("generalized: eigenvalues match the naively-reduced standard problem") {
    // Independent cross-check: form C = L⁻¹ A L⁻ᵀ by explicitly inverting L (a
    // different formula from the production solve-LY=A / solve-LC=Yᵀ path) and run
    // the standard symmetric solver on it.
    val n = 7
    val a = randomSymmetric(n, 314L)
    val b = randomSPD(n, 271L)
    val aSym = mirror(a)
    val l = b.cholesky.toOption.get.lower
    val linvCols = (0 until n).map: j =>
      TriangularSolve.lower(l, DVec.tabulate(n)(i => if i == j then 1.0 else 0.0)).toOption.get
    val linv = Matrix.tabulate(n, n)((i, j) => linvCols(j)(i))
    val cNaive = linv * aSym * linv.t
    val standard = Eigen.eigSymmetric(cNaive, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get

    val generalized = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    assertEquals(generalized.size, n)
    var i = 0
    while i < n do
      assert(math.abs(generalized.eigenvalues(i) - standard.eigenvalues(i)) < 1e-9, s"λ$i mismatch")
      i += 1
  }

  // --- selection slicing consistency vs All ----------------------------------

  test("generalized: selections are slices of the full ascending spectrum") {
    val n = 8
    val a = randomSymmetric(n, 55L)
    val b = randomSPD(n, 66L)
    val full = values(Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).toOption.get)
    val k = 3

    val smallest = values(Eigen.eigSymmetricGeneralized(a, b, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic)).toOption.get)
    smallest.zip(full.take(k)).foreach { case (x, y) => assert(math.abs(x - y) < 1e-10, s"$x != $y") }

    val largest = values(Eigen.eigSymmetricGeneralized(a, b, EigenSelection.Count(k, EigenOrder.LargestAlgebraic)).toOption.get)
    largest.zip(full.takeRight(k)).foreach { case (x, y) => assert(math.abs(x - y) < 1e-10, s"$x != $y") }

    val range = values(Eigen.eigSymmetricGeneralized(a, b, EigenSelection.IndexRange(2, 5)).toOption.get)
    range.zip(full.slice(2, 6)).foreach { case (x, y) => assert(math.abs(x - y) < 1e-10, s"$x != $y") }

    // ValueInterval: every eigenvalue in (lo, hi].
    val lo = full(1)
    val hi = full(n - 2)
    val interval = values(Eigen.eigSymmetricGeneralized(a, b, EigenSelection.ValueInterval(lo, hi)).toOption.get)
    val expected = full.filter(v => v > lo && v <= hi)
    assertEquals(interval.length, expected.length)
    interval.zip(expected).foreach { case (x, y) => assert(math.abs(x - y) < 1e-10, s"$x != $y") }
  }

  // --- lower-triangle read convention ----------------------------------------

  test("generalized: only the lower triangles of A and B are read") {
    val n = 6
    val cleanA = randomSymmetric(n, 700L)
    val cleanB = randomSPD(n, 800L)
    // Same lower triangles, garbage strict-upper triangles.
    val garbageA = Matrix.tabulate(n, n)((i, j) => if i >= j then cleanA(i, j) else 987.0)
    val garbageB = Matrix.tabulate(n, n)((i, j) => if i >= j then cleanB(i, j) else -654.0)

    val clean = Eigen.eigSymmetricGeneralized(cleanA, cleanB, EigenSelection.All).toOption.get
    val garbage = Eigen.eigSymmetricGeneralized(garbageA, garbageB, EigenSelection.All).toOption.get
    // Identical lower triangles ⇒ identical results, bit for bit.
    var i = 0
    while i < n do
      assertEquals(garbage.eigenvalues(i), clean.eigenvalues(i), s"λ$i")
      i += 1
  }

  // --- values-only -----------------------------------------------------------

  test("generalized: values-only omits eigenvectors and quality metrics") {
    val n = 8
    val a = randomSymmetric(n, 11L)
    val b = randomSPD(n, 22L)
    val full = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).toOption.get
    val vo = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    assertEquals(vo.eigenvectors.cols, 0)
    assertEquals(vo.diagnostics.orthogonalityError, 0.0)
    var i = 0
    while i < n do
      assert(math.abs(vo.eigenvalues(i) - full.eigenvalues(i)) < 1e-10, s"λ$i")
      i += 1
  }

  // --- structural Lefts ------------------------------------------------------

  test("generalized: B not positive-definite returns Left(NotPositiveDefinite)") {
    val a = Matrix.tabulate(3, 3)((i, j) => if i == j then 1.0 else 0.0)
    // Indefinite B (a negative diagonal entry): Cholesky fails.
    val b = Matrix.tabulate(3, 3)((i, j) => if i == j then Array(1.0, -1.0, 1.0)(i) else 0.0)
    Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All) match
      case Left(_: LinAlgError.NotPositiveDefinite) => ()
      case other                                    => fail(s"expected NotPositiveDefinite, got $other")
  }

  test("generalized: shape violations return the mapped Left") {
    val b3 = randomSPD(3, 1L)
    val a3 = randomSymmetric(3, 2L)

    // Non-square A / B.
    Eigen.eigSymmetricGeneralized(Matrix.tabulate(3, 4)((_, _) => 1.0), b3, EigenSelection.All) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected NonSquareMatrix (A), got $other")
    Eigen.eigSymmetricGeneralized(a3, Matrix.tabulate(3, 4)((_, _) => 1.0), EigenSelection.All) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected NonSquareMatrix (B), got $other")

    // A and B disagree in size.
    Eigen.eigSymmetricGeneralized(a3, randomSPD(4, 3L), EigenSelection.All) match
      case Left(_: LinAlgError.DimensionMismatch) => ()
      case other                                  => fail(s"expected DimensionMismatch, got $other")

    // Real-part orders are nonsymmetric-only.
    assert(Eigen.eigSymmetricGeneralized(a3, b3, EigenSelection.Count(2, EigenOrder.LargestRealPart)).isLeft)
    // Bad k.
    assert(Eigen.eigSymmetricGeneralized(a3, b3, EigenSelection.Count(0, EigenOrder.SmallestAlgebraic)).isLeft)
    assert(Eigen.eigSymmetricGeneralized(a3, b3, EigenSelection.Count(4, EigenOrder.SmallestAlgebraic)).isLeft)
    // Left / left-and-right vectors (left = right for this problem).
    for flag <- Seq(EigenVectors.Left, EigenVectors.LeftAndRight) do
      assert(Eigen.eigSymmetricGeneralized(a3, b3, EigenSelection.All, flag).isLeft)
  }
