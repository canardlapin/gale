package gale.spectral

import gale.linalg.DMat
import gale.linalg.DenseWorkspace
import gale.linalg.Matrix
import gale.linalg.Vec

/** Regression for the tridiagonal QL deflation test (mote
  * bd-01M4EKP8CWDHZ2W1WFBRQ64Z9X). A purely local criterion
  * `|e(m)| ≤ ε(|d(m)| + |d(m+1)|)` can never fire inside a cluster of
  * eigenvalues at 0 (the neighbouring diagonals are themselves `O(ε‖A‖)`), so
  * matrices `Q diag(λ) Qᵀ` with highly repeated zero eigenvalues exhausted the
  * 30-sweep cap for `n ≳ 96`. The solver now also deflates at `ε ‖T‖`.
  *
  * Checks are reference-free (no Breeze; runs on JVM and JS): every eigenvalue is
  * within `c n ε ‖A‖` of the prescribed spectrum (Weyl), `‖A V − V Λ‖_max ≤ c n ε ‖A‖`,
  * and `‖VᵀV − I‖_max ≤ c n ε`, with `c = 32`.
  */
class SymmetricEigenRepeatedZeroSuite extends munit.FunSuite:

  private val eps = Math.ulp(1.0)
  private val c = 32.0

  /** `Q diag(spectrum) Qᵀ`, exactly symmetric, with `Q` from a seeded QR. */
  private def withSpectrum(spectrum: Array[Double], seed: Long): DMat =
    val n = spectrum.length
    val rng = new scala.util.Random(seed)
    val q = Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0).qr.q
    val a = Array.ofDim[Double](n, n)
    for i <- 0 until n; j <- 0 to i do
      var s = 0.0
      var k = 0
      while k < n do
        s += q(i, k) * spectrum(k) * q(j, k)
        k += 1
      a(i)(j) = s
      a(j)(i) = s
    Matrix.tabulate(n, n)((i, j) => a(i)(j))

  private def twoLevel(n: Int): Array[Double] = Array.tabulate(n)(i => if i < n / 2 then 0.0 else 1.0)
  private def threeLevel(n: Int): Array[Double] = Array.tabulate(n)(i => (i % 3).toDouble)
  private def threeLevelWithTail(n: Int): Array[Double] =
    Array.tabulate(n)(i => if i < 3 * n / 4 then (i % 3).toDouble else 3.0 + i.toDouble / n)

  private val cases: List[(String, Array[Double], Long)] =
    List(
      ("two-level n=96", twoLevel(96), 5L),
      ("three-level n=96", threeLevel(96), 5L),
      ("three-level+tail n=128", threeLevelWithTail(128), 128L),
      ("three-level+tail n=128", threeLevelWithTail(128), 2L),
      ("three-level+tail n=160", threeLevelWithTail(160), 160L),
      ("three-level+tail n=160", threeLevelWithTail(160), 3L)
    )

  private def maxAbs(m: DMat): Double =
    var worst = 0.0
    for i <- 0 until m.rows; j <- 0 until m.cols do worst = math.max(worst, math.abs(m(i, j)))
    worst

  test("repeated zero eigenvalues converge with backward-stable residuals and orthonormal vectors") {
    for (label, spectrum, seed) <- cases do
      val n = spectrum.length
      val clue = s"$label seed=$seed"
      val a = withSpectrum(spectrum, seed)
      val norm = spectrum.map(math.abs).max
      val tol = c * n * eps * norm
      val d = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.Right) match
        case Right(value) => value
        case Left(error)  => fail(s"$clue: $error")
      val expected = spectrum.sorted
      for i <- 0 until n do
        assert(math.abs(d.eigenvalues(i) - expected(i)) <= tol, s"$clue λ[$i]=${d.eigenvalues(i)} expected ${expected(i)}")
      val v = d.eigenvectors
      val lambda = Matrix.tabulate(n, n)((i, j) => if i == j then d.eigenvalues(i) else 0.0)
      val residual = maxAbs(a * v - v * lambda)
      assert(residual <= tol, s"$clue residual=$residual tol=$tol")
      val orthogonality = maxAbs(v.t * v - Matrix.eye(n))
      assert(orthogonality <= c * n * eps, s"$clue orthogonality=$orthogonality")
  }

  test("values-only and workspace routes converge on the same inputs and agree exactly") {
    for (label, spectrum, seed) <- cases do
      val clue = s"$label seed=$seed"
      val a = withSpectrum(spectrum, seed)
      val ordinary = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly)
      val workspace = Eigen.eigSymmetricWith(a, EigenSelection.All, EigenVectors.ValuesOnly, DenseWorkspace.empty)
      assert(ordinary.isRight, s"$clue ordinary: $ordinary")
      assertEquals(
        workspace.map(_.eigenvalues.toSeq),
        ordinary.map(_.eigenvalues.toSeq),
        s"$clue: workspace and ordinary routes share one kernel"
      )
  }

  // ---------------------------------------------------------------------------
  // Non-finite and overflowing scales: the norm-scaled test must not misfire
  // ---------------------------------------------------------------------------

  private val bLow = 2.5 - math.sqrt(1.25) // eigenvalues of [[2, 1], [1, 3]]
  private val bHigh = 2.5 + math.sqrt(1.25)

  test("an infinite entry disables the norm-scaled test: the finite block is still iterated") {
    // A = [[B, 0], [0, +Inf]]. With a scale of ‖T‖ = Inf every finite e(m)
    // would deflate at once and return diag(B) unrotated.
    val a = Matrix.dense(3, 3)(
      2.0, 1.0, 0.0,
      1.0, 3.0, 0.0,
      0.0, 0.0, Double.PositiveInfinity
    )
    for vectors <- List(EigenVectors.ValuesOnly, EigenVectors.Right) do
      val d = Eigen.eigSymmetric(a, EigenSelection.All, vectors) match
        case Right(value) => value
        case Left(error)  => fail(s"$vectors: $error")
      assert(math.abs(d.eigenvalues(0) - bLow) <= 4 * eps * bHigh, s"$vectors λ0=${d.eigenvalues(0)}")
      assert(math.abs(d.eigenvalues(1) - bHigh) <= 4 * eps * bHigh, s"$vectors λ1=${d.eigenvalues(1)}")
      assert(d.eigenvalues(2).isPosInfinity, s"$vectors λ2=${d.eigenvalues(2)}")
      if vectors == EigenVectors.Right then
        // The B block was rotated, not returned as the identity basis.
        assert(math.abs(d.eigenvectors(0, 0)) < 0.99, s"unrotated eigenvector ${d.eigenvectors(0, 0)}")
    val t = DenseSpectralKernels.symmetricTridiagonalEigen(
      Vec(2.0, 3.0, Double.PositiveInfinity),
      Vec(1.0, 0.0),
      wantVectors = false
    )
    assertEquals(t.map(_.values.length), Right(3))
    val values = t.toOption.get.values
    assert(math.abs(values(0) - bLow) <= 4 * eps * bHigh && math.abs(values(1) - bHigh) <= 4 * eps * bHigh, s"tridiagonal $values")
  }

  test("finite tridiagonals whose row sums overflow: finite norm scale and safe-range rescaling") {
    // ‖T‖∞ overflows to +Inf for both although every entry is finite; a row-sum
    // scale would deflate everything and return the unrotated diagonal. The
    // unscaled QL recurrences also overflowed on these (DidNotConverge, or a
    // wrong eigenvalue a(1 - √3)/2 for the second) before the 2^-600 rescaling.
    val e = 9e307
    val a = 6e307
    val cases = List(
      ("zero diagonal", Vec(0.0, 0.0, 0.0), Vec(e, e), Seq(-math.sqrt(2.0) * e, 0.0, math.sqrt(2.0) * e), e),
      ("constant", Vec(a, a, a), Vec(a, a), Seq(a * (1.0 - math.sqrt(2.0)), a, a * (1.0 + math.sqrt(2.0))), a)
    )
    for (label, diagonal, off, expected, scale) <- cases do
      val t = DenseSpectralKernels.symmetricTridiagonalEigen(diagonal, off, wantVectors = true) match
        case Right(value) => value
        case Left(error)  => fail(s"$label: $error")
      for i <- 0 until 3 do
        assert(t.values(i).isFinite, s"$label λ[$i]=${t.values(i)}")
        assert(math.abs(t.values(i) - expected(i)) <= 32 * eps * scale, s"$label λ[$i]=${t.values(i)} expected ${expected(i)}")
  }
