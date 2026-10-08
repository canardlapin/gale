package gale.spectral

import gale.linalg.DMat
import gale.linalg.DenseWorkspace
import gale.linalg.Matrix

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
