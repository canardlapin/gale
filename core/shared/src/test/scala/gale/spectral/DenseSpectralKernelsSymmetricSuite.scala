package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Matrix

/** Residual/orthogonality tests for the symmetric dense spectral kernels:
  * Householder tridiagonalization and the tridiagonal QL/QR eigensolver. The
  * kernels are `private[gale]`, reachable directly from this in-package suite.
  */
class DenseSpectralKernelsSymmetricSuite extends munit.FunSuite:

  import DenseSpectralKernels.*

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

  private def tridiagonal(d: DVec, e: DVec): DMat =
    val n = d.length
    Matrix.tabulate(n, n): (i, j) =>
      if i == j then d(i)
      else if i == j + 1 then e(j)
      else if j == i + 1 then e(i)
      else 0.0

  private def frobenius(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        val v = a(i, j)
        sum += v * v
        j += 1
      i += 1
    math.sqrt(sum)

  private def diag(values: DVec): DMat =
    Matrix.tabulate(values.length, values.length)((i, j) => if i == j then values(i) else 0.0)

  private def orthogonalityError(q: DMat): Double =
    frobenius(q.t * q - Matrix.eye(q.cols))

  private def assertAscending(v: DVec): Unit =
    var i = 1
    while i < v.length do
      assert(v(i) >= v(i - 1) - 1e-12, s"eigenvalues not ascending at $i: ${v(i - 1)} then ${v(i)}")
      i += 1

  /** Eigenvector residual ‖A V − V Λ‖_F, scaled to be dimensionless. */
  private def eigenResidual(a: DMat, values: DVec, vectors: DMat): Double =
    frobenius(a * vectors - vectors * diag(values))

  // --- tridiagonalization ----------------------------------------------------

  test("tridiagonalize: Q orthogonal and Q T Qᵀ reconstructs A") {
    val n = 12
    val a = randomSymmetric(n, 424242L)
    val tri = tridiagonalize(a, wantQ = true)
    val q = tri.q.get
    val t = tridiagonal(tri.diagonal, tri.offDiagonal)

    assert(orthogonalityError(q) < 1e-10, "Q not orthogonal")
    assert(frobenius(q * t * q.t - a) < 1e-10, "Q T Qᵀ != A")
  }

  test("tridiagonalize: values-only path yields the same T without Q") {
    val n = 9
    val a = randomSymmetric(n, 7L)
    val withQ = tridiagonalize(a, wantQ = true)
    val noQ = tridiagonalize(a, wantQ = false)
    assertEquals(noQ.q, None)
    // The tridiagonal itself is the reduction's contract; both paths must agree.
    assert(frobenius(tridiagonal(withQ.diagonal, withQ.offDiagonal) - tridiagonal(noQ.diagonal, noQ.offDiagonal)) < 1e-14)
  }

  // --- standalone tridiagonal solver -----------------------------------------

  test("tridiagonal eigen: analytic Laplacian spectrum 2 - 2cos(kπ/(n+1))") {
    val n = 15
    val d = DVec.tabulate(n)(_ => 2.0)
    val e = DVec.tabulate(n - 1)(_ => -1.0)
    val result = symmetricTridiagonalEigen(d, e, wantVectors = true).toOption.get

    assertAscending(result.values)
    var k = 0
    while k < n do
      val analytic = 2.0 - 2.0 * math.cos((k + 1) * math.Pi / (n + 1))
      assert(math.abs(result.values(k) - analytic) < 1e-10, s"eigenvalue $k: ${result.values(k)} != $analytic")
      k += 1

    val t = tridiagonal(d, e)
    assert(eigenResidual(t, result.values, result.vectors.get) < 1e-10)
    assert(orthogonalityError(result.vectors.get) < 1e-10)
  }

  test("tridiagonal eigen: random (d, e) reconstructs the tridiagonal") {
    val rng = new scala.util.Random(99L)
    val n = 18
    val d = DVec.tabulate(n)(_ => rng.nextDouble() * 4.0 - 2.0)
    val e = DVec.tabulate(n - 1)(_ => rng.nextDouble() * 2.0 - 1.0)
    val result = symmetricTridiagonalEigen(d, e, wantVectors = true).toOption.get
    val t = tridiagonal(d, e)

    assertAscending(result.values)
    // T = V Λ Vᵀ.
    assert(frobenius(result.vectors.get * diag(result.values) * result.vectors.get.t - t) < 1e-10)
    assert(orthogonalityError(result.vectors.get) < 1e-10)
  }

  test("tridiagonal eigen: values-only matches the full solve") {
    val rng = new scala.util.Random(13L)
    val n = 10
    val d = DVec.tabulate(n)(_ => rng.nextDouble())
    val e = DVec.tabulate(n - 1)(_ => rng.nextDouble())
    val full = symmetricTridiagonalEigen(d, e, wantVectors = true).toOption.get
    val valuesOnly = symmetricTridiagonalEigen(d, e, wantVectors = false).toOption.get
    assertEquals(valuesOnly.vectors, None)
    var i = 0
    while i < n do
      assert(math.abs(full.values(i) - valuesOnly.values(i)) < 1e-14)
      i += 1
  }

  // --- dense symmetric eigen -------------------------------------------------

  test("symmetric eigen: random matrix satisfies A V = V Λ, ascending, orthonormal") {
    val n = 20
    val a = randomSymmetric(n, 20260709L)
    val result = symmetricEigen(a, wantVectors = true).toOption.get

    assertAscending(result.values)
    assert(orthogonalityError(result.vectors.get) < 1e-9, "V not orthonormal")
    assert(eigenResidual(a, result.values, result.vectors.get) < 1e-9, "A V != V Λ")
    // A = V Λ Vᵀ.
    assert(frobenius(result.vectors.get * diag(result.values) * result.vectors.get.t - a) < 1e-9)
  }

  test("symmetric eigen: 2x2 with hand-computed eigenvalues") {
    // [[2,1],[1,2]] has eigenvalues 1 and 3, eigenvectors (1,-1)/√2 and (1,1)/√2.
    val a = Matrix.dense(2, 2)(2.0, 1.0, 1.0, 2.0)
    val result = symmetricEigen(a, wantVectors = true).toOption.get
    assert(math.abs(result.values(0) - 1.0) < 1e-12)
    assert(math.abs(result.values(1) - 3.0) < 1e-12)
    assert(eigenResidual(a, result.values, result.vectors.get) < 1e-12)
  }

  test("symmetric eigen: 3x3 with hand-computed spectrum") {
    // Diagonal-plus-rank structure with known spectrum {0, 3, 3}? Use a clean one:
    // [[2,-1,0],[-1,2,-1],[0,-1,2]] (n=3 Laplacian): eigenvalues 2-√2, 2, 2+√2.
    val a = Matrix.dense(3, 3)(2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0)
    val result = symmetricEigen(a, wantVectors = true).toOption.get
    val expected = Seq(2.0 - math.sqrt(2.0), 2.0, 2.0 + math.sqrt(2.0))
    var i = 0
    while i < 3 do
      assert(math.abs(result.values(i) - expected(i)) < 1e-12, s"eig $i: ${result.values(i)}")
      i += 1
    assert(eigenResidual(a, result.values, result.vectors.get) < 1e-12)
  }

  test("symmetric eigen: already-diagonal matrix returns sorted diagonal") {
    val a = Matrix.dense(3, 3)(5.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 1.0)
    val result = symmetricEigen(a, wantVectors = true).toOption.get
    assert(math.abs(result.values(0) - -2.0) < 1e-14)
    assert(math.abs(result.values(1) - 1.0) < 1e-14)
    assert(math.abs(result.values(2) - 5.0) < 1e-14)
    assert(eigenResidual(a, result.values, result.vectors.get) < 1e-12)
  }

  test("symmetric eigen: zero matrix") {
    val a = Matrix.dense(4, 4)(Seq.fill(16)(0.0)*)
    val result = symmetricEigen(a, wantVectors = true).toOption.get
    var i = 0
    while i < 4 do
      assert(math.abs(result.values(i)) < 1e-14)
      i += 1
    assert(orthogonalityError(result.vectors.get) < 1e-12)
  }

  test("symmetric eigen: repeated eigenvalues keep an orthonormal basis") {
    // 2·I on a rotated basis stays 2·I; use [[3,0,0],[0,3,0],[0,0,3]] plus a
    // rank-one bump to force a repeated eigenvalue with a nontrivial eigenspace.
    val a = Matrix.dense(3, 3)(3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 7.0)
    val result = symmetricEigen(a, wantVectors = true).toOption.get
    assert(math.abs(result.values(0) - 3.0) < 1e-12)
    assert(math.abs(result.values(1) - 3.0) < 1e-12)
    assert(math.abs(result.values(2) - 7.0) < 1e-12)
    assert(orthogonalityError(result.vectors.get) < 1e-12, "degenerate eigenspace not orthonormal")
    assert(eigenResidual(a, result.values, result.vectors.get) < 1e-12)
  }

  test("symmetric eigen: n = 1") {
    val a = Matrix.dense(1, 1)(42.0)
    val result = symmetricEigen(a, wantVectors = true).toOption.get
    assertEquals(result.values.length, 1)
    assert(math.abs(result.values(0) - 42.0) < 1e-14)
    assert(math.abs(result.vectors.get(0, 0) - 1.0) < 1e-14)
  }

  test("symmetric eigen: n = 2 general") {
    val a = Matrix.dense(2, 2)(4.0, 3.0, 3.0, -4.0)
    val result = symmetricEigen(a, wantVectors = true).toOption.get
    // eigenvalues ±5.
    assert(math.abs(result.values(0) - -5.0) < 1e-12)
    assert(math.abs(result.values(1) - 5.0) < 1e-12)
    assert(eigenResidual(a, result.values, result.vectors.get) < 1e-12)
  }

  test("tridiagonal eigen: iteration guard reports DidNotConverge") {
    // A Laplacian block genuinely needs QL sweeps; a zero budget must fail typed
    // rather than loop forever.
    val n = 6
    val d = DVec.tabulate(n)(_ => 2.0)
    val e = DVec.tabulate(n - 1)(_ => -1.0)
    val result = symmetricTridiagonalEigen(d, e, wantVectors = false, maxSweepsPerValue = 0)
    assertEquals(result, Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(0)))
  }

  test("symmetric eigen: values-only agrees with full solve") {
    val n = 14
    val a = randomSymmetric(n, 555L)
    val full = symmetricEigen(a, wantVectors = true).toOption.get
    val valuesOnly = symmetricEigen(a, wantVectors = false).toOption.get
    assertEquals(valuesOnly.vectors, None)
    var i = 0
    while i < n do
      assert(math.abs(full.values(i) - valuesOnly.values(i)) < 1e-12)
      i += 1
  }
