package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Matrix
import gale.linalg.orThrow

/** Residual tests for the nonsymmetric dense spectral kernels: Householder
  * Hessenberg reduction and Francis double-shift QR. Eigenpairs are verified in
  * complex arithmetic against the decoded real-Schur packed vectors, and the
  * kernel output is fed through the `NonsymmetricEigenDecomposition` constructor
  * so its exact-conjugate-symmetry invariant is exercised on real output.
  */
class DenseSpectralKernelsNonsymmetricSuite extends munit.FunSuite:

  import DenseSpectralKernels.*

  // --- helpers ---------------------------------------------------------------

  private def randomReal(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

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

  private def orthogonalityError(q: DMat): Double =
    frobenius(q.t * q - Matrix.eye(q.cols))

  /** Wrap kernel output in the public result so the packing invariant is checked
    * and the typed accessors decode the vectors.
    */
  private def wrap(result: NonsymmetricEigen): NonsymmetricEigenDecomposition =
    val n = result.re.length
    NonsymmetricEigenDecomposition(
      result.re,
      result.im,
      result.vectors.getOrElse(DMat.zeros(n, 0)),
      None,
      SpectralDiagnostics(n, n, DVec.zeros(n), 0.0, 0)
    )

  /** ‖A v − λ v‖ for the i-th eigenpair in complex arithmetic, using the decoded
    * (realPart, imagPart) vectors. `A` is real, so `A v = (A vr, A vi)`.
    */
  private def eigenpairResidual(a: DMat, decomp: NonsymmetricEigenDecomposition, i: Int): Double =
    val lambda = decomp.eigenvalue(i)
    val (vr, vi) = decomp.eigenvector(i)
    val avr = a * vr
    val avi = a * vi
    // λ v = (lr·vr − li·vi) + i (lr·vi + li·vr)
    val realPart = avr - (vr * lambda.re - vi * lambda.im)
    val imagPart = avi - (vi * lambda.re + vr * lambda.im)
    math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))

  private def assertAllEigenpairs(a: DMat, result: NonsymmetricEigen, tol: Double): Unit =
    val decomp = wrap(result)
    val scale = math.max(frobenius(a), 1.0)
    var i = 0
    while i < decomp.size do
      val res = eigenpairResidual(a, decomp, i)
      assert(res < tol * scale, s"eigenpair $i residual $res exceeds ${tol * scale}")
      i += 1

  private def trace(a: DMat): Double =
    var s = 0.0
    var i = 0
    while i < a.rows do
      s += a(i, i)
      i += 1
    s

  // --- Hessenberg reduction --------------------------------------------------

  test("hessenberg: Q H Qᵀ = A, H upper Hessenberg, Q orthogonal") {
    val n = 12
    val a = randomReal(n, 314159L)
    val red = hessenberg(a, wantQ = true)
    val q = red.q.get
    val h = red.h

    // H has zeros strictly below the subdiagonal.
    var i = 2
    while i < n do
      var j = 0
      while j <= i - 2 do
        assertEquals(h(i, j), 0.0, s"H($i,$j) should be zero")
        j += 1
      i += 1

    assert(orthogonalityError(q) < 1e-10, "Q not orthogonal")
    assert(frobenius(q * h * q.t - a) < 1e-9, "Q H Qᵀ != A")
  }

  test("hessenberg: values-only path leaves Q absent but H identical") {
    val n = 8
    val a = randomReal(n, 2L)
    val withQ = hessenberg(a, wantQ = true)
    val noQ = hessenberg(a, wantQ = false)
    assertEquals(noQ.q, None)
    assert(frobenius(withQ.h - noQ.h) < 1e-14)
  }

  // --- Francis QR: fixed spectra ---------------------------------------------

  test("nonsymmetric eigen: diagonal matrix returns its diagonal, real vectors") {
    val a = Matrix.dense(3, 3)(2.0, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 7.0)
    val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
    val decomp = wrap(result)
    // Eigenvalues are {2,-5,7} in some order, all real.
    val values = (0 until 3).map(i => result.re(i)).sorted
    assertEquals(values, Seq(-5.0, 2.0, 7.0))
    var i = 0
    while i < 3 do
      assert(decomp.isRealPair(i))
      i += 1
    assertAllEigenpairs(a, result, 1e-12)
  }

  test("nonsymmetric eigen: scaled rotation block has conjugate-pair spectrum") {
    // [[1,-√3],[√3,1]] has eigenvalues 1 ± i√3 (magnitude 2).
    val s3 = math.sqrt(3.0)
    val a = Matrix.dense(2, 2)(1.0, -s3, s3, 1.0)
    val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
    // Positive-imaginary member first, exact conjugate symmetry.
    assert(result.im(0) > 0.0)
    assertEquals(result.re(0), result.re(1))
    assertEquals(result.im(1), -result.im(0))
    assert(math.abs(result.re(0) - 1.0) < 1e-12)
    assert(math.abs(result.im(0) - s3) < 1e-12)
    assertAllEigenpairs(a, result, 1e-12)
  }

  test("nonsymmetric eigen: companion matrix of (x-1)(x-2)(x-3)") {
    // Companion of x³ - 6x² + 11x - 6; eigenvalues 1, 2, 3.
    val a = Matrix.dense(3, 3)(
      0.0, 0.0, 6.0,
      1.0, 0.0, -11.0,
      0.0, 1.0, 6.0
    )
    val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
    val values = (0 until 3).map(i => result.re(i)).sorted
    var k = 0
    while k < 3 do
      assert(math.abs(values(k) - (k + 1)) < 1e-9, s"root ${values(k)}")
      assert(result.im(k) == 0.0)
      k += 1
    assertAllEigenpairs(a, result, 1e-10)
  }

  test("nonsymmetric eigen: mixed real + complex spectrum with packing") {
    // block diag: real 5, and a 2x2 with eigenvalues 1 ± 2i.
    val a = Matrix.dense(3, 3)(
      5.0, 0.0, 0.0,
      0.0, 1.0, -2.0,
      0.0, 2.0, 1.0
    )
    val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
    val decomp = wrap(result)
    // exactly one real eigenvalue (5) and one conjugate pair (1 ± 2i).
    val reals = (0 until 3).filter(i => decomp.isRealPair(i)).map(i => result.re(i))
    assertEquals(reals.length, 1)
    assert(math.abs(reals.head - 5.0) < 1e-12)
    assertAllEigenpairs(a, result, 1e-11)
  }

  // --- Francis QR: random matrices, laws -------------------------------------

  test("nonsymmetric eigen: random matrices satisfy A v = λ v and pairing laws") {
    val seeds = Seq(1L, 2L, 3L, 20260709L)
    val sizes = Seq(4, 7, 12, 18)
    for (seed, n) <- seeds.zip(sizes) do
      val a = randomReal(n, seed)
      val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
      assertAllEigenpairs(a, result, 1e-8)

      // Conjugate-pair convention: positive-imag first, exact symmetry.
      var i = 0
      while i < n do
        if result.im(i) > 0.0 then
          assert(i + 1 < n, s"positive-imag eigenvalue $i has no successor")
          assertEquals(result.im(i + 1), -result.im(i), s"pair $i imag not exact negation")
          assertEquals(result.re(i + 1), result.re(i), s"pair $i real parts differ")
        i += 1

      // trace(A) == Σ re(λ).
      var sumRe = 0.0
      i = 0
      while i < n do
        sumRe += result.re(i)
        i += 1
      assert(math.abs(sumRe - trace(a)) < 1e-8, s"trace mismatch for seed $seed")
  }

  test("nonsymmetric eigen: det(A) equals the product of eigenvalues") {
    val n = 6
    val a = randomReal(n, 77L)
    val result = nonsymmetricEigen(a, wantVectors = false).toOption.get
    // Product of complex eigenvalues (imaginary part cancels to ~0).
    var pr = 1.0
    var pi = 0.0
    var i = 0
    while i < n do
      val nr = pr * result.re(i) - pi * result.im(i)
      val ni = pr * result.im(i) + pi * result.re(i)
      pr = nr
      pi = ni
      i += 1
    val det = a.det.orThrow
    assert(math.abs(pi) < 1e-7, s"eigenvalue product not real: $pi")
    assert(math.abs(pr - det) < 1e-7 * math.max(math.abs(det), 1.0), s"det $det vs product $pr")
  }

  test("nonsymmetric eigen: values-only agrees with the vector solve") {
    val n = 10
    val a = randomReal(n, 909L)
    val full = nonsymmetricEigen(a, wantVectors = true).toOption.get
    val valuesOnly = nonsymmetricEigen(a, wantVectors = false).toOption.get
    assertEquals(valuesOnly.vectors, None)
    var i = 0
    while i < n do
      assertEquals(valuesOnly.re(i), full.re(i), s"re($i)")
      assertEquals(valuesOnly.im(i), full.im(i), s"im($i)")
      i += 1
  }

  test("nonsymmetric eigen: unit-norm eigenvectors") {
    val n = 9
    val a = randomReal(n, 4242L)
    val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
    val decomp = wrap(result)
    var i = 0
    while i < n do
      val (vr, vi) = decomp.eigenvector(i)
      val nrm = math.sqrt(vr.dot(vr) + vi.dot(vi))
      assert(math.abs(nrm - 1.0) < 1e-9, s"eigenvector $i norm $nrm")
      i += 1
  }

  test("nonsymmetric eigen: iteration guard reports DidNotConverge") {
    // A random matrix needs Francis QR sweeps; a zero budget must fail typed.
    val a = randomReal(6, 123L)
    val result = nonsymmetricEigen(a, wantVectors = false, maxSweeps = 0)
    result match
      case Left(SpectralKernelFailure.DidNotConverge(iters)) => assertEquals(iters, 0)
      case other                                             => fail(s"expected DidNotConverge, got $other")
  }

  test("nonsymmetric eigen: n = 1") {
    val a = Matrix.dense(1, 1)(-3.5)
    val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
    assertEquals(result.re(0), -3.5)
    assertEquals(result.im(0), 0.0)
    assert(math.abs(math.abs(result.vectors.get(0, 0)) - 1.0) < 1e-12)
  }

  test("nonsymmetric eigen: zero matrix has zero spectrum and identity vectors") {
    val a = Matrix.dense(3, 3)(Seq.fill(9)(0.0)*)
    val result = nonsymmetricEigen(a, wantVectors = true).toOption.get
    var i = 0
    while i < 3 do
      assertEquals(result.re(i), 0.0)
      assertEquals(result.im(i), 0.0)
      i += 1
    assert(orthogonalityError(result.vectors.get) < 1e-12)
  }
