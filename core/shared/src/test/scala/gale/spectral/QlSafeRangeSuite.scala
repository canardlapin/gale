package gale.spectral

import gale.linalg.*

/** QL at the edges of the floating-point range. Tiny inputs (every entry
  * below `2^-600`) are solved on the exact `2^600` lift of the matrix and of
  * the tridiagonal; before the lift a 1e-13 cluster at 1e-300 lost
  * orthogonality (about 1e9·ε) and 1e-310 entries exhausted the sweep cap.
  */
class QlSafeRangeSuite extends munit.FunSuite:
  private val Eps = 2.220446049250313e-16

  private def symmetric(n: Int, entry: (Int, Int) => Double): DMat =
    Matrix.tabulate(n, n)((r, c) => if r >= c then entry(r, c) else entry(c, r))

  private def gaussian(n: Int, seed: Long, scale: Double): DMat =
    val rng = new scala.util.Random(seed)
    val raw = Array.fill(n * n)(rng.nextGaussian() * scale)
    symmetric(n, (r, c) => raw(r * n + c))

  private def withSpectrum(values: Array[Double], seed: Long): DMat =
    val n = values.length
    val q = gaussian(n, seed, 1.0).qr.q
    val scaled = Matrix.tabulate(n, n)((r, c) => q(r, c) * values(c))
    val a = scaled * q.t
    symmetric(n, (r, c) => 0.5 * (a(r, c) + a(c, r)))

  /** Relative residual and orthogonality, measured on `A / max|A|`. */
  private def quality(a: DMat, values: DVec, vectors: DMat): (Double, Double) =
    val n = a.rows
    var scale = 0.0
    for r <- 0 until n; c <- 0 until n do scale = math.max(scale, math.abs(a(r, c)))
    val unit = Matrix.tabulate(n, n)((r, c) => a(r, c) / scale)
    val av = unit * vectors
    var residual = 0.0
    var norm = 0.0
    for r <- 0 until n; c <- 0 until n do
      val x = av(r, c) - vectors(r, c) * (values(c) / scale)
      residual += x * x
      norm += unit(r, c) * unit(r, c)
    val gram = vectors.t * vectors
    var orth = 0.0
    for r <- 0 until n; c <- 0 until n do
      val x = gram(r, c) - (if r == c then 1.0 else 0.0)
      orth += x * x
    (math.sqrt(residual / norm), math.sqrt(orth))

  private def assertAccurate(name: String, a: DMat, result: DenseSpectralKernels.SymmetricEigen): Unit =
    val n = a.rows
    val (residual, orth) = quality(a, result.values, result.vectors.get)
    assert(residual <= n * Eps, s"$name: residual $residual")
    assert(orth <= 4 * n * Eps, s"$name: orthogonality $orth")

  test("QL keeps orthogonality on a 1e-13 cluster at 1e-300") {
    val a = withSpectrum(Array.tabulate(200)(i => 1e-300 * (1.0 + 1e-13 * (i % 7))), 5L)
    val ql = DenseSpectralKernels.symmetricEigen(a, wantVectors = true, divideAndConquer = false)
    assertAccurate("cluster", a, ql.toOption.getOrElse(fail(s"QL failed: $ql")))
  }

  test("large finite Hadamard inputs are scaled before reduction on every dense route") {
    // Symmetric Sylvester H obeys H^2 = n I. At n=64 its spectrum is +/-8e307,
    // all representable, but summing 63 row magnitudes of 1e307 overflowed.
    // n=32 also covers the ordinary QL vector route below the D&C threshold.
    for n <- Seq(32, 64) do
      val unit = 1e307
      def entry(r: Int, c: Int): Double =
        if (java.lang.Integer.bitCount(r & c) & 1) == 0 then unit else -unit
      val a = Matrix.tabulate(n, n)(entry)
      val workspace = DenseWorkspace.empty
      for wantVectors <- Seq(false, true) do
        val ordinary = DenseSpectralKernels.symmetricEigen(a, wantVectors).toOption.get
        val reusable = DenseSpectralKernels.symmetricEigenWith(a, wantVectors, workspace).toOption.get
        val ql = DenseSpectralKernels.symmetricEigen(a, wantVectors, divideAndConquer = false).toOption.get
        assertEquals(reusable.values.toSeq, ordinary.values.toSeq)
        for result <- Seq(ordinary, reusable, ql) do
          for i <- 0 until n do
            val expected = if i < n / 2 then -math.sqrt(n.toDouble) else math.sqrt(n.toDouble)
            assert(result.values(i).isFinite)
            // Normwise eigenvalue bound n*epsilon*||H||_2, where ||H||_2=sqrt(n).
            assertEqualsDouble(result.values(i) / unit, expected, n * Eps * math.abs(expected))
          if wantVectors then assertAccurate(s"Hadamard n=$n", a, result)
        if wantVectors then
          for r <- 0 until n; c <- 0 until n do
            assertEquals(reusable.vectors.get(r, c), ordinary.vectors.get(r, c))
      // Neither route modifies the caller's matrix.
      for r <- 0 until n; c <- 0 until n do assertEquals(a(r, c), entry(r, c))
  }

  test("large two-by-two spectrum is restored after input scaling") {
    val a = Matrix.tabulate(2, 2)((r, c) => if r == c then 8e307 else 2e307)
    for wantVectors <- Seq(false, true) do
      val result = DenseSpectralKernels.symmetricEigen(a, wantVectors).toOption.get
      assertEqualsDouble(result.values(0) / 1e307, 6.0, 8 * Eps)
      assertEqualsDouble(result.values(1) / 1e307, 10.0, 16 * Eps)
      if wantVectors then assertAccurate("large 2x2", a, result)
  }

  test("QL converges and stays accurate on subnormal-scale (1e-310) entries") {
    val a = gaussian(90, 6L, 1e-310)
    val ql = DenseSpectralKernels.symmetricEigen(a, wantVectors = true, divideAndConquer = false)
    assertAccurate("1e-310 QL", a, ql.toOption.getOrElse(fail(s"QL failed: $ql")))
    val dc = DenseSpectralKernels.symmetricEigen(a, wantVectors = true)
    assertAccurate("1e-310 D&C", a, dc.toOption.getOrElse(fail(s"D&C failed: $dc")))
    val small = gaussian(40, 7L, 1e-310)
    val public = Eigen.eigSymmetric(small, EigenSelection.All)
    assert(public.isRight, s"eigSymmetric n=40 at 1e-310: $public")
    val values = DenseSpectralKernels.symmetricEigen(a, wantVectors = false).toOption.get.values
    var i = 0
    while i < a.rows do
      assert(math.abs(values(i) - ql.toOption.get.values(i)) <= a.rows * Eps * 1e-309, s"value $i")
      i += 1
  }

  test("the tridiagonal QL solver lifts a tiny tridiagonal") {
    // (-1, 2, -1) scaled to 2^-1040: eigenvalues 2^-1040 (2 - 2cos(kπ/(n+1))).
    val n = 60
    val scale = java.lang.Double.longBitsToDouble(1L << 12) // 2^-1062 (subnormal)
    val tinyUnit = scale * 4194304.0 // 2^-1040
    val result = DenseSpectralKernels
      .symmetricTridiagonalEigen(DVec.tabulate(n)(_ => 2.0 * tinyUnit), DVec.tabulate(n - 1)(_ => -tinyUnit), wantVectors = true)
      .toOption
      .get
    var k = 0
    while k < n do
      val exact = (2.0 - 2.0 * math.cos((k + 1) * math.Pi / (n + 1))) * tinyUnit
      // Absolute error limited by the subnormal grid of the returned values.
      assert(math.abs(result.values(k) - exact) <= 64 * java.lang.Double.MIN_VALUE, s"eigenvalue $k")
      k += 1
    val v = result.vectors.get
    val gram = v.t * v
    var orth = 0.0
    for r <- 0 until n; c <- 0 until n do
      val x = gram(r, c) - (if r == c then 1.0 else 0.0)
      orth += x * x
    assert(math.sqrt(orth) <= 4 * n * Eps, s"orthogonality ${math.sqrt(orth)}")
  }
