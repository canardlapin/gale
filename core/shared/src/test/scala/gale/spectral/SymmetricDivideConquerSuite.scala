package gale.spectral

import gale.TestAccess
import gale.linalg.*

/** Dense symmetric eigenvectors at divide-and-conquer orders: scratch contract,
  * exact route agreement, and accuracy against the QL solver on structured
  * spectra (clusters, repeated zeros, grading, extreme scaling).
  */
class SymmetricDivideConquerSuite extends munit.FunSuite:
  private val Eps = 2.220446049250313e-16

  private def symmetric(n: Int, entry: (Int, Int) => Double): DMat =
    Matrix.tabulate(n, n)((r, c) => if r >= c then entry(r, c) else entry(c, r))

  private def gaussian(n: Int, seed: Long, scale: Double = 1.0): DMat =
    val rng = new scala.util.Random(seed)
    val raw = Array.fill(n * n)(rng.nextGaussian() * scale)
    symmetric(n, (r, c) => raw(r * n + c))

  private def withSpectrum(values: Array[Double], seed: Long): DMat =
    val n = values.length
    val q = gaussian(n, seed).qr.q
    val scaled = Matrix.tabulate(n, n)((r, c) => q(r, c) * values(c))
    val a = scaled * q.t
    symmetric(n, (r, c) => 0.5 * (a(r, c) + a(c, r)))

  private def entries(a: DMat): Seq[Double] =
    for r <- 0 until a.rows; c <- 0 until a.cols yield a(r, c)

  private def maxAbs(a: DMat): Double =
    var m = 0.0
    for r <- 0 until a.rows; c <- 0 until a.cols do m = math.max(m, math.abs(a(r, c)))
    m

  /** `‖A V − V Λ‖_F / ‖A‖_F` and `‖VᵀV − I‖_F`, measured on `A / max|A|`. */
  private def quality(a: DMat, values: DVec, vectors: DMat): (Double, Double) =
    val n = a.rows
    val scale = math.max(maxAbs(a), Double.MinPositiveValue)
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
    (if norm == 0.0 then math.sqrt(residual) else math.sqrt(residual / norm), math.sqrt(orth))

  test("vector-route scratch grows to the divide-and-conquer size exactly at MinOrder") {
    val below = TridiagonalDivideConquer.MinOrder - 1
    val at = TridiagonalDivideConquer.MinOrder
    assertEquals(
      Eigen.symmetricScratchRequirement(below, EigenVectors.Right).map(r => (r.doubleElements, r.indexElements)),
      Right((below, 0))
    )
    val n = at.toLong
    assertEquals(
      Eigen.symmetricScratchRequirement(at, EigenVectors.Right).map(r => (r.doubleElements.toLong, r.indexElements.toLong)),
      Right((n + 2 * n * n + 105 * n, 7 * n))
    )
    assertEquals(
      Eigen.symmetricScratchRequirement(at, EigenVectors.ValuesOnly).map(r => (r.doubleElements, r.indexElements)),
      Right((at * at + at, 0))
    )
  }

  test("vector route falls back to QL where divide-and-conquer scratch is not addressable") {
    // Largest order whose n + 2n² + 105n doubles still fit an Int (no allocation here).
    def dcDoubles(n: Long): Long = n + 2 * n * n + 105 * n
    var last = 32000L
    while dcDoubles(last + 1) <= Int.MaxValue.toLong do last += 1
    val n = last.toInt
    assert(n > 32000 && n < 46340, s"boundary $n")
    assert(DenseSpectralKernels.usesDivideAndConquer(n))
    assert(!DenseSpectralKernels.usesDivideAndConquer(n + 1))
    assertEquals(
      Eigen.symmetricScratchRequirement(n, EigenVectors.Right).map(r => (r.doubleElements.toLong, r.indexElements.toLong)),
      Right((dcDoubles(last), 7 * last))
    )
    // Above it the QL requirement (the off-diagonal alone) is reported, up to
    // the n² result limit that QL itself has.
    assertEquals(
      Eigen.symmetricScratchRequirement(n + 1, EigenVectors.Right).map(r => (r.doubleElements, r.indexElements)),
      Right((n + 1, 0))
    )
    assertEquals(
      Eigen.symmetricScratchRequirement(46340, EigenVectors.Right).map(r => (r.doubleElements, r.indexElements)),
      Right((46340, 0))
    )
  }

  test("ordinary and workspace routes agree exactly and reuse exact-sized scratch") {
    for n <- Seq(TridiagonalDivideConquer.MinOrder, 101, 150) do
      val a = gaussian(n, 9000L + n)
      val ordinary = Eigen.eigSymmetric(a, EigenSelection.All).toOption.get
      val requirement = Eigen.symmetricScratchRequirement(n, EigenVectors.Right).toOption.get
      val workspace = DenseWorkspace.empty
      val first = Eigen.eigSymmetricWith(a, EigenSelection.All, workspace).toOption.get
      assertEquals(workspace.doubleCapacity, requirement.doubleElements)
      assertEquals(workspace.indexCapacity, requirement.indexElements)
      assertEquals(first.eigenvalues.toSeq, ordinary.eigenvalues.toSeq, s"n=$n values")
      assertEquals(entries(first.eigenvectors), entries(ordinary.eigenvectors), s"n=$n vectors")
      val backing = TestAccess.workBacking(workspace)
      val second = Eigen.eigSymmetricWith(a, EigenSelection.All, workspace).toOption.get
      assert(TestAccess.sameStorage(backing, TestAccess.workBacking(workspace)))
      assertEquals(entries(second.eigenvectors), entries(first.eigenvectors), s"n=$n reuse")
  }

  private val probes: Seq[(String, () => DMat)] = Seq(
    "gaussian n=64" -> (() => gaussian(64, 1L)),
    "gaussian n=100" -> (() => gaussian(100, 2L)),
    "gaussian n=257" -> (() => gaussian(257, 3L)),
    "repeated zero n=96" -> (() => withSpectrum(Array.tabulate(96)(i => if i < 48 then 0.0 else 1.0 + i), 4L)),
    "repeated zero n=128" -> (() => withSpectrum(Array.tabulate(128)(i => if i < 64 then 0.0 else 1.0 + i), 5L)),
    "repeated zero n=160" -> (() => withSpectrum(Array.tabulate(160)(i => if i < 80 then 0.0 else 1.0 + i), 6L)),
    "cluster 1e-12 n=200" -> (() => withSpectrum(Array.tabulate(200)(i => 1.0 + 1e-12 * (i % 7)), 7L)),
    "two clusters n=180" -> (() =>
      withSpectrum(Array.tabulate(180)(i => if i % 2 == 0 then -1.0 else 1.0 + 1e-15 * i), 8L)
    ),
    "graded n=120" -> (() => withSpectrum(Array.tabulate(120)(i => math.pow(10.0, -i / 8.0)), 9L)),
    "toeplitz (-1,2,-1) n=200" -> (() =>
      symmetric(200, (r, c) => if r == c then 2.0 else if r - c == 1 then -1.0 else 0.0)
    ),
    "wilkinson n=129" -> (() =>
      symmetric(129, (r, c) => if r == c then math.abs(r - 64).toDouble else if r - c == 1 then 1.0 else 0.0)
    ),
    "identity n=90" -> (() => symmetric(90, (r, c) => if r == c then 1.0 else 0.0)),
    "zero n=80" -> (() => symmetric(80, (_, _) => 0.0)),
    "huge n=90" -> (() => gaussian(90, 10L, 1e300)),
    "tiny n=90" -> (() => gaussian(90, 11L, 1e-300)),
    "rank one n=150" -> (() =>
      val rng = new scala.util.Random(12L)
      val v = Array.fill(150)(rng.nextGaussian())
      symmetric(150, (r, c) => v(r) * v(c))
    ),
    "block diagonal n=100" -> (() =>
      symmetric(100, (r, c) => if r / 25 == c / 25 then 1.0 / (1 + r + c) else 0.0)
    ),
    "diagonal with ties n=70" -> (() => symmetric(70, (r, c) => if r == c then (r % 5).toDouble else 0.0)),
    "gaussian n=500" -> (() => gaussian(500, 13L)),
    "glued W21 x20, glue 1e-14" -> (() => gluedWilkinson(10, 20, 1e-14)),
    "glued W21 x20, glue 1e-3" -> (() => gluedWilkinson(10, 20, 1e-3)),
    "glued W25 x16, glue 1e-8" -> (() => gluedWilkinson(12, 16, 1e-8)),
    "glued W25 x16, glue 1e-3" -> (() => gluedWilkinson(12, 16, 1e-3)),
    "graded tridiagonal 1e-8..1e8 n=200" -> (() =>
      val n = 200
      def level(x: Double): Double = math.pow(10.0, -8.0 + 16.0 * x / (n - 1))
      symmetric(n, (r, c) => if r == c then level(r) else if r - c == 1 then 0.5 * level(c + 0.5) else 0.0)
    )
  )

  /** `copies` Wilkinson matrices `W(2m+1)` on the diagonal, joined by `glue`. */
  private def gluedWilkinson(m: Int, copies: Int, glue: Double): DMat =
    val size = 2 * m + 1
    symmetric(
      size * copies,
      (r, c) =>
        if r == c then math.abs(r % size - m).toDouble
        else if r - c == 1 then (if r % size == 0 then glue else 1.0)
        else 0.0
    )

  for (name, build) <- probes do
    test(s"divide and conquer matches QL accuracy: $name") {
      val a = build()
      val n = a.rows
      assert(n >= TridiagonalDivideConquer.MinOrder)
      val dc = DenseSpectralKernels.symmetricEigen(a, wantVectors = true).toOption.get
      val ql = DenseSpectralKernels.symmetricEigen(a, wantVectors = true, divideAndConquer = false).toOption.get
      val (dcResidual, dcOrth) = quality(a, dc.values, dc.vectors.get)
      val (qlResidual, qlOrth) = quality(a, ql.values, ql.vectors.get)
      // Backward-stable bounds of order n·ε. The residual must also stay within
      // 4x of QL unless it is below 2√n·ε: a deflation or orthogonality bug
      // costs orders of magnitude, while QL is sometimes unusually accurate
      // on a given input (3x on a graded tridiagonal), so a pure ratio would
      // test QL's luck rather than divide and conquer.
      val bound = n * Eps
      assert(dcResidual <= bound, s"residual $dcResidual > n·ε = $bound (QL $qlResidual)")
      assert(dcOrth <= 4 * bound, s"orthogonality $dcOrth > 4n·ε (QL $qlOrth)")
      assert(
        dcResidual <= math.max(4 * qlResidual, 2 * math.sqrt(n.toDouble) * Eps),
        s"residual $dcResidual vs QL $qlResidual"
      )
      val scale = math.max((0 until n).map(i => math.abs(ql.values(i))).max, Double.MinPositiveValue)
      var i = 0
      while i < n do
        assert(
          math.abs(dc.values(i) - ql.values(i)) <= bound * scale,
          s"eigenvalue $i: ${dc.values(i)} vs QL ${ql.values(i)}"
        )
        if i > 0 then assert(dc.values(i - 1) <= dc.values(i), s"values not ascending at $i")
        i += 1
    }

  test("non-finite input keeps the QL failure semantics at divide-and-conquer orders") {
    val n = TridiagonalDivideConquer.MinOrder + 6
    val a = symmetric(n, (r, c) => if r == 7 && c == 3 then Double.NaN else 1.0 / (1 + r + c))
    val result = Eigen.eigSymmetric(a, EigenSelection.All)
    assert(result.left.exists(_.isInstanceOf[LinAlgError.DidNotConverge]), s"got $result")
  }
