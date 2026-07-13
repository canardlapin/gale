package gale.examples

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.linalg.Vec
import gale.spectral.Eigen
import gale.spectral.EigenSelection
import gale.spectral.EigenVectors
import gale.sparse.Sparse

/** Four worked, end-to-end examples of gale's public API, each recovering a value
  * with a known, independently-derived correct answer (not merely "it ran"). This
  * suite is the compiled, always-green source backing `docs/examples.md` — every
  * code block there corresponds to a call pattern exercised here.
  */
class WorkedExamplesSuite extends munit.FunSuite:

  /** Unwrap an `Either[LinAlgError, A]`, failing the test with the error's message
    * on `Left` rather than throwing it raw.
    */
  private def unwrap[A](result: Either[LinAlgError, A]): A =
    result.fold(error => fail(error.toString), identity)

  // ===========================================================================
  // 1. Linear regression via least squares (QR)
  // ===========================================================================

  test("linear regression: QR least-squares recovers a known line's slope and intercept") {
    // y = intercept + slope*x exactly, for intercept = 2.0, slope = 3.0.
    val intercept = 2.0
    val slope = 3.0
    val xs = Seq(0.0, 1.0, 2.0, 3.0, 4.0)
    val ys = xs.map(x => intercept + slope * x)

    // Design matrix: columns [1, x] (tall: 5 rows, 2 columns).
    val X = Matrix.tabulate(xs.length, 2)((i, j) => if j == 0 then 1.0 else xs(i))
    val y = DVec.fromSeq(ys)

    val coeffs = unwrap(X.leastSquares(y))
    assertEqualsDouble(coeffs(0), intercept, 1e-9)
    assertEqualsDouble(coeffs(1), slope, 1e-9)
  }

  // ===========================================================================
  // 2. Ordinary least squares via the normal equations (LU) — a second route to
  //    the same answer as example 1.
  // ===========================================================================

  test("OLS via normal equations (LU solve) matches the QR least-squares route") {
    val intercept = 2.0
    val slope = 3.0
    val xs = Seq(0.0, 1.0, 2.0, 3.0, 4.0)
    val ys = xs.map(x => intercept + slope * x)

    val X = Matrix.tabulate(xs.length, 2)((i, j) => if j == 0 then 1.0 else xs(i))
    val y = DVec.fromSeq(ys)

    val viaQr = unwrap(X.leastSquares(y))

    // Normal equations: (XᵀX) beta = Xᵀy, solved as a square system via LU.
    val xtx = X.t * X
    val xty = X.t * y
    val viaNormalEquations = unwrap(xtx.solve(xty))

    assertEqualsDouble(viaNormalEquations(0), viaQr(0), 1e-6)
    assertEqualsDouble(viaNormalEquations(1), viaQr(1), 1e-6)
    assertEqualsDouble(viaNormalEquations(0), intercept, 1e-6)
    assertEqualsDouble(viaNormalEquations(1), slope, 1e-6)
  }

  // ===========================================================================
  // 3. PCA: the leading eigenvector of the (centered) scatter matrix recovers a
  //    known principal direction.
  // ===========================================================================

  test("PCA: the leading eigenvector of the scatter matrix aligns with the known principal direction") {
    // Two orthonormal directions: u (the dominant spread) and v (a much smaller,
    // orthogonal spread). Coefficients along each are mean-zero and mutually
    // uncorrelated (dot product 0) by construction, so the true covariance matrix
    // is EXACTLY diagonal in the (u, v) basis: its eigenvectors are exactly ±u
    // (largest eigenvalue) and ±v (smallest), independent of the sample-covariance
    // scaling convention.
    val ux = 0.6
    val uy = 0.8
    val vx = -0.8
    val vy = 0.6
    val cu = Seq(-2.0, -1.0, 0.0, 1.0, 2.0) // variance 2.0, mean 0
    val cv = Seq(0.1, -0.2, 0.2, -0.2, 0.1) // variance 0.028, mean 0; dot(cu, cv) == 0
    val ox = 5.0 // arbitrary offset, removed by centering below
    val oy = -3.0
    val n = cu.length

    val points = Matrix.tabulate(n, 2): (i, j) =>
      if j == 0 then cu(i) * ux + cv(i) * vx + ox
      else cu(i) * uy + cv(i) * vy + oy

    val colMean0 = (0 until n).map(points(_, 0)).sum / n
    val colMean1 = (0 until n).map(points(_, 1)).sum / n
    val centered = Matrix.tabulate(n, 2): (i, j) =>
      if j == 0 then points(i, 0) - colMean0 else points(i, 1) - colMean1

    // Xᵀ X is the (unnormalized) scatter matrix — a positive scalar multiple of
    // the sample covariance matrix, hence identical eigenvectors and ordering.
    val scatter = centered.t * centered
    val eigen = unwrap(Eigen.eigSymmetric(scatter, EigenSelection.All, EigenVectors.Right))

    // Eigenvalues ascending and nonnegative (a scatter matrix is PSD).
    assert(eigen.eigenvalues(0) >= -1e-9, s"eigenvalue 0 should be nonnegative: ${eigen.eigenvalues(0)}")
    assert(eigen.eigenvalues(1) >= -1e-9, s"eigenvalue 1 should be nonnegative: ${eigen.eigenvalues(1)}")
    assert(
      eigen.eigenvalues(0) <= eigen.eigenvalues(1) + 1e-9,
      s"eigenvalues must be ascending: ${eigen.eigenvalues(0)} then ${eigen.eigenvalues(1)}"
    )

    // The leading eigenvector (largest eigenvalue, last index) is parallel to u
    // (up to sign): |leading . u| ~= 1, since both are unit-norm.
    val leading = eigen.eigenvectors.col(1)
    val direction = DVec.fromSeq(Seq(ux, uy))
    assertEqualsDouble(math.abs(leading.dot(direction)), 1.0, 1e-6)
  }

  // ===========================================================================
  // 4. Graph Laplacian / Fiedler vector: two triangles joined by a bridge edge.
  // ===========================================================================

  private val laplacianEdges: Seq[(Int, Int)] =
    Seq((0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3))
  private val laplacianSize = 6

  private def degrees(edges: Seq[(Int, Int)], n: Int): Array[Int] =
    val out = Array.fill(n)(0)
    edges.foreach { case (i, j) => out(i) += 1; out(j) += 1 }
    out

  private def denseLaplacian(edges: Seq[(Int, Int)], n: Int): DMat =
    val degree = degrees(edges, n)
    val adjacency = edges.toSet ++ edges.map { case (i, j) => (j, i) }
    Matrix.tabulate(n, n): (i, j) =>
      if i == j then degree(i).toDouble
      else if adjacency.contains((i, j)) then -1.0
      else 0.0

  test("graph Laplacian: zero/positive spectral gap and a Fiedler vector that splits the two triangles") {
    // Two triangles {0,1,2} and {3,4,5}, connected only by the bridge edge (2,3) —
    // the combinatorial Laplacian L = D - A of this graph.
    val L = denseLaplacian(laplacianEdges, laplacianSize)
    val eigen = unwrap(Eigen.eigSymmetric(L, EigenSelection.All, EigenVectors.Right))

    // Smallest eigenvalue is 0 (every row of L sums to zero: L * ones == 0).
    assertEqualsDouble(eigen.eigenvalues(0), 0.0, 1e-9)
    // Its eigenvector is the constant vector: every entry has the same sign and
    // magnitude.
    val constant = eigen.eigenvectors.col(0)
    val meanAbs = (0 until laplacianSize).map(i => math.abs(constant(i))).sum / laplacianSize
    (0 until laplacianSize).foreach { i =>
      assertEqualsDouble(math.abs(constant(i)), meanAbs, 1e-9)
      assert(constant(i) * constant(0) > 0.0, s"constant eigenvector must be single-signed at $i")
    }

    // Fiedler value: strictly positive (the graph is connected).
    val fiedlerValue = eigen.eigenvalues(1)
    assert(fiedlerValue > 1e-6, s"Fiedler value should be positive, got $fiedlerValue")

    // Fiedler vector: same sign within each triangle, opposite sign across the
    // bridge — it recovers the two-triangle community split.
    val fiedler = eigen.eigenvectors.col(1)
    assert(fiedler(0) * fiedler(1) > 0.0, "nodes 0,1 (triangle A) should share sign")
    assert(fiedler(1) * fiedler(2) > 0.0, "nodes 1,2 (triangle A) should share sign")
    assert(fiedler(3) * fiedler(4) > 0.0, "nodes 3,4 (triangle B) should share sign")
    assert(fiedler(4) * fiedler(5) > 0.0, "nodes 4,5 (triangle B) should share sign")
    assert(fiedler(0) * fiedler(3) < 0.0, "the two triangles should have opposite Fiedler sign")
  }

  // ===========================================================================
  // 5. The same Laplacian built as a SPARSE matrix: L * ones ~= 0 (sparse mat-vec).
  // ===========================================================================

  test("sparse graph Laplacian: L * ones is (numerically) the zero vector") {
    val degree = degrees(laplacianEdges, laplacianSize)
    val builder = Sparse.coo(laplacianSize, laplacianSize)
    laplacianEdges.foreach { case (i, j) =>
      builder.add(i, j, -1.0)
      builder.add(j, i, -1.0)
    }
    (0 until laplacianSize).foreach(i => builder.add(i, i, degree(i).toDouble))
    val sparseL = builder.toCSR()

    val ones = Vec.fill(laplacianSize)(1.0)
    val result = sparseL * ones
    assert(result.norm2 < 1e-9, s"L * ones should vanish, norm was ${result.norm2}")
  }
