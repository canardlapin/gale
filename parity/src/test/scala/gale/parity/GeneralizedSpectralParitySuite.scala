package gale.parity

import gale.linalg.*
import gale.parity.NumpyScipyFixtures.*
import gale.parity.NumpyScipySupport.*
import gale.spectral.*

/** Generalized eigen / GSVD / QZ against NumPy and SciPy.
  *
  * Breeze 2.1 has no public `eigh(A, B)`, `gsvd`, or `qz`. SciPy does:
  * `scipy.linalg.eigh(A, B)` (type 1) is the reference for
  * [[Eigen.eigSymmetricGeneralized]]; generalized singular values of a
  * full-column-rank pencil are the square roots of `eigh(AᵀA, BᵀB)` (SciPy
  * ships no high-level `gsvd`; that Gram-pencil identity is the portable
  * NumPy/SciPy check, matching R `geigen` on the same pencil). QZ /
  * `eig(A, B)` is a backend seam: with no capable backend Gale must return
  * `Left(UnsupportedOperation)` rather than invent a spectrum.
  */
class GeneralizedSpectralParitySuite extends munit.FunSuite:

  private val eigTol = 1e-8
  private val resTol = 1e-7
  private val gsvdTol = 1e-7

  test("eigSymmetricGeneralized eigenvalues match SciPy eigh(A, B) type 1") {
    for ref <- generalizedEigen do
      val a = galeMatrix(ref.a)
      val b = galeMatrix(ref.b)
      val d = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).orThrow
      val galeVals = Array.tabulate(d.size)(d.eigenvalues(_))
      assertArrayClose(galeVals, ref.eigenvalues, eigTol, s"${ref.name} eigenvalues")
  }

  test("eigSymmetricGeneralized vectors are B-orthonormal with small generalized residuals") {
    for ref <- generalizedEigen do
      val a = galeMatrix(ref.a)
      val b = galeMatrix(ref.b)
      val d = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All, EigenVectors.Right).orThrow
      val gram = d.eigenvectors.t * b * d.eigenvectors
      val eye = Matrix.eye(d.size)
      val gramErr = frobenius(gram - eye)
      assert(gramErr <= resTol, s"${ref.name}: ||Xᵀ B X − I||_F=$gramErr")
      var i = 0
      while i < d.size do
        val x = d.eigenvectors.col(i)
        val lambda = d.eigenvalues(i)
        val r = (a * x) - ((b * x) * lambda)
        assert(r.norm2 <= resTol, s"${ref.name}: residual[$i]=${r.norm2} λ=$lambda")
        i += 1
  }

  test("gsvd ratios match SciPy Gram-pencil eigh(AᵀA, BᵀB) on full-rank pencils") {
    for ref <- gsvd do
      val g = Svds.gsvd(galeMatrix(ref.a), galeMatrix(ref.b)).orThrow
      val galeRatios = Array.tabulate(g.size)(g.ratio)
      assertArrayClose(galeRatios, ref.ratios, gsvdTol, s"${ref.name} gsvd ratios")
      ref.ordinarySvd.foreach: svd =>
        assertArrayClose(galeRatios, svd, gsvdTol, s"${ref.name} gsvd vs numpy.linalg.svd (B = I)")
  }

  test("gsvd reconstructs A and B on well-determined columns") {
    for ref <- gsvd if ref.name != "analytic_infinite_finite_zero" do
      val a = galeMatrix(ref.a)
      val b = galeMatrix(ref.b)
      val g = Svds.gsvd(a, b).orThrow
      val cMat = DMat.tabulate(g.size, g.size)((i, j) => if i == j then g.c(i) else 0.0)
      val sMat = DMat.tabulate(g.size, g.size)((i, j) => if i == j then g.s(i) else 0.0)
      val aErr = frobenius(a - g.u * cMat * g.x.t) / math.max(1.0, frobenius(a))
      val bErr = frobenius(b - g.v * sMat * g.x.t) / math.max(1.0, frobenius(b))
      assert(aErr <= gsvdTol, s"${ref.name}: A reconstruction $aErr")
      assert(bErr <= gsvdTol, s"${ref.name}: B reconstruction $bErr")
  }

  test("gsvd rank-deficient stacked pencil is Left(RankDeficient) without a backend") {
    val a = Matrix.tabulate(2, 4)((i, j) => if i == 0 then (j + 1).toDouble else 0.0)
    val b = Matrix.tabulate(1, 4)((_, j) => (j + 1).toDouble)
    Svds.gsvd(a, b) match
      case Left(_: LinAlgError.RankDeficient) => ()
      case other => fail(s"expected RankDeficient, got $other")
  }

  test("QZ / eigGeneralizedNonsymmetric is UnsupportedOperation without a backend") {
    val a = galeMatrix(generalizedEigen.head.a)
    val b = galeMatrix(generalizedEigen.head.b)
    Eigen.eigGeneralizedNonsymmetric(a, b) match
      case Left(LinAlgError.UnsupportedOperation(op)) =>
        assert(
          op.contains("QZ") || op.contains("generalized nonsymmetric"),
          s"unexpected unsupported message: $op"
        )
      case other =>
        fail(s"expected UnsupportedOperation for QZ, got $other")
  }
