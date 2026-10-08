package gale.spectral

import gale.linalg.DenseWorkspace
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.linalg.Vec

/** Empty (`0 × 0`) input across the dense spectral entry points. The symmetric
  * routes return an empty decomposition (the `lu`/`qr`/`cholesky` precedent);
  * a selection that needs at least one eigenpair is a typed `Left`; the dense
  * SVD rejects an empty dimension by contract. Nothing may throw.
  */
class SpectralEmptyInputSuite extends munit.FunSuite:

  private val empty = Matrix.zeros(0, 0)

  private def assertEmpty(result: Either[LinAlgError, EigenDecomposition], clue: String): Unit =
    result match
      case Right(d) =>
        assertEquals(d.size, 0, clue)
        assertEquals(d.eigenvalues.length, 0, clue)
        assertEquals((d.eigenvectors.rows, d.eigenvectors.cols), (0, 0), clue)
      case Left(error) => fail(s"$clue: expected an empty decomposition, got $error")

  test("eigSymmetric on 0x0 returns an empty decomposition (vectors and values-only)") {
    assertEmpty(Eigen.eigSymmetric(empty, EigenSelection.All, EigenVectors.Right), "All/Right")
    assertEmpty(Eigen.eigSymmetric(empty, EigenSelection.All, EigenVectors.ValuesOnly), "All/ValuesOnly")
    assertEmpty(Eigen.eigSymmetric(empty, EigenSelection.ValueInterval(-1.0, 1.0), EigenVectors.Right), "ValueInterval")
  }

  test("eigSymmetric selections that need an eigenpair are typed Left on 0x0") {
    val top = Eigen.eigSymmetric(empty, EigenSelection.Count(1, EigenOrder.LargestAlgebraic), EigenVectors.ValuesOnly)
    assert(top.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]), s"Count: $top")
    val window = Eigen.eigSymmetric(empty, EigenSelection.IndexRange(0, 0), EigenVectors.Right)
    assert(window.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]), s"IndexRange: $window")
  }

  test("workspace eigSymmetricWith on 0x0 matches the ordinary route") {
    assertEmpty(Eigen.eigSymmetricWith(empty, EigenSelection.All, DenseWorkspace.empty), "workspace Right")
    assertEmpty(
      Eigen.eigSymmetricWith(empty, EigenSelection.All, EigenVectors.ValuesOnly, DenseWorkspace.empty),
      "workspace ValuesOnly"
    )
  }

  test("workspace ValueInterval on 0x0 is empty; an inverted interval is Left") {
    assertEmpty(
      Eigen.eigSymmetricWith(empty, EigenSelection.ValueInterval(-1.0, 1.0), DenseWorkspace.empty),
      "workspace ValueInterval"
    )
    val inverted = Eigen.eigSymmetricWith(empty, EigenSelection.ValueInterval(1.0, -1.0), DenseWorkspace.empty)
    assert(inverted.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]), s"workspace inverted: $inverted")
    val invertedOrdinary = Eigen.eigSymmetric(empty, EigenSelection.ValueInterval(1.0, -1.0), EigenVectors.Right)
    assert(invertedOrdinary.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]), s"inverted: $invertedOrdinary")
  }

  test("generalized symmetric-definite eigen on a 0x0 pencil is empty") {
    assertEmpty(Eigen.eigSymmetricGeneralized(empty, empty, EigenSelection.All), "generalized")
  }

  test("tridiagonal kernels accept order 0") {
    val t = DenseSpectralKernels.tridiagonalize(empty, wantQ = true)
    assertEquals((t.diagonal.length, t.offDiagonal.length), (0, 0))
    assertEquals(t.q.map(q => (q.rows, q.cols)), Some((0, 0)))
    val solved = DenseSpectralKernels.symmetricTridiagonalEigen(Vec.zeros(0), Vec.zeros(0), wantVectors = true)
    assertEquals(solved.map(_.values.length), Right(0))
  }

  test("nonsymmetric eigen on 0x0 is empty; dense SVD rejects the empty dimension") {
    val nonsym = Eigen.eigNonsymmetric(empty, EigenSelection.All, EigenVectors.Right)
    assertEquals(nonsym.map(_.size), Right(0))
    assert(empty.svd.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]), "svd 0x0")
    assert(empty.pinv.left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]), "pinv 0x0")
  }
