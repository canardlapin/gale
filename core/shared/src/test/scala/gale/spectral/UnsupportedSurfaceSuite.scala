package gale.spectral

import gale.linalg.DMat
import gale.linalg.LinAlgError
import gale.linalg.LinearOperator
import gale.linalg.Matrix
import gale.linalg.Vec

/** Lock tests for the public surfaces that currently return a typed
  * `Left(UnsupportedOperation)` (or the documented `RankDeficient` stand-in
  * for rank-deficient GSVD). The table in `docs/shipped-vs-deferred.md` is the
  * narrative; this suite fails if a deferred row starts succeeding silently or
  * changes error class.
  */
class UnsupportedSurfaceSuite extends munit.FunSuite:

  private def eye(n: Int): DMat =
    Matrix.eye(n)

  private def assertUnsupported(result: Either[LinAlgError, ?], clue: String): Unit =
    result match
      case Left(_: LinAlgError.UnsupportedOperation) => ()
      case other                                     => fail(s"$clue: expected UnsupportedOperation, got $other")

  test("eigGeneralizedNonsymmetric without a QZ backend is UnsupportedOperation") {
    val a = Matrix.dense(2, 2)(1.0, 2.0, 0.0, 3.0)
    val b = eye(2)
    assertUnsupported(Eigen.eigGeneralizedNonsymmetric(a, b), "QZ default backend")
  }

  test("iterative shift-invert and Around remain UnsupportedOperation") {
    val n = 8
    val a = eye(n)
    val count = EigenSelection.Count(2, EigenOrder.LargestAlgebraic)
    val mag = EigenSelection.Count(2, EigenOrder.LargestMagnitude)
    val shift = Some(SpectralTarget.ShiftInvert(0.5, LinearSolvePlan.Backend))
    val around = Some(SpectralTarget.Around(0.25))
    assertUnsupported(Eigen.eigSymmetric(a, n, count, SpectralOptions(), shift), "symmetric ShiftInvert")
    assertUnsupported(Eigen.eigSymmetric(a, n, count, SpectralOptions(), around), "symmetric Around")
    assertUnsupported(Eigen.eigNonsymmetric(a, n, mag, SpectralOptions(), shift), "nonsymmetric ShiftInvert")
    assertUnsupported(Eigen.eigNonsymmetric(a, n, mag, SpectralOptions(), around), "nonsymmetric Around")
  }

  test("Arnoldi left and LeftAndRight vectors remain UnsupportedOperation") {
    val n = 8
    val op = LinearOperator.fromFunction(n, n): (x, into) =>
      var i = 0
      while i < n do
        into(i) = (i + 1).toDouble * x(i)
        i += 1
    val sel = EigenSelection.Count(2, EigenOrder.LargestMagnitude)
    for flag <- Seq(EigenVectors.Left, EigenVectors.LeftAndRight) do
      assertUnsupported(
        Eigen.eigNonsymmetric(op, n, sel, SpectralOptions(returnVectors = flag)),
        s"Arnoldi $flag"
      )
  }

  test("underdetermined least squares remains UnsupportedOperation") {
    val wide = Matrix.dense(2, 3)(
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0
    )
    assertUnsupported(wide.leastSquares(Vec(1.0, 2.0)), "underdetermined vector RHS")
    assertUnsupported(wide.leastSquares(Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)), "underdetermined matrix RHS")
  }

  test("rank-deficient GSVD stays Left(RankDeficient) without a backend") {
    val a = Matrix.dense(3, 3)(
      1.0, 2.0, 1.0,
      3.0, 1.0, 3.0,
      0.0, 4.0, 0.0
    )
    val b = Matrix.dense(3, 3)(
      2.0, 1.0, 2.0,
      0.0, 5.0, 0.0,
      1.0, 1.0, 1.0
    )
    Svds.gsvd(a, b) match
      case Left(_: LinAlgError.RankDeficient) => ()
      case other                              => fail(s"expected RankDeficient, got $other")
  }
