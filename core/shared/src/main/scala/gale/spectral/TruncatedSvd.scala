package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Rows
import gale.linalg.Shape

final case class MinimumNormSolution(solution: DVec, residual: DVec, rank: Int, cutoff: Double)
final case class MinimumNormSolutions(solution: DMat, residual: DMat, rank: Int, cutoff: Double)

/** Retained economy-SVD factors, selected once under a caller's cutoff. Solves minimize ||A_k x - b|| and then ||x||,
  * where A_k discards every singular value at or below cutoff. This need not minimize ||A x - b|| when the caller
  * deliberately discards a nonzero singular value.
  */
final class TruncatedSvd private[spectral] (
    val u: DMat,
    val vt: DMat,
    val singularValues: DVec,
    val sigmaMax: Double,
    val cutoff: Double
):
  def rank: Int = singularValues.length
  def observationCount: Int = u.rows
  def coefficientCount: Int = vt.cols

  def rowSpace: Subspace = new Subspace(vt.t, singularValues, sigmaMax, cutoff)
  def columnSpace: Subspace = new Subspace(u, singularValues, sigmaMax, cutoff)

  /** Materialize A_k^+, without changing the retained factors. */
  def pinv: DMat =
    DMat.tabulate(coefficientCount, observationCount): (i, j) =>
      var sum = 0.0
      var k = 0
      while k < rank do
        sum += (vt(k, i) / singularValues(k)) * u(j, k)
        k += 1
      sum

  def solve(b: DVec): Either[LinAlgError, MinimumNormSolution] =
    if b.length != observationCount then Left(LinAlgError.VectorLengthMismatch(observationCount, b.length))
    else
      FiniteInput
        .vector(b)
        .map: _ =>
          val z = DVec.tabulate(rank): k =>
            var sum = 0.0
            var i = 0
            while i < observationCount do
              sum += u(i, k) * b(i)
              i += 1
            sum
          val x = DVec.tabulate(coefficientCount): i =>
            var sum = 0.0
            var k = 0
            while k < rank do
              sum += vt(k, i) * (z(k) / singularValues(k))
              k += 1
            sum
          val residual = DVec.tabulate(observationCount): i =>
            var value = b(i)
            var k = 0
            while k < rank do
              value -= u(i, k) * z(k)
              k += 1
            value
          MinimumNormSolution(x, residual, rank, cutoff)

  def solve(b: DMat): Either[LinAlgError, MinimumNormSolutions] =
    if b.rows != observationCount then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(observationCount), Cols(b.cols)), b.shape))
    else
      FiniteInput
        .matrix(b)
        .map: _ =>
          val z = DMat.tabulate(rank, b.cols): (k, j) =>
            var sum = 0.0
            var i = 0
            while i < observationCount do
              sum += u(i, k) * b(i, j)
              i += 1
            sum
          val x = DMat.tabulate(coefficientCount, b.cols): (i, j) =>
            var sum = 0.0
            var k = 0
            while k < rank do
              sum += vt(k, i) * (z(k, j) / singularValues(k))
              k += 1
            sum
          val residual = DMat.tabulate(observationCount, b.cols): (i, j) =>
            var value = b(i, j)
            var k = 0
            while k < rank do
              value -= u(i, k) * z(k, j)
              k += 1
            value
          MinimumNormSolutions(x, residual, rank, cutoff)

object TruncatedSvd:
  def factor(a: DMat, cutoff: SvdCutoff = SvdCutoff.Default)(using SpectralBackend): Either[LinAlgError, TruncatedSvd] =
    for
      _ <- cutoff.threshold(a.rows, a.cols, 0.0)
      _ <- FiniteInput.matrix(a)
      s <- Svds.svd(a, SingularSelection.All)
      _ <- FiniteInput.vector(s.singularValues)
      _ <- FiniteInput.matrix(s.u)
      _ <- FiniteInput.matrix(s.vt)
      threshold <- cutoff.threshold(a.rows, a.cols, s.singularValues(0))
    yield
      val kept = (0 until s.size).count(i => s.singularValues(i) > threshold)
      new TruncatedSvd(
        DMat.tabulate(a.rows, kept)((i, k) => s.u(i, k)),
        DMat.tabulate(kept, a.cols)((k, j) => s.vt(k, j)),
        DVec.tabulate(kept)(s.singularValues(_)),
        s.singularValues(0),
        threshold
      )
