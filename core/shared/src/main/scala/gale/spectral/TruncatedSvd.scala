package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Rows
import gale.linalg.Shape
import gale.platform.DoubleArray

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
    val n = coefficientCount
    val m = observationCount
    val r = rank
    // A_k^+ = V Σ⁻¹ Uᵀ: entry (i, j) is the row-i dot of W = V Σ⁻¹ (n×r) with
    // row j of U (m×r), both contiguous, summed over k ascending.
    val w = DoubleArray.alloc(n * r)
    var k = 0
    while k < r do
      val sigma = singularValues(k)
      var i = 0
      while i < n do
        w(i * r + k) = vt(k, i) / sigma
        i += 1
      k += 1
    val uData = u.toDoubleArrayCopyRowMajor
    val out = DoubleArray.alloc(n * m)
    var i = 0
    while i < n do
      val rowW = i * r
      var j = 0
      // Four independent dots per pass (each still summed in k order).
      while j + 3 < m do
        val u0 = j * r
        val u1 = u0 + r
        val u2 = u1 + r
        val u3 = u2 + r
        var s0 = 0.0
        var s1 = 0.0
        var s2 = 0.0
        var s3 = 0.0
        k = 0
        while k < r do
          val wk = w(rowW + k)
          s0 += wk * uData(u0 + k)
          s1 += wk * uData(u1 + k)
          s2 += wk * uData(u2 + k)
          s3 += wk * uData(u3 + k)
          k += 1
        out(i * m + j) = s0
        out(i * m + j + 1) = s1
        out(i * m + j + 2) = s2
        out(i * m + j + 3) = s3
        j += 4
      while j < m do
        val rowU = j * r
        var sum = 0.0
        k = 0
        while k < r do
          sum += w(rowW + k) * uData(rowU + k)
          k += 1
        out(i * m + j) = sum
        j += 1
      i += 1
    DMat.fromDoubleArrayOwned(n, m, out)

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
      // Full rank keeps the (immutable) factors as they are; only a truncation
      // slices them.
      val full = kept == s.size
      new TruncatedSvd(
        if full then s.u else DMat.tabulate(a.rows, kept)((i, k) => s.u(i, k)),
        if full then s.vt else DMat.tabulate(kept, a.cols)((k, j) => s.vt(k, j)),
        if full then s.singularValues else DVec.tabulate(kept)(s.singularValues(_)),
        s.singularValues(0),
        threshold
      )
