package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Rows
import gale.linalg.Shape

/** An orthonormal basis stored as dimension x rank columns. Projection of q vectors costs O(dimension * rank * q),
  * without a dimension-square projector. Singular values and cutoff describe the source matrix, not the unit basis.
  */
final class Subspace private[spectral] (
    val basis: DMat,
    val singularValues: DVec,
    val sigmaMax: Double,
    val cutoff: Double
):
  def dimension: Int = basis.rows
  def rank: Int = basis.cols

  def project(b: DVec): Either[LinAlgError, DVec] =
    if b.length != dimension then Left(LinAlgError.VectorLengthMismatch(dimension, b.length))
    else
      FiniteInput
        .vector(b)
        .map: _ =>
          val coordinates = DVec.tabulate(rank): k =>
            var sum = 0.0
            var i = 0
            while i < dimension do
              sum += basis(i, k) * b(i)
              i += 1
            sum
          DVec.tabulate(dimension): i =>
            var sum = 0.0
            var k = 0
            while k < rank do
              sum += basis(i, k) * coordinates(k)
              k += 1
            sum

  /** Each column of b is a vector in this subspace's ambient space. */
  def project(b: DMat): Either[LinAlgError, DMat] =
    if b.rows != dimension then Left(LinAlgError.DimensionMismatch(Shape(Rows(dimension), Cols(b.cols)), b.shape))
    else
      FiniteInput
        .matrix(b)
        .map: _ =>
          val coordinates = DMat.tabulate(rank, b.cols): (k, j) =>
            var sum = 0.0
            var i = 0
            while i < dimension do
              sum += basis(i, k) * b(i, j)
              i += 1
            sum
          DMat.tabulate(dimension, b.cols): (i, j) =>
            var sum = 0.0
            var k = 0
            while k < rank do
              sum += basis(i, k) * coordinates(k, j)
              k += 1
            sum

  def residual(b: DVec): Either[LinAlgError, DVec] = project(b).map(b - _)
  def residual(b: DMat): Either[LinAlgError, DMat] = project(b).map(b - _)

  /** Euclidean distance of one vector; no implicit membership decision. */
  def distance(b: DVec): Either[LinAlgError, Double] = residual(b).map(_.norm2)

object Subspace:
  def rowSpace(a: DMat, cutoff: SvdCutoff = SvdCutoff.Default)(using SpectralBackend): Either[LinAlgError, Subspace] =
    TruncatedSvd.factor(a, cutoff).map(_.rowSpace)

  def columnSpace(a: DMat, cutoff: SvdCutoff = SvdCutoff.Default)(using
      SpectralBackend
  ): Either[LinAlgError, Subspace] =
    TruncatedSvd.factor(a, cutoff).map(_.columnSpace)
