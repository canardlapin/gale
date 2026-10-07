package gale.optim

import gale.linalg.DMat

/** Owned elementwise bounds. Infinite open sides and finite fixed coordinates are supported. */
final class BoxBounds private (val lower: DMat, val upper: DMat):
  def contains(point: DMat): Boolean =
    if point.rows != lower.rows || point.cols != lower.cols then false
    else
      var r = 0
      while r < point.rows do
        var c = 0
        while c < point.cols do
          if !point(r, c).isFinite || point(r, c) < lower(r, c) || point(r, c) > upper(r, c) then return false
          c += 1
        r += 1
      true

object BoxBounds:
  def from(lower: DMat, upper: DMat): Either[FirstOrderError, BoxBounds] =
    if lower.rows <= 0 || lower.cols <= 0 || lower.rows != upper.rows || lower.cols != upper.cols then
      Left(FirstOrderError.InvalidConfiguration("box bounds must have equal nonempty shapes"))
    else
      var r = 0
      while r < lower.rows do
        var c = 0
        while c < lower.cols do
          val l = lower(r, c)
          val u = upper(r, c)
          if l.isNaN || u.isNaN || l > u || l == Double.PositiveInfinity || u == Double.NegativeInfinity then
            return Left(FirstOrderError.InvalidConfiguration("box bounds must admit a finite point"))
          c += 1
        r += 1
      Right(
        new BoxBounds(
          DMat.tabulate(lower.rows, lower.cols)((r, c) => lower(r, c)),
          DMat.tabulate(upper.rows, upper.cols)((r, c) => upper(r, c))
        )
      )

  def uniform(rows: Int, lower: Double, upper: Double, columns: Int = 1): Either[FirstOrderError, BoxBounds] =
    if rows <= 0 || columns <= 0 then Left(FirstOrderError.InvalidConfiguration("box dimensions must be positive"))
    else from(DMat.tabulate(rows, columns)((_, _) => lower), DMat.tabulate(rows, columns)((_, _) => upper))
