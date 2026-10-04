package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError

/** Keep exactly singular values strictly above the resolved threshold. Explicit cutoffs must be finite and
  * non-negative; zero keeps every positive singular value. Default preserves Gale's existing pinv convention.
  */
enum SvdCutoff:
  case Default
  case Relative(tolerance: Double)
  case Absolute(threshold: Double)

  def threshold(rows: Int, columns: Int, sigmaMax: Double): Either[LinAlgError, Double] =
    val value = this match
      case Default             => math.max(rows, columns).toDouble * SvdCutoff.MachineEpsilon
      case Relative(tolerance) => tolerance
      case Absolute(threshold) => threshold
    if rows < 0 || columns < 0 || !sigmaMax.isFinite || sigmaMax < 0.0 then
      Left(LinAlgError.InvalidArgument("invalid cutoff dimensions or largest singular value"))
    else if !value.isFinite || value < 0.0 then
      Left(LinAlgError.InvalidArgument("SVD cutoff must be finite and non-negative"))
    else
      val resolved = this match
        case Absolute(_) => value
        case _           => value * sigmaMax
      if resolved.isFinite then Right(resolved)
      else Left(LinAlgError.InvalidArgument("resolved SVD cutoff overflows"))

object SvdCutoff:
  private[gale] val MachineEpsilon: Double = 2.220446049250313e-16

/** Validation for the policy-driven dense APIs, before invoking a kernel. */
private[gale] object FiniteInput:
  def matrix(a: DMat): Either[LinAlgError, Unit] =
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        if !a(i, j).isFinite then return Left(LinAlgError.InvalidArgument(s"non-finite matrix entry ($i,$j)"))
        j += 1
      i += 1
    Right(())

  def vector(b: DVec): Either[LinAlgError, Unit] =
    var i = 0
    while i < b.length do
      if !b(i).isFinite then return Left(LinAlgError.InvalidArgument(s"non-finite vector entry $i"))
      i += 1
    Right(())
