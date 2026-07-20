package scalafim.compat

import gale.linalg.LinAlgError
import gale.spectral.EigenDecomposition

/** Compile-only probe for the public API shape consumed by Scalafim's Scala
  * 3.4.2 modules. Keeping both policies visible demonstrates that the strict
  * addition does not alter the historical residual-only helper.
  */
object ExtremeCertificationConsumer:
  def acceptResidualConvergence(
      result: EigenDecomposition
  ): Either[LinAlgError, EigenDecomposition] =
    result.requireConverged

  def requireCertifiedExtreme(
      result: EigenDecomposition
  ): Either[LinAlgError, EigenDecomposition] =
    result.requireExtremeCertified

  def failureKind(error: LinAlgError): String =
    error match
      case _: LinAlgError.DidNotConverge               => "not-converged"
      case _: LinAlgError.SpectralExtremeNotCertified => "extreme-not-certified"
      case _                                           => "other"
