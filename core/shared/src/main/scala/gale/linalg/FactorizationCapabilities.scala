package gale.linalg

/** Common read-only metadata carried by factorization capabilities.
  *
  * This deliberately contains no determinant, inverse, orthogonal-factor, or
  * solve operation. Generic code opts into only the mathematical capability it
  * needs through [[ExactSolveFactor]] or [[LeastSquaresFactor]].
  */
trait FactorizationCapability:
  def diagnostics: FactorizationDiagnostics

/** A factorization of a square system that supports exact solves.
  *
  * "Exact" distinguishes the problem contract from least squares; numerical
  * answers still have floating-point error. Both right-hand-side forms return
  * typed failures and leave inputs unchanged.
  */
trait ExactSolveFactor extends FactorizationCapability:
  def size: Int
  def solve(b: DVec): Either[LinAlgError, DVec]
  def solve(b: DMat): Either[LinAlgError, DMat]

/** A factorization of an observation-by-coefficient system that supports
  * overdetermined least-squares solves.
  */
trait LeastSquaresFactor extends FactorizationCapability:
  def observationCount: Int
  def coefficientCount: Int
  def solveLeastSquares(b: DVec): Either[LinAlgError, DVec]
  def solveLeastSquares(b: DMat): Either[LinAlgError, DMat]
