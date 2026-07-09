package gale.linalg

def norm(x: DVec): Double =
  x.norm2

extension [A](either: Either[LinAlgError, A])
  /** Unwrap a total result, throwing the `LinAlgError` on failure. */
  def orThrow: A =
    either match
      case Right(value) => value
      case Left(error)  => throw error
