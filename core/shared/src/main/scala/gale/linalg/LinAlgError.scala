package gale.linalg

sealed abstract class LinAlgError(message: String) extends RuntimeException(message)

object LinAlgError:
  final case class DimensionMismatch(expected: Shape, actual: Shape)
      extends LinAlgError(s"dimension mismatch: expected ${expected.rows.value}x${expected.cols.value}, got ${actual.rows.value}x${actual.cols.value}")

  final case class VectorLengthMismatch(expected: Int, actual: Int)
      extends LinAlgError(s"vector length mismatch: expected $expected, got $actual")

  final case class IndexOutOfBounds(index: Int, bound: Int)
      extends LinAlgError(s"index $index is outside bound $bound")

  final case class NonSquareMatrix(shape: Shape)
      extends LinAlgError(s"matrix must be square, got ${shape.rows.value}x${shape.cols.value}")

  final case class SingularMatrix(index: Int)
      extends LinAlgError(s"matrix is singular at pivot $index")

  final case class NotPositiveDefinite(index: Int)
      extends LinAlgError(s"matrix is not positive definite at leading minor $index")

  final case class RankDeficient(rank: Int, cols: Int)
      extends LinAlgError(s"matrix is rank deficient: rank $rank for $cols columns")

  final case class UnsupportedOperation(operation: String)
      extends LinAlgError(s"unsupported linear algebra operation: $operation")

  final case class InvalidArgument(message: String)
      extends LinAlgError(message)

  final case class UnsupportedRepresentation(message: String)
      extends LinAlgError(message)

  final case class DidNotConverge(iterations: Int, residual: Double)
      extends LinAlgError(s"solver did not converge after $iterations iterations; residual=$residual")

  /** Every requested spectral pair passed its residual test, but the solver did
    * not certify that those pairs belong to the requested global spectral
    * extreme.
    */
  final case class SpectralExtremeNotCertified(iterations: Int, residual: Double)
      extends LinAlgError(
        s"spectral residuals converged after $iterations iterations, but membership in the requested extreme is not certified; residual=$residual"
      )
