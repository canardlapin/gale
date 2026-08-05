package gale.laws

import gale.linalg.*
import gale.spectral.EigenDecomposition
import munit.Assertions

/** Reusable invariants for partial generalized symmetric-definite operator results.
  */
object GeneralizedOperatorLaws extends Assertions:

  /** Every returned pair satisfies the true ambient residual and the diagnostic records that same norm.
    */
  def trueResiduals(
      a: DoubleLinearOperator,
      b: DoubleLinearOperator,
      result: EigenDecomposition,
      tolerance: Double
  ): Unit =
    var column = 0
    while column < result.size do
      val vector = result.eigenvectors.col(column)
      val residual =
        (a(vector) - b(vector) * result.eigenvalues(column)).norm2
      assert(
        residual <= tolerance,
        s"generalized residual $column = $residual exceeds $tolerance"
      )
      assert(
        math.abs(residual - result.diagnostics.residuals(column)) <=
          math.max(1e-12, tolerance * 1e-2),
        s"diagnostic residual $column ${result.diagnostics.residuals(column)} != $residual"
      )
      column += 1

  /** Returned vectors form a B-orthonormal block and diagnostics report the same Frobenius Gram error.
    */
  def bOrthonormal(
      b: DoubleLinearOperator,
      result: EigenDecomposition,
      tolerance: Double
  ): Unit =
    val bx = b.applyTo(result.eigenvectors).toOption.get
    val gram = result.eigenvectors.t * bx
    val error = identityError(gram)
    assert(error <= tolerance, s"B-orthogonality error $error exceeds $tolerance")
    assert(
      math.abs(error - result.diagnostics.orthogonalityError) <=
        math.max(1e-14, tolerance * 1e-4),
      s"diagnostic B-orthogonality ${result.diagnostics.orthogonalityError} != $error"
    )

  /** `actual` equals `scale * reference` under a combined absolute/relative tolerance.
    */
  def scaledSpectrum(
      actual: EigenDecomposition,
      reference: EigenDecomposition,
      scale: Double,
      absoluteTolerance: Double,
      relativeTolerance: Double
  ): Unit =
    assertEquals(actual.size, reference.size)
    var i = 0
    while i < actual.size do
      val expected = reference.eigenvalues(i) * scale
      val bound =
        absoluteTolerance + relativeTolerance * math.max(1.0, math.abs(expected))
      assert(
        math.abs(actual.eigenvalues(i) - expected) <= bound,
        s"eigenvalue $i ${actual.eigenvalues(i)} != $expected within $bound"
      )
      i += 1

  /** B-orthogonal projector `X Xᵀ B` onto the returned invariant subspace. */
  def metricProjector(
      vectors: DMat,
      b: DoubleLinearOperator
  ): DMat =
    val bx = b.applyTo(vectors).toOption.get
    vectors * bx.t

  /** Two invariant subspaces agree through their projectors, avoiding eigenvector sign and repeated-eigenspace basis
    * choices.
    */
  def sameMetricProjector(
      actualVectors: DMat,
      expectedVectors: DMat,
      b: DoubleLinearOperator,
      tolerance: Double
  ): Unit =
    val error =
      frobenius(metricProjector(actualVectors, b) - metricProjector(expectedVectors, b))
    assert(error <= tolerance, s"metric projector error $error exceeds $tolerance")

  private def identityError(matrix: DMat): Double =
    var sum = 0.0
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        val delta = matrix(row, col) - (if row == col then 1.0 else 0.0)
        sum += delta * delta
        col += 1
      row += 1
    math.sqrt(sum)

  private def frobenius(matrix: DMat): Double =
    var sum = 0.0
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        sum += matrix(row, col) * matrix(row, col)
        col += 1
      row += 1
    math.sqrt(sum)
