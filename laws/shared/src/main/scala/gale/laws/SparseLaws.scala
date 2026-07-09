package gale.laws

import gale.linalg.*
import munit.Assertions

/** Reusable laws for the sparse matrix types, expressed against the public API.
  *
  * The bundle deliberately takes a [[gale.linalg.DoubleLinearOperator]] for the
  * matvec and a dense [[gale.linalg.DMat]] as the reference, so it applies to any
  * sparse representation (CSR, CSC, COO, Banded, Diagonal, ...) without depending
  * on a shared `toDense` on the sparse trait.
  */
object SparseLaws extends Assertions:
  /** A sparse operator's matvec matches the matvec of its dense equivalent. */
  def matvecMatchesDense(sparse: DoubleLinearOperator, dense: DMat, x: DVec, tolerance: Double = 1e-12): Unit =
    VecLaws.assertClose(sparse * x, dense * x, tolerance)

  /** A sparse operator's transpose-matvec matches the dense transpose's matvec. */
  def transposeMatvecMatchesDense(
      sparse: DoubleLinearOperator,
      dense: DMat,
      x: DVec,
      tolerance: Double = 1e-12
  ): Unit =
    val out = MutableDVec.zeros(sparse.cols)
    sparse.transposeApplyTo(x, out)
    VecLaws.assertClose(out.asVec, dense.t * x, tolerance)

  /** A densified sparse matrix equals the expected dense matrix elementwise. */
  def densifiesTo(actualDense: DMat, expectedDense: DMat, tolerance: Double = 1e-12): Unit =
    MatrixLaws.assertClose(actualDense, expectedDense, tolerance)

  /** Canonicalization is idempotent at the dense level: re-canonicalizing a
    * matrix that already advertises canonical format leaves its dense image
    * unchanged. `densify` maps each stage to its dense form.
    */
  def canonicalizeIsIdempotentDense(first: DMat, second: DMat): Unit =
    MatrixLaws.assertExact(second, first)
