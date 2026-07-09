package gale.linalg

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Reusable scratch for the one dense factorization that genuinely needs it.
  *
  * After P3 the sole owner is [[DenseDecompositions.qr]], whose Householder
  * reflector vector is length `rows`. LU and Cholesky own their result buffers
  * outright and require no external scratch, so they have no workspace-aware
  * entry point; that is why this holds only a single `work` buffer.
  */
final class DenseWorkspace private (
    private var workValue: DoubleArray
):
  def workCapacity: Int =
    workValue.length

  private[gale] def work(minLength: Int): DoubleArray =
    if workValue.length < minLength then
      workValue = DoubleArray.alloc(minLength)
    workValue

  /** The backing scratch array, exposed so tests can assert genuine reuse
    * (identity) across successive factorizations.
    */
  private[gale] def workBacking: DoubleArray =
    workValue

object DenseWorkspace:
  def empty: DenseWorkspace =
    new DenseWorkspace(DoubleArray.alloc(0))

  /** Pre-size the scratch for a QR of the given shape. The reflector vector is
    * length `rows`; `cols` names the problem but does not affect the buffer.
    */
  def forQR(rows: Int, cols: Int): DenseWorkspace =
    require(rows >= 0 && cols >= 0, "rows and cols must be non-negative")
    val workspace = empty
    workspace.work(rows)
    workspace
