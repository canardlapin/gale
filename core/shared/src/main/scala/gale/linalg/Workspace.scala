package gale.linalg

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Reusable scratch for the one dense factorization that genuinely needs it.
  *
  * The QR implementation partitions this single buffer differently by path: a
  * small unblocked factorization needs one row-update vector, while blocked QR
  * reuses packed `Vᵀ`, compact-WY `T`, and trailing-update `W` regions. LU and
  * Cholesky own their result buffers outright and require no external scratch.
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
  private[gale] inline val QrBlockSize = 32
  private[gale] inline val QrBlockedMin = 96

  private[gale] def usesBlockedQR(rows: Int, cols: Int): Boolean =
    math.min(rows, cols) >= QrBlockedMin

  private[gale] def qrWorkSize(rows: Int, cols: Int): Int =
    if usesBlockedQR(rows, cols) then
      // Packed V^T (block*rows), W (block*cols), and compact-WY T (block^2).
      val size = QrBlockSize.toLong * rows + QrBlockSize.toLong * cols + QrBlockSize.toLong * QrBlockSize
      require(size <= Int.MaxValue.toLong, s"QR workspace exceeds ${Int.MaxValue} elements")
      size.toInt
    else
      math.max(rows, cols)

  def empty: DenseWorkspace =
    new DenseWorkspace(DoubleArray.alloc(0))

  /** Pre-size the reusable scratch for the QR path selected by this shape. */
  def forQR(rows: Int, cols: Int): DenseWorkspace =
    require(rows >= 0 && cols >= 0, "rows and cols must be non-negative")
    val workspace = empty
    workspace.work(qrWorkSize(rows, cols))
    workspace
