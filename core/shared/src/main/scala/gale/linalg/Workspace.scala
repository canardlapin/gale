package gale.linalg

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*

/** Checked primitive scratch counts for one algorithm path.
  *
  * Counts describe simultaneously live storage, not result ownership. Use
  * [[simultaneous]] when two requirements must coexist (component-wise sum) and
  * [[alternative]] when only one branch executes (component-wise maximum).
  */
final class ScratchRequirement private (
    val doubleElements: Int,
    val indexElements: Int
):
  def isEmpty: Boolean = doubleElements == 0 && indexElements == 0

  def simultaneous(that: ScratchRequirement): Either[LinAlgError, ScratchRequirement] =
    ScratchRequirement.checked(
      doubleElements.toLong + that.doubleElements.toLong,
      indexElements.toLong + that.indexElements.toLong
    )

  def alternative(that: ScratchRequirement): ScratchRequirement =
    new ScratchRequirement(
      math.max(doubleElements, that.doubleElements),
      math.max(indexElements, that.indexElements)
    )

  override def equals(other: Any): Boolean =
    other match
      case that: ScratchRequirement =>
        doubleElements == that.doubleElements && indexElements == that.indexElements
      case _ => false

  override def hashCode(): Int = 31 * doubleElements + indexElements

  override def toString: String =
    s"ScratchRequirement(doubles=$doubleElements,indices=$indexElements)"

object ScratchRequirement:
  val empty: ScratchRequirement = new ScratchRequirement(0, 0)

  /** Total constructor that rejects negative or non-addressable primitive counts. */
  def checked(
      doubleElements: Long,
      indexElements: Long
  ): Either[LinAlgError, ScratchRequirement] =
    if doubleElements < 0L || indexElements < 0L then
      Left(
        LinAlgError.InvalidArgument(
          s"scratch counts must be non-negative, got doubles=$doubleElements indices=$indexElements"
        )
      )
    else if doubleElements > Int.MaxValue.toLong || indexElements > Int.MaxValue.toLong then
      Left(
        LinAlgError.InvalidArgument(
          s"scratch counts exceed ${Int.MaxValue}: doubles=$doubleElements indices=$indexElements"
        )
      )
    else Right(new ScratchRequirement(doubleElements.toInt, indexElements.toInt))

  def doubles(elements: Long): Either[LinAlgError, ScratchRequirement] =
    checked(elements, 0L)

  def indices(elements: Long): Either[LinAlgError, ScratchRequirement] =
    checked(0L, elements)

/** Grow-only reusable primitive scratch for portable dense and sparse algorithms.
  *
  * The workspace owns one Double region and one index region. [[reserve]] grows
  * either region only when required and never exposes their backing arrays.
  * Instances are mutable single-owner execution resources and are not safe for
  * concurrent use.
  */
final class DenseWorkspace private (
    private var doubleWorkValue: DoubleArray,
    private var indexWorkValue: IndexArray
):
  def doubleCapacity: Int = doubleWorkValue.length

  /** Compatibility name for the original QR-only Double capacity. */
  def workCapacity: Int = doubleCapacity

  def indexCapacity: Int = indexWorkValue.length

  /** Ensure both primitive regions satisfy `requirement`. Existing larger
    * regions retain identity and contents.
    */
  def reserve(requirement: ScratchRequirement): Unit =
    if doubleWorkValue.length < requirement.doubleElements then
      doubleWorkValue = DoubleArray.alloc(requirement.doubleElements)
    if indexWorkValue.length < requirement.indexElements then
      indexWorkValue = IndexArray.alloc(requirement.indexElements)

  private[gale] def doubles(requirement: ScratchRequirement): DoubleArray =
    reserve(requirement)
    doubleWorkValue

  private[gale] def indices(requirement: ScratchRequirement): IndexArray =
    reserve(requirement)
    indexWorkValue

  /** Compatibility doorway for existing QR internals. */
  private[gale] def work(minLength: Int): DoubleArray =
    ScratchRequirement.checked(minLength.toLong, 0L) match
      case Left(error)        => throw error
      case Right(requirement) => doubles(requirement)

  private[gale] def workBacking: DoubleArray = doubleWorkValue

  private[gale] def indexBacking: IndexArray = indexWorkValue

object DenseWorkspace:
  private[gale] inline val QrBlockSize = 32
  private[gale] inline val QrBlockedMin = 96

  private[gale] def usesBlockedQR(rows: Int, cols: Int): Boolean =
    math.min(rows, cols) >= QrBlockedMin

  /** Report the QR scratch requirement without allocating storage. */
  def qrRequirement(rows: Int, cols: Int): Either[LinAlgError, ScratchRequirement] =
    if rows < 0 || cols < 0 then
      Left(LinAlgError.InvalidArgument(s"QR shape must be non-negative, got ${rows}x${cols}"))
    else
      val doubles =
        if usesBlockedQR(rows, cols) then
          // Packed V^T (block*rows), W (block*cols), compact-WY T (block^2).
          QrBlockSize.toLong * rows + QrBlockSize.toLong * cols + QrBlockSize.toLong * QrBlockSize
        else math.max(rows, cols).toLong
      ScratchRequirement.checked(doubles, 0L)

  private[gale] def qrWorkSize(rows: Int, cols: Int): Int =
    qrRequirement(rows, cols) match
      case Left(error)        => throw error
      case Right(requirement) => requirement.doubleElements

  def empty: DenseWorkspace =
    new DenseWorkspace(DoubleArray.alloc(0), IndexArray.alloc(0))

  def forRequirement(requirement: ScratchRequirement): DenseWorkspace =
    val workspace = empty
    workspace.reserve(requirement)
    workspace

  /** Pre-size the reusable scratch for the QR path selected by this shape. */
  def forQR(rows: Int, cols: Int): DenseWorkspace =
    qrRequirement(rows, cols) match
      case Left(error) => throw error
      case Right(requirement) =>
        forRequirement(requirement)
