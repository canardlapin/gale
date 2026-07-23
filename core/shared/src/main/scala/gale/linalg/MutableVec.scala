package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** A vector that may be written in place.
  *
  * `MutableVec` is the explicit destination tier: pure APIs return immutable
  * [[Vec]] values, while `mulInto`-style APIs write into a `MutableVec`
  * provided by the caller. Element access is a convenience API, not a
  * performance contract; kernels operate on validated storage directly.
  */
trait MutableVec[A] extends Vec[A]:
  def update(index: Int, value: A): Unit

/** The primitive `Double` mutable vector.
  *
  * Shares the same storage discipline as [[DVec]] (platform array, offset,
  * stride). Obtain one from `MutableVec.zeros`, `MutableVec.from`, or
  * `DVec.mutableCopy`; [[toVec]] is the public conversion to an independently
  * owned immutable value. Gale internals may also borrow an aliasing read-only
  * view while they retain the mutable owner, but that view is deliberately not
  * part of the public API.
  */
final class MutableDVec private[gale] (
    private[gale] val data: DoubleArray,
    private[gale] val offsetValue: Offset,
    private[gale] val lengthValue: Length,
    private[gale] val strideValue: Stride
) extends MutableVec[Double]:
  def length: Int = lengthValue.value
  def offset: Offset = offsetValue
  def stride: Stride = strideValue

  def apply(index: Int): Double =
    checkIndex(index)
    data(offset.value + index * stride.value)

  def update(index: Int, value: Double): Unit =
    checkIndex(index)
    data(offset.value + index * stride.value) = value

  /** Internal borrowed view sharing this vector's storage. Mutations made
    * through this `MutableDVec` remain visible through the returned value.
    * Public callers must use [[toVec]], which returns an independent snapshot.
    */
  private[gale] def asVec: DVec =
    new DVec(data, offsetValue, lengthValue, strideValue)

  /** Independent immutable snapshot of the current contents. */
  def toVec: DVec =
    asVec.copy

  def :=(that: DVec): Unit =
    requireSameLength(that.length)
    DoubleKernels.dcopy(
      length,
      that.data,
      that.offset.value,
      that.stride.value,
      data,
      offset.value,
      stride.value
    )

  def +=(that: DVec): Unit =
    axpyInPlace(1.0, that)

  def -=(that: DVec): Unit =
    axpyInPlace(-1.0, that)

  def *=(alpha: Double): Unit =
    DoubleKernels.dscal(length, alpha, data, offset.value, stride.value)

  /** `this += alpha * x`. */
  def axpyInPlace(alpha: Double, x: DVec): Unit =
    requireSameLength(x.length)
    DoubleKernels.daxpy(
      length,
      alpha,
      x.data,
      x.offset.value,
      x.stride.value,
      data,
      offset.value,
      stride.value
    )

  /** Set every element to zero. */
  def clear(): Unit =
    var i = 0
    var xi = offset.value
    val step = stride.value
    while i < length do
      data(xi) = 0.0
      xi += step
      i += 1

  private def requireSameLength(thatLength: Int): Unit =
    if length != thatLength then
      throw LinAlgError.VectorLengthMismatch(length, thatLength)

  private def checkIndex(index: Int): Unit =
    if index < 0 || index >= length then
      throw LinAlgError.IndexOutOfBounds(index, length)

object MutableDVec:
  def zeros(length: Int): MutableDVec =
    require(length >= 0, "length must be non-negative")
    new MutableDVec(DoubleArray.alloc(length), Offset.unsafe(0), Length.unsafe(length), Stride.unsafe(1))

  /** Mutable copy of `x`; the original is unaffected by later writes. */
  def from(x: DVec): MutableDVec =
    val out = zeros(x.length)
    out := x
    out

object MutableVec:
  def zeros(length: Int): MutableDVec =
    MutableDVec.zeros(length)

  def from(x: DVec): MutableDVec =
    MutableDVec.from(x)
