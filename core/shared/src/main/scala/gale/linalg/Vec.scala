package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

trait Vec[A]:
  def length: Int
  def apply(index: Int): A

/** An immutable-facing dense `Double` vector: a length-`n` window over platform
  * storage described by an offset and a (strictly positive) stride.
  *
  *   - '''Views vs copies.''' [[slice]] returns an aliasing view — it shares the
  *     backing storage, so it is `O(1)` and reflects later writes to that storage.
  *     [[copy]] (and the `Array`/`Seq` exporters) return independent data. `DVec`
  *     itself exposes no mutators; obtain a writable view via [[mutableCopy]] or
  *     `MutableVec`, and remember that a [[MutableDVec]] built as a view (e.g. a
  *     matrix row) writes through to the shared storage.
  *   - '''Positive stride only.''' `gale` supports positive strides; there is no
  *     negative-stride (reversed) view. Element `i` lives at `offset + i*stride`.
  */
final class DVec private[gale] (
    private[gale] val data: DoubleArray,
    private[gale] val offsetValue: Offset,
    private[gale] val lengthValue: Length,
    private[gale] val strideValue: Stride
) extends Vec[Double]:
  def length: Int = lengthValue.value
  def offset: Offset = offsetValue
  def stride: Stride = strideValue

  def apply(index: Int): Double =
    checkIndex(index)
    data(offset.value + index * stride.value)

  def updated(index: Int, value: Double): DVec =
    checkIndex(index)
    val out = toDoubleArrayOwnedCopy
    out(index) = value
    DVec.fromDoubleArrayOwned(out)

  def slice(from: Int, until: Int): DVec =
    require(from >= 0 && until >= from && until <= length, s"invalid slice [$from, $until) for length $length")
    new DVec(
      data,
      Offset.unsafe(offset.value + from * stride.value),
      Length.unsafe(until - from),
      stride
    )

  def copy: DVec =
    DVec.fromDoubleArrayOwned(toDoubleArrayOwnedCopy)

  /** Immutable copy of this vector's elements in logical order. */
  def toSeq: Seq[Double] =
    val builder = Vector.newBuilder[Double]
    builder.sizeHint(length)
    var i = 0
    var xi = offset.value
    while i < length do
      builder += data(xi)
      xi += stride.value
      i += 1
    builder.result()

  private[gale] def toArray: Array[Double] =
    val out = new Array[Double](length)
    var i = 0
    var xi = offset.value
    while i < length do
      out(i) = data(xi)
      xi += stride.value
      i += 1
    out

  /** Contiguous, freshly-allocated platform copy of the logical elements; the
    * returned storage is owned by the caller (stride 1, offset 0).
    */
  private[gale] def toDoubleArrayOwnedCopy: DoubleArray =
    val out = DoubleArray.alloc(length)
    var i = 0
    var xi = offset.value
    while i < length do
      out(i) = data(xi)
      xi += stride.value
      i += 1
    out

  def +(that: DVec): DVec =
    requireSameLength(that)
    val out = DVec.zeros(length)
    DoubleKernels.dadd(
      length,
      data,
      offset.value,
      stride.value,
      that.data,
      that.offset.value,
      that.stride.value,
      out.data,
      out.offset.value,
      out.stride.value
    )
    out

  def -(that: DVec): DVec =
    requireSameLength(that)
    val out = DVec.zeros(length)
    DoubleKernels.dsub(
      length,
      data,
      offset.value,
      stride.value,
      that.data,
      that.offset.value,
      that.stride.value,
      out.data,
      out.offset.value,
      out.stride.value
    )
    out

  def *(alpha: Double): DVec =
    val out = copy
    DoubleKernels.dscal(out.length, alpha, out.data, out.offset.value, out.stride.value)
    out

  def dot(that: DVec): Double =
    requireSameLength(that)
    DoubleKernels.ddot(
      length,
      data,
      offset.value,
      stride.value,
      that.data,
      that.offset.value,
      that.stride.value
    )

  def norm2: Double =
    DoubleKernels.dnrm2(length, data, offset.value, stride.value)

  /** Mutable copy of this vector; writes to it never affect this value. */
  def mutableCopy: MutableDVec =
    MutableDVec.from(this)

  private def requireSameLength(that: DVec): Unit =
    if length != that.length then
      throw LinAlgError.VectorLengthMismatch(length, that.length)

  private def checkIndex(index: Int): Unit =
    if index < 0 || index >= length then
      throw LinAlgError.IndexOutOfBounds(index, length)

object DVec:
  def zeros(length: Int): DVec =
    require(length >= 0, "length must be non-negative")
    new DVec(DoubleArray.alloc(length), Offset.unsafe(0), Length.unsafe(length), Stride.unsafe(1))

  def fill(length: Int)(value: Double): DVec =
    val out = zeros(length)
    var i = 0
    while i < length do
      out.data(i) = value
      i += 1
    out

  def tabulate(length: Int)(f: Int => Double): DVec =
    val out = zeros(length)
    var i = 0
    while i < length do
      out.data(i) = f(i)
      i += 1
    out

  def fromSeq(values: Seq[Double]): DVec =
    val out = zeros(values.length)
    var i = 0
    values.foreach { value =>
      out.data(i) = value
      i += 1
    }
    out

  private[gale] def fromArray(values: Array[Double]): DVec =
    new DVec(DoubleArray.fromArray(values), Offset.unsafe(0), Length.unsafe(values.length), Stride.unsafe(1))

  /** Wrap a platform array as a contiguous vector without copying; the caller
    * transfers ownership of `data` and must not mutate it afterwards.
    */
  private[gale] def fromDoubleArrayOwned(data: DoubleArray): DVec =
    new DVec(data, Offset.unsafe(0), Length.unsafe(data.length), Stride.unsafe(1))

object Vec:
  def apply(values: Double*): DVec =
    DVec.fromSeq(values)

  def zeros(length: Int): DVec =
    DVec.zeros(length)

  def fill(length: Int)(value: Double): DVec =
    DVec.fill(length)(value)

  def tabulate(length: Int)(f: Int => Double): DVec =
    DVec.tabulate(length)(f)

extension (alpha: Double)
  def *(x: DVec): DVec =
    x * alpha
