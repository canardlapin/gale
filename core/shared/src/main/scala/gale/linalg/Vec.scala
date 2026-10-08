package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

trait Vec[A]:
  def length: Int
  def apply(index: Int): A

/** An immutable dense `Double` vector: a length-`n` window over platform
  * storage described by an offset and a (strictly positive) stride.
  *
  *   - '''Views vs copies.''' [[slice]] returns an aliasing view — it shares the
  *     immutable backing storage, so it is `O(1)`. [[copy]] (and the `Array`/`Seq`
  *     exporters) return independent data. No unqualified public core API can
  *     create a mutable alias of a `DVec`; [[mutableCopy]] always owns an
  *     independent copy. Optional interop modules must prefix any deliberate
  *     external-storage alias with `unsafe`.
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

  /** Gather elements in caller-specified order into one independently owned
    * contiguous vector. Repeated indices repeat values.
    */
  def gather(indices: IndexedSeq[Int]): DVec =
    var i = 0
    while i < indices.length do
      checkIndex(indices(i))
      i += 1
    val out = DVecBuilder.zeros(indices.length)
    i = 0
    while i < indices.length do
      out(i) = apply(indices(i))
      i += 1
    out.result()

  def copy: DVec =
    DVec.fromDoubleArrayOwned(toDoubleArrayOwnedCopy)

  /** One owned logical copy, ready for in-place construction edits. */
  def toBuilder: DVecBuilder =
    DVecBuilder.from(this)

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

  /** Visit logical elements without allocating an intermediate collection. */
  def foreachValue(f: Double => Unit): Unit =
    var i = 0
    var xi = offset.value
    while i < length do
      f(data(xi))
      xi += stride.value
      i += 1

  /** Copy logical elements into caller-owned primitive storage. No Gale backing
    * storage is exposed or adopted.
    */
  def copyTo(destination: Array[Double], destinationOffset: Int = 0): Unit =
    val end = destinationOffset.toLong + length.toLong
    if destinationOffset < 0 || end > destination.length.toLong then
      throw LinAlgError.InvalidArgument(
        s"destination range [$destinationOffset, $end) does not fit array length ${destination.length}"
      )
    var write = destinationOffset
    var i = 0
    var xi = offset.value
    while i < length do
      destination(write) = data(xi)
      write += 1
      xi += stride.value
      i += 1

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

  /** Sum of absolute values; `0.0` for an empty vector. */
  def norm1: Double =
    DoubleKernels.dasum(length, data, offset.value, stride.value)

  /** Largest absolute value; `0.0` for an empty vector, NaN if any entry is NaN. */
  def normInf: Double =
    DoubleKernels.damax(length, data, offset.value, stride.value)

  /** Fast sum; `0.0` for an empty vector.
    *
    * The kernel may reassociate (it keeps several partial sums), so the result
    * is deterministic for one build and platform but can differ from a
    * left-to-right loop, and between JVM and Scala.js, in the final bits. NaN
    * and infinities follow IEEE addition. Use [[sumExact]] when the result must
    * be correctly rounded and independent of order and platform.
    */
  def sum: Double =
    DoubleKernels.dsum(length, data, offset.value, stride.value)

  /** The exact sum of the entries, rounded once to the nearest `Double`
    * (ties to even), via [[gale.numeric.ExactSum]]. Identical on every platform
    * and for every permutation of the entries; NaN and infinities follow IEEE
    * addition. Slower than [[sum]].
    */
  def sumExact: Double =
    DVec.exactSum(data, offset.value, length, stride.value, 1, 0)

  /** Arithmetic mean, `sum / length`. As in NumPy, the intermediate sum can
    * overflow to `±Inf` even when the mean is representable (entries near
    * `Double.MaxValue`). Throws [[LinAlgError.EmptyInput]] when empty.
    */
  def mean: Double =
    requireNonEmpty("mean")
    sum / length

  /** Largest entry; NaN if any entry is NaN. Throws [[LinAlgError.EmptyInput]]
    * when empty.
    */
  def max: Double =
    requireNonEmpty("max")
    data(offset.value + DoubleKernels.dmaxIndex(length, data, offset.value, stride.value) * stride.value)

  /** Smallest entry; NaN if any entry is NaN. Throws [[LinAlgError.EmptyInput]]
    * when empty.
    */
  def min: Double =
    requireNonEmpty("min")
    data(offset.value + DoubleKernels.dminIndex(length, data, offset.value, stride.value) * stride.value)

  /** Index of the first largest entry, or of the first NaN if any entry is NaN.
    * Throws [[LinAlgError.EmptyInput]] when empty.
    */
  def argmax: Int =
    requireNonEmpty("argmax")
    DoubleKernels.dmaxIndex(length, data, offset.value, stride.value)

  /** Index of the first smallest entry, or of the first NaN if any entry is NaN.
    * Throws [[LinAlgError.EmptyInput]] when empty.
    */
  def argmin: Int =
    requireNonEmpty("argmin")
    DoubleKernels.dminIndex(length, data, offset.value, stride.value)

  private def requireNonEmpty(operation: String): Unit =
    if length == 0 then throw LinAlgError.EmptyInput(s"DVec.$operation")

  /** Mutable copy of this vector; writes to it never affect this value. */
  def mutableCopy: MutableDVec =
    MutableDVec.from(this)

  private def requireSameLength(that: DVec): Unit =
    if length != that.length then
      throw LinAlgError.VectorLengthMismatch(length, that.length)

  private def checkIndex(index: Int): Unit =
    if index < 0 || index >= length then
      throw LinAlgError.IndexOutOfBounds(index, length)

/** Single-owner construction buffer for an immutable [[DVec]].
  *
  * The builder is the primitive-loop construction seam for callers that already
  * know the vector length. [[result]] transfers the backing storage without a
  * copy and permanently closes the builder, so no mutable alias survives through
  * the public API.
  */
final class DVecBuilder private (val length: Int, private[gale] val data: DoubleArray):
  private var open = true

  def apply(index: Int): Double =
    requireOpen()
    checkIndex(index)
    data(index)

  def update(index: Int, value: Double): Unit =
    requireOpen()
    checkIndex(index)
    data(index) = value

  def fill(value: Double): Unit =
    requireOpen()
    var i = 0
    while i < length do
      data(i) = value
      i += 1

  /** Transfer this builder's storage to an immutable vector without copying. */
  def result(): DVec =
    requireOpen()
    open = false
    DVec.fromDoubleArrayOwned(data)

  private def requireOpen(): Unit =
    if !open then
      throw LinAlgError.UnsupportedOperation("DVecBuilder is closed after result()")

  private def checkIndex(index: Int): Unit =
    if index < 0 || index >= length then
      throw LinAlgError.IndexOutOfBounds(index, length)

object DVecBuilder:
  def zeros(length: Int): DVecBuilder =
    require(length >= 0, "length must be non-negative")
    new DVecBuilder(length, DoubleArray.alloc(length))

  /** Initialize from one owned logical copy of a contiguous or strided vector. */
  def from(vector: DVec): DVecBuilder =
    new DVecBuilder(vector.length, vector.toDoubleArrayOwnedCopy)

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

  /** Allocate a single-owner builder whose [[DVecBuilder.result]] transfers its
    * storage into an immutable vector without copying.
    */
  def newBuilder(length: Int): DVecBuilder =
    DVecBuilder.zeros(length)

  def builderFrom(vector: DVec): DVecBuilder =
    DVecBuilder.from(vector)

  /** Correctly rounded sum of `lines` strided runs of `lineLength` elements. */
  private[gale] def exactSum(
      data: DoubleArray,
      offset: Int,
      lineLength: Int,
      elementStride: Int,
      lines: Int,
      lineStride: Int
  ): Double =
    val accumulator = gale.numeric.ExactSum.zero()
    var line = 0
    while line < lines do
      var i = 0
      var xi = offset + line * lineStride
      while i < lineLength do
        // Gale storage holds at most Int.MaxValue entries, far below
        // `ExactSum.MaxTerms` (2^60), so this refusal is unreachable; it is
        // still mapped to a typed error rather than discarded.
        accumulator.add(data(xi)) match
          case Left(error) => throw LinAlgError.InvalidArgument(error.message)
          case Right(())   => ()
        xi += elementStride
        i += 1
      line += 1
    accumulator.value

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

  def newBuilder(length: Int): DVecBuilder =
    DVec.newBuilder(length)

  def builderFrom(vector: DVec): DVecBuilder =
    DVec.builderFrom(vector)

extension (alpha: Double)
  def *(x: DVec): DVec =
    x * alpha

  def *(a: DMat): DMat =
    a * alpha
