package gale.sparse

import gale.kernel.DoubleKernels
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*

/** Allocation-free callback for the active entries of a [[SparseVector]]. Like
  * [[SparseEntryConsumer]], it has a primitive JVM method signature and never
  * constructs a tuple.
  */
trait SparseVectorEntryConsumer:
  def apply(index: Int, value: Double): Unit

/** An immutable sparse `Double` vector of fixed `length`.
  *
  * Storage is a strictly increasing index array and a parallel value array;
  * every other position is an implicit zero. The stored ("active") entries
  * follow these rules:
  *
  *   - '''Explicit zeros are kept.''' A stored `0.0` remains active through
  *     construction from entries, scaling, [[mapActive]], `+` and `-`, as in
  *     Breeze. [[compact]] is the only operation that drops them, and
  *     [[SparseVector.fromDense]] never stores them in the first place.
  *   - '''Structure-preserving arithmetic.''' `+`/`-` store the union of both
  *     patterns and `*` a scalar keeps this pattern, even when a result
  *     cancels to zero. [[mapActive]] applies to active entries only.
  *   - '''Implicit zeros do not participate in products.''' [[dot]] and
  *     [[axpyInto]] visit active entries only, so a non-finite value facing an
  *     implicit zero contributes nothing (a dense product would form
  *     `NaN * 0.0`). Reductions ([[sum]], [[max]], [[min]], the norms) are the
  *     dense reductions of [[toDense]]: `max`/`min` include `0.0` whenever an
  *     implicit zero exists, and any active `NaN` propagates.
  */
final class SparseVector private[gale] (
    val length: Int,
    private[gale] val indices: IndexArray,
    private[gale] val values: DoubleArray
):
  /** Number of stored entries, including explicit zeros. */
  def activeSize: Int =
    values.length

  /** Value at `index`: the stored value, or `0.0` for an implicit zero.
    * `O(log activeSize)` by binary search.
    */
  def apply(index: Int): Double =
    if index < 0 || index >= length then throw LinAlgError.IndexOutOfBounds(index, length)
    val p = find(index)
    if p >= 0 then values(p) else 0.0

  /** Visit every active entry once in increasing index order. */
  def foreachActive(consumer: SparseVectorEntryConsumer): Unit =
    var p = 0
    while p < values.length do
      consumer(indices(p), values(p))
      p += 1

  /** Independent copy of the active indices, strictly increasing. */
  def activeIndices: IndexedSeq[Int] =
    indices.toArray.toIndexedSeq

  /** Independent copy of the active values, in [[activeIndices]] order. */
  def activeValues: IndexedSeq[Double] =
    values.toArray.toIndexedSeq

  def toDense: DVec =
    val out = DoubleArray.alloc(length)
    var p = 0
    while p < values.length do
      out(indices(p)) = values(p)
      p += 1
    DVec.fromDoubleArrayOwned(out)

  /** Sparse inner product by a two-pointer merge over both active patterns. */
  def dot(that: SparseVector): Double =
    requireSameLength(that.length)
    val aIdx = indices
    val aVal = values
    val bIdx = that.indices
    val bVal = that.values
    val aEnd = aVal.length
    val bEnd = bVal.length
    var pa = 0
    var pb = 0
    var acc = 0.0
    while pa < aEnd && pb < bEnd do
      val ia = aIdx(pa)
      val ib = bIdx(pb)
      if ia < ib then pa += 1
      else if ib < ia then pb += 1
      else
        acc += aVal(pa) * bVal(pb)
        pa += 1
        pb += 1
    acc

  /** Inner product with a dense vector, gathering only the active positions. */
  def dot(that: DVec): Double =
    requireSameLength(that.length)
    val xData = that.data
    val xOff = that.offset.value
    val xStep = that.stride.value
    var acc = 0.0
    var p = 0
    while p < values.length do
      acc += values(p) * xData(xOff + indices(p) * xStep)
      p += 1
    acc

  /** Elementwise sum over the union of both patterns; cancellations stay stored. */
  def +(that: SparseVector): SparseVector =
    merge(that, subtract = false)

  /** Elementwise difference over the union of both patterns; cancellations stay stored. */
  def -(that: SparseVector): SparseVector =
    merge(that, subtract = true)

  /** Scale every active value; the pattern is shared, even for `alpha == 0.0`. */
  def *(alpha: Double): SparseVector =
    mapActive(_ * alpha)

  /** `y += alpha * this`, writing only the active positions of `y`. */
  def axpyInto(alpha: Double, y: MutableDVec): Unit =
    requireSameLength(y.length)
    val yData = y.data
    val yOff = y.offset.value
    val yStep = y.stride.value
    var p = 0
    while p < values.length do
      val at = yOff + indices(p) * yStep
      yData(at) = yData(at) + alpha * values(p)
      p += 1

  /** Apply `f` to the active values only. Implicit zeros are not visited, so
    * the result is the dense map only when `f(0.0) == 0.0`. The pattern,
    * including explicit zeros, is shared.
    */
  def mapActive(f: Double => Double): SparseVector =
    val n = values.length
    val out = DoubleArray.alloc(n)
    var p = 0
    while p < n do
      out(p) = f(values(p))
      p += 1
    new SparseVector(length, indices, out)

  /** Drop explicit zeros (`0.0` and `-0.0`); every other entry is kept. */
  def compact: SparseVector =
    var nonzero = 0
    var p = 0
    while p < values.length do
      if values(p) != 0.0 then nonzero += 1
      p += 1
    if nonzero == values.length then this
    else
      val outIdx = IndexArray.alloc(nonzero)
      val outVal = DoubleArray.alloc(nonzero)
      var write = 0
      p = 0
      while p < values.length do
        if values(p) != 0.0 then
          outIdx(write) = indices(p)
          outVal(write) = values(p)
          write += 1
        p += 1
      new SparseVector(length, outIdx, outVal)

  /** Sum of the active values, accumulated in index order. */
  def sum: Double =
    var acc = 0.0
    var p = 0
    while p < values.length do
      acc += values(p)
      p += 1
    acc

  /** Largest element of the dense vector: includes `0.0` when an implicit zero
    * exists. `NaN` propagates.
    *
    * @throws LinAlgError.EmptyInput for a length-0 vector
    */
  def max: Double =
    extreme("SparseVector.max", largest = true)

  /** Smallest element of the dense vector: includes `0.0` when an implicit
    * zero exists. `NaN` propagates.
    *
    * @throws LinAlgError.EmptyInput for a length-0 vector
    */
  def min: Double =
    extreme("SparseVector.min", largest = false)

  /** Sum of absolute values; `0.0` for an empty vector. */
  def norm1: Double =
    var acc = 0.0
    var p = 0
    while p < values.length do
      acc += math.abs(values(p))
      p += 1
    acc

  /** Euclidean norm with scaled accumulation, so it neither overflows for
    * large entries nor underflows for tiny ones. `NaN` propagates; otherwise
    * any infinite entry gives `+Infinity`.
    */
  def norm2: Double =
    var infinite = false
    var p = 0
    while p < values.length do
      val v = values(p)
      if v.isNaN then return Double.NaN
      if v.isInfinite then infinite = true
      p += 1
    if infinite then Double.PositiveInfinity
    else DoubleKernels.dnrm2(values.length, values, 0, 1)

  /** Largest absolute value; `0.0` for an empty vector. `NaN` propagates. */
  def normInf: Double =
    var out = 0.0
    var p = 0
    while p < values.length do
      val a = math.abs(values(p))
      if a.isNaN then return Double.NaN
      if a > out then out = a
      p += 1
    out

  override def toString: String =
    val shown = math.min(activeSize, 16)
    val entries = (0 until shown).map(p => s"${indices(p)} -> ${values(p)}")
    val more = if activeSize > shown then Seq(s"... ${activeSize - shown} more") else Nil
    (entries ++ more).mkString(s"SparseVector($length)(", ", ", ")")

  private def extreme(operation: String, largest: Boolean): Double =
    if length == 0 then throw LinAlgError.EmptyInput(operation)
    var out = if values.length < length then 0.0 else values(0)
    var p = 0
    while p < values.length do
      val v = values(p)
      if v.isNaN then return Double.NaN
      if (largest && v > out) || (!largest && v < out) then out = v
      p += 1
    out

  private def merge(that: SparseVector, subtract: Boolean): SparseVector =
    requireSameLength(that.length)
    val aIdx = indices
    val aVal = values
    val bIdx = that.indices
    val bVal = that.values
    val aEnd = aVal.length
    val bEnd = bVal.length
    // One merge pass into arrays sized for a disjoint union; a pattern overlap
    // leaves slack that one bulk prefix copy trims. Sizing exactly would need
    // the overlap count, which costs a second branchy merge pass.
    val capacity = math.min(aEnd.toLong + bEnd, length.toLong).toInt
    var outIdx = IndexArray.alloc(capacity)
    var outVal = DoubleArray.alloc(capacity)
    var pa = 0
    var pb = 0
    var write = 0
    while pa < aEnd && pb < bEnd do
      val ia = aIdx(pa)
      val ib = bIdx(pb)
      if ia < ib then
        outIdx(write) = ia
        outVal(write) = aVal(pa)
        pa += 1
      else if ib < ia then
        outIdx(write) = ib
        outVal(write) = if subtract then 0.0 - bVal(pb) else bVal(pb)
        pb += 1
      else
        outIdx(write) = ia
        outVal(write) = if subtract then aVal(pa) - bVal(pb) else aVal(pa) + bVal(pb)
        pa += 1
        pb += 1
      write += 1
    while pa < aEnd do
      outIdx(write) = aIdx(pa)
      outVal(write) = aVal(pa)
      pa += 1
      write += 1
    while pb < bEnd do
      outIdx(write) = bIdx(pb)
      outVal(write) = if subtract then 0.0 - bVal(pb) else bVal(pb)
      pb += 1
      write += 1
    if write < capacity then
      outIdx = IndexArray.copyPrefix(outIdx, write)
      outVal = DoubleArray.copyPrefix(outVal, write)
    new SparseVector(length, outIdx, outVal)

  /** Position of `index` among the active entries, or `-1`. */
  private[sparse] def find(index: Int): Int =
    var lo = 0
    var hi = values.length - 1
    while lo <= hi do
      val mid = (lo + hi) >>> 1
      val at = indices(mid)
      if at < index then lo = mid + 1
      else if at > index then hi = mid - 1
      else return mid
    -1

  private def requireSameLength(thatLength: Int): Unit =
    if length != thatLength then throw LinAlgError.VectorLengthMismatch(length, thatLength)

object SparseVector:
  /** The all-implicit-zero vector of `length`. */
  def zeros(length: Int): SparseVector =
    require(length >= 0, "length must be non-negative")
    new SparseVector(length, IndexArray.alloc(0), DoubleArray.alloc(0))

  /** Store every element of `x` that is not `0.0` (or `-0.0`). Non-finite
    * values are stored verbatim.
    */
  def fromDense(x: DVec): SparseVector =
    val n = x.length
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    var nonzero = 0
    var i = 0
    var at = xOff
    while i < n do
      if xData(at) != 0.0 then nonzero += 1
      at += xStep
      i += 1
    val outIdx = IndexArray.alloc(nonzero)
    val outVal = DoubleArray.alloc(nonzero)
    var write = 0
    i = 0
    at = xOff
    while i < n do
      val v = xData(at)
      if v != 0.0 then
        outIdx(write) = i
        outVal(write) = v
        write += 1
      at += xStep
      i += 1
    new SparseVector(n, outIdx, outVal)

  /** Throwing convenience over [[tryFromEntries]]. */
  def fromEntries(
      length: Int,
      entries: Seq[(Int, Double)],
      duplicates: DuplicatePolicy = DuplicatePolicy.Sum,
      valuePolicy: SparseValuePolicy = SparseValuePolicy.AllowNonFinite
  ): SparseVector =
    tryFromEntries(length, entries, duplicates, valuePolicy) match
      case Left(error)   => throw error
      case Right(vector) => vector

  /** Total construction from `(index, value)` pairs in any order.
    *
    * Entries are stably sorted by index, then repeated indices are resolved by
    * `duplicates`: `Sum` adds them, `Last` keeps the last in input order, and
    * `Error` rejects them. Explicit zeros are stored; use
    * [[SparseVector.compact]] to drop them. `valuePolicy` decides whether
    * `NaN` and infinities are stored or rejected.
    *
    * Failures are values: `InvalidArgument` for a negative length, a rejected
    * duplicate, or a rejected non-finite value, and `IndexOutOfBounds` for an
    * index outside `[0, length)`.
    */
  def tryFromEntries(
      length: Int,
      entries: Seq[(Int, Double)],
      duplicates: DuplicatePolicy = DuplicatePolicy.Sum,
      valuePolicy: SparseValuePolicy = SparseValuePolicy.AllowNonFinite
  ): Either[LinAlgError, SparseVector] =
    val n = entries.length
    val idx = new Array[Int](n)
    val vals = new Array[Double](n)
    var i = 0
    entries.foreach { case (index, value) =>
      idx(i) = index
      vals(i) = value
      i += 1
    }
    build(length, idx, vals, duplicates, valuePolicy)

  /** Array form of [[tryFromEntries]]; the inputs are read, never retained. */
  def tryFromArrays(
      length: Int,
      indices: Array[Int],
      values: Array[Double],
      duplicates: DuplicatePolicy = DuplicatePolicy.Sum,
      valuePolicy: SparseValuePolicy = SparseValuePolicy.AllowNonFinite
  ): Either[LinAlgError, SparseVector] =
    if indices.length != values.length then
      Left(
        LinAlgError.InvalidArgument(
          s"sparse vector index and value arrays differ in length: ${indices.length} vs ${values.length}"
        )
      )
    else build(length, indices, values, duplicates, valuePolicy)

  private def build(
      length: Int,
      idx: Array[Int],
      vals: Array[Double],
      duplicates: DuplicatePolicy,
      valuePolicy: SparseValuePolicy
  ): Either[LinAlgError, SparseVector] =
    if length < 0 then
      return Left(LinAlgError.InvalidArgument(s"sparse vector length must be non-negative, got $length"))
    val n = idx.length
    var i = 0
    while i < n do
      val index = idx(i)
      if index < 0 || index >= length then return Left(LinAlgError.IndexOutOfBounds(index, length))
      if valuePolicy == SparseValuePolicy.RequireFinite && !vals(i).isFinite then
        return Left(LinAlgError.InvalidArgument(s"non-finite sparse value ${vals(i)} at index $index"))
      i += 1
    val order = stableIndexOrder(idx)
    // First pass validates the duplicate policy and counts distinct indices,
    // so the second pass writes exact-size storage once.
    var distinct = 0
    var read = 0
    while read < n do
      val index = idx(order(read))
      var next = read + 1
      while next < n && idx(order(next)) == index do
        if duplicates == DuplicatePolicy.Error then
          return Left(LinAlgError.InvalidArgument(s"duplicate sparse entry at index $index"))
        next += 1
      distinct += 1
      read = next
    val outIdx = IndexArray.alloc(distinct)
    val outVal = DoubleArray.alloc(distinct)
    read = 0
    var write = 0
    while read < n do
      val index = idx(order(read))
      var value = vals(order(read))
      var next = read + 1
      while next < n && idx(order(next)) == index do
        duplicates match
          case DuplicatePolicy.Sum  => value += vals(order(next))
          case DuplicatePolicy.Last => value = vals(order(next))
          case DuplicatePolicy.Error =>
            throw IllegalStateException("duplicate indices were rejected in the first pass")
        next += 1
      outIdx(write) = index
      outVal(write) = value
      write += 1
      read = next
    Right(new SparseVector(length, outIdx, outVal))

  /** Stable ascending order of positions in `keys` by bottom-up primitive
    * mergesort, so `Last` keeps input-order semantics. Already-sorted input is
    * detected in one pass.
    */
  private def stableIndexOrder(keys: Array[Int]): Array[Int] =
    val n = keys.length
    var order = new Array[Int](n)
    var i = 0
    while i < n do
      order(i) = i
      i += 1
    var sorted = true
    i = 1
    while i < n && sorted do
      if keys(i) < keys(i - 1) then sorted = false
      i += 1
    if !sorted then
      var scratch = new Array[Int](n)
      var width = 1
      while width < n do
        var from = 0
        while from < n do
          val middle = math.min(from + width, n)
          val until = math.min(from + 2 * width, n)
          var left = from
          var right = middle
          var out = from
          while out < until do
            if right >= until || (left < middle && keys(order(left)) <= keys(order(right))) then
              scratch(out) = order(left)
              left += 1
            else
              scratch(out) = order(right)
              right += 1
            out += 1
          from = until
        val swap = order
        order = scratch
        scratch = swap
        width *= 2
    order

  /** Sorted, duplicate-summed sparse vector from one compressed slice
    * `[start, end)` of a CSR row or CSC column. Explicit stored zeros are kept;
    * duplicates (possible only in a non-canonical matrix) are summed, matching
    * the matrix-vector sum semantics. A strictly increasing slice, the
    * canonical case, is copied once into exact-size storage.
    */
  private[sparse] def fromCompressedSlice(
      length: Int,
      minor: IndexArray,
      data: DoubleArray,
      start: Int,
      end: Int
  ): SparseVector =
    val width = end - start
    var increasing = true
    var p = start + 1
    while p < end && increasing do
      if minor(p) <= minor(p - 1) then increasing = false
      p += 1
    if increasing then
      val outIdx = IndexArray.alloc(width)
      val outVal = DoubleArray.alloc(width)
      var k = 0
      while k < width do
        outIdx(k) = minor(start + k)
        outVal(k) = data(start + k)
        k += 1
      new SparseVector(length, outIdx, outVal)
    else
      val keys = new Array[Int](width)
      val vals = new Array[Double](width)
      var k = 0
      while k < width do
        keys(k) = minor(start + k)
        vals(k) = data(start + k)
        k += 1
      insertionSortRange(keys, vals, 0, width)
      var distinct = 0
      k = 0
      while k < width do
        if k == 0 || keys(k) != keys(k - 1) then distinct += 1
        k += 1
      val outIdx = IndexArray.alloc(distinct)
      val outVal = DoubleArray.alloc(distinct)
      var write = -1
      k = 0
      while k < width do
        if k == 0 || keys(k) != keys(k - 1) then
          write += 1
          outIdx(write) = keys(k)
          outVal(write) = vals(k)
        else outVal(write) = outVal(write) + vals(k)
        k += 1
      new SparseVector(length, outIdx, outVal)
