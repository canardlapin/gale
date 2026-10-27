package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** One value per row or per column. The invariant axis parameter preserves orientation through arithmetic; the stored
  * axis also supports dynamic choices. This is deliberately not a `DVec` subtype: use [[toDVec]] to discard the axis
  * explicitly. Tagging and unwrapping share immutable storage without copying.
  *
  * @since 0.2
  */
final class Margin[A <: Axis] private[linalg] (val axis: A, private val values: DVec):
  def length: Int = values.length
  def apply(index: Int): Double = values(index)
  def toDVec: DVec = values

  /** Map values while retaining the axis and length, into independent storage. */
  def map(f: Double => Double): Margin[A] =
    Margin.tag(Vec.tabulate(length)(i => f(values(i))), axis)

  /** Dynamic axes must agree as well as statically known axes. */
  def zipMap(that: Margin[A])(f: (Double, Double) => Double): Margin[A] =
    requireSameAxis(that)
    if length != that.length then throw LinAlgError.VectorLengthMismatch(length, that.length)
    Margin.tag(Vec.tabulate(length)(i => f(values(i), that(i))), axis)

  def +(that: Margin[A]): Margin[A] = zipMap(that)(_ + _)
  def -(that: Margin[A]): Margin[A] = zipMap(that)(_ - _)
  def +(scalar: Double): Margin[A] = map(_ + scalar)
  def -(scalar: Double): Margin[A] = map(_ - scalar)
  def *(scalar: Double): Margin[A] = map(_ * scalar)
  def /(scalar: Double): Margin[A] = map(_ / scalar)

  private[linalg] def requireSameAxis(that: Margin[A]): Unit =
    if axis != that.axis then throw LinAlgError.InvalidArgument(s"margin axes differ: $axis and ${that.axis}")

object Margin:
  /** One value per row. */
  type ByRow = Margin[Axis.Rows.type]

  /** One value per column. */
  type ByCol = Margin[Axis.Cols.type]

  private[linalg] def tag[A <: Axis](values: DVec, axis: A): Margin[A] = new Margin(axis, values)

/** Additive axis-aware facade over the existing reduction algorithms. Empty and nonfinite behavior is exactly that of
  * the corresponding `DMat` reduction.
  *
  * @since 0.2
  */
final class AxisReductions[A <: Axis] private[linalg] (matrix: DMat, val axis: A):
  def sum: Margin[A] = Margin.tag(matrix.sum(axis), axis)
  def mean: Margin[A] = Margin.tag(matrix.mean(axis), axis)
  def min: Margin[A] = Margin.tag(matrix.min(axis), axis)
  def max: Margin[A] = Margin.tag(matrix.max(axis), axis)

/** Broadcast traversal stays internal: zero read increments never escape as public DMat/DVec strides. Every destination
  * is an open, nonaliased builder.
  */
private[gale] object MarginOps:
  private def requireLength(a: DMat, margin: Margin[?]): Unit =
    val expected = if margin.axis == Axis.Rows then a.rows else a.cols
    if margin.length != expected then throw LinAlgError.VectorLengthMismatch(expected, margin.length)

  private def check(a: DMat, margin: Margin[?], destination: DMatBuilder): DoubleArray =
    val out = destination.writableData
    requireLength(a, margin)
    if destination.rows != a.rows || destination.cols != a.cols then
      throw LinAlgError.DimensionMismatch(a.shape, Shape(Rows(destination.rows), Cols(destination.cols)))
    if DoubleArray.sameStorage(a.data, out) || DoubleArray.sameStorage(margin.toDVec.data, out) then
      throw LinAlgError.UnsupportedOperation("aliased margin destination")
    out

  inline def binary(a: DMat, margin: Margin[?])(inline f: (Double, Double) => Double): DMat =
    requireLength(a, margin)
    val out = DMat.newBuilder(a.rows, a.cols)
    binaryInto(a, margin, out)(f)
    out.result()

  inline def binaryInto(a: DMat, margin: Margin[?], destination: DMatBuilder)(
      inline f: (Double, Double) => Double
  ): Unit =
    val out = check(a, margin, destination)
    val v = margin.toDVec
    val perRow = margin.axis == Axis.Rows
    val readStep = if perRow then 0 else v.stride.value
    var row = 0
    while row < a.rows do
      val start = v.offset.value + (if perRow then row * v.stride.value else 0)
      if perRow && a.cols > 0 then
        val value = v.data(start)
        DoubleKernels.dmapInto(
          a.cols,
          a.data,
          a.offset.value + row * a.rowStride.value,
          a.colStride.value,
          out,
          row * a.cols,
          1
        )(x => f(x, value))
      else
        DoubleKernels.dzipInto(
          a.cols,
          a.data,
          a.offset.value + row * a.rowStride.value,
          a.colStride.value,
          v.data,
          start,
          readStep,
          out,
          row * a.cols,
          1
        )(f)
      row += 1

  inline def fused[A <: Axis](a: DMat, first: Margin[A], second: Margin[A])(
      inline f: (Double, Double, Double) => Double
  ): DMat =
    first.requireSameAxis(second)
    requireLength(a, first)
    requireLength(a, second)
    val out = DMat.newBuilder(a.rows, a.cols)
    fusedInto(a, first, second, out)(f)
    out.result()

  inline def fusedInto[A <: Axis](a: DMat, first: Margin[A], second: Margin[A], destination: DMatBuilder)(
      inline f: (Double, Double, Double) => Double
  ): Unit =
    first.requireSameAxis(second)
    val out = check(a, first, destination)
    requireLength(a, second)
    val x = first.toDVec
    val y = second.toDVec
    if DoubleArray.sameStorage(y.data, out) then throw LinAlgError.UnsupportedOperation("aliased margin destination")
    val perRow = first.axis == Axis.Rows
    val xStep = if perRow then 0 else x.stride.value
    val yStep = if perRow then 0 else y.stride.value
    var row = 0
    while row < a.rows do
      var ai = a.offset.value + row * a.rowStride.value
      var xi = x.offset.value + (if perRow then row * x.stride.value else 0)
      var yi = y.offset.value + (if perRow then row * y.stride.value else 0)
      var oi = row * a.cols
      val end = oi + a.cols
      if perRow && a.cols > 0 then
        val firstValue = x.data(xi)
        val secondValue = y.data(yi)
        DoubleKernels.dmapInto(a.cols, a.data, ai, a.colStride.value, out, oi, 1)(value =>
          f(value, firstValue, secondValue)
        )
      else
        while oi < end do
          out(oi) = f(a.data(ai), x.data(xi), y.data(yi))
          ai += a.colStride.value
          xi += xStep
          yi += yStep
          oi += 1
      row += 1
