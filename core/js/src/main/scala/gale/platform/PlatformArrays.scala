package gale.platform

import scala.scalajs.js.typedarray.Float32Array
import scala.scalajs.js.typedarray.Float64Array
import scala.scalajs.js.typedarray.Int32Array

opaque type DoubleArray = Float64Array
object DoubleArray:
  private[gale] inline def alloc(length: Int): DoubleArray =
    new Float64Array(length)

  private[gale] def fromArray(values: Array[Double]): DoubleArray =
    val out = alloc(values.length)
    var i = 0
    while i < values.length do
      out(i) = values(i)
      i += 1
    out

  /** Adopt a `Float64Array` as backing storage without copying. The caller
    * hands over ownership: it must not mutate `values` afterwards. Used by the
    * implementation only to avoid a redundant copy of freshly allocated data.
    */
  private[gale] def adopt(values: Float64Array): DoubleArray =
    values

  /** Independent copy of `array`. */
  private[gale] def copy(array: DoubleArray): DoubleArray =
    val out = alloc(array.length)
    var i = 0
    while i < array.length do
      out(i) = array(i)
      i += 1
    out

  /** View an owned `DoubleArray` as its underlying `Float64Array` without
    * copying, for copy-only JS export. The argument must be a
    * freshly-allocated array the caller owns (e.g. from a `*OwnedCopy`).
    */
  private[gale] def asFloat64Array(owned: DoubleArray): Float64Array =
    owned

  /** True when both handles refer to the same underlying storage. Used to
    * reject aliased in-place destinations before a kernel corrupts them.
    */
  private[gale] def sameStorage(x: DoubleArray, y: DoubleArray): Boolean =
    x.asInstanceOf[AnyRef] eq y.asInstanceOf[AnyRef]

  extension (array: DoubleArray)
    private[gale] inline def apply(index: Int): Double =
      array(index)

    private[gale] inline def update(index: Int, value: Double): Unit =
      array(index) = value

    private[gale] inline def length: Int =
      array.length

    private[gale] def toArray: Array[Double] =
      val out = new Array[Double](array.length)
      var i = 0
      while i < array.length do
        out(i) = array(i)
        i += 1
      out

opaque type FloatArray = Float32Array
object FloatArray:
  private[gale] inline def alloc(length: Int): FloatArray =
    new Float32Array(length)

  extension (array: FloatArray)
    private[gale] inline def apply(index: Int): Float =
      array(index).toFloat

    private[gale] inline def update(index: Int, value: Float): Unit =
      array(index) = value

    private[gale] inline def length: Int =
      array.length

opaque type IndexArray = Int32Array
object IndexArray:
  private[gale] inline def alloc(length: Int): IndexArray =
    new Int32Array(length)

  private[gale] def fromArray(values: Array[Int]): IndexArray =
    val out = alloc(values.length)
    var i = 0
    while i < values.length do
      out(i) = values(i)
      i += 1
    out

  /** Adopt an `Int32Array` as backing storage without copying. */
  private[gale] def adopt(values: Int32Array): IndexArray =
    values

  /** Independent copy of `array`. */
  private[gale] def copy(array: IndexArray): IndexArray =
    val out = alloc(array.length)
    var i = 0
    while i < array.length do
      out(i) = array(i)
      i += 1
    out

  extension (array: IndexArray)
    private[gale] inline def apply(index: Int): Int =
      array(index)

    private[gale] inline def update(index: Int, value: Int): Unit =
      array(index) = value

    private[gale] inline def length: Int =
      array.length

    private[gale] def toArray: Array[Int] =
      val out = new Array[Int](array.length)
      var i = 0
      while i < array.length do
        out(i) = array(i)
        i += 1
      out
