package gale.platform

opaque type DoubleArray = Array[Double]
object DoubleArray:
  private[gale] inline def alloc(length: Int): DoubleArray =
    new Array[Double](length)

  private[gale] def fromArray(values: Array[Double]): DoubleArray =
    values.clone()

  /** Adopt a raw `Array[Double]` as backing storage without copying. The caller
    * hands over ownership: it must not mutate `values` afterwards. Used by the
    * JVM interop doorway and internally to avoid a redundant clone.
    */
  private[gale] def adopt(values: Array[Double]): DoubleArray =
    values

  /** Independent copy of `array`. */
  private[gale] def copy(array: DoubleArray): DoubleArray =
    array.clone()

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
      array.clone()

opaque type FloatArray = Array[Float]
object FloatArray:
  private[gale] inline def alloc(length: Int): FloatArray =
    new Array[Float](length)

  extension (array: FloatArray)
    private[gale] inline def apply(index: Int): Float =
      array(index)

    private[gale] inline def update(index: Int, value: Float): Unit =
      array(index) = value

    private[gale] inline def length: Int =
      array.length

opaque type IndexArray = Array[Int]
object IndexArray:
  private[gale] inline def alloc(length: Int): IndexArray =
    new Array[Int](length)

  private[gale] def fromArray(values: Array[Int]): IndexArray =
    values.clone()

  /** Adopt a raw `Array[Int]` as backing storage without copying. */
  private[gale] def adopt(values: Array[Int]): IndexArray =
    values

  /** Independent copy of `array`. */
  private[gale] def copy(array: IndexArray): IndexArray =
    array.clone()

  extension (array: IndexArray)
    private[gale] inline def apply(index: Int): Int =
      array(index)

    private[gale] inline def update(index: Int, value: Int): Unit =
      array(index) = value

    private[gale] inline def length: Int =
      array.length

    private[gale] def toArray: Array[Int] =
      array.clone()
