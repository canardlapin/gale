package gale.linalg

import gale.platform.DoubleArray
import scala.scalajs.js.typedarray.Float64Array

/** JS-only interop between gale's dense types and `Float64Array` storage.
  *
  * This is the single public doorway to gale's backing typed arrays: shared code
  * keeps raw storage off its public surface. The `*Unsafe` constructors adopt the
  * caller's `Float64Array` without copying; the exporters return fresh copies.
  *
  * SUBTLE CONTRACT: unlike the JVM exporters, these are named `toFloat64Array`, so
  * they do not shadow a like-named member. They do, however, call the
  * `private[gale]` members `DVec.toDoubleArrayOwnedCopy` /
  * `DMat.toDoubleArrayCopyRowMajor` by name; renaming those would break these
  * exporters at compile time. The `userland.InteropSuite` (in a non-`gale`
  * package) exercises the public surface and is the canary for both platforms.
  */

extension (companion: Vec.type)
  /** Wrap `values` as a vector without copying. The vector adopts the typed
    * array, so the caller must not mutate `values` afterwards. Use [[Vec.apply]]
    * / `DVec.fromSeq` when a copy is wanted.
    */
  def fromFloat64ArrayUnsafe(values: Float64Array): DVec =
    DVec.fromDoubleArrayOwned(DoubleArray.adopt(values))

extension (companion: Matrix.type)
  /** Wrap `values` (row-major) as a matrix without copying. The matrix adopts
    * the typed array, so the caller must not mutate `values` afterwards. Use
    * [[Matrix.dense]] when a copy is wanted.
    */
  def fromFloat64ArrayUnsafe(rows: Int, cols: Int, values: Float64Array): DMat =
    DMat.fromDoubleArrayOwned(rows, cols, DoubleArray.adopt(values))

extension (v: DVec)
  /** Fresh `Float64Array` copy of this vector's elements in logical order. */
  def toFloat64Array: Float64Array =
    DoubleArray.asFloat64Array(v.toDoubleArrayOwnedCopy)

extension (m: DMat)
  /** Fresh `Float64Array` copy of this matrix's elements in row-major order,
    * materialising strided or transposed views.
    */
  def toFloat64Array: Float64Array =
    DoubleArray.asFloat64Array(m.toDoubleArrayCopyRowMajor)
