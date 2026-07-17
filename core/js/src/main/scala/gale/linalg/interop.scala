package gale.linalg

import gale.platform.DoubleArray
import scala.scalajs.js.typedarray.Float64Array

/** JS-only copy interop between gale's dense types and `Float64Array` values.
  *
  * Shared code keeps raw storage off its public surface. These helpers copy in
  * and copy out, so Gale does not expose an owned typed-array representation.
  *
  * SUBTLE CONTRACT: unlike the JVM exporters, these are named `toFloat64Array`, so
  * they do not shadow a like-named member. They do, however, call the
  * `private[gale]` members `DVec.toDoubleArrayOwnedCopy` /
  * `DMat.toDoubleArrayCopyRowMajor` by name; renaming those would break these
  * exporters at compile time. The `userland.InteropSuite` (in a non-`gale`
  * package) exercises the public surface and is the canary for both platforms.
  */

extension (companion: Vec.type)
  /** Copy `values` into an independently owned vector. */
  def fromFloat64ArrayCopy(values: Float64Array): DVec =
    DVec.fromDoubleArrayOwned(DoubleArray.copy(DoubleArray.adopt(values)))

extension (companion: Matrix.type)
  /** Copy row-major `values` into an independently owned matrix. */
  def fromFloat64ArrayCopy(rows: Int, cols: Int, values: Float64Array): DMat =
    DMat.fromDoubleArrayOwned(rows, cols, DoubleArray.copy(DoubleArray.adopt(values)))

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
