package gale.linalg

import gale.platform.DoubleArray

/** JVM-only copy interop between gale's dense types and `Array[Double]` values.
  *
  * Shared code keeps `Array[Double]` off its public surface. These platform
  * helpers copy in and copy out, so the public API does not promise or expose
  * gale's owned storage representation.
  *
  * SUBTLE CONTRACT: the `toArray` / `toArrayRowMajor` exporters below share their
  * names with the `private[gale]` members `DVec.toArray` / `DMat.toArrayRowMajor`.
  * Inside these extension bodies, `v.toArray` / `m.toArrayRowMajor` resolve to the
  * private members (a member wins over a same-named extension, and the private
  * members are visible from this `gale.linalg` file), while external callers —
  * outside `gale` — see only the public extension. Renaming those private members
  * would silently break this resolution: `v.toArray` inside the extension would
  * then bind to the extension itself and recurse. The `userland.InteropSuite`
  * (deliberately in a non-`gale` package) is the canary that keeps this honest.
  */

extension (companion: Vec.type)
  /** Copy `values` into an independently owned vector. */
  def fromArrayCopy(values: Array[Double]): DVec =
    DVec.fromDoubleArrayOwned(DoubleArray.fromArray(values))

extension (companion: Matrix.type)
  /** Copy row-major `values` into an independently owned matrix. */
  def fromArrayCopy(rows: Int, cols: Int, values: Array[Double]): DMat =
    DMat.fromDoubleArrayOwned(rows, cols, DoubleArray.fromArray(values))

extension (v: DVec)
  /** Fresh `Array[Double]` copy of this vector's elements in logical order. */
  def toArray: Array[Double] =
    v.toArray

extension (m: DMat)
  /** Fresh `Array[Double]` copy of this matrix's elements in row-major order,
    * materialising strided or transposed views.
    */
  def toArrayRowMajor: Array[Double] =
    m.toArrayRowMajor
