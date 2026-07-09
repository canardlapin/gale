package gale.linalg

import gale.platform.DoubleArray

/** JVM-only interop between gale's dense types and raw `Array[Double]` storage.
  *
  * This is the single public doorway to gale's backing arrays: shared code keeps
  * `Array[Double]` off its public surface, and callers that genuinely need raw
  * interop reach for these platform-specific helpers. The `*Unsafe` constructors
  * adopt the caller's array without copying; the exporters return fresh copies.
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
  /** Wrap `values` as a vector without copying. The vector adopts the array, so
    * the caller must not mutate `values` afterwards or the vector's contents
    * change with it. Use [[Vec.apply]] / `DVec.fromSeq` when a copy is wanted.
    */
  def fromArrayUnsafe(values: Array[Double]): DVec =
    DVec.fromDoubleArrayOwned(DoubleArray.adopt(values))

extension (companion: Matrix.type)
  /** Wrap `values` (row-major) as a matrix without copying. The matrix adopts
    * the array, so the caller must not mutate `values` afterwards. Use
    * [[Matrix.dense]] when a copy is wanted.
    */
  def fromArrayUnsafe(rows: Int, cols: Int, values: Array[Double]): DMat =
    DMat.fromDoubleArrayOwned(rows, cols, DoubleArray.adopt(values))

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
