package gale.interop.breeze

import _root_.breeze.linalg.CSCMatrix
import _root_.breeze.linalg.DenseMatrix
import _root_.breeze.linalg.DenseVector
import gale.linalg.*
import gale.platform.DoubleArray
import gale.sparse.CSC
import gale.sparse.CSR
import gale.sparse.Sparse

/** JVM-only conversions between gale's dense/sparse types and Scala Breeze 2.1.0.
  *
  * This module is the only place gale meets Breeze; `gale-core` stays Breeze-free.
  * Designed for `import gale.interop.breeze.*`, exposing plain functions (not
  * extension methods) so the copy-vs-view distinction is unmissable at every call
  * site.
  *
  * ==Copy vs view contract==
  *
  * The direction of a possible zero-copy view is '''asymmetric''', and this is a
  * deliberate consequence of the two storage models:
  *
  *   - '''gale → Breeze is copy-only.''' gale's backing store is an `opaque type
  *     DoubleArray` whose raw `Array[Double]` never escapes the `gale` package
  *     (the only exporters clone). There is therefore no way to alias gale storage
  *     into a Breeze matrix without a copy — the encapsulation boundary forbids it
  *     by design. [[toBreezeCopy]] always allocates.
  *   - '''Breeze → gale can be a zero-copy view.''' gale's dense types carry an
  *     arbitrary `(offset, rowStride, colStride)` layout, which is a strict
  *     superset of Breeze's `(offset, majorStride, isTranspose)` layout, so '''any'''
  *     Breeze matrix or vector maps onto a gale view that shares its `Array[Double]`.
  *     [[fromBreezeView]] performs no copy and '''aliases''': later writes through
  *     either handle are visible through the other. [[fromBreezeCopy]] is the safe,
  *     non-aliasing alternative.
  *
  * Sparse conversions are '''always copies''' (both directions): gale `CSR` is
  * row-compressed while Breeze only has column-compressed `CSCMatrix`, and the two
  * differ in canonical-ordering/dedup bookkeeping, so aliasing sparse storage would
  * be unsound.
  *
  * All copies are '''bit-exact''': only `Double` values are read and written, never
  * recomputed, so a round trip reproduces every element exactly (no tolerance).
  *
  * | Direction | Function | Copy or view | Complexity |
  * |---|---|---|---|
  * | `DMat`  → `DenseMatrix` | [[toBreezeCopy]]   | copy | O(rows·cols) |
  * | `DVec`  → `DenseVector` | [[toBreezeCopy]]   | copy | O(n) |
  * | `DenseMatrix` → `DMat`  | [[fromBreezeCopy]] | copy | O(rows·cols) |
  * | `DenseMatrix` → `DMat`  | [[fromBreezeView]] | '''view''' (aliases) | O(1) |
  * | `DenseVector` → `DVec`  | [[fromBreezeCopy]] | copy | O(n) |
  * | `DenseVector` → `DVec`  | [[fromBreezeView]] | '''view''' (aliases) | O(1) |
  * | `CSR`/`CSC` → `CSCMatrix` | [[toBreezeCopy]] | copy | O(nnz·log nnz) |
  * | `CSCMatrix` → `CSC` | [[fromBreezeToCsc]] | copy | O(nnz·log nnz) |
  * | `CSCMatrix` → `CSR` | [[fromBreezeToCsr]] | copy | O(nnz·log nnz) |
  */

// ===========================================================================
// Dense matrix
// ===========================================================================

/** Copy `m` into a fresh column-major Breeze `DenseMatrix` (the conventional
  * Breeze layout, `isTranspose = false`). Reads through `m(i, j)`, so any gale
  * layout — strided or transposed view — is materialised correctly.
  */
def toBreezeCopy(m: DMat): DenseMatrix[Double] =
  val rows = m.rows
  val cols = m.cols
  val out  = new Array[Double](rows * cols)
  var j    = 0
  while j < cols do
    var i = 0
    while i < rows do
      out(j * rows + i) = m(i, j) // column-major packing
      i += 1
    j += 1
  new DenseMatrix(rows, cols, out)

/** Copy a Breeze `DenseMatrix` (any layout: column-major, transposed, or a slice)
  * into a fresh contiguous row-major gale `DMat`. Reads through `bm(i, j)`.
  */
def fromBreezeCopy(bm: DenseMatrix[Double]): DMat =
  val rows = bm.rows
  val cols = bm.cols
  val out  = new Array[Double](rows * cols)
  var i    = 0
  while i < rows do
    var j = 0
    while j < cols do
      out(i * cols + j) = bm(i, j) // row-major packing
      j += 1
    i += 1
  Matrix.fromArrayCopy(rows, cols, out)

/** A zero-copy gale `DMat` '''view''' over a Breeze matrix's storage. The result
  * shares `bm.data`: writing through either handle is visible through the other.
  *
  * Breeze's `linearIndex` is `offset + col + row·majorStride` when transposed and
  * `offset + row + col·majorStride` otherwise, which is exactly a gale
  * `(offset, rowStride, colStride)` layout with one stride equal to 1. Use
  * [[fromBreezeCopy]] when independent storage is wanted.
  */
def fromBreezeView(bm: DenseMatrix[Double]): DMat =
  val ms = bm.majorStride
  val (rowStride, colStride) = if bm.isTranspose then (ms, 1) else (1, ms)
  require(
    rowStride > 0 && colStride > 0,
    "cannot view a Breeze matrix with a non-positive stride; use fromBreezeCopy"
  )
  new DMat(
    DoubleArray.adopt(bm.data),
    Offset.unsafe(bm.offset),
    Rows.unsafe(bm.rows),
    Cols.unsafe(bm.cols),
    Stride.unsafe(rowStride),
    Stride.unsafe(colStride)
  )

// ===========================================================================
// Dense vector
// ===========================================================================

/** Copy `v` into a fresh contiguous Breeze `DenseVector`. Reads through `v(i)`, so
  * a strided gale view is materialised correctly.
  */
def toBreezeCopy(v: DVec): DenseVector[Double] =
  val n   = v.length
  val out = new Array[Double](n)
  var i   = 0
  while i < n do
    out(i) = v(i)
    i += 1
  new DenseVector(out)

/** Copy a Breeze `DenseVector` (any offset/stride) into a fresh contiguous gale
  * `DVec`.
  */
def fromBreezeCopy(bv: DenseVector[Double]): DVec =
  val n   = bv.length
  val out = new Array[Double](n)
  var i   = 0
  while i < n do
    out(i) = bv(i)
    i += 1
  Vec.fromArrayCopy(out)

/** A zero-copy gale `DVec` '''view''' over a Breeze vector's storage (shares
  * `bv.data`, aliases). Breeze and gale share the same `offset + i·stride` element
  * map. A reversed (negative-stride) Breeze vector cannot be viewed — gale requires
  * a positive stride — so use [[fromBreezeCopy]] there.
  */
def fromBreezeView(bv: DenseVector[Double]): DVec =
  require(
    bv.stride > 0,
    "cannot view a Breeze vector with a non-positive stride; use fromBreezeCopy"
  )
  new DVec(
    DoubleArray.adopt(bv.data),
    Offset.unsafe(bv.offset),
    Length.unsafe(bv.length),
    Stride.unsafe(bv.stride)
  )

// ===========================================================================
// Sparse (copy-only)
// ===========================================================================

/** Copy a gale `CSR` into a Breeze `CSCMatrix` (same matrix, column-compressed).
  * The gale side is canonicalized first, so the Breeze result is sorted and
  * duplicate-free.
  */
def toBreezeCopy(m: CSR): CSCMatrix[Double] =
  buildBreezeCsc(m.rows, m.cols, m.canonicalize.toCOO)

/** Copy a gale `CSC` into a Breeze `CSCMatrix` (same matrix). */
def toBreezeCopy(m: CSC): CSCMatrix[Double] =
  buildBreezeCsc(m.rows, m.cols, m.canonicalize.toCSR.toCOO)

/** Copy a Breeze `CSCMatrix` into a gale `CSC`. */
def fromBreezeToCsc(bm: CSCMatrix[Double]): CSC =
  breezeToCooBuilder(bm).toCSC()

/** Copy a Breeze `CSCMatrix` into a gale `CSR`. */
def fromBreezeToCsr(bm: CSCMatrix[Double]): CSR =
  breezeToCooBuilder(bm).toCSR()

private def buildBreezeCsc(rows: Int, cols: Int, coo: gale.sparse.COO): CSCMatrix[Double] =
  val builder = new CSCMatrix.Builder[Double](rows, cols)
  coo.entries.foreach(e => builder.add(e.row, e.col, e.value))
  builder.result()

private def breezeToCooBuilder(bm: CSCMatrix[Double]): gale.sparse.COOBuilder =
  val builder = Sparse.coo(bm.rows, bm.cols)
  bm.activeIterator.foreach { case ((row, col), value) => builder.add(row, col, value) }
  builder
