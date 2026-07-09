package gale

import gale.kernel.DoubleKernels
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*
import gale.sparse.COO
import gale.sparse.CSC
import gale.sparse.CSR

/** Controlled access to `private[gale]` storage for the few suites that
  * legitimately need to build strided/aliased layouts or drive a kernel
  * directly. Centralising it here keeps the rest of the test suite on the public
  * API, so any test poking at raw storage does so through this one auditable door.
  */
object TestAccess:
  /** A contiguous platform array holding `values`. */
  def doubleArray(values: Double*): DoubleArray =
    DoubleArray.fromArray(values.toArray)

  /** A platform array of `length` cells all set to `value`. Useful for building
    * a backing buffer whose untouched cells act as sentinels.
    */
  def filled(length: Int, value: Double): DoubleArray =
    val out = DoubleArray.alloc(length)
    var i = 0
    while i < length do
      out(i) = value
      i += 1
    out

  /** A stride-`stride` view of `x`'s values scattered across a fresh backing
    * buffer (the intervening cells are zero). Exercises strided kernels.
    */
  def stridedCopy(x: DVec, stride: Int): DVec =
    require(stride >= 1, "stride must be positive")
    val n = x.length
    val backing = DoubleArray.alloc(stride * math.max(n, 1))
    var i = 0
    while i < n do
      backing(stride * i) = x(i)
      i += 1
    new DVec(backing, Offset.unsafe(0), Length.unsafe(n), Stride.unsafe(stride))

  /** Snapshot of the entire backing storage behind `v` (not just the view),
    * so a test can assert that a strided write left neighbouring cells intact.
    */
  def backingSnapshot(v: MutableDVec): Seq[Double] =
    val data = v.data
    val builder = Vector.newBuilder[Double]
    builder.sizeHint(data.length)
    var i = 0
    while i < data.length do
      builder += data(i)
      i += 1
    builder.result()

  /** The reusable scratch buffer behind a workspace, for asserting genuine
    * reuse across factorizations.
    */
  def workBacking(workspace: DenseWorkspace): DoubleArray =
    workspace.workBacking

  /** True when both handles refer to the same underlying storage. */
  def sameStorage(a: DoubleArray, b: DoubleArray): Boolean =
    DoubleArray.sameStorage(a, b)

  /** A CSR matrix's row-pointer array, as an immutable sequence. */
  def rowPtr(csr: CSR): Seq[Int] =
    indexArraySeq(csr.rowPtr)

  /** A CSR matrix's column-index array, as an immutable sequence. */
  def colIdx(csr: CSR): Seq[Int] =
    indexArraySeq(csr.colIdx)

  /** The raw value storage behind a CSR matrix (for zero-copy assertions). */
  def csrValues(csr: CSR): DoubleArray =
    csr.values

  /** The raw value storage behind a CSC matrix (for zero-copy assertions). */
  def cscValues(csc: CSC): DoubleArray =
    csc.values

  /** Build a CSR directly from raw row-pointer/column-index/value arrays, with
    * no canonicalization. Lets a suite construct deliberately non-canonical input
    * (unsorted columns, duplicate columns, explicit zeros) that the public API
    * never produces, to exercise `canonicalize`/`sortedIndices`/`prune`.
    */
  def csr(rows: Int, cols: Int, rowPtr: Seq[Int], colIdx: Seq[Int], values: Seq[Double]): CSR =
    new CSR(
      rows,
      cols,
      IndexArray.fromArray(rowPtr.toArray),
      IndexArray.fromArray(colIdx.toArray),
      DoubleArray.fromArray(values.toArray)
    )

  /** Build a CSC directly from raw column-pointer/row-index/value arrays, with no
    * canonicalization. The CSC dual of [[csr]].
    */
  def csc(rows: Int, cols: Int, colPtr: Seq[Int], rowIdx: Seq[Int], values: Seq[Double]): CSC =
    new CSC(
      rows,
      cols,
      IndexArray.fromArray(colPtr.toArray),
      IndexArray.fromArray(rowIdx.toArray),
      DoubleArray.fromArray(values.toArray)
    )

  /** Build a COO directly from raw triples, letting a suite include duplicate
    * coordinates that the public builder always sums away. `canonical` sets the
    * advertised flag without validation.
    */
  def coo(rows: Int, cols: Int, entries: Seq[(Int, Int, Double)], canonical: Boolean): COO =
    COO(
      rows,
      cols,
      IndexArray.fromArray(entries.map(_._1).toArray),
      IndexArray.fromArray(entries.map(_._2).toArray),
      DoubleArray.fromArray(entries.map(_._3).toArray),
      canonical
    )

  private def indexArraySeq(indices: IndexArray): Seq[Int] =
    val builder = Vector.newBuilder[Int]
    builder.sizeHint(indices.length)
    var i = 0
    while i < indices.length do
      builder += indices(i)
      i += 1
    builder.result()

  /** A vector view over `data` with an explicit offset/length/stride. */
  def vec(data: DoubleArray, offset: Int, length: Int, stride: Int): DVec =
    new DVec(data, Offset.unsafe(offset), Length.unsafe(length), Stride.unsafe(stride))

  /** A mutable vector view over `data` with an explicit offset/length/stride. */
  def mutableVec(data: DoubleArray, offset: Int, length: Int, stride: Int): MutableDVec =
    new MutableDVec(data, Offset.unsafe(offset), Length.unsafe(length), Stride.unsafe(stride))

  /** A mutable vector view onto `base`'s backing storage, offset into it. Used
    * to write a kernel result into the middle of a larger buffer.
    */
  def mutableViewInto(base: MutableDVec, offset: Int, length: Int, stride: Int): MutableDVec =
    new MutableDVec(base.data, Offset.unsafe(offset), Length.unsafe(length), Stride.unsafe(stride))

  /** A matrix view over `data` with an explicit offset and row/col strides. */
  def mat(data: DoubleArray, offset: Int, rows: Int, cols: Int, rowStride: Int, colStride: Int): DMat =
    new DMat(
      data,
      Offset.unsafe(offset),
      Rows.unsafe(rows),
      Cols.unsafe(cols),
      Stride.unsafe(rowStride),
      Stride.unsafe(colStride)
    )

  /** Overwrite every backing element of `m` with NaN, to prove that kernels
    * claiming `beta == 0` assign their destination rather than read-and-scale it.
    */
  def poisonWithNaN(m: DMat): Unit =
    val data = m.data
    var i = 0
    while i < data.length do
      data(i) = Double.NaN
      i += 1

  /** Drive the raw dense matrix-matrix kernel `C := A * B` with `beta == 0`,
    * exercising whatever storage layouts `a`, `b`, and `c` present.
    */
  def gemm(a: DMat, b: DMat, c: DMat): Unit =
    DoubleKernels.dgemm(
      a.rows,
      b.cols,
      a.cols,
      1.0,
      a.data,
      a.offset.value,
      a.rowStride.value,
      a.colStride.value,
      b.data,
      b.offset.value,
      b.rowStride.value,
      b.colStride.value,
      0.0,
      c.data,
      c.offset.value,
      c.rowStride.value,
      c.colStride.value
    )
