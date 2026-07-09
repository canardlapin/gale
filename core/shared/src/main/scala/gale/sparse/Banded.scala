package gale.sparse

import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** A banded matrix with `kl` subdiagonals and `ku` superdiagonals, stored in the
  * LAPACK general-band layout.
  *
  * The band is packed into a `(kl + ku + 1) x cols` dense array `band` (row-major
  * here), with `A(i, j)` living at `band((ku + i - j) * cols + j)` whenever
  * `-ku <= i - j <= kl`; entries outside the band are structural zeros and are
  * not stored. `A(i, i)` sits on band row `ku`. Matvec and transpose-matvec walk
  * only the stored band, so cost is `O(nnz)` rather than `O(rows * cols)`.
  */
final class Banded private[gale] (
    val rows: Int,
    val cols: Int,
    val kl: Int,
    val ku: Int,
    private[gale] val band: DoubleArray
) extends SparseMatrix[Double]
    with DoubleLinearOperator:

  def nnz: Int =
    var count = 0
    var i = 0
    while i < rows do
      val jStart = math.max(0, i - kl)
      val jEnd = math.min(cols - 1, i + ku)
      var bandIdx = (ku + i - jStart) * cols + jStart
      var j = jStart
      while j <= jEnd do
        if band(bandIdx) != 0.0 then count += 1
        bandIdx -= (cols - 1)
        j += 1
      i += 1
    count

  def apply(row: Int, col: Int): Double =
    if row < 0 || row >= rows then throw LinAlgError.IndexOutOfBounds(row, rows)
    if col < 0 || col >= cols then throw LinAlgError.IndexOutOfBounds(col, cols)
    val d = row - col
    if d > kl || d < -ku then 0.0
    else band((ku + d) * cols + col)

  def row(index: Int): DVec =
    if index < 0 || index >= rows then throw LinAlgError.IndexOutOfBounds(index, rows)
    val out = DVec.zeros(cols)
    val jStart = math.max(0, index - kl)
    val jEnd = math.min(cols - 1, index + ku)
    var bandIdx = (ku + index - jStart) * cols + jStart
    var j = jStart
    while j <= jEnd do
      out.data(j) = band(bandIdx)
      bandIdx -= (cols - 1)
      j += 1
    out

  def col(index: Int): DVec =
    if index < 0 || index >= cols then throw LinAlgError.IndexOutOfBounds(index, cols)
    val out = DVec.zeros(rows)
    val iStart = math.max(0, index - ku)
    val iEnd = math.min(rows - 1, index + kl)
    var bandIdx = (ku + iStart - index) * cols + index
    var i = iStart
    while i <= iEnd do
      out.data(i) = band(bandIdx)
      bandIdx += cols
      i += 1
    out

  override def applyTo(x: DVec, into: MutableDVec): Unit =
    mulInto(x, into)

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    tMulInto(x, into)

  override def *(x: DVec): DVec =
    val out = MutableDVec.zeros(rows)
    mulInto(x, out)
    out.asVec

  /** `y := A x`, summing each row's stored band against `x`. */
  def mulInto(x: DVec, y: MutableDVec): Unit =
    if x.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if y.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(y.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, y.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val bandData = band
    val n = cols
    val yData = y.data
    val yOff = y.offset.value
    val yStep = y.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    var i = 0
    while i < rows do
      val jStart = math.max(0, i - kl)
      val jEnd = math.min(n - 1, i + ku)
      // Fixed i: band index (ku + i - j) * n + j decreases by (n - 1) as j grows.
      var bandIdx = (ku + i - jStart) * n + jStart
      var xIdx = xOff + jStart * xStep
      var acc = 0.0
      var j = jStart
      while j <= jEnd do
        acc += bandData(bandIdx) * xData(xIdx)
        bandIdx -= (n - 1)
        xIdx += xStep
        j += 1
      yData(yOff + i * yStep) = acc
      i += 1

  /** `y := Aᵀ x`, gathering each column's stored band against `x`. */
  def tMulInto(x: DVec, y: MutableDVec): Unit =
    if x.length != rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(rows), Cols(1)), Shape(Rows(x.length), Cols(1)))
    if y.length != cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(cols), Cols(1)), Shape(Rows(y.length), Cols(1)))
    if DoubleArray.sameStorage(x.data, y.data) then
      throw LinAlgError.UnsupportedOperation("aliased mulInto destination")
    val bandData = band
    val n = cols
    val yData = y.data
    val yOff = y.offset.value
    val yStep = y.stride.value
    val xData = x.data
    val xOff = x.offset.value
    val xStep = x.stride.value
    var j = 0
    while j < n do
      val iStart = math.max(0, j - ku)
      val iEnd = math.min(rows - 1, j + kl)
      // Fixed j: band index (ku + i - j) * n + j grows by n as i grows.
      var bandIdx = (ku + iStart - j) * n + j
      var xIdx = xOff + iStart * xStep
      var acc = 0.0
      var i = iStart
      while i <= iEnd do
        acc += bandData(bandIdx) * xData(xIdx)
        bandIdx += n
        xIdx += xStep
        i += 1
      yData(yOff + j * yStep) = acc
      j += 1

  def toDense(maxEntries: Int = Int.MaxValue): DMat =
    val entries = rows.toLong * cols.toLong
    if entries > maxEntries.toLong then
      throw LinAlgError.UnsupportedOperation(s"dense conversion would allocate $entries entries")
    val out = DMat.zeros(rows, cols)
    val bandData = band
    val n = cols
    var i = 0
    while i < rows do
      val jStart = math.max(0, i - kl)
      val jEnd = math.min(n - 1, i + ku)
      var bandIdx = (ku + i - jStart) * n + jStart
      var j = jStart
      while j <= jEnd do
        out.data(i * n + j) = bandData(bandIdx)
        bandIdx -= (n - 1)
        j += 1
      i += 1
    out

  def diagonal: DVec =
    val k = math.min(rows, cols)
    val out = DVec.zeros(k)
    val bandData = band
    val n = cols
    var i = 0
    while i < k do
      out.data(i) = bandData(ku * n + i)
      i += 1
    out

  def trace: Double =
    val k = math.min(rows, cols)
    val bandData = band
    val n = cols
    var out = 0.0
    var i = 0
    while i < k do
      out += bandData(ku * n + i)
      i += 1
    out

  /** Transpose: an `n x m` banded matrix with `kl` and `ku` swapped. The stored
    * band is remapped (not shared, since the packed layout differs).
    */
  def t: Banded =
    val m = rows
    val n = cols
    val bandData = band
    val out = DoubleArray.alloc((kl + ku + 1) * m)
    var i = 0
    while i < m do
      val jStart = math.max(0, i - kl)
      val jEnd = math.min(n - 1, i + ku)
      var srcIdx = (ku + i - jStart) * n + jStart
      var j = jStart
      while j <= jEnd do
        // Aᵀ(j, i) = A(i, j); in the transpose (kuᵀ = kl, colsᵀ = m) it lands at
        // (kl + j - i) * m + i.
        out((kl + j - i) * m + i) = bandData(srcIdx)
        srcIdx -= (n - 1)
        j += 1
      i += 1
    new Banded(n, m, ku, kl, out)

object Banded:
  /** Build from diagonals keyed by offset: `0` is the main diagonal, `d > 0` the
    * `d`-th superdiagonal, `d < 0` the `|d|`-th subdiagonal. Every entry is
    * `A(i, i + d)` as `i` runs over the diagonal, and each sequence must match its
    * diagonal's length. `kl`/`ku` are inferred from the most extreme offsets.
    */
  def fromDiagonals(rows: Int, cols: Int, diagonals: Map[Int, Seq[Double]]): Banded =
    require(rows >= 0 && cols >= 0, "banded shape must be non-negative")
    var kl = 0
    var ku = 0
    diagonals.keysIterator.foreach { d =>
      if d > 0 then ku = math.max(ku, d)
      else if d < 0 then kl = math.max(kl, -d)
    }
    val n = cols
    val out = DoubleArray.alloc((kl + ku + 1) * n)
    diagonals.foreach { case (d, values) =>
      val iStart = math.max(0, -d)
      val iEnd = math.min(rows - 1, cols - 1 - d)
      val len = math.max(0, iEnd - iStart + 1)
      require(
        values.length == len,
        s"diagonal $d must have $len entries for a ${rows}x$cols matrix, got ${values.length}"
      )
      var k = 0
      val it = values.iterator
      while k < len do
        val i = iStart + k
        val j = i + d
        out((ku + i - j) * n + j) = it.next()
        k += 1
    }
    new Banded(rows, cols, kl, ku, out)

  /** Build from a dense matrix, detecting the tightest lower/upper bandwidth that
    * contains every nonzero, then packing that band (in-band zeros are stored).
    */
  def fromDense(A: DMat): Banded =
    val (kl, ku) = detectBandwidth(A)
    val rows = A.rows
    val cols = A.cols
    val out = DoubleArray.alloc((kl + ku + 1) * cols)
    var i = 0
    while i < rows do
      val jStart = math.max(0, i - kl)
      val jEnd = math.min(cols - 1, i + ku)
      var j = jStart
      while j <= jEnd do
        out((ku + i - j) * cols + j) = A(i, j)
        j += 1
      i += 1
    new Banded(rows, cols, kl, ku, out)

  /** The tightest `(kl, ku)` whose band contains every nonzero of `A`. An all-zero
    * matrix yields `(0, 0)`.
    */
  def detectBandwidth(A: DMat): (Int, Int) =
    var kl = 0
    var ku = 0
    var i = 0
    while i < A.rows do
      var j = 0
      while j < A.cols do
        if A(i, j) != 0.0 then
          val d = i - j
          if d > 0 then kl = math.max(kl, d)
          else if d < 0 then ku = math.max(ku, -d)
        j += 1
      i += 1
    (kl, ku)
