package gale.sparse.direct.pure

import gale.linalg.LinAlgError
import gale.sparse.CSR
import gale.sparse.CSRPattern

import scala.collection.mutable
import scala.util.Sorting

/** Symbolic and numeric kernels for the portable sparse Cholesky provider.
  *
  * The graph is the undirected lower triangle of the input pattern (upper
  * stored entries are ignored, matching dense Cholesky). `ProviderDefault`
  * is minimum degree with lowest-index tie-break, not the full AMD
  * (Amestoy–Davis–Duff) algorithm.
  */
private[pure] object SparseCholeskyKernels:
  val FillFactor: Long = 64L
  val AbsoluteNnzCap: Long = 8_000_000L

  def fillCap(inputNnz: Int): Long =
    math.max(inputNnz.toLong * FillFactor, AbsoluteNnzCap)

  def exceedsFill(inputNnz: Int, predicted: Long): Boolean =
    predicted > fillCap(inputNnz)

  def identityPerm(n: Int): Array[Int] =
    val perm = new Array[Int](n)
    var i = 0
    while i < n do
      perm(i) = i
      i += 1
    perm

  def invertPerm(perm: Array[Int]): Array[Int] =
    val inv = new Array[Int](perm.length)
    var i = 0
    while i < perm.length do
      inv(perm(i)) = i
      i += 1
    inv

  /** Undirected adjacency from the strict lower triangle. */
  def undirectedLowerAdj(pattern: CSRPattern): Array[Array[Int]] =
    val n = pattern.rows
    val buf = Array.fill(n)(new mutable.ArrayBuffer[Int](4))
    pattern.foreachStoredPosition: (row, col) =>
      if row > col then
        buf(row) += col
        buf(col) += row
    val out = new Array[Array[Int]](n)
    var i = 0
    while i < n do
      val uniq = buf(i).distinct.toArray
      Sorting.quickSort(uniq)
      out(i) = uniq
      i += 1
    out

  /** Minimum degree on an undirected graph. `perm(k)` is the original vertex
    * placed at elimination position `k`. Ties take the lowest index.
    */
  def minimumDegree(adj: Array[Array[Int]]): Array[Int] =
    val n = adj.length
    val nbr = Array.fill(n)(mutable.HashSet.empty[Int])
    var i = 0
    while i < n do
      var t = 0
      val row = adj(i)
      while t < row.length do
        nbr(i) += row(t)
        t += 1
      i += 1
    val active = Array.fill(n)(true)
    val perm = new Array[Int](n)
    var k = 0
    while k < n do
      var best = -1
      var bestDeg = Int.MaxValue
      var v = 0
      while v < n do
        if active(v) then
          val deg = nbr(v).size
          if deg < bestDeg then
            bestDeg = deg
            best = v
        v += 1
      perm(k) = best
      active(best) = false
      val survivors = nbr(best).iterator.filter(active).toArray
      var a = 0
      while a < survivors.length do
        nbr(survivors(a)) -= best
        var b = a + 1
        while b < survivors.length do
          nbr(survivors(a)) += survivors(b)
          nbr(survivors(b)) += survivors(a)
          b += 1
        a += 1
      nbr(best).clear()
      k += 1
    perm

  def permuteAdj(adj: Array[Array[Int]], perm: Array[Int], inv: Array[Int]): Array[Array[Int]] =
    val n = adj.length
    val out = new Array[Array[Int]](n)
    var newI = 0
    while newI < n do
      val oldNbrs = adj(perm(newI))
      val mapped = new Array[Int](oldNbrs.length)
      var t = 0
      while t < oldNbrs.length do
        mapped(t) = inv(oldNbrs(t))
        t += 1
      Sorting.quickSort(mapped)
      out(newI) = mapped
      newI += 1
    out

  /** Elimination tree: `parent(k) = -1` for a root, otherwise a later column. */
  def eliminationTree(n: Int, adjPerm: Array[Array[Int]]): Array[Int] =
    val parent = Array.fill(n)(-1)
    val ancestor = Array.fill(n)(-1)
    var j = 0
    while j < n do
      val nbrs = adjPerm(j)
      var t = 0
      while t < nbrs.length do
        val i = nbrs(t)
        if i < j then
          var r = i
          while ancestor(r) != -1 && ancestor(r) != j do
            val next = ancestor(r)
            ancestor(r) = j
            r = next
          if ancestor(r) == -1 then
            ancestor(r) = j
            parent(r) = j
        t += 1
      j += 1
    parent

  def symbolicColumns(
      n: Int,
      adjPerm: Array[Array[Int]],
      parent: Array[Int],
      inputNnz: Int
  ): Either[LinAlgError, Array[Array[Int]]] =
    val cap = fillCap(inputNnz)
    val children = Array.fill(n)(new mutable.ArrayBuffer[Int](2))
    var k = 0
    while k < n do
      val p = parent(k)
      if p >= 0 then children(p) += k
      k += 1
    val cols = new Array[Array[Int]](n)
    val mark = Array.fill(n)(-1)
    var predicted = 0L
    var j = 0
    while j < n do
      val buf = new mutable.ArrayBuffer[Int](8)
      buf += j
      mark(j) = j
      val nbrs = adjPerm(j)
      var t = 0
      while t < nbrs.length do
        val i = nbrs(t)
        if i > j && mark(i) != j then
          mark(i) = j
          buf += i
        t += 1
      val ch = children(j)
      var c = 0
      while c < ch.length do
        val childCol = cols(ch(c))
        var p = 0
        while p < childCol.length do
          val i = childCol(p)
          if i > j && mark(i) != j then
            mark(i) = j
            buf += i
          p += 1
        c += 1
      val col = buf.toArray
      Sorting.quickSort(col)
      predicted += col.length.toLong
      if exceedsFill(inputNnz, predicted) then
        return Left(
          LinAlgError.InvalidArgument(
            s"sparse Cholesky fill estimate $predicted exceeds guard max(${inputNnz} * $FillFactor, $AbsoluteNnzCap) = $cap"
          )
        )
      cols(j) = col
      j += 1
    Right(cols)

  def packCsc(columns: Array[Array[Int]]): (Array[Int], Array[Int]) =
    val n = columns.length
    val colPtr = new Array[Int](n + 1)
    var j = 0
    while j < n do
      colPtr(j + 1) = colPtr(j) + columns(j).length
      j += 1
    val rowIdx = new Array[Int](colPtr(n))
    j = 0
    while j < n do
      val col = columns(j)
      var p = 0
      while p < col.length do
        rowIdx(colPtr(j) + p) = col(p)
        p += 1
      j += 1
    (colPtr, rowIdx)

  /** For each row `i`, the columns `k < i` that store `L(i,k)`, with CSC value offsets. */
  def rowHits(n: Int, colPtr: Array[Int], rowIdx: Array[Int]): (Array[Int], Array[Int], Array[Int]) =
    val counts = new Array[Int](n)
    var j = 0
    var p = 0
    while j < n do
      p = colPtr(j)
      val end = colPtr(j + 1)
      while p < end do
        val i = rowIdx(p)
        if i > j then counts(i) += 1
        p += 1
      j += 1
    val rowPtr = new Array[Int](n + 1)
    var i = 0
    while i < n do
      rowPtr(i + 1) = rowPtr(i) + counts(i)
      i += 1
    val hitCol = new Array[Int](rowPtr(n))
    val hitOff = new Array[Int](rowPtr(n))
    val next = rowPtr.clone()
    j = 0
    while j < n do
      p = colPtr(j)
      val end = colPtr(j + 1)
      while p < end do
        val row = rowIdx(p)
        if row > j then
          val dest = next(row)
          hitCol(dest) = j
          hitOff(dest) = p
          next(row) += 1
        p += 1
      j += 1
    (rowPtr, hitCol, hitOff)

  /** Lower-triangle symmetric value: original lower triangle only. */
  def lowerSym(a: CSR, row: Int, col: Int): Double =
    if row >= col then a(row, col) else a(col, row)

  def factorValues(
      a: CSR,
      perm: Array[Int],
      colPtr: Array[Int],
      rowIdx: Array[Int],
      rowPtr: Array[Int],
      hitCol: Array[Int],
      hitOff: Array[Int],
      work: Array[Double]
  ): Either[LinAlgError, Array[Double]] =
    val n = perm.length
    val values = new Array[Double](colPtr(n))
    var j = 0
    while j < n do
      val start = colPtr(j)
      val end = colPtr(j + 1)
      var p = start
      while p < end do
        work(rowIdx(p)) = 0.0
        p += 1
      val pj = perm(j)
      p = start
      while p < end do
        val i = rowIdx(p)
        work(i) = lowerSym(a, perm(i), pj)
        p += 1
      var h = rowPtr(j)
      val hEnd = rowPtr(j + 1)
      while h < hEnd do
        val k = hitCol(h)
        val ljk = values(hitOff(h))
        var q = colPtr(k)
        val qEnd = colPtr(k + 1)
        while q < qEnd do
          val i = rowIdx(q)
          if i >= j then work(i) -= values(q) * ljk
          q += 1
        h += 1
      val diag = work(j)
      if diag <= 0.0 || diag.isNaN then return Left(LinAlgError.NotPositiveDefinite(j))
      val ljj = math.sqrt(diag)
      p = start
      while p < end do
        val i = rowIdx(p)
        if i == j then values(p) = ljj
        else values(p) = work(i) / ljj
        work(i) = 0.0
        p += 1
      j += 1
    Right(values)

  def solveInPlace(
      colPtr: Array[Int],
      rowIdx: Array[Int],
      values: Array[Double],
      work: Array[Double]
  ): Unit =
    val n = colPtr.length - 1
    var j = 0
    while j < n do
      val diagOff = colPtr(j)
      work(j) = work(j) / values(diagOff)
      var p = diagOff + 1
      val end = colPtr(j + 1)
      val yj = work(j)
      while p < end do
        work(rowIdx(p)) -= values(p) * yj
        p += 1
      j += 1
    j = n - 1
    while j >= 0 do
      var sum = work(j)
      var p = colPtr(j) + 1
      val end = colPtr(j + 1)
      while p < end do
        sum -= values(p) * work(rowIdx(p))
        p += 1
      work(j) = sum / values(colPtr(j))
      j -= 1
