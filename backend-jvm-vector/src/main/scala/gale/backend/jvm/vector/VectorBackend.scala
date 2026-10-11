package gale.backend.jvm.vector

import gale.backend.Backend
import gale.backend.BackendConfig
import gale.backend.BackendThresholds
import gale.backend.Capability
import gale.backend.DenseDoubleKernel
import gale.backend.PureDenseDoubleKernel
import gale.platform.DoubleArray

import jdk.incubator.vector.DoubleVector
import jdk.incubator.vector.VectorSpecies

/** SIMD dense-`Double` kernels built on the JDK Vector API
  * ([[jdk.incubator.vector.DoubleVector]]). Contiguous row-major `gemv`, `gemm` and
  * `syrk` have explicit SIMD kernels; every other operation forwards to
  * [[gale.backend.PureDenseDoubleKernel]], the portable reference.
  *
  * GEMV has two SIMD layouts: a four-output row tile for row-major `A`, and an
  * eight-column axpy sweep for column-contiguous `A` (the transpose view of a
  * row-major matrix). GEMM packs `B` by columns, uses a 3×3 SIMD dot-product tile,
  * and handles only fully row-major inputs. Any other strided operand falls back
  * verbatim to the pure kernel — correctness over cleverness.
  *
  * `syrk` packs the columns of `A` and reuses the GEMM tile on the upper triangle.
  *
  * The row-major GEMV, the GEMM and the `syrk` reassociate (SIMD lane order, FMA),
  * so they are law-equivalent to the pure kernel within a small tolerance, NOT bit-identical.
  * The column-contiguous GEMV keeps the pure column-major kernel's per-element FMA
  * order and is bit-identical to it.
  */
object VectorDenseDoubleKernel extends DenseDoubleKernel:
  // The Vector API intrinsics need a compile-time-constant species, or every vector
  // is boxed. The JDK's static final species is one; a Scala object `val` (an
  // instance final) is not, hence an inline accessor. For the same reason a helper
  // reads `Species` itself rather than taking the species as a parameter: the GEMM
  // edge `vectorDot` once did, and whenever the JIT compiled it standalone it boxed
  // (about 10× the allocation and half the throughput of a 64³ product).
  private inline def Species: VectorSpecies[java.lang.Double] = DoubleVector.SPECIES_PREFERRED
  private[vector] def preferredLaneCount: Int = Species.length()

  // --- level 1: forward to the pure reference ---------------------------------------

  def dot(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Double =
    PureDenseDoubleKernel.dot(n, x, xOffset, xStride, y, yOffset, yStride)

  def nrm2(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    PureDenseDoubleKernel.nrm2(n, x, xOffset, xStride)

  def copy(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    PureDenseDoubleKernel.copy(n, x, xOffset, xStride, y, yOffset, yStride)

  def axpy(
      n: Int,
      alpha: Double,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    PureDenseDoubleKernel.axpy(n, alpha, x, xOffset, xStride, y, yOffset, yStride)

  def scal(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, xStride: Int): Unit =
    PureDenseDoubleKernel.scal(n, alpha, x, xOffset, xStride)

  def gemv(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      rowStride: Int,
      colStride: Int,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      beta: Double,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    if Species.length() >= 2 && colStride == 1 && xStride == 1 && yStride == 1 then
      gemvRowMajorSimd(
        rows,
        cols,
        alpha,
        DoubleArray.asArray(a),
        aOffset,
        rowStride,
        DoubleArray.asArray(x),
        xOffset,
        beta,
        DoubleArray.asArray(y),
        yOffset
      )
    else if Species.length() >= 2 && rowStride == 1 && yStride == 1 then
      gemvColMajorSimd(
        rows,
        cols,
        alpha,
        DoubleArray.asArray(a),
        aOffset,
        colStride,
        DoubleArray.asArray(x),
        xOffset,
        xStride,
        beta,
        DoubleArray.asArray(y),
        yOffset
      )
    else
      PureDenseDoubleKernel.gemv(
        rows, cols, alpha, a, aOffset, rowStride, colStride,
        x, xOffset, xStride, beta, y, yOffset, yStride
      )

  /** Four-output row-major matrix-vector tile. Each output row and `x` are
    * contiguous, so four independent vector accumulators hide the horizontal
    * reduction latency even on a two-double ARM species. This is the row-major
    * dual of a column-major BLAS transpose-GEMV kernel.
    */
  private def gemvRowMajorSimd(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      x: Array[Double],
      xOffset: Int,
      beta: Double,
      y: Array[Double],
      yOffset: Int
  ): Unit =
    val species = Species
    val lanes = species.length()
    val colBound = species.loopBound(cols)
    val pairBound = cols - cols % (2 * lanes)
    val rowBound = rows - (rows & 3)
    var i = 0
    while i < rowBound do
      val aRow0 = aOffset + i * aRowStride
      val aRow1 = aRow0 + aRowStride
      val aRow2 = aRow1 + aRowStride
      val aRow3 = aRow2 + aRowStride
      var s0 = DoubleVector.zero(species)
      var s1 = DoubleVector.zero(species)
      var s2 = DoubleVector.zero(species)
      var s3 = DoubleVector.zero(species)
      // Two accumulators per row: eight independent vector FMA chains.
      var t0 = DoubleVector.zero(species)
      var t1 = DoubleVector.zero(species)
      var t2 = DoubleVector.zero(species)
      var t3 = DoubleVector.zero(species)
      var j = 0
      while j < pairBound do
        val xv = DoubleVector.fromArray(species, x, xOffset + j)
        val xw = DoubleVector.fromArray(species, x, xOffset + j + lanes)
        s0 = DoubleVector.fromArray(species, a, aRow0 + j).fma(xv, s0)
        s1 = DoubleVector.fromArray(species, a, aRow1 + j).fma(xv, s1)
        s2 = DoubleVector.fromArray(species, a, aRow2 + j).fma(xv, s2)
        s3 = DoubleVector.fromArray(species, a, aRow3 + j).fma(xv, s3)
        t0 = DoubleVector.fromArray(species, a, aRow0 + j + lanes).fma(xw, t0)
        t1 = DoubleVector.fromArray(species, a, aRow1 + j + lanes).fma(xw, t1)
        t2 = DoubleVector.fromArray(species, a, aRow2 + j + lanes).fma(xw, t2)
        t3 = DoubleVector.fromArray(species, a, aRow3 + j + lanes).fma(xw, t3)
        j += 2 * lanes
      s0 = s0.add(t0)
      s1 = s1.add(t1)
      s2 = s2.add(t2)
      s3 = s3.add(t3)
      while j < colBound do
        val xv = DoubleVector.fromArray(species, x, xOffset + j)
        s0 = DoubleVector.fromArray(species, a, aRow0 + j).fma(xv, s0)
        s1 = DoubleVector.fromArray(species, a, aRow1 + j).fma(xv, s1)
        s2 = DoubleVector.fromArray(species, a, aRow2 + j).fma(xv, s2)
        s3 = DoubleVector.fromArray(species, a, aRow3 + j).fma(xv, s3)
        j += species.length()
      var sum0 = s0.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      var sum1 = s1.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      var sum2 = s2.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      var sum3 = s3.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      while j < cols do
        val xj = x(xOffset + j)
        sum0 = Math.fma(a(aRow0 + j), xj, sum0)
        sum1 = Math.fma(a(aRow1 + j), xj, sum1)
        sum2 = Math.fma(a(aRow2 + j), xj, sum2)
        sum3 = Math.fma(a(aRow3 + j), xj, sum3)
        j += 1
      y(yOffset + i) = combine(alpha, sum0, beta, y(yOffset + i))
      y(yOffset + i + 1) = combine(alpha, sum1, beta, y(yOffset + i + 1))
      y(yOffset + i + 2) = combine(alpha, sum2, beta, y(yOffset + i + 2))
      y(yOffset + i + 3) = combine(alpha, sum3, beta, y(yOffset + i + 3))
      i += 4

    while i < rows do
      val aRow = aOffset + i * aRowStride
      var sum = DoubleVector.zero(species)
      var j = 0
      while j < colBound do
        val av = DoubleVector.fromArray(species, a, aRow + j)
        val xv = DoubleVector.fromArray(species, x, xOffset + j)
        sum = av.fma(xv, sum)
        j += species.length()
      var scalar = sum.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      while j < cols do
        scalar = Math.fma(a(aRow + j), x(xOffset + j), scalar)
        j += 1
      y(yOffset + i) = combine(alpha, scalar, beta, y(yOffset + i))
      i += 1

  /** Column-contiguous (e.g. the transpose view of a row-major matrix) GEMV as a
    * sequence of SIMD axpys over `y`, eight (then four) columns per pass so each
    * `y` chunk is loaded and stored once per pass. Per element the FMA order is the
    * column order, matching the pure column-major kernel.
    */
  private def gemvColMajorSimd(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aColStride: Int,
      x: Array[Double],
      xOffset: Int,
      xStride: Int,
      beta: Double,
      y: Array[Double],
      yOffset: Int
  ): Unit =
    val species = Species
    val rowBound = species.loopBound(rows)
    var i = 0
    if beta == 0.0 then java.util.Arrays.fill(y, yOffset, yOffset + rows, 0.0)
    else if beta != 1.0 then
      while i < rows do
        y(yOffset + i) = beta * y(yOffset + i)
        i += 1

    // Eight columns per pass halves the `y` load/store traffic of the four-column
    // pass; a remaining block of four uses the four-column pass.
    val colBound8 = cols - (cols & 7)
    var j = 0
    while j < colBound8 do
      val s0 = alpha * x(xOffset + j * xStride)
      val s1 = alpha * x(xOffset + (j + 1) * xStride)
      val s2 = alpha * x(xOffset + (j + 2) * xStride)
      val s3 = alpha * x(xOffset + (j + 3) * xStride)
      val s4 = alpha * x(xOffset + (j + 4) * xStride)
      val s5 = alpha * x(xOffset + (j + 5) * xStride)
      val s6 = alpha * x(xOffset + (j + 6) * xStride)
      val s7 = alpha * x(xOffset + (j + 7) * xStride)
      val col0 = aOffset + j * aColStride
      val col1 = col0 + aColStride
      val col2 = col1 + aColStride
      val col3 = col2 + aColStride
      val col4 = col3 + aColStride
      val col5 = col4 + aColStride
      val col6 = col5 + aColStride
      val col7 = col6 + aColStride
      val v0 = DoubleVector.broadcast(species, s0)
      val v1 = DoubleVector.broadcast(species, s1)
      val v2 = DoubleVector.broadcast(species, s2)
      val v3 = DoubleVector.broadcast(species, s3)
      val v4 = DoubleVector.broadcast(species, s4)
      val v5 = DoubleVector.broadcast(species, s5)
      val v6 = DoubleVector.broadcast(species, s6)
      val v7 = DoubleVector.broadcast(species, s7)
      i = 0
      while i < rowBound do
        var yv = DoubleVector.fromArray(species, y, yOffset + i)
        yv = DoubleVector.fromArray(species, a, col0 + i).fma(v0, yv)
        yv = DoubleVector.fromArray(species, a, col1 + i).fma(v1, yv)
        yv = DoubleVector.fromArray(species, a, col2 + i).fma(v2, yv)
        yv = DoubleVector.fromArray(species, a, col3 + i).fma(v3, yv)
        yv = DoubleVector.fromArray(species, a, col4 + i).fma(v4, yv)
        yv = DoubleVector.fromArray(species, a, col5 + i).fma(v5, yv)
        yv = DoubleVector.fromArray(species, a, col6 + i).fma(v6, yv)
        yv = DoubleVector.fromArray(species, a, col7 + i).fma(v7, yv)
        yv.intoArray(y, yOffset + i)
        i += species.length()
      while i < rows do
        var yi = y(yOffset + i)
        yi = Math.fma(s0, a(col0 + i), yi)
        yi = Math.fma(s1, a(col1 + i), yi)
        yi = Math.fma(s2, a(col2 + i), yi)
        yi = Math.fma(s3, a(col3 + i), yi)
        yi = Math.fma(s4, a(col4 + i), yi)
        yi = Math.fma(s5, a(col5 + i), yi)
        yi = Math.fma(s6, a(col6 + i), yi)
        yi = Math.fma(s7, a(col7 + i), yi)
        y(yOffset + i) = yi
        i += 1
      j += 8

    val colBound = cols - (cols & 3)
    while j < colBound do
      val s0 = alpha * x(xOffset + j * xStride)
      val s1 = alpha * x(xOffset + (j + 1) * xStride)
      val s2 = alpha * x(xOffset + (j + 2) * xStride)
      val s3 = alpha * x(xOffset + (j + 3) * xStride)
      val col0 = aOffset + j * aColStride
      val col1 = col0 + aColStride
      val col2 = col1 + aColStride
      val col3 = col2 + aColStride
      val v0 = DoubleVector.broadcast(species, s0)
      val v1 = DoubleVector.broadcast(species, s1)
      val v2 = DoubleVector.broadcast(species, s2)
      val v3 = DoubleVector.broadcast(species, s3)
      i = 0
      while i < rowBound do
        var yv = DoubleVector.fromArray(species, y, yOffset + i)
        yv = DoubleVector.fromArray(species, a, col0 + i).fma(v0, yv)
        yv = DoubleVector.fromArray(species, a, col1 + i).fma(v1, yv)
        yv = DoubleVector.fromArray(species, a, col2 + i).fma(v2, yv)
        yv = DoubleVector.fromArray(species, a, col3 + i).fma(v3, yv)
        yv.intoArray(y, yOffset + i)
        i += species.length()
      while i < rows do
        var yi = y(yOffset + i)
        yi = Math.fma(s0, a(col0 + i), yi)
        yi = Math.fma(s1, a(col1 + i), yi)
        yi = Math.fma(s2, a(col2 + i), yi)
        yi = Math.fma(s3, a(col3 + i), yi)
        y(yOffset + i) = yi
        i += 1
      j += 4

    while j < cols do
      val s = alpha * x(xOffset + j * xStride)
      val col = aOffset + j * aColStride
      val v = DoubleVector.broadcast(species, s)
      i = 0
      while i < rowBound do
        DoubleVector
          .fromArray(species, a, col + i)
          .fma(v, DoubleVector.fromArray(species, y, yOffset + i))
          .intoArray(y, yOffset + i)
        i += species.length()
      while i < rows do
        y(yOffset + i) = Math.fma(s, a(col + i), y(yOffset + i))
        i += 1
      j += 1

  private inline def combine(alpha: Double, sum: Double, beta: Double, prior: Double): Double =
    if beta == 0.0 then alpha * sum
    else Math.fma(beta, prior, alpha * sum)

  /** `C := AᵀA` for a row-major `A`. Packs the columns of `A` contiguously (the
    * same transpose-pack as the GEMM's `B`), then runs the packed 3×3 dot tile over
    * the upper triangle only and mirrors it, so the result is exactly symmetric.
    * The shared (row) dimension is processed in [[SyrkChunk]]-row slices, each
    * packed into one reused panel; each later slice accumulates into `C`.
    * Reassociates like the GEMM.
    */
  def syrk(
      m: Int,
      k: Int,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    if syrkUsesSimd(m, k) then
      syrkPackedSimd(m, k, DoubleArray.asArray(a), aOffset, aRowStride, DoubleArray.asArray(c), cOffset, cRowStride)
    else PureDenseDoubleKernel.syrk(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)

  /** The `syrk`'s own crossover, inside the kernel because the shared gemm gate
    * (`m·k·k ≥ 64³`) does not fit it. The pure kernel already does half the work
    * of a GEMM, so the SIMD path needs twice the GEMM's element-work (`AᵀA` of
    * `64×64` stays pure, `128×64` goes SIMD). Each dot product runs over the `m`
    * rows of `A`, so short-wide `A` (`m < 4·lanes`, where a dot is at most a few
    * vectors or wholly scalar tail against nine horizontal reductions per tile)
    * also stays pure. Measured on ARM64 (two-double species), JDK 25.
    */
  private[vector] def syrkUsesSimd(m: Int, k: Int): Boolean =
    val lanes = Species.length()
    lanes >= 2 && m >= 4 * lanes && m.toLong * k.toLong * k.toLong >= SyrkMinWork

  private final val SyrkMinWork = 2L * 64L * 64L * 64L

  /** Shared-dimension slice length for the SIMD `syrk`. A slice of `A` is packed
    * into a reused contiguous `k × 256` panel (column `j` at `j·256`), so transient
    * memory is bounded by `256·k` doubles whatever `m` is; the six live packed rows
    * of a tile (12 KB) stay in L1, and for `k ≤ 256` the whole panel (≤ 512 KB) fits
    * a 1 MB L2.
    */
  private final val SyrkChunk = 256

  private[vector] def syrkPackedSimd(
      m: Int,
      k: Int,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      c: Array[Double],
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    val panel = new Array[Double](k * math.min(SyrkChunk, m))
    // At least one slice, so m = 0 still assigns the zero Gram matrix.
    var start = 0
    while start == 0 || start < m do
      val length = math.min(SyrkChunk, m - start)
      packColumnsInto(length, k, a, aOffset + start * aRowStride, aRowStride, panel)
      packedDotTiles(
        k, k, length, 1.0,
        panel, 0, length,
        panel, 0, length,
        if start == 0 then 0.0 else 1.0, c, cOffset, cRowStride,
        upperOnly = true
      )
      start += math.max(length, 1)
    var i = 1
    while i < k do
      val cRowI = cOffset + i * cRowStride
      var j = 0
      while j < i do
        c(cRowI + j) = c(cOffset + j * cRowStride + i)
        j += 1
      i += 1

  // --- level 3: the accelerated lever -----------------------------------------------

  def gemm(
      rows: Int,
      cols: Int,
      shared: Int,
      alpha: Double,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      aColStride: Int,
      b: DoubleArray,
      bOffset: Int,
      bRowStride: Int,
      bColStride: Int,
      beta: Double,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int,
      cColStride: Int
  ): Unit =
    // SIMD only for the fully row-major layout (unit column stride everywhere), where a
    // row of B and a row of C are contiguous. Any strided/transposed operand routes to
    // the pure kernel unchanged — it already honours arbitrary strides correctly.
    if Species.length() >= 2 && aColStride == 1 && bColStride == 1 && cColStride == 1 then
      gemmRowMajorDotSimd(
        rows, cols, shared, alpha,
        DoubleArray.asArray(a), aOffset, aRowStride,
        DoubleArray.asArray(b), bOffset, bRowStride,
        beta, DoubleArray.asArray(c), cOffset, cRowStride
      )
    else
      PureDenseDoubleKernel.gemm(
        rows, cols, shared, alpha,
        a, aOffset, aRowStride, aColStride,
        b, bOffset, bRowStride, bColStride,
        beta, c, cOffset, cRowStride, cColStride
      )

  /** Pack `B` by columns, then compute 3×3 output tiles as nine independent SIMD
    * dot products over the shared dimension. The old across-column kernel pays too
    * much broadcast/load overhead on a two-double ARM species; this shape instead
    * keeps nine accumulators live and hides horizontal-reduction latency. Packing
    * is included in the public-facade benchmark and is amortized only above the
    * measured threshold.
    */
  private def gemmRowMajorDotSimd(
      rows: Int,
      cols: Int,
      shared: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      b: Array[Double],
      bOffset: Int,
      bRowStride: Int,
      beta: Double,
      c: Array[Double],
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    val packedB = new Array[Double](cols * shared)
    packColumnsInto(shared, cols, b, bOffset, bRowStride, packedB)
    packedDotTiles(
      rows, cols, shared, alpha,
      a, aOffset, aRowStride,
      packedB, 0, shared,
      beta, c, cOffset, cRowStride,
      upperOnly = false
    )

  /** Copy the columns of the row-major `rows × cols` matrix at `b` into `packed`,
    * column `j` contiguous at `j * rows`, in 32×32 blocks.
    */
  private def packColumnsInto(
      rows: Int,
      cols: Int,
      b: Array[Double],
      bOffset: Int,
      bRowStride: Int,
      packed: Array[Double]
  ): Unit =
    val PackBlock = 32
    var kb = 0
    while kb < rows do
      val kend = math.min(kb + PackBlock, rows)
      var jb = 0
      while jb < cols do
        val jend = math.min(jb + PackBlock, cols)
        var k = kb
        while k < kend do
          val bRow = bOffset + k * bRowStride
          var j = jb
          while j < jend do
            packed(j * rows + k) = b(bRow + j)
            j += 1
          k += 1
        jb += PackBlock
      kb += PackBlock

  /** `C(i, j) := alpha · ⟨a_i, p_j⟩ + beta · C(i, j)` over `length`-long contiguous
    * rows `a_i` (at `aOffset + i·aRowStride`) and `p_j` (at `pOffset + j·pRowStride`),
    * as 3×3 tiles of nine independent SIMD dot products. With `upperOnly`, tiles and
    * edge cells strictly below the diagonal are skipped (diagonal tiles are written
    * whole); the caller owns the lower triangle.
    */
  private def packedDotTiles(
      rows: Int,
      cols: Int,
      length: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      p: Array[Double],
      pOffset: Int,
      pRowStride: Int,
      beta: Double,
      c: Array[Double],
      cOffset: Int,
      cRowStride: Int,
      upperOnly: Boolean
  ): Unit =
    val species = Species
    val sharedBound = species.loopBound(length)
    val rowBound = rows - rows % 3
    val colBound = cols - cols % 3
    var i = 0
    while i < rowBound do
      val a0 = aOffset + i * aRowStride
      val a1 = a0 + aRowStride
      val a2 = a1 + aRowStride
      var j = if upperOnly then i else 0
      while j < colBound do
        val b0 = pOffset + j * pRowStride
        val b1 = b0 + pRowStride
        val b2 = b1 + pRowStride
        var s00 = DoubleVector.zero(species)
        var s01 = DoubleVector.zero(species)
        var s02 = DoubleVector.zero(species)
        var s10 = DoubleVector.zero(species)
        var s11 = DoubleVector.zero(species)
        var s12 = DoubleVector.zero(species)
        var s20 = DoubleVector.zero(species)
        var s21 = DoubleVector.zero(species)
        var s22 = DoubleVector.zero(species)
        var k = 0
        while k < sharedBound do
          val av0 = DoubleVector.fromArray(species, a, a0 + k)
          val av1 = DoubleVector.fromArray(species, a, a1 + k)
          val av2 = DoubleVector.fromArray(species, a, a2 + k)
          val bv0 = DoubleVector.fromArray(species, p, b0 + k)
          val bv1 = DoubleVector.fromArray(species, p, b1 + k)
          val bv2 = DoubleVector.fromArray(species, p, b2 + k)
          s00 = av0.fma(bv0, s00)
          s01 = av0.fma(bv1, s01)
          s02 = av0.fma(bv2, s02)
          s10 = av1.fma(bv0, s10)
          s11 = av1.fma(bv1, s11)
          s12 = av1.fma(bv2, s12)
          s20 = av2.fma(bv0, s20)
          s21 = av2.fma(bv1, s21)
          s22 = av2.fma(bv2, s22)
          k += species.length()
        var x00 = s00.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x01 = s01.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x02 = s02.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x10 = s10.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x11 = s11.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x12 = s12.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x20 = s20.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x21 = s21.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        var x22 = s22.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
        while k < length do
          val av0 = a(a0 + k)
          val av1 = a(a1 + k)
          val av2 = a(a2 + k)
          val bv0 = p(b0 + k)
          val bv1 = p(b1 + k)
          val bv2 = p(b2 + k)
          x00 = Math.fma(av0, bv0, x00)
          x01 = Math.fma(av0, bv1, x01)
          x02 = Math.fma(av0, bv2, x02)
          x10 = Math.fma(av1, bv0, x10)
          x11 = Math.fma(av1, bv1, x11)
          x12 = Math.fma(av1, bv2, x12)
          x20 = Math.fma(av2, bv0, x20)
          x21 = Math.fma(av2, bv1, x21)
          x22 = Math.fma(av2, bv2, x22)
          k += 1
        val c0 = cOffset + i * cRowStride + j
        val c1 = c0 + cRowStride
        val c2 = c1 + cRowStride
        c(c0) = combine(alpha, x00, beta, c(c0))
        c(c0 + 1) = combine(alpha, x01, beta, c(c0 + 1))
        c(c0 + 2) = combine(alpha, x02, beta, c(c0 + 2))
        c(c1) = combine(alpha, x10, beta, c(c1))
        c(c1 + 1) = combine(alpha, x11, beta, c(c1 + 1))
        c(c1 + 2) = combine(alpha, x12, beta, c(c1 + 2))
        c(c2) = combine(alpha, x20, beta, c(c2))
        c(c2 + 1) = combine(alpha, x21, beta, c(c2 + 1))
        c(c2 + 2) = combine(alpha, x22, beta, c(c2 + 2))
        j += 3
      i += 3

    // Edge rows/columns are a small fraction of large products. Keep them on a
    // vector dot against the already-packed B rather than complicating the 3×3 tile.
    i = 0
    while i < rows do
      val tileEnd = if i < rowBound then colBound else 0
      val jStart = if upperOnly then math.max(tileEnd, i) else tileEnd
      val aRow = aOffset + i * aRowStride
      var j = jStart
      while j < cols do
        val sum = vectorDot(length, a, aRow, p, pOffset + j * pRowStride)
        val cij = cOffset + i * cRowStride + j
        c(cij) = combine(alpha, sum, beta, c(cij))
        j += 1
      i += 1

  private def vectorDot(
      length: Int,
      x: Array[Double],
      xOffset: Int,
      y: Array[Double],
      yOffset: Int
  ): Double =
    val species = Species
    val bound = species.loopBound(length)
    var sum = DoubleVector.zero(species)
    var k = 0
    while k < bound do
      val xv = DoubleVector.fromArray(species, x, xOffset + k)
      val yv = DoubleVector.fromArray(species, y, yOffset + k)
      sum = xv.fma(yv, sum)
      k += species.length()
    var scalar = sum.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
    while k < length do
      scalar = Math.fma(x(xOffset + k), y(yOffset + k), scalar)
      k += 1
    scalar

end VectorDenseDoubleKernel

/** Measured Vector-API crossovers. On ARM64 (two-double species, JDK 25) the packed
  * GEMM and `syrk` beat the pure kernels from 64³ element-work (`64×64·64×64`,
  * `256×64·64×64`, `AᵀA` of `256×64`); the four-row GEMV wins from 128 square.
  * `syrk` additionally applies its own shape guard (`VectorDenseDoubleKernel.syrkUsesSimd`).
  * Narrower-than-two-double runtimes retain the pure implementation. Factorizations
  * themselves stay pure, but blocked QR's trailing updates go through the gemm
  * seam (`Factorizations.dispatchGemm`), so above this threshold they run on the
  * SIMD GEMM.
  */
object VectorThresholds extends BackendThresholds:
  // The preferred species is fixed for the life of the JVM, so the crossovers are
  // values, not per-call re-evaluations on the hot dispatch path.
  private val simdCapable: Boolean = VectorDenseDoubleKernel.preferredLaneCount >= 2
  val nativeGemmMinFlops: Long = if simdCapable then 64L * 64L * 64L else Long.MaxValue
  val nativeGemvMinWork: Long = if simdCapable then 128L * 128L else Long.MaxValue
  val nativeFactorizationMinSize: Int = Int.MaxValue

/** The JVM-only SIMD acceleration backend. Advertises [[Capability.Vectorized]], so the
  * coarse gemm seam (`Backend.acceleratesGemm`) routes level-3 products to
  * [[VectorDenseDoubleKernel]] above [[VectorThresholds.nativeGemmMinFlops]]. On a
  * runtime whose preferred species is narrower than two doubles, that threshold is
  * infinite and the backend deliberately retains Gale's faster pure kernel. Bring it
  * into scope with `import gale.backend.jvm.vector.given` (or `using VectorBackend`).
  */
object VectorBackend extends Backend:
  val name: String = "jvm-vector"
  val capabilities: Set[Capability] = Set(Capability.Vectorized)
  val denseDouble: DenseDoubleKernel = VectorDenseDoubleKernel
  val thresholds: BackendThresholds = VectorThresholds
  val config: BackendConfig = BackendConfig.singleThreaded

/** The `given` a user brings into scope with `import gale.backend.jvm.vector.given`.
  * Top-level in the package (not nested in [[VectorBackend]]) so that wildcard-given
  * import resolves it; being in lexical scope it outranks the always-present companion
  * fallback `Backend.pure`. The module also exposes the plain value [[VectorBackend]], so
  * a multi-backend user can `Backend.compose` the values and declare one composite given.
  */
given vectorBackend: Backend = VectorBackend
