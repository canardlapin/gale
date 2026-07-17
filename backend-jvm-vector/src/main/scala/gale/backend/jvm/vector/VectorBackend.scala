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
  * ([[jdk.incubator.vector.DoubleVector]]). Contiguous row-major `gemv` and `gemm`
  * have explicit SIMD kernels; every other operation forwards to
  * [[gale.backend.PureDenseDoubleKernel]], the portable reference.
  *
  * The SIMD path handles only fully row-major inputs. GEMV uses a four-output
  * row tile; GEMM packs `B` by columns and uses a 3×3 SIMD dot-product tile.
  * Any strided or transposed operand falls back verbatim to the pure kernel — correctness over
  * cleverness. Reassociation (SIMD lane order, FMA) makes the result law-equivalent to
  * the pure kernel within a small tolerance, NOT bit-identical.
  */
object VectorDenseDoubleKernel extends DenseDoubleKernel:
  private final val Species: VectorSpecies[java.lang.Double] = DoubleVector.SPECIES_PREFERRED
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
    val colBound = species.loopBound(cols)
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
      var j = 0
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

  private inline def combine(alpha: Double, sum: Double, beta: Double, prior: Double): Double =
    if beta == 0.0 then alpha * sum
    else Math.fma(beta, prior, alpha * sum)

  /** The coarse `AᵀA` seam can route here, but this first Vector backend has no
    * independently faster symmetric kernel. Preserve the half-work pure routine.
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
    PureDenseDoubleKernel.syrk(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)

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
    val PackBlock = 32
    var kb = 0
    while kb < shared do
      val kend = math.min(kb + PackBlock, shared)
      var jb = 0
      while jb < cols do
        val jend = math.min(jb + PackBlock, cols)
        var k = kb
        while k < kend do
          val bRow = bOffset + k * bRowStride
          var j = jb
          while j < jend do
            packedB(j * shared + k) = b(bRow + j)
            j += 1
          k += 1
        jb += PackBlock
      kb += PackBlock

    val species = Species
    val sharedBound = species.loopBound(shared)
    val rowBound = rows - rows % 3
    val colBound = cols - cols % 3
    var i = 0
    while i < rowBound do
      val a0 = aOffset + i * aRowStride
      val a1 = a0 + aRowStride
      val a2 = a1 + aRowStride
      var j = 0
      while j < colBound do
        val b0 = j * shared
        val b1 = b0 + shared
        val b2 = b1 + shared
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
          val bv0 = DoubleVector.fromArray(species, packedB, b0 + k)
          val bv1 = DoubleVector.fromArray(species, packedB, b1 + k)
          val bv2 = DoubleVector.fromArray(species, packedB, b2 + k)
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
        while k < shared do
          val av0 = a(a0 + k)
          val av1 = a(a1 + k)
          val av2 = a(a2 + k)
          val bv0 = packedB(b0 + k)
          val bv1 = packedB(b1 + k)
          val bv2 = packedB(b2 + k)
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
      val jStart = if i < rowBound then colBound else 0
      val aRow = aOffset + i * aRowStride
      var j = jStart
      while j < cols do
        val sum = vectorDot(species, shared, a, aRow, packedB, j * shared)
        val cij = cOffset + i * cRowStride + j
        c(cij) = combine(alpha, sum, beta, c(cij))
        j += 1
      i += 1

  private def vectorDot(
      species: VectorSpecies[java.lang.Double],
      length: Int,
      x: Array[Double],
      xOffset: Int,
      y: Array[Double],
      yOffset: Int
  ): Double =
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

/** Measured Vector-API crossovers. On ARM64/JDK 22 the packed GEMM and four-row
  * GEMV kernels both win from 128 square; narrower-than-two-double runtimes retain
  * the pure implementation. Factorizations remain on the pure backend.
  */
object VectorThresholds extends BackendThresholds:
  // The preferred species is fixed for the life of the JVM, so the crossovers are
  // values, not per-call re-evaluations on the hot dispatch path.
  private val simdCapable: Boolean = VectorDenseDoubleKernel.preferredLaneCount >= 2
  val nativeGemmMinFlops: Long = if simdCapable then 128L * 128L * 128L else Long.MaxValue
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
