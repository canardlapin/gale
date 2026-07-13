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
  * ([[jdk.incubator.vector.DoubleVector]]). Only `gemm` is accelerated — the design's
  * lever is the level-3 product (doc §A.2); everything else forwards to
  * [[gale.backend.PureDenseDoubleKernel]], the portable reference.
  *
  * The SIMD path handles ONLY the fully row-major case (unit column stride on all
  * three operands): it streams whole rows of `B` and `C` with an i-k-j loop, issuing
  * one broadcast-`A` × vector-`B` fused-multiply-add per lane group. Any strided or
  * transposed operand falls back verbatim to the pure kernel — correctness over
  * cleverness. Reassociation (SIMD lane order, FMA) makes the result law-equivalent to
  * the pure kernel within a small tolerance, NOT bit-identical.
  */
object VectorDenseDoubleKernel extends DenseDoubleKernel:
  private final val Species: VectorSpecies[java.lang.Double] = DoubleVector.SPECIES_PREFERRED

  // --- level 1 / level 2: forward to the pure reference (not the lever) -------------

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
    PureDenseDoubleKernel.gemv(rows, cols, alpha, a, aOffset, rowStride, colStride, x, xOffset, xStride, beta, y, yOffset, yStride)

  /** `syrk` isn't routed by the coarse gemm seam yet (the `AᵀA` fast-path stays on the
    * dedicated pure symmetric kernel), so delegate to the pure reference rather than
    * carry a second SIMD kernel that nothing exercises.
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
    if aColStride == 1 && bColStride == 1 && cColStride == 1 then
      gemmRowMajorSimd(
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

  /** Row-major `C := beta·C + alpha·A·B` over raw arrays. Scale `C` by `beta` once,
    * then an i-k-j accumulation: each k-step broadcasts `alpha·A[i,k]` and fuses it
    * against a contiguous lane group of `B`'s k-th row into `C`'s i-th row.
    */
  private def gemmRowMajorSimd(
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
    scaleRowMajor(rows, cols, beta, c, cOffset, cRowStride)
    val species = Species
    val bound = species.loopBound(cols) // largest multiple of the lane count ≤ cols
    val laneCount = species.length()
    var i = 0
    while i < rows do
      val cRow = cOffset + i * cRowStride
      val aRow = aOffset + i * aRowStride
      var k = 0
      while k < shared do
        val aik = alpha * a(aRow + k)
        val aVec = DoubleVector.broadcast(species, aik)
        val bRow = bOffset + k * bRowStride
        var j = 0
        while j < bound do
          val bVec = DoubleVector.fromArray(species, b, bRow + j)
          val cVec = DoubleVector.fromArray(species, c, cRow + j)
          aVec.fma(bVec, cVec).intoArray(c, cRow + j) // aik·B[k,j:] + C[i,j:]
          j += laneCount
        while j < cols do
          c(cRow + j) = Math.fma(aik, b(bRow + j), c(cRow + j))
          j += 1
        k += 1
      i += 1

  /** Scale a row-major `C` in place by `beta` (zeroing when `beta == 0`, a no-op when
    * `beta == 1`) before the accumulation — mirrors the pure kernel's `scaleRowMajor`.
    */
  private def scaleRowMajor(
      rows: Int,
      cols: Int,
      beta: Double,
      c: Array[Double],
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    if beta == 1.0 then ()
    else
      val zero = beta == 0.0
      var i = 0
      while i < rows do
        val cRow = cOffset + i * cRowStride
        var j = 0
        while j < cols do
          c(cRow + j) = if zero then 0.0 else beta * c(cRow + j)
          j += 1
        i += 1
end VectorDenseDoubleKernel

/** Crossovers for the Vector-API backend. `nativeGemmMinFlops` is a MODEST placeholder
  * (`32·32·32 = 32768` on the element-count scale the gemm seam compares against, i.e.
  * `rows·cols·shared`) so mid-size products route to SIMD; the real value is set by the
  * later measurement bead. `gemv`/factorization are not routed here, hence effectively
  * infinite.
  */
object VectorThresholds extends BackendThresholds:
  // Placeholder pending the measurement bead: routes ~32×32×32 and larger to SIMD.
  def nativeGemmMinFlops: Long = 32L * 32L * 32L
  def nativeGemvMinWork: Long = Long.MaxValue
  def nativeFactorizationMinSize: Int = Int.MaxValue

/** The JVM-only SIMD acceleration backend. Advertises [[Capability.Vectorized]], so the
  * coarse gemm seam (`Backend.acceleratesGemm`) routes level-3 products to
  * [[VectorDenseDoubleKernel]] above [[VectorThresholds.nativeGemmMinFlops]]. Bring it
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
