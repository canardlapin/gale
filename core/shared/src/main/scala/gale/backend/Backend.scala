package gale.backend

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.spectral.SpectralBackend

/** A backend '''acceleration class''' — the yes/no question "does this backend make
  * an operation fast?", carried as a `Set[Capability]` on every [[Backend]]. It is
  * the dense-kernel counterpart of [[gale.spectral.SpectralCapability]] (which
  * classifies spectral operations); the two stay distinct enums but read as one
  * family. See `docs/backend-architecture.md` §A.3.
  */
enum Capability:
  case NativeBlas, NativeLapack, NativeSparse   // native providers (FFM)
  case Vectorized                               // JDK Vector API SIMD kernels (pure JVM)
  case Multithreaded                            // uses a compute pool (§B)
  case WasmFriendly, Deterministic              // pure-path properties

/** The accelerable dense-`Double` kernel surface (BLAS levels 1–3), declared with the
  * exact strided signatures of [[gale.kernel.DoubleKernels]] so an implementation is a
  * drop-in for the pure kernels. The FULL surface is declared so a native-BLAS backend
  * can implement all of it, but the '''dispatch policy''' (the facade + [[BackendThresholds]])
  * decides what is actually ROUTED: only the coarse ops — `gemm`, `syrk`, and standalone
  * `gemv` — dispatch here; loop-called L1/L2 (`dot`/`axpy`/… as inner-loop steps) never
  * route through a backend (doc §A.2, §A.2.1). All methods assume validated dimensions
  * (unsafe). `trsm` and the whole-op factorization hooks belong to a future LAPACK backend
  * and are intentionally not part of this first surface.
  */
trait DenseDoubleKernel:
  def dot(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Double
  def nrm2(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double
  def copy(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit
  def axpy(
      n: Int,
      alpha: Double,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit
  def scal(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, xStride: Int): Unit
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
  ): Unit
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
  ): Unit
  def syrk(
      m: Int,
      k: Int,
      a: DoubleArray,
      aOffset: Int,
      aRowStride: Int,
      c: DoubleArray,
      cOffset: Int,
      cRowStride: Int
  ): Unit

/** The always-present kernel set: forwards verbatim to [[gale.kernel.DoubleKernels]], so
  * it is the pure, portable, deterministic reference — and the only kernel set on JS.
  */
object PureDenseDoubleKernel extends DenseDoubleKernel:
  def dot(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Double =
    DoubleKernels.ddot(n, x, xOffset, xStride, y, yOffset, yStride)

  def nrm2(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    DoubleKernels.dnrm2(n, x, xOffset, xStride)

  def copy(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    DoubleKernels.dcopy(n, x, xOffset, xStride, y, yOffset, yStride)

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
    DoubleKernels.daxpy(n, alpha, x, xOffset, xStride, y, yOffset, yStride)

  def scal(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, xStride: Int): Unit =
    DoubleKernels.dscal(n, alpha, x, xOffset, xStride)

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
    DoubleKernels.dgemv(rows, cols, alpha, a, aOffset, rowStride, colStride, x, xOffset, xStride, beta, y, yOffset, yStride)

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
    DoubleKernels.dgemm(
      rows, cols, shared, alpha,
      a, aOffset, aRowStride, aColStride,
      b, bOffset, bRowStride, bColStride,
      beta, c, cOffset, cRowStride, cColStride
    )

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
    DoubleKernels.dsyrkRowMajor(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)

/** Measured, per-backend-family crossovers (doc §D.4). Above a threshold the coarse op
  * routes to the backend; below it (and on JS always) the pure kernel runs.
  */
trait BackendThresholds:
  def nativeGemmMinFlops: Long           // the primary crossover; heap vs NativeDMat variants later
  def nativeGemvMinWork: Long            // standalone `A*x` only
  def nativeFactorizationMinSize: Int    // per-routine (LU/Cholesky/QR differ)

/** The pure path is always chosen, so every native crossover is effectively infinite. */
object PureThresholds extends BackendThresholds:
  def nativeGemmMinFlops: Long = Long.MaxValue
  def nativeGemvMinWork: Long = Long.MaxValue
  def nativeFactorizationMinSize: Int = Int.MaxValue

/** Backend threading configuration, read once at construction, never per call (doc §B).
  * `jvmThreads` sizes gale's one fixed compute pool (1 = single-threaded); `nativeThreads`
  * is pinned on the native library so a JVM pool and a native pool never oversubscribe.
  */
final case class BackendConfig(jvmThreads: Int, nativeThreads: Int, allowNestedParallelism: Boolean = false)
object BackendConfig:
  /** The default: gale adds no parallelism, so it is safe embedded in a parallel host. */
  val singleThreaded: BackendConfig = BackendConfig(1, 1)
  /** Standalone batch use: gale owns `cores` threads. */
  def dedicated(cores: Int): BackendConfig = BackendConfig(cores, cores)

/** A dense/kernel acceleration backend, selected by an explicit `given` import
  * (PRD Doctrine 8). The always-present, lowest-priority [[Backend.pure]] computes with
  * the portable kernels; an acceleration module provides its own `given Backend`.
  * See `docs/backend-architecture.md` §A.
  */
trait Backend:
  def name: String
  def capabilities: Set[Capability]
  def denseDouble: DenseDoubleKernel
  def thresholds: BackendThresholds
  def config: BackendConfig
  /** The spectral half, if this backend also provides one (the bridge, §A.5). A backend
    * that advertises [[Capability.NativeLapack]] MUST provide it (A-R3).
    */
  def spectral: Option[SpectralBackend] = None

  /** Does this backend provide an accelerated `gemm`/`syrk`? The coarse gemm seam routes
    * on THIS, not on [[Capability.NativeBlas]] alone, so a `Vectorized` backend triggers
    * the accelerated path.
    */
  final def acceleratesGemm: Boolean =
    capabilities.contains(Capability.NativeBlas) || capabilities.contains(Capability.Vectorized)

/** The always-present pure backend: portable, deterministic, no acceleration. With no
  * acceleration import this is the only `given Backend`, so every coarse facade call takes
  * the pure branch — today's behaviour, byte-identical.
  */
object PureBackend extends Backend:
  val name: String = "pure"
  val capabilities: Set[Capability] = Set(Capability.WasmFriendly, Capability.Deterministic)
  val denseDouble: DenseDoubleKernel = PureDenseDoubleKernel
  val thresholds: BackendThresholds = PureThresholds
  val config: BackendConfig = BackendConfig.singleThreaded

object Backend:
  /** The always-in-scope, lowest-priority fallback (mirrors `SpectralBackend.none`, but
    * `pure` COMPUTES rather than returning `Left(Unsupported)`).
    */
  given pure: Backend = PureBackend

  /** Union the parts' capabilities and route the coarse dense-kernel surface to the first
    * part that accelerates `gemm` (earlier parts win), falling back to `primary`. Because
    * loop-called L1/L2 never routes at the facade, this whole-surface choice only affects
    * the coarse ops that do. A single-part compose is that part. Mirrors
    * `SpectralBackend.compose`; two `given Backend`s in scope are a compile-time ambiguity,
    * so a multi-backend user imports the VALUES and declares one composite given.
    */
  def compose(primary: Backend, rest: Backend*): Backend =
    if rest.isEmpty then primary
    else
      val parts: Seq[Backend] = primary +: rest
      val gemmProvider: Backend = parts.find(_.acceleratesGemm).getOrElse(primary)
      new Backend:
        val name: String = parts.map(_.name).mkString("compose(", " + ", ")")
        val capabilities: Set[Capability] = parts.foldLeft(Set.empty[Capability])(_ ++ _.capabilities)
        val denseDouble: DenseDoubleKernel = gemmProvider.denseDouble
        val thresholds: BackendThresholds = gemmProvider.thresholds
        // Threading policy follows the PRIMARY deliberately: a composite is "primary,
        // augmented with the rest's capabilities/kernels", and §B makes threading a
        // process-level, construction-time concern the primary declares — not a
        // per-op property that should flip with which part serves gemm.
        val config: BackendConfig = primary.config
        override val spectral: Option[SpectralBackend] =
          val specs = parts.flatMap(_.spectral)
          if specs.isEmpty then None
          else Some(SpectralBackend.compose(specs.head, specs.drop(1)*))
