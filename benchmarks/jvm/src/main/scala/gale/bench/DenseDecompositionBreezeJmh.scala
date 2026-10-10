package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.det
import breeze.linalg.eig
import breeze.linalg.inv
import breeze.linalg.pinv
import breeze.linalg.svd
import gale.backend.Backend
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import gale.spectral.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams
import org.openjdk.jmh.infra.Blackhole

/** Dense inverse, determinant, pseudo-inverse, full SVD and nonsymmetric eigen,
  * gale vs Breeze (`n` in {16, 64, 256}):
  *
  *   - `inv`: gale `A.inverse` (LU, then `U⁻¹ L⁻¹ P` formed from the factors) vs
  *     breeze `inv(A)` (`dgetrf` + `dgetri`).
  *   - `det`: gale `A.det` (LU) vs breeze `det(A)` (LU).
  *   - `pinv`: gale `A.pinv` (economy SVD, cutoff `max(m,n)·ε·σ_max`) vs breeze
  *     `pinv(A)` (also SVD-based; Breeze's cutoff differs, which does not change the
  *     work on the full-rank input used here).
  *   - `svd`: gale `A.svd` (economy, values and both factors) vs breeze `svd(A)`
  *     (`dgesdd`, full factors). The input is square, so economy and full shapes agree.
  *   - `eig`: gale `Eigen.eigNonsymmetric(A, All, Right)` vs breeze `eig(A)` (`dgeev`
  *     with right vectors) on a general nonsymmetric matrix.
  *
  * `inv` and `det` take the gale backend (they factor through LU, which reaches the
  * backend only where it routes a product through gemm); the spectral ops resolve a
  * `SpectralBackend` the Vector backend does not supply, so they are backend-insensitive.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class DenseDecompositionBreezeJmh:
  @Param(Array("16", "64", "256"))
  var n: Int = 0

  private var gA: DMat        = uninitialized
  private var gGeneral: DMat  = uninitialized
  private var bA: BDM[Double] = uninitialized
  private var bGeneral: BDM[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val aData       = diagonallyDominant(n, 2100L)
    val generalData = matrixData(n, n, 2200L)
    gA = galeMatrix(aData)
    gGeneral = galeMatrix(generalData)
    bA = breezeMatrix(aData)
    bGeneral = breezeMatrix(generalData)

  @Benchmark def galeInv(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.inverse)
  @Benchmark def breezeInv(bh: Blackhole): Unit = bh.consume(inv(bA))

  @Benchmark def galeDet(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.det)
  @Benchmark def breezeDet(): Double = det(bA)

  @Benchmark def galePinv(bh: Blackhole): Unit   = bh.consume(gGeneral.pinv)
  @Benchmark def breezePinv(bh: Blackhole): Unit = bh.consume(pinv(bGeneral))

  @Benchmark def galeSvd(bh: Blackhole): Unit   = bh.consume(gGeneral.svd)
  @Benchmark def breezeSvd(bh: Blackhole): Unit = bh.consume(svd(bGeneral))

  @Benchmark def galeEig(bh: Blackhole): Unit =
    bh.consume(Eigen.eigNonsymmetric(gGeneral, EigenSelection.All, EigenVectors.Right))
  @Benchmark def breezeEig(bh: Blackhole): Unit = bh.consume(eig(bGeneral))

/** Small fixed-size dense work (`n` in {3, 4}), gale vs Breeze — report-only: per-call
  * overhead dominates, so the plan tracks these but sets no pass/fail target.
  *
  *   - `gemm`: `A * B` both ways (gale through the selected backend).
  *   - `solve`: gale `A.solve(b)` vs breeze `A \ b`.
  *   - `inv`: gale `A.inverse` vs breeze `inv(A)`.
  *   - `det`: gale `A.det` (LU) vs breeze `det(A)` (LU).
  *
  * gale's `Mat3`/`Mat4` value types (`gale.linalg.Tiny`) only offer matrix-vector
  * products and `det`, so both sides use the general dense types here.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SmallDenseBreezeJmh:
  @Param(Array("3", "4"))
  var n: Int = 0

  private var gA: DMat        = uninitialized
  private var gB: DMat        = uninitialized
  private var gb: DVec        = uninitialized
  private var bA: BDM[Double] = uninitialized
  private var bB: BDM[Double] = uninitialized
  private var bb: BDV[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val aData = diagonallyDominant(n, 2300L)
    val bData = matrixData(n, n, 2400L)
    val vData = vectorData(n, 2500L)
    gA = galeMatrix(aData)
    gB = galeMatrix(bData)
    gb = galeVector(vData)
    bA = breezeMatrix(aData)
    bB = breezeMatrix(bData)
    bb = breezeVector(vData)

  @Benchmark def galeGemm(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA * gB)
  @Benchmark def breezeGemm(bh: Blackhole): Unit = bh.consume(bA * bB)

  @Benchmark def galeSolve(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.solve(gb))
  @Benchmark def breezeSolve(bh: Blackhole): Unit = bh.consume(bA \ bb)

  @Benchmark def galeInv(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.inverse)
  @Benchmark def breezeInv(bh: Blackhole): Unit = bh.consume(inv(bA))

  @Benchmark def galeDet(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.det)
  @Benchmark def breezeDet(): Double = det(bA)
