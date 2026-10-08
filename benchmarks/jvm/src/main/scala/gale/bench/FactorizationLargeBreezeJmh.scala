package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.LU as BreezeLU
import breeze.linalg.cholesky
import breeze.linalg.eigSym
import breeze.linalg.qr
import dev.ludovic.netlib.lapack.LAPACK
import gale.backend.Backend
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import gale.spectral.*
import java.util.concurrent.TimeUnit
import org.netlib.util.intW
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams
import org.openjdk.jmh.infra.Blackhole

/** Large dense factorizations (`n` in {512, 1024}), gale vs Breeze — the sizes at
  * which blocked LAPACK pays off. Same op names and pairings as
  * [[FactorizationBreezeJmh]], [[LeastSquaresBreezeJmh]] and [[SymEigenBreezeJmh]],
  * kept in a separate class so the default small-size sweeps stay bounded:
  *
  *   - `solve`: gale `A.solve(b)` vs breeze `A \ b` (`dgesv`).
  *   - `lu`: gale `A.lu` vs breeze `LU.primitive` (`dgetrf` only).
  *   - `chol`: gale `S.cholesky` vs breeze `cholesky(S)` on an SPD matrix.
  *   - `qr`: gale `A.qr` (Q lazy) vs breeze `qr.justR(A)`.
  *   - `lstsq`: tall `m = 2n` (not the small class's `4n`, to bound run time) — gale
  *     `A.leastSquares(b)` vs breeze `A \ b` (`dgels`).
  *   - `eigSym`: values and vectors — gale `Eigen.eigSymmetric(All, Right)` vs
  *     breeze `eigSym` (`dsyev`); backend-insensitive (no `GaleBackendState`).
  *
  * A single operation takes from tens of milliseconds to seconds here, so the class
  * reports average time with one-second iterations.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3, time = 1, timeUnit = TimeUnit.SECONDS)
@Measurement(iterations = 5, time = 1, timeUnit = TimeUnit.SECONDS)
@Fork(2)
@State(Scope.Thread)
class FactorizationLargeBreezeJmh:
  @Param(Array("512", "1024"))
  var n: Int = 0

  private var gA: DMat        = uninitialized
  private var gS: DMat        = uninitialized
  private var gSym: DMat      = uninitialized
  private var gb: DVec        = uninitialized
  private var gTall: DMat     = uninitialized
  private var gTallB: DVec    = uninitialized
  private var bA: BDM[Double] = uninitialized
  private var bS: BDM[Double] = uninitialized
  private var bSym: BDM[Double] = uninitialized
  private var bb: BDV[Double] = uninitialized
  private var bTall: BDM[Double] = uninitialized
  private var bTallB: BDV[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val aData    = diagonallyDominant(n, 1100L)
    val sData    = spd(n, 1200L)
    val symData  = symmetric(n, 1300L)
    val bData    = vectorData(n, 1400L)
    val tallData = matrixData(2 * n, n, 1500L)
    val tallB    = vectorData(2 * n, 1600L)
    gA = galeMatrix(aData)
    gS = galeMatrix(sData)
    gSym = galeMatrix(symData)
    gb = galeVector(bData)
    gTall = galeMatrix(tallData)
    gTallB = galeVector(tallB)
    bA = breezeMatrix(aData)
    bS = breezeMatrix(sData)
    bSym = breezeMatrix(symData)
    bb = breezeVector(bData)
    bTall = breezeMatrix(tallData)
    bTallB = breezeVector(tallB)

  @Benchmark def galeSolve(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.solve(gb))
  @Benchmark def breezeSolve(bh: Blackhole): Unit = bh.consume(bA \ bb)

  @Benchmark def galeLu(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.lu)
  @Benchmark def breezeLu(bh: Blackhole): Unit = bh.consume(BreezeLU.primitive(bA))

  @Benchmark def galeChol(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gS.cholesky)
  @Benchmark def breezeChol(bh: Blackhole): Unit = bh.consume(cholesky(bS))

  @Benchmark def galeQr(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.qr)
  @Benchmark def breezeQr(bh: Blackhole): Unit = bh.consume(qr.justR(bA))

  @Benchmark def galeLstsq(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gTall.leastSquares(gTallB))
  @Benchmark def breezeLstsq(bh: Blackhole): Unit = bh.consume(bTall \ bTallB)

  @Benchmark def galeEigSym(bh: Blackhole): Unit =
    bh.consume(Eigen.eigSymmetric(gSym, EigenSelection.All, EigenVectors.Right))
  @Benchmark def breezeEigSym(bh: Blackhole): Unit = bh.consume(eigSym(bSym))

/** Multi-right-hand-side solves, `n = 256` with `k = 64` columns, gale vs Breeze:
  *
  *   - `luSolve`: gale `A.solve(B)` (one LU, all columns) vs breeze `A \ B` (`dgesv`).
  *   - `cholSolve`: gale `S.cholesky` then `solve(B)` vs Breeze's LAPACK stack — Breeze
  *     has no Cholesky-solve entry point, so its twin calls the same netlib `LAPACK`
  *     instance Breeze resolves: `dpotrf` + `dpotrs` on copies of `S` and `B`
  *     (column-major, Breeze's layout). Both sides copy their inputs once, as
  *     `cholesky(S)` would.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class MultiRhsBreezeJmh:
  @Param(Array("256"))
  var n: Int = 0

  @Param(Array("64"))
  var k: Int = 0

  private var gA: DMat        = uninitialized
  private var gS: DMat        = uninitialized
  private var gB: DMat        = uninitialized
  private var bA: BDM[Double] = uninitialized
  private var bS: BDM[Double] = uninitialized
  private var bB: BDM[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val aData = diagonallyDominant(n, 1700L)
    val sData = spd(n, 1800L)
    val bData = matrixData(n, k, 1900L)
    gA = galeMatrix(aData)
    gS = galeMatrix(sData)
    gB = galeMatrix(bData)
    bA = breezeMatrix(aData)
    bS = breezeMatrix(sData)
    bB = breezeMatrix(bData)

  @Benchmark def galeLuSolve(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gA.solve(gB))
  @Benchmark def breezeLuSolve(bh: Blackhole): Unit = bh.consume(bA \ bB)

  @Benchmark def galeCholSolve(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(gS.cholesky.flatMap(_.solve(gB)))
  @Benchmark def breezeCholSolve(bh: Blackhole): Unit =
    val lapack = LAPACK.getInstance()
    val a      = bS.data.clone()
    val x      = bB.data.clone()
    val info   = new intW(0)
    lapack.dpotrf("L", n, a, n, info)
    if info.`val` != 0 then throw new IllegalStateException(s"dpotrf info=${info.`val`}")
    lapack.dpotrs("L", n, k, a, n, x, n, info)
    if info.`val` != 0 then throw new IllegalStateException(s"dpotrs info=${info.`val`}")
    bh.consume(x)
