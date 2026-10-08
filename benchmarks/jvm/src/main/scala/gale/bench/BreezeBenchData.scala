package gale.bench

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.backend.Backend
import gale.backend.PureBackend
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Matrix
import gale.linalg.Vec
import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.StandardOpenOption
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams

/** Seeded input generators shared by the paired gale-vs-Breeze JMH benchmarks.
  *
  * Every input is generated once as a plain array and then handed to '''both'''
  * libraries, so a gale benchmark and its Breeze twin run on identical data at the
  * same `@Param` size. Generation happens in `@Setup` (never in a timed method), so
  * these allocations do not count against per-invocation cost.
  *
  * Self-contained on purpose (it does not reach into the parity module): the two
  * harnesses evolve independently.
  */
object BreezeBenchData:

  /** Deterministic pseudo-random `[-1, 1)` sequence — a tiny LCG so a given seed
    * yields the same data on every fork.
    */
  private def fill(out: Array[Double], seed: Long): Unit =
    var state = seed * 6364136223846793005L + 1442695040888963407L
    var i = 0
    while i < out.length do
      state = state * 6364136223846793005L + 1442695040888963407L
      out(i) = ((state >>> 11).toDouble / (1L << 53).toDouble) * 2.0 - 1.0
      i += 1

  def vectorData(n: Int, seed: Long): Array[Double] =
    val out = new Array[Double](n)
    fill(out, seed)
    out

  def matrixData(rows: Int, cols: Int, seed: Long): Array[Array[Double]] =
    Array.tabulate(rows)(i => vectorData(cols, seed + i * 0x9e3779b9L))

  /** Strictly diagonally dominant (well-conditioned, nonsingular): each diagonal is
    * set to its row's absolute off-diagonal sum plus one.
    */
  def diagonallyDominant(n: Int, seed: Long): Array[Array[Double]] =
    val a = matrixData(n, n, seed)
    var i = 0
    while i < n do
      var sum = 0.0
      var j = 0
      while j < n do
        if j != i then sum += math.abs(a(i)(j))
        j += 1
      a(i)(i) = sum + 1.0
      i += 1
    a

  /** Symmetric positive-definite `B Bᵀ + n·I`, built exactly symmetric. */
  def spd(n: Int, seed: Long): Array[Array[Double]] =
    val b = matrixData(n, n, seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        var s = 0.0
        var k = 0
        while k < n do
          s += b(i)(k) * b(j)(k)
          k += 1
        if i == j then s += n.toDouble
        a(i)(j) = s
        a(j)(i) = s
        j += 1
      i += 1
    a

  /** Symmetric matrix with entries in `[-1, 1)` (lower triangle mirrored). */
  def symmetric(n: Int, seed: Long): Array[Array[Double]] =
    val src = matrixData(n, n, seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        a(i)(j) = src(i)(j)
        a(j)(i) = src(i)(j)
        j += 1
      i += 1
    a

  def galeMatrix(data: Array[Array[Double]]): DMat =
    Matrix.tabulate(data.length, if data.isEmpty then 0 else data(0).length)((i, j) => data(i)(j))

  def breezeMatrix(data: Array[Array[Double]]): BDM[Double] =
    BDM.tabulate(data.length, if data.isEmpty then 0 else data(0).length)((i, j) => data(i)(j))

  def galeVector(data: Array[Double]): DVec =
    Vec.tabulate(data.length)(i => data(i))

  def breezeVector(data: Array[Double]): BDV[Double] =
    BDV(data.clone())

  /** System property naming the netlib sidecar file (JSON Lines, appended). */
  final val NetlibSidecarProperty = "gale.bench.netlibSidecar"

  /** Default sidecar path, relative to the forked JVM's working directory. */
  final val DefaultNetlibSidecar = "target/breeze-netlib.jsonl"

  /** Records which `dev.ludovic.netlib` BLAS/LAPACK implementation Breeze resolved in
    * this fork. Called from every Breeze bench's trial setup: one line to stderr and one
    * JSON object appended to the sidecar (`-Dgale.bench.netlibSidecar=...`, default
    * [[DefaultNetlibSidecar]]), which `tools/bench/breeze_scoreboard.py` reads to
    * label — and, for lane A, reject — a receipt. Netlib picks native JNI first, then
    * VectorBLAS when `jdk.incubator.vector` is resolvable, then the scalar Java BLAS.
    */
  def recordNetlib(params: BenchmarkParams): Unit =
    val blas   = dev.ludovic.netlib.blas.BLAS.getInstance().getClass.getName
    val lapack = dev.ludovic.netlib.lapack.LAPACK.getInstance().getClass.getName
    val lane   = System.getProperty("gale.bench.lane", "unset")
    val vectorModule = ModuleLayer.boot().findModule("jdk.incubator.vector").isPresent
    val paramText =
      params.getParamsKeys.toArray.map(k => s"$k=${params.getParam(k.toString)}").mkString(",")
    System.err.println(
      s"[breeze-netlib] lane=$lane benchmark=${params.getBenchmark} params=$paramText blas=$blas lapack=$lapack"
    )
    def q(v: String): String = "\"" + v.replace("\\", "\\\\").replace("\"", "\\\"") + "\""
    val line =
      s"""{"lane":${q(lane)},"benchmark":${q(params.getBenchmark)},"params":${q(paramText)},""" +
        s""""jdk":${q(System.getProperty("java.version"))},"vectorModule":$vectorModule,""" +
        s""""blas":${q(blas)},"lapack":${q(lapack)}}\n"""
    val path = Paths.get(System.getProperty(NetlibSidecarProperty, DefaultNetlibSidecar))
    Option(path.toAbsolutePath.getParent).foreach(Files.createDirectories(_))
    Files.write(
      path,
      line.getBytes(StandardCharsets.UTF_8),
      StandardOpenOption.CREATE,
      StandardOpenOption.APPEND
    )

/** The gale-side backend switch for the paired Breeze benchmarks. Only the
  * backend-sensitive `gale*` methods take this state, so JMH expands `backend` for them
  * alone and each Breeze twin runs once per size. Gale methods whose operation takes no
  * `Backend` (L1 `dot`/`norm2`/`axpyInPlace`, symmetric eigen) omit it and run once.
  *
  *   - `pure` — [[gale.backend.PureBackend]]: lane A (out-of-box) and the lane-B control.
  *   - `vector` — the `backend-jvm-vector` SIMD backend: lane B only. It needs
  *     `--add-modules=jdk.incubator.vector`; without the module the trial fails fast with
  *     an explicit message instead of a `NoClassDefFoundError` (run lane A with
  *     `-p backend=pure`, which the `breezeLaneA` alias does).
  */
@State(Scope.Thread)
class GaleBackendState:
  @Param(Array("pure", "vector"))
  var backend: String = "pure"

  var selected: Backend = PureBackend

  @Setup(Level.Trial)
  def selectBackend(): Unit =
    selected = backend match
      case "pure"   => PureBackend
      case "vector" =>
        if !ModuleLayer.boot().findModule("jdk.incubator.vector").isPresent then
          throw new IllegalStateException(
            "backend=vector needs --add-modules=jdk.incubator.vector (lane B); lane A must run with -p backend=pure"
          )
        VectorBackendLoader.load()
      case other    => throw new IllegalArgumentException(s"unknown gale backend '$other' (expected pure|vector)")

/** Isolates the only reference to the Vector API backend so that loading
  * [[GaleBackendState]] never resolves `jdk.incubator.vector` classes.
  */
private object VectorBackendLoader:
  def load(): Backend = gale.backend.jvm.vector.VectorBackend
