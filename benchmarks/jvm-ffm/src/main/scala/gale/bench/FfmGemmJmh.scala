package gale.bench

import gale.backend.{Backend, PureBackend}
import gale.backend.jvm.blas.{FfmBlasBackend, FfmBlasThresholds}
import gale.backend.jvm.`native`.{Layout, NativeDMat}
import gale.linalg.{DMat, Matrix}

import org.openjdk.jmh.annotations.*
import java.util.concurrent.TimeUnit

/** Heap-copy-inclusive crossover: both methods enter through the public `DMat.*`
  * facade and allocate their result; the FFM case additionally includes all
  * heap/native copies and the native call.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 1)
@Measurement(iterations = 8, time = 1)
@Fork(2)
@State(Scope.Benchmark)
class FfmGemmJmh:
  @Param(Array("64", "128", "256", "512", "1024"))
  var n: Int = _

  private var a: DMat = _
  private var b: DMat = _
  private var nativeBackend: FfmBlasBackend = _
  private var nativeA: NativeDMat = _
  private var nativeB: NativeDMat = _
  private var nativeC: NativeDMat = _

  @Setup(Level.Trial)
  def setup(): Unit =
    a = Matrix.tabulate(n, n)((i, j) => ((i * 17 + j * 11) % 31 - 15).toDouble / 16.0)
    b = Matrix.tabulate(n, n)((i, j) => ((i * 7 - j * 13) % 29 - 14).toDouble / 15.0)
    nativeBackend = FfmBlasBackend
      .load(thresholds = Some(FfmBlasThresholds(nativeGemmMinFlops = 0L)))
      .fold(throw _, identity)
    nativeA = NativeDMat.allocate(n, n, Layout.RowMajor)
    nativeB = NativeDMat.allocate(n, n, Layout.RowMajor)
    nativeC = NativeDMat.allocate(n, n, Layout.RowMajor)
    var i = 0
    while i < n do
      var j = 0
      while j < n do
        nativeA(i, j) = a(i, j)
        nativeB(i, j) = b(i, j)
        j += 1
      i += 1

  @TearDown(Level.Trial)
  def tearDown(): Unit =
    nativeC.close()
    nativeB.close()
    nativeA.close()
    nativeBackend.close()

  @Benchmark
  def pureGemm(): DMat =
    a.*(b)(using PureBackend)

  @Benchmark
  def ffmGemm(): DMat =
    a.*(b)(using nativeBackend)

  @Benchmark
  def ffmNativeGemm(): NativeDMat =
    nativeBackend.gemm(nativeA, nativeB, nativeC)
    nativeC
