package gale.bench

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.backend.PureBackend
import gale.backend.jvm.vector.VectorBackend
import gale.linalg.{DMat, DVec, Matrix, Vec}
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

import java.util.concurrent.TimeUnit

/** Copy-inclusive public-facade comparison used to derive Vector backend thresholds. */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class VectorGemmJmh:
  @Param(Array("64", "128", "256", "512"))
  var n: Int = _

  private var a: DMat = _
  private var b: DMat = _
  private var breezeA: BDM[Double] = _
  private var breezeB: BDM[Double] = _

  @Setup(Level.Trial)
  def setup(): Unit =
    a = Matrix.tabulate(n, n)((i, j) => (((i.toLong * 131 + j.toLong * 17 + 3) % 97).toDouble - 48) / 49)
    b = Matrix.tabulate(n, n)((i, j) => (((i.toLong * 113 + j.toLong * 19 + 11) % 89).toDouble - 44) / 47)
    breezeA = BDM.tabulate(n, n)((i, j) => a(i, j))
    breezeB = BDM.tabulate(n, n)((i, j) => b(i, j))

  @Benchmark
  def pure(bh: Blackhole): Unit =
    bh.consume(a.*(b)(using PureBackend))

  @Benchmark
  def vector(bh: Blackhole): Unit =
    bh.consume(a.*(b)(using VectorBackend))

  /** Same-process control: with the Vector module enabled, Breeze 2.1 delegates
    * this operation to dev.ludovic's pure-Java VectorBLAS on JDK 22.
    */
  @Benchmark
  def breezeVectorBlas(bh: Blackhole): Unit =
    bh.consume(breezeA * breezeB)

/** Public-facade GEMV crossover for the optional Vector backend. */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 8, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Benchmark)
class VectorGemvJmh:
  @Param(Array("64", "128", "256", "1024", "2048"))
  var n: Int = _

  private var a: DMat = _
  private var x: DVec = _
  private var breezeA: BDM[Double] = _
  private var breezeX: BDV[Double] = _

  @Setup(Level.Trial)
  def setup(): Unit =
    a = Matrix.tabulate(n, n)((i, j) => (((i.toLong * 131 + j.toLong * 17 + 3) % 97).toDouble - 48) / 49)
    x = Vec.tabulate(n)(i => (((i.toLong * 29 + 5) % 53).toDouble - 26) / 27)
    breezeA = BDM.tabulate(n, n)((i, j) => a(i, j))
    breezeX = BDV.tabulate(n)(i => x(i))

  @Benchmark
  def pure(bh: Blackhole): Unit =
    bh.consume(a.*(x)(using PureBackend))

  @Benchmark
  def vector(bh: Blackhole): Unit =
    bh.consume(a.*(x)(using VectorBackend))

  @Benchmark
  def breezeVectorBlas(bh: Blackhole): Unit =
    bh.consume(breezeA * breezeX)
