package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseMatrix as BDM
import gale.bench.BreezeBenchData.*
import gale.linalg.DMat
import gale.syntax.all.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams
import org.openjdk.jmh.infra.Blackhole

/** Elementwise matrix arithmetic, gale vs Breeze, on square `n × n` operands; every
  * operation allocates its result on both sides:
  *
  *   - `add` / `sub`: `A + B`, `A - B`.
  *   - `hadamard`: gale `A.pointwise * B` vs breeze `A *:* B`.
  *
  * Backend-insensitive: gale's elementwise ops take no `Backend`.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class ElementwiseBreezeJmh:
  @Param(Array("256", "1024"))
  var n: Int = 0

  private var gA: DMat        = uninitialized
  private var gB: DMat        = uninitialized
  private var bA: BDM[Double] = uninitialized
  private var bB: BDM[Double] = uninitialized

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    val aData = matrixData(n, n, 5100L)
    val bData = matrixData(n, n, 5200L)
    gA = galeMatrix(aData)
    gB = galeMatrix(bData)
    bA = breezeMatrix(aData)
    bB = breezeMatrix(bData)

  @Benchmark def galeAdd(bh: Blackhole): Unit   = bh.consume(gA + gB)
  @Benchmark def breezeAdd(bh: Blackhole): Unit = bh.consume(bA + bB)

  @Benchmark def galeSub(bh: Blackhole): Unit   = bh.consume(gA - gB)
  @Benchmark def breezeSub(bh: Blackhole): Unit = bh.consume(bA - bB)

  @Benchmark def galeHadamard(bh: Blackhole): Unit   = bh.consume(gA.pointwise * gB)
  @Benchmark def breezeHadamard(bh: Blackhole): Unit = bh.consume(bA *:* bB)
