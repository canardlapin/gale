package gale.bench

import gale.linalg.*
import gale.sparse.*
import java.util.Random
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*

/** Sparse matrix-vector product (CSR SpMV) at 1% density over `n x n` matrices. */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SpmvJmh:
  @Param(Array("1024", "16384"))
  var n: Int = 0

  private val density = 0.01

  private var a: CSR = _
  private var x: DVec = _
  private var y: MutableDVec = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    val rng = new Random(20260709L)
    val perRow = math.max(1, (density * n).toInt)
    val builder = Sparse.coo(n, n)
    var row = 0
    while row < n do
      var k = 0
      while k < perRow do
        builder.add(row, rng.nextInt(n), rng.nextDouble() * 2.0 - 1.0)
        k += 1
      row += 1
    a = builder.toCSR()
    x = Vec.tabulate(n)(i => ((i % 13) - 6).toDouble * 0.1)
    y = MutableVec.zeros(n)

  @Benchmark
  def spmv(): Double =
    a.mulInto(x, y)
    y(0)
