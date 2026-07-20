package gale.bench

import gale.linalg.DMat
import gale.linalg.DVec
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Allocation baseline for the owned dense transformation APIs. The primitive
  * copy methods write into caller-owned arrays through direct loops; run with
  * `-prof gc` to contrast them with the boxed compatibility exports.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 3, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class DenseTransformJmh:
  @Param(Array("1024", "16384"))
  var entries: Int = 0

  private var matrix: DMat = _
  private var vector: DVec = _
  private var matrixDestination: Array[Double] = _
  private var vectorDestination: Array[Double] = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    val columns = 16
    val rows = entries / columns
    // Both inputs are non-contiguous logical views: a transpose and a strided
    // column. The timed copies therefore exercise the general traversal paths.
    matrix = DMat.tabulate(columns, rows)((row, col) => row.toDouble + col.toDouble).t
    vector = DMat.tabulate(entries, 2)((row, col) => row.toDouble + col.toDouble).col(1)
    matrixDestination = new Array[Double](entries)
    vectorDestination = new Array[Double](entries)

  @Benchmark
  def boxedMatrixValues(blackhole: Blackhole): Unit =
    blackhole.consume(matrix.valuesRowMajor)

  @Benchmark
  def primitiveMatrixCopy(): Double =
    matrix.copyRowMajorTo(matrixDestination)
    matrixDestination(entries - 1)

  @Benchmark
  def boxedVectorValues(blackhole: Blackhole): Unit =
    blackhole.consume(vector.toSeq)

  @Benchmark
  def primitiveVectorCopy(): Double =
    vector.copyTo(vectorDestination)
    vectorDestination(entries - 1)
