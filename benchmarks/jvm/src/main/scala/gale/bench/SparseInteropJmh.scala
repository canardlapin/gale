package gale.bench

import gale.sparse.COO
import gale.sparse.CSR
import gale.sparse.Sparse
import gale.sparse.SparseEntryConsumer
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Allocation baseline contrasting the compatibility `Seq[COOEntry]` export
  * with reusable primitive COO/CSR consumers. Run with `-prof gc` to record
  * bytes/op; the primitive paths allocate no per-entry objects or collections.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 3, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SparseInteropJmh:
  @Param(Array("1024", "16384"))
  var entries: Int = 0

  private var coo: COO = _
  private var csr: CSR = _
  private val consumer = new SumConsumer

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    val builder = Sparse.coo(entries, 1)
    var row = 0
    while row < entries do
      builder.add(row, 0, row.toDouble + 1.0)
      row += 1
    coo = builder.toCOO()
    csr = coo.toCSR

  @Benchmark
  def boxedCooEntries(blackhole: Blackhole): Unit =
    blackhole.consume(coo.entries)

  @Benchmark
  def primitiveCooTraversal(): Double =
    consumer.sum = 0.0
    coo.foreachStoredEntry(consumer)
    consumer.sum

  @Benchmark
  def primitiveCsrTraversal(): Double =
    consumer.sum = 0.0
    csr.foreachStoredEntry(consumer)
    consumer.sum

  private final class SumConsumer extends SparseEntryConsumer:
    var sum = 0.0

    def apply(row: Int, col: Int, value: Double): Unit =
      sum += row.toDouble + col.toDouble + value
