package gale.bench

import gale.linalg.*
import gale.solvers.Preconditioner
import gale.spectral.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Matrix-free generalized symmetric-definite LOBPCG scenarios.
  *
  * The parameter product covers three ambient dimensions, three block sizes,
  * two pencils, and three preconditioners. Run with `-prof gc` for allocation
  * data. [[GeneralizedLobpcgWorkReceipt]] executes the same fixed fixtures once
  * each and records convergence plus exact operator/preconditioner work.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 1, time = 200, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 2, time = 200, timeUnit = TimeUnit.MILLISECONDS)
@Fork(1)
@State(Scope.Thread)
class GeneralizedLobpcgJmh:
  @Param(Array("128", "512", "2048"))
  var n: Int = 0

  @Param(Array("4", "8", "16"))
  var k: Int = 0

  @Param(Array("clustered-diagonal", "stiffness-mass"))
  var pencil: String = ""

  @Param(Array("identity", "jacobi", "block-jacobi"))
  var preconditioner: String = ""

  private var scenario: GeneralizedLobpcgScenario = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    scenario = GeneralizedLobpcgScenario(n, k, pencil, preconditioner)

  @Benchmark
  def solve(blackhole: Blackhole): Unit =
    blackhole.consume(scenario.run())

/** One untimed work-accounting pass over every JMH parameter combination.
  *
  * Run with:
  * `sbt "benchmarksJVM/Jmh/runMain gale.bench.GeneralizedLobpcgWorkReceipt"`.
  */
object GeneralizedLobpcgWorkReceipt:
  def main(args: Array[String]): Unit =
    println(
      "n,k,pencil,preconditioner,converged,requested,iterations," +
        "a_applications,b_applications,preconditioner_applications,worst_residual"
    )
    for
      n <- Seq(128, 512, 2048)
      k <- Seq(4, 8, 16)
      pencil <- Seq("clustered-diagonal", "stiffness-mass")
      preconditioner <- Seq("identity", "jacobi", "block-jacobi")
    do
      val scenario = GeneralizedLobpcgScenario(n, k, pencil, preconditioner)
      val result = scenario.run()
      val diagnostics = result.diagnostics
      println(
        Seq(
          n,
          k,
          pencil,
          preconditioner,
          diagnostics.converged,
          diagnostics.requested,
          diagnostics.iterations,
          scenario.aApplications,
          scenario.bApplications,
          scenario.preconditionerApplications,
          diagnostics.worstResidual
        ).mkString(",")
      )

private final class GeneralizedLobpcgScenario private (
    val n: Int,
    val k: Int,
    operator: DoubleLinearOperator,
    metric: DoubleLinearOperator,
    preconditioner: Preconditioner,
    counters: GeneralizedLobpcgCounters
):
  private val options =
    GeneralizedSpectralOptions(tolerance = 1e-6, maxIterations = 40)

  def run(): EigenDecomposition =
    counters.reset()
    Eigen
      .eigSymmetricGeneralized(
        operator.assumeSymmetricOperator,
        metric.assumePositiveDefiniteOperator,
        n,
        EigenSelection.Count(k, EigenOrder.SmallestAlgebraic),
        options,
        preconditioner
      )(using SpectralBackend.none)
      .fold(throw _, identity)

  def aApplications: Long = counters.aApplications
  def bApplications: Long = counters.bApplications
  def preconditionerApplications: Long = counters.preconditionerApplications

private object GeneralizedLobpcgScenario:
  def apply(
      n: Int,
      k: Int,
      pencil: String,
      preconditionerKind: String
  ): GeneralizedLobpcgScenario =
    require(n > 0, s"n must be positive, got $n")
    require(k > 0 && k < n, s"k must be in [1, n), got k=$k and n=$n")
    val counters = new GeneralizedLobpcgCounters
    val data =
      pencil match
        case "clustered-diagonal" => clusteredDiagonal(n, k, counters)
        case "stiffness-mass"     => stiffnessMass(n, counters)
        case other                => throw new IllegalArgumentException(s"unknown pencil '$other'")
    val preconditioner =
      preconditionerKind match
        case "identity" =>
          new CountingPreconditioner(counters, IndexedSeq.fill(n)(1.0), IndexedSeq.fill(math.max(0, n - 1))(0.0))
        case "jacobi" =>
          new CountingPreconditioner(counters, data.operatorDiagonal, IndexedSeq.fill(math.max(0, n - 1))(0.0))
        case "block-jacobi" =>
          new CountingPreconditioner(counters, data.operatorDiagonal, data.operatorOffDiagonal)
        case other =>
          throw new IllegalArgumentException(s"unknown preconditioner '$other'")
    new GeneralizedLobpcgScenario(
      n,
      k,
      data.operator,
      data.metric,
      preconditioner,
      counters
    )

  private final case class PencilData(
      operator: DoubleLinearOperator,
      metric: DoubleLinearOperator,
      operatorDiagonal: IndexedSeq[Double],
      operatorOffDiagonal: IndexedSeq[Double]
  )

  private def clusteredDiagonal(
      n: Int,
      k: Int,
      counters: GeneralizedLobpcgCounters
  ): PencilData =
    val metricDiagonal =
      IndexedSeq.tabulate(n)(i => 0.75 + ((i * 17 + 11) % 29).toDouble / 29.0)
    val generalizedValues =
      IndexedSeq.tabulate(n): i =>
        if i < 2 * k then 1.0 + (i % math.max(1, k)).toDouble * 1e-7
        else 2.0 + (i - 2 * k).toDouble / math.max(1, n - 2 * k).toDouble
    val operatorDiagonal =
      generalizedValues.indices.map(i => generalizedValues(i) * metricDiagonal(i))
    PencilData(
      new DiagonalOperator(operatorDiagonal, () => counters.aApplications += 1),
      new DiagonalOperator(metricDiagonal, () => counters.bApplications += 1),
      operatorDiagonal,
      IndexedSeq.fill(math.max(0, n - 1))(0.0)
    )

  private def stiffnessMass(
      n: Int,
      counters: GeneralizedLobpcgCounters
  ): PencilData =
    val operatorDiagonal =
      IndexedSeq.tabulate(n)(i => 2.5 + 0.2 * math.sin((i + 1).toDouble * 0.17))
    val operatorOffDiagonal = IndexedSeq.fill(math.max(0, n - 1))(-1.0)
    val metricDiagonal = IndexedSeq.fill(n)(4.0 / 6.0)
    val metricOffDiagonal = IndexedSeq.fill(math.max(0, n - 1))(1.0 / 6.0)
    PencilData(
      new TridiagonalOperator(
        operatorDiagonal,
        operatorOffDiagonal,
        () => counters.aApplications += 1
      ),
      new TridiagonalOperator(
        metricDiagonal,
        metricOffDiagonal,
        () => counters.bApplications += 1
      ),
      operatorDiagonal,
      operatorOffDiagonal
    )

private final class GeneralizedLobpcgCounters:
  var aApplications: Long = 0L
  var bApplications: Long = 0L
  var preconditionerApplications: Long = 0L

  def reset(): Unit =
    aApplications = 0L
    bApplications = 0L
    preconditionerApplications = 0L

private final class DiagonalOperator(
    diagonal: IndexedSeq[Double],
    count: () => Unit
) extends DoubleLinearOperator:
  def rows: Int = diagonal.length
  def cols: Int = diagonal.length

  def applyTo(x: DVec, into: MutableDVec): Unit =
    count()
    var i = 0
    while i < diagonal.length do
      into(i) = diagonal(i) * x(i)
      i += 1

private final class TridiagonalOperator(
    diagonal: IndexedSeq[Double],
    offDiagonal: IndexedSeq[Double],
    count: () => Unit
) extends DoubleLinearOperator:
  def rows: Int = diagonal.length
  def cols: Int = diagonal.length

  def applyTo(x: DVec, into: MutableDVec): Unit =
    count()
    var i = 0
    while i < diagonal.length do
      var value = diagonal(i) * x(i)
      if i > 0 then value += offDiagonal(i - 1) * x(i - 1)
      if i + 1 < diagonal.length then value += offDiagonal(i) * x(i + 1)
      into(i) = value
      i += 1

/** Pairwise block-Jacobi for the stiffness matrix; a zero off-diagonal reduces
  * exactly to diagonal Jacobi for the clustered pencil.
  */
private final class CountingPreconditioner(
    counters: GeneralizedLobpcgCounters,
    diagonal: IndexedSeq[Double],
    offDiagonal: IndexedSeq[Double]
) extends Preconditioner:
  def solve(r: DVec, into: MutableVec[Double]): Unit =
    counters.preconditionerApplications += 1
    var i = 0
    while i + 1 < diagonal.length do
      val off = offDiagonal(i)
      val determinant = diagonal(i) * diagonal(i + 1) - off * off
      into(i) = (diagonal(i + 1) * r(i) - off * r(i + 1)) / determinant
      into(i + 1) = (diagonal(i) * r(i + 1) - off * r(i)) / determinant
      i += 2
    if i < diagonal.length then into(i) = r(i) / diagonal(i)
