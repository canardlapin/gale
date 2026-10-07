package gale.bench

import gale.linalg.DMat
import gale.optim.FirstOrderSolution
import org.openjdk.jmh.annotations.*
import java.util.concurrent.TimeUnit

@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@State(Scope.Thread)
class OptimizationJmhBench:
  @Param(Array("diag-quadratic-n16-k1e3", "rosenbrock-n32", "lasso-m256-n32"))
  var caseId: String = ""
  private var run: () => FirstOrderSolution = () => throw new IllegalStateException("setup required")
  @Setup(Level.Trial)
  def setup(): Unit =
    val fixture = OptimizationBench.fixtures(OptimizationInputs.read()).find(_.id == caseId).get
    run = OptimizationBench.solver(fixture)
    require(run().certificate.primalResidual <= 1e-6)
  @Benchmark
  def solve(): FirstOrderSolution = run()

/** Isolates the removed full-matrix temporary from the original affine update. */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@State(Scope.Thread)
class OptimizationAffineBench:
  @Param(Array("16", "512", "4096"))
  var size: Int = 0
  private var x = DMat.zeros(1, 1)
  private var g = DMat.zeros(1, 1)
  @Setup(Level.Trial)
  def setup(): Unit =
    x = DMat.tabulate(size, 1)((r, _) => r.toDouble)
    g = DMat.tabulate(size, 1)((r, _) => Math.sin(r.toDouble))
  @Benchmark
  def previousTwoPass(): DMat =
    val scaled = DMat.tabulate(size, 1)((r, _) => -.01 * g(r, 0))
    DMat.tabulate(size, 1)((r, _) => x(r, 0) + scaled(r, 0))
  @Benchmark
  def fused(): DMat = DMat.tabulate(size, 1)((r, _) => x(r, 0) - .01 * g(r, 0))

@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@State(Scope.Thread)
class OptimizationProximalBench:
  @Setup(Level.Trial)
  def check(): Unit =
    ProximalComparison.verify(ProximalComparison.plain())
    ProximalComparison.verify(ProximalComparison.accelerated())
    ()
  @Benchmark
  def fixedStep(): FirstOrderSolution = ProximalComparison.plain()
  @Benchmark
  def accelerated(): FirstOrderSolution = ProximalComparison.accelerated()

/** Same prepared full solves as the portable comparisons, with forked allocation/stack profiling. */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@State(Scope.Thread)
class OptimizationPortfolioBench:
  @Param(
    Array(
      "lm-exp-offset-256",
      "lm-scaled-linear-8",
      "box-diagonal-512",
      "box-rotated-32",
      "box-diag-quadratic-n16-k1e3",
      "lasso-m256-n32",
      "lasso-m1024-n128"
    )
  )
  var caseId: String = ""
  private var run: () => FirstOrderSolution = () => throw new IllegalStateException("setup required")
  @Setup(Level.Trial)
  def setup(): Unit =
    val extensions =
      OptimizationBench.fixtures(OptimizationInputs.read("docs/verification/optim-extensions/fixtures.tsv"))
    run = extensions.find(_.id == caseId) match
      case Some(fixture) => OptimizationExtensionsBench.runner(fixture)
      case None          =>
        OptimizationBench.solver(OptimizationBench.fixtures(OptimizationInputs.read()).find(_.id == caseId).get)
    require(run().certificate.primalResidual <= 1e-7)
  @Benchmark
  def solve(): FirstOrderSolution = run()
