package gale.bench

import scala.compiletime.uninitialized

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.optimize.DiffFunction
import breeze.optimize.LBFGS as BreezeLBFGS
import gale.backend.Backend
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import gale.optim.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.BenchmarkParams
import org.openjdk.jmh.infra.Blackhole

/** Full L-BFGS solves, gale `LBFGS.minimize` vs breeze `LBFGS.minimizeAndReturnState`,
  * timing one complete solve per operation from the same start point:
  *
  *   - `rosenbrock`: the chained Rosenbrock function in 100 variables,
  *     `Σ 100(x_{i+1} − x_i²)² + (1 − x_i)²`, from `(−1.2, 1, −1.2, 1, …)`.
  *   - `logistic`: L2-regularized logistic regression,
  *     `Σ log(1 + exp(−y_i x_iᵀw)) + (λ/2)‖w‖²` with `λ = 1`, on fixed synthetic data
  *     (`1000 × 50` design, labels `±1` from a seeded linear model plus noise). The
  *     objective uses each library's own matrix-vector products (`X w`, `Xᵀ r`), so
  *     only this gale method takes the backend.
  *
  * `budget` selects the stopping rule:
  *
  *   - `fixed`: exactly [[LbfgsBreezeJmh.FixedIterations]] iterations on both sides,
  *     tolerances disabled — the fair per-iteration comparison. Setup verifies both
  *     libraries actually run the full budget.
  *   - `tolerance`: each library's own convergence test at `1e-8` (gale: gradient
  *     ∞-norm against `1e-8 + 1e-8·scale`; Breeze: relative function-value or gradient
  *     convergence). The criteria differ, so this is time-to-own-convergence.
  *
  * Both use history 10 and a strong-Wolfe line search (`c1 = 1e-4`, `c2 = 0.9`).
  * Evaluation counts are deterministic per solve, so rather than JMH `@AuxCounters`
  * (whose per-iteration sums would scale with the iteration length) the trial setup
  * runs one solve and logs iterations, objective/gradient evaluations and the final
  * gradient ∞-norm to stderr as `[breeze-lbfgs] ...`.
  *
  * Allocation: gale's callback receives an `n × 1` `DMat` and must return an owned
  * gradient `DMat`, built with one `Matrix.fromArrayCopy` (one extra `O(n)` copy that
  * Breeze's `DenseVector` wrap avoids).
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class LbfgsBreezeJmh:
  import LbfgsBreezeJmh.*

  @Param(Array("fixed", "tolerance"))
  var budget: String = "fixed"

  private var rosenbrockStart: Array[Double] = uninitialized
  private var gRosenbrockStart: DMat         = uninitialized
  private var gX: DMat                       = uninitialized
  private var gLabels: Array[Double]         = uninitialized
  private var bX: BDM[Double]                = uninitialized
  private var bLabels: Array[Double]         = uninitialized
  private var galeConfig: LBFGSConfig        = uninitialized
  private var breezeMaxIter: Int             = 0
  private var breezeTolerance: Double        = 0.0

  @Setup(Level.Trial)
  def setupTrial(params: BenchmarkParams): Unit =
    recordNetlib(params)
    rosenbrockStart = Array.tabulate(RosenbrockSize)(i => if i % 2 == 0 then -1.2 else 1.0)
    gRosenbrockStart = Matrix.fromArrayCopy(RosenbrockSize, 1, rosenbrockStart)
    val design = matrixData(LogisticRows, LogisticCols, 6100L)
    val truth  = vectorData(LogisticCols, 6200L)
    val noise  = vectorData(LogisticRows, 6300L)
    val labels = Array.tabulate(LogisticRows): i =>
      var z = 0.0
      var j = 0
      while j < LogisticCols do
        z += design(i)(j) * truth(j)
        j += 1
      if z + 0.5 * noise(i) >= 0.0 then 1.0 else -1.0
    gX = galeMatrix(design)
    gLabels = labels.clone()
    bX = breezeMatrix(design)
    bLabels = labels.clone()
    budget match
      case "fixed" =>
        galeConfig = LBFGSConfig(
          maxIterations = FixedIterations,
          tolerance = FirstOrderTolerance.from(Double.MinPositiveValue, 0.0).fold(e => throw new IllegalStateException(e.message), identity)
        )
        breezeMaxIter = FixedIterations
        breezeTolerance = 0.0
      case "tolerance" =>
        galeConfig = LBFGSConfig(maxIterations = 10000)
        breezeMaxIter = 10000
        breezeTolerance = 1e-8
      case other => throw new IllegalArgumentException(s"unknown budget '$other' (expected fixed|tolerance)")
    report("rosenbrock", runGaleRosenbrock(), runBreezeRosenbrock())
    report("logistic", runGaleLogistic(using gale.backend.PureBackend), runBreezeLogistic())

  private def report(
      problem: String,
      galeRun: (FirstOrderSolution, Int),
      breezeRun: (BreezeLBFGS[BDV[Double]]#State, Int)
  ): Unit =
    val (gs, gEvals)  = galeRun
    val (bs, bEvals)  = breezeRun
    val gIter         = gs.certificate.iterations
    val bGrad         = breeze.linalg.norm(bs.grad, Double.PositiveInfinity)
    System.err.println(
      s"[breeze-lbfgs] budget=$budget problem=$problem " +
        s"gale(iterations=$gIter evaluations=$gEvals status=${gs.status} f=${gs.objective} gradInf=${gs.certificate.primalResidual}) " +
        s"breeze(iterations=${bs.iter} evaluations=$bEvals reason=${bs.convergenceReason.map(_.reason).getOrElse("-")} f=${bs.value} gradInf=$bGrad)"
    )
    if budget == "fixed" && (gIter != FixedIterations || bs.iter != FixedIterations) then
      throw new IllegalStateException(
        s"$problem fixed budget not honoured: gale ran $gIter, breeze ${bs.iter} iterations (expected $FixedIterations)"
      )

  // gale ---------------------------------------------------------------------

  private def galeRosenbrockObjective: DifferentiableObjective =
    new DifferentiableObjective:
      def variableRows: Int = RosenbrockSize
      def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
        val g = new Array[Double](RosenbrockSize)
        var f = 0.0
        var i = 0
        while i < RosenbrockSize - 1 do
          val xi = at(i, 0)
          val t  = at(i + 1, 0) - xi * xi
          val u  = 1.0 - xi
          f += 100.0 * t * t + u * u
          g(i) += -400.0 * t * xi - 2.0 * u
          g(i + 1) += 200.0 * t
          i += 1
        Right(ObjectiveEvaluation(f, Matrix.fromArrayCopy(RosenbrockSize, 1, g)))

  private def galeLogisticObjective(using Backend): DifferentiableObjective =
    new DifferentiableObjective:
      def variableRows: Int = LogisticCols
      def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
        val w = at.col(0)
        val z = gX * w
        val r = new Array[Double](LogisticRows)
        var f = 0.0
        var i = 0
        while i < LogisticRows do
          val m = gLabels(i) * z(i)
          f += logOnePlusExp(-m)
          r(i) = -gLabels(i) * sigmoid(-m)
          i += 1
        val xtr = gX.t * Vec.fromArrayCopy(r)
        val g   = new Array[Double](LogisticCols)
        var ww  = 0.0
        var j   = 0
        while j < LogisticCols do
          val wj = w(j)
          ww += wj * wj
          g(j) = xtr(j) + Lambda * wj
          j += 1
        Right(ObjectiveEvaluation(f + 0.5 * Lambda * ww, Matrix.fromArrayCopy(LogisticCols, 1, g)))

  private def runGaleRosenbrock(): (FirstOrderSolution, Int) =
    val solution = LBFGS.minimize(galeRosenbrockObjective, gRosenbrockStart, galeConfig).fold(e => throw new IllegalStateException(e.message), identity)
    (solution, solution.evaluations.callbacks)

  private def runGaleLogistic(using Backend): (FirstOrderSolution, Int) =
    val start    = DMat.zeros(LogisticCols, 1)
    val solution = LBFGS.minimize(galeLogisticObjective, start, galeConfig).fold(e => throw new IllegalStateException(e.message), identity)
    (solution, solution.evaluations.callbacks)

  // breeze -------------------------------------------------------------------

  private final class Counting(f: BDV[Double] => (Double, BDV[Double])) extends DiffFunction[BDV[Double]]:
    var evaluations = 0
    def calculate(x: BDV[Double]): (Double, BDV[Double]) =
      evaluations += 1
      f(x)

  private def rosenbrockBreeze(x: BDV[Double]): (Double, BDV[Double]) =
    val g = new Array[Double](RosenbrockSize)
    var f = 0.0
    var i = 0
    while i < RosenbrockSize - 1 do
      val xi = x(i)
      val t  = x(i + 1) - xi * xi
      val u  = 1.0 - xi
      f += 100.0 * t * t + u * u
      g(i) += -400.0 * t * xi - 2.0 * u
      g(i + 1) += 200.0 * t
      i += 1
    (f, BDV(g))

  private def logisticBreeze(w: BDV[Double]): (Double, BDV[Double]) =
    val z = bX * w
    val r = new Array[Double](LogisticRows)
    var f = 0.0
    var i = 0
    while i < LogisticRows do
      val m = bLabels(i) * z(i)
      f += logOnePlusExp(-m)
      r(i) = -bLabels(i) * sigmoid(-m)
      i += 1
    val xtr = bX.t * BDV(r)
    val g   = new Array[Double](LogisticCols)
    var ww  = 0.0
    var j   = 0
    while j < LogisticCols do
      val wj = w(j)
      ww += wj * wj
      g(j) = xtr(j) + Lambda * wj
      j += 1
    (f + 0.5 * Lambda * ww, BDV(g))

  private def breezeSolver = new BreezeLBFGS[BDV[Double]](breezeMaxIter, History, breezeTolerance)

  private def runBreezeRosenbrock(): (BreezeLBFGS[BDV[Double]]#State, Int) =
    val f     = new Counting(rosenbrockBreeze)
    val state = breezeSolver.minimizeAndReturnState(f, BDV(rosenbrockStart.clone()))
    (state, f.evaluations)

  private def runBreezeLogistic(): (BreezeLBFGS[BDV[Double]]#State, Int) =
    val f     = new Counting(logisticBreeze)
    val state = breezeSolver.minimizeAndReturnState(f, BDV.zeros[Double](LogisticCols))
    (state, f.evaluations)

  // timed --------------------------------------------------------------------

  @Benchmark def galeRosenbrock(bh: Blackhole): Unit   = bh.consume(runGaleRosenbrock())
  @Benchmark def breezeRosenbrock(bh: Blackhole): Unit = bh.consume(runBreezeRosenbrock())

  @Benchmark def galeLogistic(g: GaleBackendState, bh: Blackhole): Unit =
    val backend = g.selected
    given Backend = backend
    bh.consume(runGaleLogistic)
  @Benchmark def breezeLogistic(bh: Blackhole): Unit = bh.consume(runBreezeLogistic())

object LbfgsBreezeJmh:
  final val RosenbrockSize  = 100
  final val LogisticRows    = 1000
  final val LogisticCols    = 50
  final val Lambda          = 1.0
  final val History         = 10
  final val FixedIterations = 20

  /** `log(1 + exp(t))` without overflow. */
  def logOnePlusExp(t: Double): Double =
    if t > 0.0 then t + math.log1p(math.exp(-t)) else math.log1p(math.exp(t))

  /** Stable logistic `1 / (1 + exp(-t))`. */
  def sigmoid(t: Double): Double =
    if t >= 0.0 then 1.0 / (1.0 + math.exp(-t))
    else
      val e = math.exp(t)
      e / (1.0 + e)
