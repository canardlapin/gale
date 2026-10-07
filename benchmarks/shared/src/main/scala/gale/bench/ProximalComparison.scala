package gale.bench

import gale.linalg.DMat
import gale.optim.*

/** Equal-accuracy comparison against the retained fixed-step solver, with analytic KKT checks. */
object ProximalComparison:
  private val n = 64
  private val diagonal = Array.tabulate(n)(i => Math.pow(1000.0, i.toDouble / (n - 1)))
  private val center = Array.tabulate(n)(i => Math.sin(i + .5))
  private val tolerance = FirstOrderTolerance.from(1e-7, 0.0).toOption.get
  private val term = ProximalTerms.l1(n, .1).toOption.get
  private val initial = DMat.zeros(n, 1)
  private val smooth = new SmoothObjective:
    val variableRows = n
    val lipschitz = 1000.0
    def value(x: DMat): Either[FirstOrderError, Double] =
      var value = 0.0
      var i = 0
      while i < n do
        val d = x(i, 0) - center(i)
        value += .5 * diagonal(i) * d * d
        i += 1
      Right(value)
    def gradient(x: DMat): Either[FirstOrderError, DMat] = Right(
      DMat.tabulate(n, 1)((r, _) => diagonal(r) * (x(r, 0) - center(r)))
    )
  private val fused = DifferentiableObjective.fromSmooth(smooth)
  def plain(): FirstOrderSolution = FirstOrderSolvers
    .proximalGradient(smooth, term, initial, FirstOrderConfig.from(100000, tolerance).toOption.get)
    .toOption
    .get
  def accelerated(): FirstOrderSolution = AcceleratedProximal
    .minimize(fused, term, initial, AcceleratedConfig(tolerance = tolerance, initialLipschitz = 1000.0))
    .toOption
    .get
  def verify(result: FirstOrderSolution): Double =
    var residual = 0.0
    var i = 0
    while i < n do
      val x = result.primal(i, 0)
      val gradient = diagonal(i) * (x - center(i))
      val kkt = if x == 0.0 then Math.max(0.0, Math.abs(gradient) - .1) else Math.abs(gradient + .1 * Math.signum(x))
      residual = Math.max(residual, kkt)
      i += 1
    require(residual <= 1e-6 && result.status == FirstOrderStoppingStatus.Converged)
    residual
  def main(args: Array[String]): Unit =
    println("method\tmedian_ms\tkkt\titerations\tcallbacks")
    Vector("fixed-step" -> (() => plain()), "accelerated" -> (() => accelerated())).foreach: (name, solve) =>
      var i = 0
      while i < 20 do
        solve()
        i += 1
      val samples = new Array[Double](11)
      var last = solve()
      i = 0
      while i < samples.length do
        val start = System.nanoTime()
        last = solve()
        samples(i) = (System.nanoTime() - start) / 1e6
        i += 1
      println(
        s"$name\t${samples.sorted.apply(5)}\t${verify(last)}\t${last.certificate.iterations}\t${last.evaluations.callbacks}"
      )
