package gale.bench

import gale.linalg.{DMat, DVec}
import gale.optim.*

/** Prepared analytic callbacks; timing includes the entire solve and its final diagnostics. */
object OptimizationALMBench:
  def runner(
      fixture: OptimizationBench.Fixture,
      config: AugmentedLagrangianConfig = AugmentedLagrangianConfig(control = SolverControl(maxEvaluations = 50000))
  ): () => AugmentedLagrangianResult =
    val d = fixture.fields
    val n = d("X0").length
    val initial = DMat.fromArrayRowMajor(n, 1, d("X0"))
    val q = DMat.fromArrayRowMajor(n, n, d("Q"))
    val center = d("CENTER")
    val objective = DifferentiableObjective(n): x =>
      val delta = DVec.tabulate(n)(i => x(i, 0) - center(i))
      val gradient = q * delta
      var value = 0.0
      var i = 0
      while i < n do
        value += delta(i) * gradient(i)
        i += 1
      Right(ObjectiveEvaluation(0.5 * value, DMat.tabulate(n, 1)((i, _) => gradient(i))))
    val eq = d("EQUALITIES")(0).toInt
    val ineq = d("INEQUALITIES")(0).toInt
    val kind = d("KIND")(0).toInt
    val a = d.get("A").map(v => DMat.fromArrayRowMajor(eq + ineq, n, v))
    val constraints = NonlinearConstraints(n, eq, ineq): x =>
      if kind == 0 then
        val image = a.get * x.col(0)
        Right(ConstraintEvaluation(DVec.tabulate(eq + ineq)(i => image(i) - d("B")(i)), a.get))
      else if kind == 1 then
        var squared = 0.0
        var i = 0
        while i < n do
          squared += x(i, 0) * x(i, 0)
          i += 1
        Right(ConstraintEvaluation(DVec.fromSeq(Seq(squared - 1.0)), DMat.tabulate(1, n)((_, j) => 2 * x(j, 0))))
      else
        Right(
          ConstraintEvaluation(
            DVec.tabulate(eq)(i => x(2 * i + 1, 0) - x(2 * i, 0) * x(2 * i, 0)),
            DMat.tabulate(eq, n)((i, j) => if j == 2 * i then -2 * x(j, 0) else if j == 2 * i + 1 then 1.0 else 0.0)
          )
        )
    val bounds = d
      .get("LOWER")
      .map(l =>
        MatrixBoxBounds.from(DMat.fromArrayRowMajor(n, 1, l), DMat.fromArrayRowMajor(n, 1, d("UPPER"))).toOption.get
      )
    val f = objective.toOption.get
    val c = constraints.toOption.get
    () => AugmentedLagrangian.minimize(f, c, initial, bounds, config).fold(e => sys.error(e.message), identity)

  def main(args: Array[String]): Unit =
    val warmups = args.headOption.map(_.toInt).getOrElse(20)
    val repeats = args.lift(1).map(_.toInt).getOrElse(11)
    val path = args.lift(2).getOrElse("docs/verification/optim-alm/fixtures.tsv")
    val runtime = args.lift(3).getOrElse("unknown")
    require(warmups >= 0 && repeats > 0)
    val variant = args.lift(4).getOrElse("default")
    require(Set("default", "fixed", "adaptive", "scale", "combined").contains(variant))
    val config = AugmentedLagrangianConfig(
      control = SolverControl(maxEvaluations = 50000, traceCapacity = if args.contains("trace") then 100 else 0),
      adaptiveInnerTolerance = variant == "default" || variant == "adaptive" || variant == "combined",
      reuseCurvatureScale = variant == "default" || variant == "scale" || variant == "combined"
    )
    val prepared = OptimizationBench.fixtures(OptimizationInputs.read(path)).map(f => f -> runner(f, config))
    // Warm every family before measuring any family, so early cases also see compiled shared solver paths.
    for (_, run) <- prepared do for _ <- 0 until warmups do run()
    for (fixture, run) <- prepared do
      for repeat <- 0 until repeats do
        val before = System.nanoTime()
        val r = run()
        val milliseconds = (System.nanoTime() - before) / 1e6
        println(
          Vector(
            "ALM",
            runtime,
            fixture.id,
            repeat.toString,
            milliseconds.toString,
            r.status.toString,
            r.objective.toString,
            r.iterations.toString,
            r.innerIterations.toString,
            r.evaluations.callbacks.toString,
            r.evaluations.values.toString,
            r.evaluations.jacobians.toString,
            r.diagnostics.feasibility.toString,
            r.diagnostics.scaledFeasibility.toString,
            r.diagnostics.stationarity.toString,
            r.diagnostics.complementarity.toString,
            r.penalty.toString,
            r.primal.col(0).toSeq.mkString(","),
            r.multipliers.toSeq.mkString(",")
          ).mkString("\t")
        )
        r.innerTrace.foreach(t => println(s"INNER\t${fixture.id}\t$variant\t$t"))
