package gale.bench

import gale.linalg.{DMat, DVec}
import gale.optim.*

/** Shared fixtures and prepared callbacks: parsing and data/model construction are outside timed solves. */
object OptimizationExtensionsBench:
  def runner(fixture: OptimizationBench.Fixture): () => FirstOrderSolution =
    val data = fixture.fields
    val start = data("X0")
    val n = start.length
    val initial = DVec.fromSeq(start.toIndexedSeq)
    val tolerance = FirstOrderTolerance.from(1e-8, 0.0).toOption.get
    if fixture.id.startsWith("lm-") then
      val objective = if fixture.id.startsWith("lm-exp") then
        val t = data("T")
        val y = data("Y")
        LeastSquaresObjective(n, t.length)(
          x =>
            Right(DVec.tabulate(t.length)(i => x(0) * Math.exp(-x(1) * t(i)) + (if n == 3 then x(2) else 0.0) - y(i))),
          x =>
            Right(
              DMat.tabulate(t.length, n)((r, c) =>
                if c == 0 then Math.exp(-x(1) * t(r))
                else if c == 1 then -x(0) * t(r) * Math.exp(-x(1) * t(r))
                else 1.0
              )
            )
        ).toOption.get
      else if fixture.id == "lm-rosenbrock" then
        LeastSquaresObjective(2, 2)(
          x => Right(DVec.fromSeq(Seq(10.0 * (x(1) - x(0) * x(0)), 1.0 - x(0)))),
          x => Right(DMat.dense(2, 2, Seq(-20.0 * x(0), 10.0, -1.0, 0.0)))
        ).toOption.get
      else
        val design = DMat.fromArrayRowMajor(data("Y").length, n, data("X"))
        val response = DVec.fromSeq(data("Y").toIndexedSeq)
        LeastSquaresObjective(n, response.length)(
          x =>
            val fitted = design * x
            Right(DVec.tabulate(response.length)(i => fitted(i) - response(i)))
          ,
          _ => Right(design)
        ).toOption.get
      () =>
        LevenbergMarquardt
          .minimize(
            objective,
            initial,
            LevenbergMarquardtConfig(
              maxIterations = 1000,
              tolerance = tolerance,
              control = SolverControl(maxEvaluations = 20000)
            )
          )
          .fold(e => throw new IllegalStateException(e.message), _.solution)
    else
      val objective =
        if fixture.id.startsWith("box-logistic") then
          DenseObjectives
            .logistic(
              DMat.fromArrayRowMajor(data("Y").length, n, data("X")),
              DVec.fromSeq(data("Y").toIndexedSeq),
              data("LAMBDA")(0)
            )
            .toOption
            .get
        else
          // The existing harness recognizes Rosenbrock by its id prefix.
          OptimizationBench.objective(
            if fixture.id.startsWith("box-rosenbrock") then fixture.copy(id = "rosenbrock-n2") else fixture
          )
      val bounds = MatrixBoxBounds
        .from(DMat.fromArrayRowMajor(n, 1, data("LOWER")), DMat.fromArrayRowMajor(n, 1, data("UPPER")))
        .toOption
        .get
      val at = DMat.fromArrayRowMajor(n, 1, start)
      () =>
        LBFGSB
          .minimize(
            objective,
            bounds,
            at,
            LBFGSBConfig(maxIterations = 10000, tolerance = tolerance, control = SolverControl(maxEvaluations = 50000))
          )
          .fold(e => throw new IllegalStateException(e.message), identity)

  def main(args: Array[String]): Unit =
    val warmups = args.headOption.map(_.toInt).getOrElse(100)
    val repeats = args.lift(1).map(_.toInt).getOrElse(11)
    require(warmups >= 0 && repeats > 0)
    println(
      "case_id\tstatus\tobjective\tresidual\titerations\tcallbacks\tresidual_calls\tjacobian_calls\tmedian_ms\ttimes_ms\tx"
    )
    OptimizationBench
      .fixtures(OptimizationInputs.read("docs/verification/optim-extensions/fixtures.tsv"))
      .foreach: fixture =>
        val run = runner(fixture)
        var i = 0
        while i < warmups do
          run()
          i += 1
        val times = new Array[Double](repeats)
        var result: Option[FirstOrderSolution] = None
        i = 0
        while i < repeats do
          val before = System.nanoTime()
          result = Some(run())
          times(i) = (System.nanoTime() - before) / 1e6
          i += 1
        val solution = result.get
        val sorted = times.sorted
        val median =
          if repeats % 2 == 0 then .5 * (sorted(repeats / 2 - 1) + sorted(repeats / 2)) else sorted(repeats / 2)
        val x = Array.tabulate(solution.primal.rows)(r => solution.primal(r, 0)).mkString(",")
        println(
          s"${fixture.id}\t${solution.status}\t${solution.objective}\t${solution.certificate.primalResidual}\t${solution.certificate.iterations}\t${solution.evaluations.callbacks}\t${solution.evaluations.residuals}\t${solution.evaluations.jacobians}\t$median\t${times.mkString(",")}\t$x"
        )
