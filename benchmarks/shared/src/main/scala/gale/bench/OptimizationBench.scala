package gale.bench

import gale.linalg.{DMat, DVec}
import gale.optim.*

/** Paired full-solve benchmark. Fixture I/O and oracle construction are outside timing. */
object OptimizationBench:
  final case class Fixture(id: String, fields: Map[String, Array[Double]])
  def fixtures(text: String): Vector[Fixture] =
    val out = Vector.newBuilder[Fixture]
    var id = ""
    var fields = Map.empty[String, Array[Double]]
    text.linesIterator.foreach: line =>
      val cells = line.split("\t")
      cells.headOption match
        case Some("CASE")                  => id = cells(1); fields = Map.empty
        case Some("END") if id.nonEmpty    => out += Fixture(id, fields); id = ""
        case Some("META")                  => ()
        case Some(key) if cells.length > 1 => fields += key -> cells.drop(1).map(_.toDouble)
        case _                             => ()
    out.result()

  def objective(fixture: Fixture): DifferentiableObjective =
    val data = fixture.fields
    val n = data("X0").length
    if data.contains("X") then
      val response = DVec.fromSeq(data("Y").toIndexedSeq)
      val design = DMat.fromArrayRowMajor(response.length, n, data("X"))
      return (if fixture.id.startsWith("logistic") then DenseObjectives.logistic(design, response, data("LAMBDA")(0))
              else DenseObjectives.leastSquares(design, response)).toOption.get
    DifferentiableObjective(n): point =>
      val x = new Array[Double](n)
      point.copyRowMajorTo(x)
      val g = new Array[Double](n)
      var f = 0.0
      var i = 0
      if data.contains("DIAGONAL") then
        val d = data("DIAGONAL"); val center = data("CENTER")
        while i < n do
          val delta = x(i) - center(i)
          g(i) = d(i) * delta
          f += .5 * delta * g(i)
          i += 1
      else if data.contains("HESSIAN") then
        val h = data("HESSIAN"); val center = data("CENTER")
        while i < n do
          var j = 0
          while j < n do
            g(i) += h(i * n + j) * (x(j) - center(j))
            j += 1
          f += .5 * (x(i) - center(i)) * g(i)
          i += 1
      else if fixture.id.startsWith("rosenbrock") then
        while i < n - 1 do
          val r = x(i + 1) - x(i) * x(i)
          f += 100 * r * r + (1 - x(i)) * (1 - x(i))
          g(i) += -400 * x(i) * r - 2 * (1 - x(i))
          g(i + 1) += 200 * r
          i += 1
      Right(ObjectiveEvaluation(f, DMat.tabulate(n, 1)((r, _) => g(r))))
    .toOption.get

  def solver(fixture: Fixture): () => FirstOrderSolution =
    val obj = objective(fixture)
    val x0 = fixture.fields("X0")
    val initial = DMat.tabulate(x0.length, 1)((r, _) => x0(r))
    val tolerance = FirstOrderTolerance.from(1e-7, 0.0).toOption.get
    val control = SolverControl(maxEvaluations = 250000)
    val run: () => Either[FirstOrderError, FirstOrderSolution] =
      if fixture.id.startsWith("lasso") then
        val term = ProximalTerms.l1(x0.length, fixture.fields("LAMBDA")(0)).toOption.get
        () =>
          AcceleratedProximal.minimize(
            obj,
            term,
            initial,
            AcceleratedConfig(maxIterations = 100000, tolerance = tolerance, control = control)
          )
      else if fixture.fields.contains("LOWER") then
        val lower = fixture.fields("LOWER"); val upper = fixture.fields("UPPER")
        val box = ProjectionSets
          .box(DMat.tabulate(x0.length, 1)((r, _) => lower(r)), DMat.tabulate(x0.length, 1)((r, _) => upper(r)))
          .toOption
          .get
        () =>
          AcceleratedProximal.minimize(
            obj,
            box.indicator,
            initial,
            AcceleratedConfig(maxIterations = 20000, tolerance = tolerance, control = control)
          )
      else
        () => LBFGS.minimize(obj, initial, LBFGSConfig(maxIterations = 50000, tolerance = tolerance, control = control))
    () => run().fold(error => throw new IllegalStateException(s"${fixture.id}: $error"), identity)

  def main(args: Array[String]): Unit =
    val warmups = args.headOption.map(_.toInt).getOrElse(100)
    val repeats = args.lift(1).map(_.toInt).getOrElse(11)
    require(warmups >= 0 && repeats > 0)
    println("case_id\tmethod\tstatus\tobjective\tresidual\titerations\tcallbacks\tmedian_ms\ttimes_ms\tx")
    fixtures(OptimizationInputs.read()).foreach: fixture =>
      val solve = solver(fixture)
      var i = 0
      while i < warmups do
        solve()
        i += 1
      val times = new Array[Double](repeats)
      var result: Option[FirstOrderSolution] = None
      i = 0
      while i < repeats do
        val start = System.nanoTime()
        val solved = solve()
        times(i) = (System.nanoTime() - start) / 1e6
        result = Some(solved)
        i += 1
      val s = result.get
      val ordered = times.sorted
      val median =
        if repeats % 2 == 1 then ordered(repeats / 2) else .5 * (ordered(repeats / 2 - 1) + ordered(repeats / 2))
      val coordinates = Array.tabulate(s.primal.rows)(r => s.primal(r, 0)).mkString(",")
      println(
        s"${fixture.id}\t${s.certificate.settings.method}\t${s.status}\t${s.objective}\t${s.certificate.primalResidual}\t${s.certificate.iterations}\t${s.evaluations.callbacks}\t$median\t${times.mkString(",")}\t$coordinates"
      )
