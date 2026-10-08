package gale.bench

import gale.linalg.DVec
import gale.optim.*

/** Full solves, including solver construction and callbacks; fixture preparation is outside timing. */
object CMAESBench:
  final case class Fixture(id: String, n: Int, sigma: Double, budget: Int, target: Double, initial: DVec)
  def fixtures: Vector[Fixture] = OptimizationInputs
    .read("docs/verification/optim-cmaes/fixtures.tsv")
    .linesIterator
    .drop(1)
    .filter(_.nonEmpty)
    .map { line =>
      val p = line.split("\t")
      Fixture(
        p(0),
        p(1).toInt,
        p(2).toDouble,
        p(3).toInt,
        p(4).toDouble,
        DVec.fromSeq(p(5).split(",").toSeq.map(_.toDouble))
      )
    }
    .toVector

  def clean(f: Fixture, x: DVec): Double =
    if f.id == "rosenbrock-4" then
      var sum = 0.0
      var i = 0
      while i < f.n - 1 do
        val residual = x(i + 1) - x(i) * x(i)
        sum += 100.0 * residual * residual + (1.0 - x(i)) * (1.0 - x(i))
        i += 1
      sum
    else
      var projection = 0.0
      var denominator = 0.0
      if f.id == "rotated-ellipsoid-8" then
        var i = 0
        while i < f.n do
          projection += (i + 1) * x(i)
          denominator += (i + 1.0) * (i + 1.0)
          i += 1
      var sum = 0.0
      var i = 0
      while i < f.n do
        val v = x(i)
        sum += (f.id match
          case "rotated-ellipsoid-8" =>
            val y = v - 2.0 * (i + 1) * projection / denominator
            math.pow(1e4, i.toDouble / (f.n - 1)) * y * y
          case "nonsmooth-8"  => (i + 1.0) * math.abs(v)
          case "rastrigin-4"  => v * v + 10.0 * (1.0 - math.cos(2.0 * math.Pi * v))
          case "box-sphere-8" =>
            val d = v - (if i % 2 == 0 then 0.7 else -0.6)
            d * d
          case "box-boundary-4" => (v - 2.0) * (v - 2.0)
          case _                => v * v)
        i += 1
      sum

  def run(f: Fixture, seed: Int): (CMAESResult, Double) =
    val rng = new scala.util.Random(seed.toLong ^ 0x123456789abcdefL)
    val objective: DVec => Double =
      if f.id == "noisy-sphere-8" then x => clean(f, x) * math.exp(0.1 * rng.nextGaussian() - 0.005)
      else x => clean(f, x)
    val bounds = if f.id.startsWith("box-") then Some(MatrixBoxBounds.uniform(f.n, -1.0, 1.0).toOption.get) else None
    val config = CMAESConfig(
      seed = seed.toLong,
      initialStepSize = f.sigma,
      maxEvaluations = f.budget,
      maxGenerations = 10000,
      targetValue = Some(f.target),
      stepTolerance = 0.0,
      stagnationGenerations = 0
    )
    val begin = System.nanoTime()
    val result = CMAES.minimize(objective, f.initial, config, bounds).fold(e => sys.error(e.toString), identity)
    (result, (System.nanoTime() - begin) / 1e6)

  def references(): Unit =
    for seed <- Vector(7L, 42L) do
      val s = CMAES
        .start(
          DVec.fromSeq(Seq(1.0, -2.0, 0.5)),
          CMAESConfig(
            seed = seed,
            initialStepSize = 0.3,
            populationSize = 6,
            maxGenerations = 6,
            stepTolerance = 0.0,
            stagnationGenerations = 0,
            eigenUpdatePeriod = 1
          )
        )
        .toOption
        .get
      for generation <- 1 to 6 do
        val batch = s.ask().toOption.get.get
        val values = batch.points.map(x => x(0) * x(0) + 2.0 * x(1) * x(1) + 4.0 * x(2) * x(2))
        s.tell(batch, values).toOption.get
        val d = s.snapshot.distribution
        val c = (0 until 3).flatMap(i => (0 until 3).map(j => d.covariance(i, j)))
        println(
          Vector(
            "REF",
            seed.toString,
            generation.toString,
            batch.points.map(_.toSeq.mkString(",")).mkString(";"),
            values.mkString(","),
            d.mean.toSeq.mkString(","),
            d.stepSize.toString,
            c.mkString(",")
          ).mkString("\t")
        )

  def main(args: Array[String]): Unit =
    if args.headOption.contains("reference") then references()
    else
      val warmups = args.headOption.map(_.toInt).getOrElse(5)
      val seeds = args.lift(1).map(_.toInt).getOrElse(20)
      val repeats = args.lift(2).map(_.toInt).getOrElse(3)
      println("case_id\tseed\trepeat\tobserved\tclean\tstatus\tevaluations\tgenerations\tsampled\trejected\tms\tx")
      for f <- fixtures do
        for i <- 0 until warmups do run(f, 10000 + i)
        for seed <- 1 to seeds; repeat <- 0 until repeats do
          val (result, milliseconds) = run(f, seed)
          val best = result.best.getOrElse(sys.error("fixture produced no finite observation"))
          println(
            Vector(
              f.id,
              seed.toString,
              repeat.toString,
              best.value.toString,
              clean(f, best.coordinates).toString,
              result.stopReason.toString,
              result.work.evaluations.toString,
              result.work.generations.toString,
              result.work.sampled.toString,
              result.work.rejected.toString,
              milliseconds.toString,
              best.coordinates.toSeq.mkString(",")
            ).mkString("\t")
          )
      references()
