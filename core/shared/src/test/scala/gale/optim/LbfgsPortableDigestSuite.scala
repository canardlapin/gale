package gale.optim

import gale.linalg.*

/** Cross-platform bit guard for [[LBFGS]] and [[LBFGSB]]: the same digests must hold on the JVM and on Scala.js.
  *
  * The corpus uses only `+ − × ÷` in its objectives and exactly representable starting points (a linear congruential
  * sequence scaled by a power of two), so the solvers' own arithmetic, which uses no transcendental functions, is the
  * only thing hashed. Each case is reduced, as in the JVM-only `LbfgsDigestSuite`, to a 64-bit FNV-1a digest of every
  * point handed to the objective, the returned point and fields, the evaluation counts and the progress trace. The
  * expected digests were generated on the JVM from breeze/integration 5d3950e (`LbfgsPortableDigests.main`). On
  * Scala.js this also covers the `Float64Array`-backed storage the solvers allocate.
  */
class LbfgsPortableDigestSuite extends munit.FunSuite:
  import LbfgsPortableDigests.*

  test("LBFGS and LBFGSB solves match the stored per-case bit digests on every platform") {
    val actual = digests
    assertEquals(actual.map(_._1).distinct.size, actual.size, "duplicate case keys")
    assertEquals(Expected.size, ExpectedEntries.size, "duplicate expected keys")
    assert(actual.size > 100, s"corpus shrank to ${actual.size} cases")
    assertEquals(actual.map(_._1).toSet, Expected.keySet)
    val mismatches = actual.filter((key, digest) => !Expected.get(key).contains(digest)).map(_._1).sorted
    assert(mismatches.isEmpty, s"${mismatches.size} cases changed bits, e.g. ${mismatches.take(5).mkString("; ")}")
  }

object LbfgsPortableDigests:
  /** Print the expected-digest entries for this file: `coreJVM/Test/runMain gale.optim.LbfgsPortableDigests`. */
  def main(args: Array[String]): Unit =
    digests.sortBy(_._1).foreach((key, digest) => println(s"""    "$key" -> "$digest","""))

  /** FNV-1a over 64-bit words, each word's bytes folded high to low. */
  final class Digest:
    private var state = 0xcbf29ce484222325L
    def long(x: Long): Digest =
      var shift = 56
      while shift >= 0 do
        state = (state ^ ((x >>> shift) & 0xffL)) * 0x100000001b3L
        shift -= 8
      this
    def double(x: Double): Digest = long(java.lang.Double.doubleToRawLongBits(x))
    def string(s: String): Digest =
      s.foreach(ch => long(ch.toLong))
      this
    def mat(m: DMat): Digest =
      long(m.rows.toLong).long(m.cols.toLong)
      var i = 0
      while i < m.rows do
        var j = 0
        while j < m.cols do
          double(m(i, j))
          j += 1
        i += 1
      this
    def counts(c: EvaluationCounts): Digest =
      c.productIterator.foreach(value => long(value.asInstanceOf[Int].toLong))
      this
    def trace(events: Vector[OptimizationProgress]): Digest =
      long(events.size.toLong)
      events.foreach(e => long(e.iteration.toLong).double(e.objective).double(e.residual).counts(e.evaluations))
      this
    def hex: String =
      val digits = java.lang.Long.toHexString(state)
      "0" * (16 - digits.length) + digits

  /** Chained Rosenbrock over all `x.length` variables. */
  def rosenbrock(x: Array[Double], g: Array[Double]): Double =
    var f = 0.0
    var i = 0
    while i < x.length - 1 do
      val t = x(i + 1) - x(i) * x(i)
      val u = 1.0 - x(i)
      f += 100.0 * t * t + u * u
      g(i) += -400.0 * t * x(i) - 2.0 * u
      g(i + 1) += 200.0 * t
      i += 1
    f

  /** Separable quadratic with exact curvatures `(1 + i mod 3) · 10^(i mod 4)` and centres `i mod 3 − 1`. */
  def quadratic(x: Array[Double], g: Array[Double]): Double =
    var f = 0.0
    var i = 0
    while i < x.length do
      var d = 1.0 + i % 3
      var k = 0
      while k < i % 4 do
        d *= 10.0
        k += 1
      val r = x(i) - (i % 3 - 1).toDouble
      f += 0.5 * d * r * r
      g(i) = d * r
      i += 1
    f

  /** Two-variable Rosenbrock plus a third variable that is inert until `x(0) >= gate` and then follows
    * `slope · v + v²/20`. Once a correction pair has seen the third variable's small curvature, `gᵀHg` overflows
    * although `gᵀg` does not, so the two-loop direction is rejected, the history is cleared (after the ring has wrapped
    * when the history is short) and the solve restarts from steepest descent.
    */
  def gatedRestart(gate: Double, slope: Double)(x: Array[Double], g: Array[Double]): Double =
    val t = x(1) - x(0) * x(0)
    val u = 1.0 - x(0)
    g(0) = -400.0 * t * x(0) - 2.0 * u
    g(1) = 200.0 * t
    val open = x(0) >= gate
    g(2) = if open then slope + 0.1 * x(2) else 0.0
    100.0 * t * t + u * u + (if open then slope * x(2) + 0.05 * x(2) * x(2) else 0.0)

  /** Exactly representable pseudo-random values in `[-scale/2, scale/2)` (`scale` a power of two). */
  def lcg(n: Int, seed: Long, scale: Double): Array[Double] =
    var state = seed
    Array.fill(n):
      state = (state * 1103515245L + 12345L) & 0x7fffffffL
      (state.toDouble / 2147483648.0 - 0.5) * scale

  /** Variable layouts: plain owned columns and strided transposed views. */
  enum Layout:
    case Column, Transposed

  final case class Problem(
      name: String,
      rows: Int,
      cols: Int,
      start: Array[Double],
      value: (Array[Double], Array[Double]) => Double
  )

  def matrixOf(values: Array[Double], rows: Int, cols: Int, layout: Layout): DMat =
    layout match
      case Layout.Transposed => Matrix.tabulate(cols, rows)((c, r) => values(r * cols + c)).t
      case Layout.Column     => Matrix.tabulate(rows, cols)((r, c) => values(r * cols + c))

  /** The objective hashes every point it is handed, in call order, into `trials`. */
  def objective(problem: Problem, layout: Layout, trials: Digest): DifferentiableObjective =
    new DifferentiableObjective:
      def variableRows: Int = problem.rows
      def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
        trials.mat(at)
        val x = new Array[Double](at.rows * at.cols)
        at.copyRowMajorTo(x)
        val g = new Array[Double](x.length)
        val f = problem.value(x, g)
        Right(ObjectiveEvaluation(f, matrixOf(g, at.rows, at.cols, layout)))

  /** Errors contribute only their case name: `Double.toString` differs between the JVM and Scala.js. */
  def solutionDigest(trials: Digest, result: Either[FirstOrderError, FirstOrderSolution]): String =
    result match
      case Left(error) => trials.long(2L).string(error.productPrefix).hex
      case Right(s)    =>
        trials
          .long(1L)
          .mat(s.primal)
          .double(s.objective)
          .double(s.certificate.primalResidual)
          .double(s.primalStep)
          .string(s.status.toString)
          .long(s.certificate.iterations.toLong)
          .counts(s.evaluations)
          .trace(s.trace)
          .hex

  private val traced = SolverControl(traceCapacity = 100000)

  def problems: Seq[Problem] =
    Seq(2, 10, 100).flatMap(n =>
      Seq(
        Problem(s"rosenbrock n=$n classic", n, 1, Array.tabulate(n)(i => if i % 2 == 0 then -1.2 else 1.0), rosenbrock),
        Problem(s"rosenbrock n=$n lcg", n, 1, lcg(n, 11L + n, 2.0), rosenbrock)
      )
    ) ++ Seq(
      Problem("rosenbrock 4x3", 4, 3, lcg(12, 5L, 1.0), rosenbrock),
      Problem("quadratic n=7", 7, 1, lcg(7, 21L, 8.0), quadratic),
      Problem("quadratic n=40", 40, 1, lcg(40, 22L, 8.0), quadratic),
      Problem("quadratic 6x4", 6, 4, lcg(24, 23L, 4.0), quadratic)
    )

  private def lbfgsConfigs: Seq[(String, LBFGSConfig)] =
    Seq(
      "default" -> LBFGSConfig(control = traced),
      "iterations=3 m=1" -> LBFGSConfig(maxIterations = 3, historySize = 1, control = traced),
      "iterations=25 m=3" -> LBFGSConfig(maxIterations = 25, historySize = 3, control = traced),
      "iterations=40 m=10 untolerant" -> LBFGSConfig(
        maxIterations = 40,
        tolerance = FirstOrderTolerance.from(Double.MinPositiveValue, 0.0).toOption.get,
        control = traced
      ),
      "c2=0.1 linesearch=3" -> LBFGSConfig(c2 = 0.1, maxLineSearch = 3, control = traced),
      "maxstep=0.5" -> LBFGSConfig(maxIterations = 200, maximumStep = 0.5, control = traced),
      "evaluations=17" -> LBFGSConfig(control = SolverControl(maxEvaluations = 17, traceCapacity = 1000))
    )

  /** Restart cases: `(name, gate, slope, history)`, each reaching the two-loop fallback (see [[gatedRestart]]). */
  val restartCases: Seq[(String, Double, Double, Int)] = Seq(
    ("gate=0 slope=6e153 m=3", 0.0, 6e153, 3),
    ("gate=0 slope=5.5e153 m=3", 0.0, 5.5e153, 3),
    ("gate=0.5 slope=5.5e153 m=2", 0.5, 5.5e153, 2)
  )

  def restartProblem(name: String, gate: Double, slope: Double): Problem =
    Problem(s"restart $name", 3, 1, Array(-1.2, 1.0, 0.0), gatedRestart(gate, slope))

  def digests: Seq[(String, String)] =
    val lbfgs =
      for
        problem <- problems
        layout <- Layout.values.toSeq
        (configName, config) <- lbfgsConfigs
      yield
        val trials = new Digest
        val start = matrixOf(problem.start, problem.rows, problem.cols, layout)
        s"lbfgs ${problem.name} layout=$layout $configName" ->
          solutionDigest(trials, LBFGS.minimize(objective(problem, layout, trials), start, config))
    val restarts =
      for
        (name, gate, slope, history) <- restartCases
        layout <- Layout.values.toSeq
      yield
        val trials = new Digest
        val problem = restartProblem(name, gate, slope)
        val config = LBFGSConfig(maxIterations = 60, historySize = history, control = traced)
        val start = matrixOf(problem.start, problem.rows, problem.cols, layout)
        s"lbfgs ${problem.name} layout=$layout" ->
          solutionDigest(trials, LBFGS.minimize(objective(problem, layout, trials), start, config))
    val lbfgsb =
      for
        problem <- problems if problem.rows * problem.cols <= 24
        layout <- Layout.values.toSeq
        (configName, config) <- Seq(
          "default" -> LBFGSBConfig(control = traced),
          "iterations=15 m=2" -> LBFGSBConfig(maxIterations = 15, historySize = 2, control = traced)
        )
      yield
        val trials = new Digest
        val bounds = MatrixBoxBounds.uniform(problem.rows, -0.5, 0.75, problem.cols).toOption.get
        val start =
          matrixOf(problem.start.map(v => Math.max(-0.5, Math.min(0.75, v))), problem.rows, problem.cols, layout)
        s"lbfgsb ${problem.name} box=[-0.5,0.75] layout=$layout $configName" ->
          solutionDigest(trials, LBFGSB.minimize(objective(problem, layout, trials), bounds, start, config))
    lbfgs ++ restarts ++ lbfgsb

  /** Generated on the JVM from breeze/integration 5d3950e. */
  val ExpectedEntries: Seq[(String, String)] = Seq(
    "lbfgs quadratic 6x4 layout=Column c2=0.1 linesearch=3" -> "d883f2fcf2eb0e93",
    "lbfgs quadratic 6x4 layout=Column default" -> "ee3ec5384a4ac37e",
    "lbfgs quadratic 6x4 layout=Column evaluations=17" -> "571b0449b670648c",
    "lbfgs quadratic 6x4 layout=Column iterations=25 m=3" -> "9ea58fdc36a37781",
    "lbfgs quadratic 6x4 layout=Column iterations=3 m=1" -> "a74320bc6384c502",
    "lbfgs quadratic 6x4 layout=Column iterations=40 m=10 untolerant" -> "4dcb28a4509aff63",
    "lbfgs quadratic 6x4 layout=Column maxstep=0.5" -> "18e414adec684671",
    "lbfgs quadratic 6x4 layout=Transposed c2=0.1 linesearch=3" -> "d883f2fcf2eb0e93",
    "lbfgs quadratic 6x4 layout=Transposed default" -> "ee3ec5384a4ac37e",
    "lbfgs quadratic 6x4 layout=Transposed evaluations=17" -> "571b0449b670648c",
    "lbfgs quadratic 6x4 layout=Transposed iterations=25 m=3" -> "9ea58fdc36a37781",
    "lbfgs quadratic 6x4 layout=Transposed iterations=3 m=1" -> "a74320bc6384c502",
    "lbfgs quadratic 6x4 layout=Transposed iterations=40 m=10 untolerant" -> "4dcb28a4509aff63",
    "lbfgs quadratic 6x4 layout=Transposed maxstep=0.5" -> "18e414adec684671",
    "lbfgs quadratic n=40 layout=Column c2=0.1 linesearch=3" -> "5748a9a401304e0f",
    "lbfgs quadratic n=40 layout=Column default" -> "78864f7ea0a7e166",
    "lbfgs quadratic n=40 layout=Column evaluations=17" -> "12e10bd834cf73ab",
    "lbfgs quadratic n=40 layout=Column iterations=25 m=3" -> "ee9f5cf38bc84868",
    "lbfgs quadratic n=40 layout=Column iterations=3 m=1" -> "17c7d852ba8b608e",
    "lbfgs quadratic n=40 layout=Column iterations=40 m=10 untolerant" -> "d6f2d8e795e41114",
    "lbfgs quadratic n=40 layout=Column maxstep=0.5" -> "54e086049f71bade",
    "lbfgs quadratic n=40 layout=Transposed c2=0.1 linesearch=3" -> "5748a9a401304e0f",
    "lbfgs quadratic n=40 layout=Transposed default" -> "78864f7ea0a7e166",
    "lbfgs quadratic n=40 layout=Transposed evaluations=17" -> "12e10bd834cf73ab",
    "lbfgs quadratic n=40 layout=Transposed iterations=25 m=3" -> "ee9f5cf38bc84868",
    "lbfgs quadratic n=40 layout=Transposed iterations=3 m=1" -> "17c7d852ba8b608e",
    "lbfgs quadratic n=40 layout=Transposed iterations=40 m=10 untolerant" -> "d6f2d8e795e41114",
    "lbfgs quadratic n=40 layout=Transposed maxstep=0.5" -> "54e086049f71bade",
    "lbfgs quadratic n=7 layout=Column c2=0.1 linesearch=3" -> "32af3782794ae024",
    "lbfgs quadratic n=7 layout=Column default" -> "8c8283bbe00936d2",
    "lbfgs quadratic n=7 layout=Column evaluations=17" -> "cd33069c90f59243",
    "lbfgs quadratic n=7 layout=Column iterations=25 m=3" -> "0ee884a1400bf9e4",
    "lbfgs quadratic n=7 layout=Column iterations=3 m=1" -> "35c7a2b10ab8a45d",
    "lbfgs quadratic n=7 layout=Column iterations=40 m=10 untolerant" -> "380cc60901feac31",
    "lbfgs quadratic n=7 layout=Column maxstep=0.5" -> "91cd0286eb40ca06",
    "lbfgs quadratic n=7 layout=Transposed c2=0.1 linesearch=3" -> "32af3782794ae024",
    "lbfgs quadratic n=7 layout=Transposed default" -> "8c8283bbe00936d2",
    "lbfgs quadratic n=7 layout=Transposed evaluations=17" -> "cd33069c90f59243",
    "lbfgs quadratic n=7 layout=Transposed iterations=25 m=3" -> "0ee884a1400bf9e4",
    "lbfgs quadratic n=7 layout=Transposed iterations=3 m=1" -> "35c7a2b10ab8a45d",
    "lbfgs quadratic n=7 layout=Transposed iterations=40 m=10 untolerant" -> "380cc60901feac31",
    "lbfgs quadratic n=7 layout=Transposed maxstep=0.5" -> "91cd0286eb40ca06",
    "lbfgs restart gate=0 slope=5.5e153 m=3 layout=Column" -> "2883f872c5976f74",
    "lbfgs restart gate=0 slope=5.5e153 m=3 layout=Transposed" -> "2883f872c5976f74",
    "lbfgs restart gate=0 slope=6e153 m=3 layout=Column" -> "d2e88bb41df4aac5",
    "lbfgs restart gate=0 slope=6e153 m=3 layout=Transposed" -> "d2e88bb41df4aac5",
    "lbfgs restart gate=0.5 slope=5.5e153 m=2 layout=Column" -> "ab88abf88aaecd38",
    "lbfgs restart gate=0.5 slope=5.5e153 m=2 layout=Transposed" -> "ab88abf88aaecd38",
    "lbfgs rosenbrock 4x3 layout=Column c2=0.1 linesearch=3" -> "d2b52a5238d90d3e",
    "lbfgs rosenbrock 4x3 layout=Column default" -> "868a6e4e65530f55",
    "lbfgs rosenbrock 4x3 layout=Column evaluations=17" -> "777bac528ef04ba1",
    "lbfgs rosenbrock 4x3 layout=Column iterations=25 m=3" -> "387f39baba4ac3c4",
    "lbfgs rosenbrock 4x3 layout=Column iterations=3 m=1" -> "3f1087e67d778de4",
    "lbfgs rosenbrock 4x3 layout=Column iterations=40 m=10 untolerant" -> "8efb4bcb989b5012",
    "lbfgs rosenbrock 4x3 layout=Column maxstep=0.5" -> "0b747cf5c7a85cb5",
    "lbfgs rosenbrock 4x3 layout=Transposed c2=0.1 linesearch=3" -> "d2b52a5238d90d3e",
    "lbfgs rosenbrock 4x3 layout=Transposed default" -> "868a6e4e65530f55",
    "lbfgs rosenbrock 4x3 layout=Transposed evaluations=17" -> "777bac528ef04ba1",
    "lbfgs rosenbrock 4x3 layout=Transposed iterations=25 m=3" -> "387f39baba4ac3c4",
    "lbfgs rosenbrock 4x3 layout=Transposed iterations=3 m=1" -> "3f1087e67d778de4",
    "lbfgs rosenbrock 4x3 layout=Transposed iterations=40 m=10 untolerant" -> "8efb4bcb989b5012",
    "lbfgs rosenbrock 4x3 layout=Transposed maxstep=0.5" -> "0b747cf5c7a85cb5",
    "lbfgs rosenbrock n=10 classic layout=Column c2=0.1 linesearch=3" -> "6a0bb09cac28b0b0",
    "lbfgs rosenbrock n=10 classic layout=Column default" -> "e5c228c5c9b7fa68",
    "lbfgs rosenbrock n=10 classic layout=Column evaluations=17" -> "e0b96a5cc6eebb82",
    "lbfgs rosenbrock n=10 classic layout=Column iterations=25 m=3" -> "3e80535ec00e6481",
    "lbfgs rosenbrock n=10 classic layout=Column iterations=3 m=1" -> "836bbcefcae05b1f",
    "lbfgs rosenbrock n=10 classic layout=Column iterations=40 m=10 untolerant" -> "4596d00c127d5000",
    "lbfgs rosenbrock n=10 classic layout=Column maxstep=0.5" -> "eae77c9baf912390",
    "lbfgs rosenbrock n=10 classic layout=Transposed c2=0.1 linesearch=3" -> "6a0bb09cac28b0b0",
    "lbfgs rosenbrock n=10 classic layout=Transposed default" -> "e5c228c5c9b7fa68",
    "lbfgs rosenbrock n=10 classic layout=Transposed evaluations=17" -> "e0b96a5cc6eebb82",
    "lbfgs rosenbrock n=10 classic layout=Transposed iterations=25 m=3" -> "3e80535ec00e6481",
    "lbfgs rosenbrock n=10 classic layout=Transposed iterations=3 m=1" -> "836bbcefcae05b1f",
    "lbfgs rosenbrock n=10 classic layout=Transposed iterations=40 m=10 untolerant" -> "4596d00c127d5000",
    "lbfgs rosenbrock n=10 classic layout=Transposed maxstep=0.5" -> "eae77c9baf912390",
    "lbfgs rosenbrock n=10 lcg layout=Column c2=0.1 linesearch=3" -> "c1ae6f28007b0c2f",
    "lbfgs rosenbrock n=10 lcg layout=Column default" -> "051f63950b3a9bf7",
    "lbfgs rosenbrock n=10 lcg layout=Column evaluations=17" -> "dd0780e96be8fae5",
    "lbfgs rosenbrock n=10 lcg layout=Column iterations=25 m=3" -> "bcea73586baa57dc",
    "lbfgs rosenbrock n=10 lcg layout=Column iterations=3 m=1" -> "e98f168b069a28af",
    "lbfgs rosenbrock n=10 lcg layout=Column iterations=40 m=10 untolerant" -> "63d4308e9992e45d",
    "lbfgs rosenbrock n=10 lcg layout=Column maxstep=0.5" -> "37553aa04e713a9a",
    "lbfgs rosenbrock n=10 lcg layout=Transposed c2=0.1 linesearch=3" -> "c1ae6f28007b0c2f",
    "lbfgs rosenbrock n=10 lcg layout=Transposed default" -> "051f63950b3a9bf7",
    "lbfgs rosenbrock n=10 lcg layout=Transposed evaluations=17" -> "dd0780e96be8fae5",
    "lbfgs rosenbrock n=10 lcg layout=Transposed iterations=25 m=3" -> "bcea73586baa57dc",
    "lbfgs rosenbrock n=10 lcg layout=Transposed iterations=3 m=1" -> "e98f168b069a28af",
    "lbfgs rosenbrock n=10 lcg layout=Transposed iterations=40 m=10 untolerant" -> "63d4308e9992e45d",
    "lbfgs rosenbrock n=10 lcg layout=Transposed maxstep=0.5" -> "37553aa04e713a9a",
    "lbfgs rosenbrock n=100 classic layout=Column c2=0.1 linesearch=3" -> "1ff797c5bfdbd023",
    "lbfgs rosenbrock n=100 classic layout=Column default" -> "2fb720f77fd55d6e",
    "lbfgs rosenbrock n=100 classic layout=Column evaluations=17" -> "f6b6ac93c06933c3",
    "lbfgs rosenbrock n=100 classic layout=Column iterations=25 m=3" -> "3e2da70778403b65",
    "lbfgs rosenbrock n=100 classic layout=Column iterations=3 m=1" -> "59696f339478d56e",
    "lbfgs rosenbrock n=100 classic layout=Column iterations=40 m=10 untolerant" -> "625d48824e04fd8c",
    "lbfgs rosenbrock n=100 classic layout=Column maxstep=0.5" -> "5570cc608f594fc3",
    "lbfgs rosenbrock n=100 classic layout=Transposed c2=0.1 linesearch=3" -> "1ff797c5bfdbd023",
    "lbfgs rosenbrock n=100 classic layout=Transposed default" -> "2fb720f77fd55d6e",
    "lbfgs rosenbrock n=100 classic layout=Transposed evaluations=17" -> "f6b6ac93c06933c3",
    "lbfgs rosenbrock n=100 classic layout=Transposed iterations=25 m=3" -> "3e2da70778403b65",
    "lbfgs rosenbrock n=100 classic layout=Transposed iterations=3 m=1" -> "59696f339478d56e",
    "lbfgs rosenbrock n=100 classic layout=Transposed iterations=40 m=10 untolerant" -> "625d48824e04fd8c",
    "lbfgs rosenbrock n=100 classic layout=Transposed maxstep=0.5" -> "5570cc608f594fc3",
    "lbfgs rosenbrock n=100 lcg layout=Column c2=0.1 linesearch=3" -> "8ed1266d3ead9e26",
    "lbfgs rosenbrock n=100 lcg layout=Column default" -> "1c9d0354d767a0f5",
    "lbfgs rosenbrock n=100 lcg layout=Column evaluations=17" -> "a8b43d1e6e346e58",
    "lbfgs rosenbrock n=100 lcg layout=Column iterations=25 m=3" -> "398040328b54d2b0",
    "lbfgs rosenbrock n=100 lcg layout=Column iterations=3 m=1" -> "4a4c299ca63cfada",
    "lbfgs rosenbrock n=100 lcg layout=Column iterations=40 m=10 untolerant" -> "f91d067604e65701",
    "lbfgs rosenbrock n=100 lcg layout=Column maxstep=0.5" -> "91e736ad42668394",
    "lbfgs rosenbrock n=100 lcg layout=Transposed c2=0.1 linesearch=3" -> "8ed1266d3ead9e26",
    "lbfgs rosenbrock n=100 lcg layout=Transposed default" -> "1c9d0354d767a0f5",
    "lbfgs rosenbrock n=100 lcg layout=Transposed evaluations=17" -> "a8b43d1e6e346e58",
    "lbfgs rosenbrock n=100 lcg layout=Transposed iterations=25 m=3" -> "398040328b54d2b0",
    "lbfgs rosenbrock n=100 lcg layout=Transposed iterations=3 m=1" -> "4a4c299ca63cfada",
    "lbfgs rosenbrock n=100 lcg layout=Transposed iterations=40 m=10 untolerant" -> "f91d067604e65701",
    "lbfgs rosenbrock n=100 lcg layout=Transposed maxstep=0.5" -> "91e736ad42668394",
    "lbfgs rosenbrock n=2 classic layout=Column c2=0.1 linesearch=3" -> "b5fffc88e8c94440",
    "lbfgs rosenbrock n=2 classic layout=Column default" -> "f3020b2be3d4e6b6",
    "lbfgs rosenbrock n=2 classic layout=Column evaluations=17" -> "06dadc34202d33dc",
    "lbfgs rosenbrock n=2 classic layout=Column iterations=25 m=3" -> "c497674fe683ca5f",
    "lbfgs rosenbrock n=2 classic layout=Column iterations=3 m=1" -> "c8b0c173247c250e",
    "lbfgs rosenbrock n=2 classic layout=Column iterations=40 m=10 untolerant" -> "da8b6fcbad67075c",
    "lbfgs rosenbrock n=2 classic layout=Column maxstep=0.5" -> "22bd4636d0198388",
    "lbfgs rosenbrock n=2 classic layout=Transposed c2=0.1 linesearch=3" -> "b5fffc88e8c94440",
    "lbfgs rosenbrock n=2 classic layout=Transposed default" -> "f3020b2be3d4e6b6",
    "lbfgs rosenbrock n=2 classic layout=Transposed evaluations=17" -> "06dadc34202d33dc",
    "lbfgs rosenbrock n=2 classic layout=Transposed iterations=25 m=3" -> "c497674fe683ca5f",
    "lbfgs rosenbrock n=2 classic layout=Transposed iterations=3 m=1" -> "c8b0c173247c250e",
    "lbfgs rosenbrock n=2 classic layout=Transposed iterations=40 m=10 untolerant" -> "da8b6fcbad67075c",
    "lbfgs rosenbrock n=2 classic layout=Transposed maxstep=0.5" -> "22bd4636d0198388",
    "lbfgs rosenbrock n=2 lcg layout=Column c2=0.1 linesearch=3" -> "002f7cf15e9d7a52",
    "lbfgs rosenbrock n=2 lcg layout=Column default" -> "8da138d89e1928e2",
    "lbfgs rosenbrock n=2 lcg layout=Column evaluations=17" -> "5262e02feb1aeb4c",
    "lbfgs rosenbrock n=2 lcg layout=Column iterations=25 m=3" -> "3e0a613e883910f8",
    "lbfgs rosenbrock n=2 lcg layout=Column iterations=3 m=1" -> "8c1dda4b304abf45",
    "lbfgs rosenbrock n=2 lcg layout=Column iterations=40 m=10 untolerant" -> "d1ceca5c44344542",
    "lbfgs rosenbrock n=2 lcg layout=Column maxstep=0.5" -> "087902383d44fc8c",
    "lbfgs rosenbrock n=2 lcg layout=Transposed c2=0.1 linesearch=3" -> "002f7cf15e9d7a52",
    "lbfgs rosenbrock n=2 lcg layout=Transposed default" -> "8da138d89e1928e2",
    "lbfgs rosenbrock n=2 lcg layout=Transposed evaluations=17" -> "5262e02feb1aeb4c",
    "lbfgs rosenbrock n=2 lcg layout=Transposed iterations=25 m=3" -> "3e0a613e883910f8",
    "lbfgs rosenbrock n=2 lcg layout=Transposed iterations=3 m=1" -> "8c1dda4b304abf45",
    "lbfgs rosenbrock n=2 lcg layout=Transposed iterations=40 m=10 untolerant" -> "d1ceca5c44344542",
    "lbfgs rosenbrock n=2 lcg layout=Transposed maxstep=0.5" -> "087902383d44fc8c",
    "lbfgsb quadratic 6x4 box=[-0.5,0.75] layout=Column default" -> "82ba3f735856e086",
    "lbfgsb quadratic 6x4 box=[-0.5,0.75] layout=Column iterations=15 m=2" -> "3516e880258359f4",
    "lbfgsb quadratic 6x4 box=[-0.5,0.75] layout=Transposed default" -> "82ba3f735856e086",
    "lbfgsb quadratic 6x4 box=[-0.5,0.75] layout=Transposed iterations=15 m=2" -> "3516e880258359f4",
    "lbfgsb quadratic n=7 box=[-0.5,0.75] layout=Column default" -> "b63014b2eac2bdab",
    "lbfgsb quadratic n=7 box=[-0.5,0.75] layout=Column iterations=15 m=2" -> "dfd9538f73e27c75",
    "lbfgsb quadratic n=7 box=[-0.5,0.75] layout=Transposed default" -> "b63014b2eac2bdab",
    "lbfgsb quadratic n=7 box=[-0.5,0.75] layout=Transposed iterations=15 m=2" -> "dfd9538f73e27c75",
    "lbfgsb rosenbrock 4x3 box=[-0.5,0.75] layout=Column default" -> "3fcfb66c808159d9",
    "lbfgsb rosenbrock 4x3 box=[-0.5,0.75] layout=Column iterations=15 m=2" -> "f49e1dde6bc07c5d",
    "lbfgsb rosenbrock 4x3 box=[-0.5,0.75] layout=Transposed default" -> "3fcfb66c808159d9",
    "lbfgsb rosenbrock 4x3 box=[-0.5,0.75] layout=Transposed iterations=15 m=2" -> "f49e1dde6bc07c5d",
    "lbfgsb rosenbrock n=10 classic box=[-0.5,0.75] layout=Column default" -> "711acef5851c104b",
    "lbfgsb rosenbrock n=10 classic box=[-0.5,0.75] layout=Column iterations=15 m=2" -> "f8213fb077b3e487",
    "lbfgsb rosenbrock n=10 classic box=[-0.5,0.75] layout=Transposed default" -> "711acef5851c104b",
    "lbfgsb rosenbrock n=10 classic box=[-0.5,0.75] layout=Transposed iterations=15 m=2" -> "f8213fb077b3e487",
    "lbfgsb rosenbrock n=10 lcg box=[-0.5,0.75] layout=Column default" -> "a1039a1422dd122c",
    "lbfgsb rosenbrock n=10 lcg box=[-0.5,0.75] layout=Column iterations=15 m=2" -> "82d8b8603089d24d",
    "lbfgsb rosenbrock n=10 lcg box=[-0.5,0.75] layout=Transposed default" -> "a1039a1422dd122c",
    "lbfgsb rosenbrock n=10 lcg box=[-0.5,0.75] layout=Transposed iterations=15 m=2" -> "82d8b8603089d24d",
    "lbfgsb rosenbrock n=2 classic box=[-0.5,0.75] layout=Column default" -> "1ae3a41cc96d22fa",
    "lbfgsb rosenbrock n=2 classic box=[-0.5,0.75] layout=Column iterations=15 m=2" -> "2f3cd6e3d456b092",
    "lbfgsb rosenbrock n=2 classic box=[-0.5,0.75] layout=Transposed default" -> "1ae3a41cc96d22fa",
    "lbfgsb rosenbrock n=2 classic box=[-0.5,0.75] layout=Transposed iterations=15 m=2" -> "2f3cd6e3d456b092",
    "lbfgsb rosenbrock n=2 lcg box=[-0.5,0.75] layout=Column default" -> "2da44b6ab57fdf1d",
    "lbfgsb rosenbrock n=2 lcg box=[-0.5,0.75] layout=Column iterations=15 m=2" -> "2cb26d32744251e0",
    "lbfgsb rosenbrock n=2 lcg box=[-0.5,0.75] layout=Transposed default" -> "2da44b6ab57fdf1d",
    "lbfgsb rosenbrock n=2 lcg box=[-0.5,0.75] layout=Transposed iterations=15 m=2" -> "2cb26d32744251e0"
  )

  val Expected: Map[String, String] = ExpectedEntries.toMap
