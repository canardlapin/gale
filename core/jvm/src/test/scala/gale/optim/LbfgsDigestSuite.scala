package gale.optim

import gale.linalg.*

/** Bit-exact regression guard for the smooth quasi-Newton solvers that share the strong-Wolfe search and
  * [[OptimNumerics]] on the JVM: [[LBFGS]], [[LBFGSB]] and the [[AugmentedLagrangian]] inner solves (which add
  * interpolating zoom and the approximate-Wolfe test).
  *
  * Each case of a fixed corpus is reduced to one 64-bit digest of the raw bits of every point handed to the objective
  * (so every line-search trial, in order), the returned point, objective, residual, step, status, iteration and
  * evaluation counts, and the full progress trace. The stored digests in `lbfgs-digests.txt` were generated from the
  * solvers as of breeze/integration 5d3950e, so any change to their arithmetic, summation order or line-search
  * decisions fails here; such a change must regenerate the file deliberately (`LbfgsDigestSuite.main`) and justify the
  * new bits. The corpus includes restart cases that drive the two-loop direction non-finite, so the history is cleared
  * (after its ring has wrapped) and refilled.
  *
  * JVM-only because several test objectives call `exp`, `log1p`, `log` and `pow` and the starting points come from
  * `Random.nextGaussian`, whose last bits are not specified identically on Scala.js; the solvers themselves use no
  * transcendental functions. `LbfgsPortableDigestSuite` pins a transcendental-free subset on both platforms.
  */
class LbfgsDigestSuite extends munit.FunSuite:
  import LbfgsDigestSuite.*

  private lazy val expectedLines: Seq[(String, String)] =
    val stream = getClass.getResourceAsStream(DigestResource)
    assert(stream != null, s"missing resource $DigestResource")
    val source = scala.io.Source.fromInputStream(stream, "UTF-8")
    try
      source.getLines
        .filter(line => line.nonEmpty && !line.startsWith("#"))
        .map: line =>
          val tab = line.lastIndexOf('\t')
          line.substring(0, tab) -> line.substring(tab + 1)
        .toSeq
    finally source.close()

  private lazy val expected: Map[String, String] = expectedLines.toMap

  test("the digest resource has one line per case") {
    assertEquals(expected.size, expectedLines.size, "duplicate case keys in the digest resource")
  }

  test("LBFGS, LBFGSB and augmented-Lagrangian solves match the stored per-case bit digests") {
    val actual = digests
    assert(actual.size > 400, s"corpus shrank to ${actual.size} cases")
    assertEquals(actual.keySet, expected.keySet)
    val mismatches = actual.filter((key, digest) => !expected.get(key).contains(digest)).keys.toSeq.sorted
    assert(mismatches.isEmpty, s"${mismatches.size} cases changed bits, e.g. ${mismatches.take(5).mkString("; ")}")
  }

object LbfgsDigestSuite:
  val DigestResource = "/gale/optim/lbfgs-digests.txt"

  /** Regenerate the digest file: `coreJVM/Test/runMain gale.optim.LbfgsDigestSuite <path>`. */
  def main(args: Array[String]): Unit =
    val out = new java.io.PrintWriter(args(0), "UTF-8")
    try
      out.println("# LbfgsDigestSuite: case<TAB>FNV-1a digest of trial points, result fields and trace")
      digests.toSeq.sortBy(_._1).foreach((key, digest) => out.println(s"$key\t$digest"))
    finally out.close()

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
    def vec(v: DVec): Digest =
      long(v.length.toLong)
      var i = 0
      while i < v.length do
        double(v(i))
        i += 1
      this
    def counts(c: EvaluationCounts): Digest =
      c.productIterator.foreach(value => long(value.asInstanceOf[Int].toLong))
      this
    def trace(events: Vector[OptimizationProgress]): Digest =
      long(events.size.toLong)
      events.foreach(e => long(e.iteration.toLong).double(e.objective).double(e.residual).counts(e.evaluations))
      this
    def hex: String = f"$state%016x"

  /** Variable layouts: plain owned columns, multi-column matrices, and strided transposed views. */
  enum Layout:
    case Column, Wide, Transposed

  /** A smooth test problem on `rows × cols` variables; `value(x, g)` fills `g` (row-major) and returns `f`. */
  final case class Problem(
      name: String,
      rows: Int,
      cols: Int,
      start: Array[Double],
      value: (Array[Double], Array[Double]) => Double
  )

  private def rosenbrock(n: Int): (Array[Double], Array[Double]) => Double = (x, g) =>
    var f = 0.0
    var i = 0
    while i < n - 1 do
      val t = x(i + 1) - x(i) * x(i)
      val u = 1.0 - x(i)
      f += 100.0 * t * t + u * u
      g(i) += -400.0 * t * x(i) - 2.0 * u
      g(i + 1) += 200.0 * t
      i += 1
    f

  private def illConditioned(n: Int, decades: Double): (Array[Double], Array[Double]) => Double = (x, g) =>
    var f = 0.0
    var i = 0
    while i < n do
      val d = Math.pow(10.0, decades * i / Math.max(1, n - 1))
      val r = x(i) - (i % 3 - 1).toDouble
      f += 0.5 * d * r * r
      g(i) = d * r
      i += 1
    f

  private def logistic(rows: Int, cols: Int, seed: Long): (Array[Double], Array[Double]) => Double =
    val random = new scala.util.Random(seed)
    val design = Array.fill(rows, cols)(random.nextGaussian())
    val labels = Array.tabulate(rows)(i => if design(i).sum + 0.3 * random.nextGaussian() >= 0.0 then 1.0 else -1.0)
    (w, g) =>
      var f = 0.0
      var j = 0
      while j < cols do
        f += 0.005 * w(j) * w(j)
        g(j) = 0.01 * w(j)
        j += 1
      var i = 0
      while i < rows do
        var z = 0.0
        j = 0
        while j < cols do
          z += design(i)(j) * w(j)
          j += 1
        val m = labels(i) * z
        f += (if m < 0.0 then -m + Math.log1p(Math.exp(m)) else Math.log1p(Math.exp(-m)))
        val s = -labels(i) / (1.0 + Math.exp(m))
        j = 0
        while j < cols do
          g(j) += s * design(i)(j)
          j += 1
        i += 1
      f

  /** Finite only inside the open unit box, so trials beyond it are rejected as non-finite. */
  private def barrier(n: Int): (Array[Double], Array[Double]) => Double = (x, g) =>
    var f = 0.0
    var i = 0
    while i < n do
      val q = 1.0 - x(i) * x(i)
      f += -Math.log(q) + (i + 1) * 0.1 * x(i)
      g(i) = 2.0 * x(i) / q + (i + 1) * 0.1
      i += 1
    if f.isNaN then Double.PositiveInfinity else f

  private def concave: (Array[Double], Array[Double]) => Double = (x, g) =>
    g(0) = -2.0 * x(0)
    -x(0) * x(0)

  private def gaussian(n: Int, seed: Long, scale: Double): Array[Double] =
    val random = new scala.util.Random(seed)
    Array.fill(n)(scale * random.nextGaussian())

  def problems: Seq[Problem] =
    Seq(2, 10, 100).flatMap(n =>
      Seq(
        Problem(
          s"rosenbrock n=$n classic",
          n,
          1,
          Array.tabulate(n)(i => if i % 2 == 0 then -1.2 else 1.0),
          rosenbrock(n)
        ),
        Problem(s"rosenbrock n=$n gaussian", n, 1, gaussian(n, 11L + n, 1.0), rosenbrock(n))
      )
    ) ++ Seq(
      Problem("rosenbrock 4x3", 4, 3, gaussian(12, 5L, 0.7), rosenbrock(12)),
      Problem("quadratic n=5 1e2", 5, 1, gaussian(5, 21L, 3.0), illConditioned(5, 2.0)),
      Problem("quadratic n=50 1e6", 50, 1, gaussian(50, 22L, 3.0), illConditioned(50, 6.0)),
      Problem("quadratic 6x4 1e3", 6, 4, gaussian(24, 23L, 2.0), illConditioned(24, 3.0)),
      Problem("logistic 200x10", 10, 1, new Array[Double](10), logistic(200, 10, 31L)),
      Problem("barrier n=8", 8, 1, gaussian(8, 41L, 0.3), barrier(8)),
      Problem("concave", 1, 1, Array(1.0), concave)
    )

  private def matrixOf(values: Array[Double], rows: Int, cols: Int, layout: Layout): DMat =
    layout match
      case Layout.Transposed => Matrix.tabulate(cols, rows)((c, r) => values(r * cols + c)).t
      case _                 => Matrix.tabulate(rows, cols)((r, c) => values(r * cols + c))

  /** The objective hashes every point it is handed, in call order, into `trials`. */
  private def objective(problem: Problem, layout: Layout, trials: Digest): DifferentiableObjective =
    new DifferentiableObjective:
      def variableRows: Int = problem.rows
      def evaluate(at: DMat): Either[FirstOrderError, ObjectiveEvaluation] =
        trials.mat(at)
        val x = new Array[Double](at.rows * at.cols)
        at.copyRowMajorTo(x)
        val g = new Array[Double](x.length)
        val f = problem.value(x, g)
        Right(ObjectiveEvaluation(f, matrixOf(g, at.rows, at.cols, layout)))

  private def solutionDigest(trials: Digest, result: Either[FirstOrderError, FirstOrderSolution]): String =
    result match
      case Left(error) => trials.long(2L).string(error.toString).hex
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

  private def lbfgsConfigs: Seq[(String, LBFGSConfig)] =
    val loose = FirstOrderTolerance.from(1e-10, 1e-10).toOption.get
    Seq(
      "default" -> LBFGSConfig(control = traced),
      "iterations=0" -> LBFGSConfig(maxIterations = 0, control = traced),
      "iterations=3 m=1" -> LBFGSConfig(maxIterations = 3, historySize = 1, control = traced),
      "iterations=25 m=3" -> LBFGSConfig(maxIterations = 25, historySize = 3, control = traced),
      "iterations=40 m=10 untolerant" -> LBFGSConfig(
        maxIterations = 40,
        tolerance = FirstOrderTolerance.from(Double.MinPositiveValue, 0.0).toOption.get,
        control = traced
      ),
      "loose c2=0.1 linesearch=3" -> LBFGSConfig(tolerance = loose, c2 = 0.1, maxLineSearch = 3, control = traced),
      "maxstep=0.5" -> LBFGSConfig(maxIterations = 200, maximumStep = 0.5, control = traced),
      "evaluations=17" -> LBFGSConfig(control = SolverControl(maxEvaluations = 17, traceCapacity = 1000)),
      "untraced iterations=60" -> LBFGSConfig(maxIterations = 60)
    )

  private def lbfgsDigests: Seq[(String, String)] =
    for
      problem <- problems
      layout <- Layout.values.toSeq if layout != Layout.Wide || problem.cols > 1
      (configName, config) <- lbfgsConfigs
    yield
      val trials = new Digest
      val result = LBFGS.minimize(
        objective(problem, layout, trials),
        matrixOf(problem.start, problem.rows, problem.cols, layout),
        config
      )
      s"lbfgs ${problem.name} layout=$layout $configName" -> solutionDigest(trials, result)

  private def lbfgsbDigests: Seq[(String, String)] =
    val configs = Seq(
      "default" -> LBFGSBConfig(control = traced),
      "iterations=15 m=2" -> LBFGSBConfig(maxIterations = 15, historySize = 2, control = traced),
      "linesearch=4" -> LBFGSBConfig(maxIterations = 100, maxLineSearch = 4, control = traced)
    )
    for
      problem <- problems if problem.name != "barrier n=8" && problem.name != "concave"
      (boxName, low, high) <- Seq(("box=[-0.5,0.5]", -0.5, 0.5), ("box=[-2,0.8]", -2.0, 0.8))
      layout <- Seq(Layout.Column, Layout.Transposed)
      (configName, config) <- configs
    yield
      val trials = new Digest
      val bounds = MatrixBoxBounds.uniform(problem.rows, low, high, problem.cols).toOption.get
      val start = problem.start.map(v => Math.max(low, Math.min(high, v)))
      val result = LBFGSB.minimize(
        objective(problem, layout, trials),
        bounds,
        matrixOf(start, problem.rows, problem.cols, layout),
        config
      )
      s"lbfgsb ${problem.name} $boxName layout=$layout $configName" -> solutionDigest(trials, result)

  private def almDigests: Seq[(String, String)] =
    def sumConstraint(n: Int): NonlinearConstraints =
      NonlinearConstraints(n, 1, 0)(x =>
        var s = 0.0
        var i = 0
        while i < n do
          s += x(i, 0)
          i += 1
        Right(ConstraintEvaluation(DVec.fromSeq(Seq(s - 1.0)), Matrix.tabulate(1, n)((_, _) => 1.0)))
      ).toOption.get
    def ballConstraint(n: Int): NonlinearConstraints =
      NonlinearConstraints(n, 0, 1)(x =>
        var s = 0.0
        var i = 0
        while i < n do
          s += x(i, 0) * x(i, 0)
          i += 1
        Right(ConstraintEvaluation(DVec.fromSeq(Seq(s - 0.5)), Matrix.tabulate(1, n)((_, j) => 2.0 * x(j, 0))))
      ).toOption.get
    val columnProblems = problems.filter(p => p.cols == 1 && p.name != "concave" && p.name != "barrier n=8")
    for
      problem <- columnProblems
      (constraintName, constraints) <- Seq(
        "sum=1" -> sumConstraint(problem.rows),
        "ball<=0.5" -> ballConstraint(problem.rows)
      )
      bounded <- Seq(false, true)
      roundoff <- Seq(true, false)
    yield
      val trials = new Digest
      val bounds = if bounded then Some(MatrixBoxBounds.uniform(problem.rows, -0.9, 0.9).toOption.get) else None
      val start = problem.start.map(v => Math.max(-0.9, Math.min(0.9, v)))
      val config = AugmentedLagrangianConfig(maxIterations = 12, roundoffAwareLineSearch = roundoff, control = traced)
      val result = AugmentedLagrangian.minimize(
        objective(problem, Layout.Column, trials),
        constraints,
        matrixOf(start, problem.rows, 1, Layout.Column),
        bounds,
        config
      )
      val digest = result match
        case Left(error) => trials.long(2L).string(error.toString).hex
        case Right(r)    =>
          r.innerTrace
            .foldLeft(
              trials
                .long(1L)
                .mat(r.primal)
                .double(r.objective)
                .vec(r.constraints)
                .vec(r.multipliers)
                .string(r.diagnostics.toString)
                .string(r.status.toString)
                .long(r.iterations.toLong)
                .long(r.innerIterations)
                .double(r.penalty)
                .counts(r.evaluations)
                .trace(r.trace)
                .string(r.lastInnerStatus.toString)
            )((d, inner) => d.string(inner.toString))
            .hex
      s"alm ${problem.name} $constraintName bounded=$bounded roundoff=$roundoff" -> digest

  /** Two-loop fallback: see [[LbfgsPortableDigests.gatedRestart]]. */
  private def restartDigests: Seq[(String, String)] =
    for
      (name, gate, slope, history) <- LbfgsPortableDigests.restartCases
      layout <- Seq(Layout.Column, Layout.Transposed)
    yield
      val trials = new Digest
      val problem =
        Problem(s"restart $name", 3, 1, Array(-1.2, 1.0, 0.0), LbfgsPortableDigests.gatedRestart(gate, slope))
      val config = LBFGSConfig(maxIterations = 60, historySize = history, control = traced)
      val result = LBFGS.minimize(objective(problem, layout, trials), matrixOf(problem.start, 3, 1, layout), config)
      s"lbfgs ${problem.name} layout=$layout" -> solutionDigest(trials, result)

  def digests: Map[String, String] =
    val all = lbfgsDigests ++ lbfgsbDigests ++ almDigests ++ restartDigests
    assert(all.map(_._1).distinct.size == all.size, "duplicate case keys in the corpus")
    all.toMap
