package gale.optim

import gale.linalg.{DMat, DVec}
import gale.spectral.{Eigen, EigenSelection}

final class CMAESSuite extends munit.FunSuite:
  private def vec(xs: Double*): DVec = DVec.fromSeq(xs)
  private def sphere(x: DVec): Double = x.toSeq.map(v => v * v).sum
  private def start(config: CMAESConfig = CMAESConfig(seed = 42)): CMAESSession =
    CMAES.start(vec(2, -3), config).toOption.get
  private def batch(s: CMAESSession): CMAESBatch = s.ask().toOption.get.get

  test("seeded sampling and results repeat within a runtime; distinct seeds differ"):
    val a = start()
    val b = start()
    val first = batch(a).points.map(_.toSeq)
    assertEquals(first, batch(b).points.map(_.toSeq))
    assertNotEquals(batch(start(CMAESConfig(seed = 43))).points.map(_.toSeq), first)
    val config = CMAESConfig(seed = 7, targetValue = Some(1e-12), maxEvaluations = 10000)
    val r1 = CMAES.minimize(sphere, vec(2, -3), config).toOption.get
    val r2 = CMAES.minimize(sphere, vec(2, -3), config).toOption.get
    assertEquals(r1.best.get.coordinates.toSeq, r2.best.get.coordinates.toSeq)
    assertEquals(r1.work, r2.work)
    assertEquals(r1.stopReason, CMAESStop.TargetReached)
    assert(r1.best.get.value <= 1e-12)

  test("ask/tell rejects outstanding, foreign, stale and malformed batches atomically"):
    val s = start()
    val first = batch(s)
    assert(s.ask().isLeft)
    assert(s.tell(batch(start()), Vector.fill(first.points.size)(1.0)).isLeft)
    assert(s.tell(first, Vector(1.0)).isLeft)
    assert(s.tell(first, Vector.fill(first.points.size)(Double.NaN)).isLeft)
    assertEquals(s.work.evaluations, 0)
    assertEquals(s.snapshot.distribution.mean.toSeq, Seq(2.0, -3.0))
    assert(s.tell(first, first.points.map(sphere)).isRight)
    assertEquals(s.work.evaluations, first.points.size)
    assert(s.tell(first, first.points.map(sphere)).isLeft)
    assert(s.cancel(first).isLeft)

  test("partial evaluation budget records the best without adapting"):
    val s = start(CMAESConfig(maxEvaluations = 3, populationSize = 6))
    val initial = s.snapshot.distribution
    val candidates = batch(s)
    assertEquals(candidates.points.size, 3)
    s.tell(candidates, Vector(3.0, 1.0, 2.0)).toOption.get
    assertEquals(s.work.evaluations, 3)
    assertEquals(s.work.generations, 0)
    assertEquals(s.snapshot.distribution.mean.toSeq, initial.mean.toSeq)
    assertEquals(s.bestPoint.get.coordinates.toSeq, candidates.points(1).toSeq)
    assertEquals(s.stopReason, Some(CMAESStop.EvaluationLimit))
    assertEquals(s.ask(), Right(None))

  test("cancellation receipts preserve a finite prefix and cannot be replayed"):
    val s = start()
    val candidates = batch(s)
    assert(s.cancel(candidates, Vector(Double.PositiveInfinity)).isLeft)
    s.cancel(candidates, Vector(4.0, 2.0)).toOption.get
    assertEquals(s.work.evaluations, 2)
    assertEquals(s.work.generations, 0)
    assertEquals(s.bestPoint.get.value, 2.0)
    assertEquals(s.snapshot.distribution.mean.toSeq, Seq(2.0, -3.0))
    assert(s.cancel(candidates).isLeft)

  test("minimize charges failed calls and preserves preceding best; cancellation is per evaluation"):
    var calls = 0
    val failed = CMAES
      .minimize(
        x => { calls += 1; if calls == 3 then throw new IllegalArgumentException("oracle") else sphere(x) },
        vec(2, -3)
      )
      .toOption
      .get
    assertEquals(calls, 3)
    assertEquals(failed.work.evaluations, 3)
    assertEquals(failed.work.generations, 0)
    assert(failed.best.nonEmpty)
    assertEquals(failed.stopReason, CMAESStop.ObjectiveFailure("oracle"))
    calls = 0
    val cancelled = CMAES
      .minimize(x => { calls += 1; sphere(x) }, vec(2, -3), control = CMAESControl(cancelled = () => calls >= 2))
      .toOption
      .get
    assertEquals(cancelled.stopReason, CMAESStop.Cancelled)
    assertEquals(cancelled.work.evaluations, 2)
    assertEquals(cancelled.work.generations, 0)
    val invalid = CMAES.minimize(_ => Double.PositiveInfinity, vec(1, 2)).toOption.get
    assertEquals(invalid.work.evaluations, 1)
    assertEquals(invalid.best, None)

  test("bounded candidates are feasible; fixed coordinates are excluded from covariance"):
    val bounds = MatrixBoxBounds
      .from(
        DMat.tabulate(3, 1)((i, _) => if i == 0 then 0.25 else -1.0),
        DMat.tabulate(3, 1)((i, _) => if i == 0 then 0.25 else 1.0)
      )
      .toOption
      .get
    var calls = 0
    val result = CMAES
      .minimize(
        x => {
          calls += 1
          assertEquals(x(0), 0.25)
          assert((1 until 3).forall(i => x(i) >= -1.0 && x(i) <= 1.0))
          (x(1) - 0.8) * (x(1) - 0.8) + (x(2) + 0.4) * (x(2) + 0.4)
        },
        vec(0.25, 0.0, 0.0),
        CMAESConfig(seed = 18, initialStepSize = 1.0, targetValue = Some(1e-10), maxEvaluations = 10000),
        Some(bounds)
      )
      .toOption
      .get
    assertEquals(result.stopReason, CMAESStop.TargetReached)
    assertEquals(result.work.evaluations, calls)
    assert(result.work.rejected > 0)
    assertEquals(result.distribution.freeCoordinates, Vector(1, 2))
    assertEquals(result.distribution.covariance.rows, 2)

  test("all-fixed boxes evaluate once; infeasible starts and invalid configurations fail before callbacks"):
    val bounds = MatrixBoxBounds.uniform(2, 1.0, 1.0).toOption.get
    val result = CMAES.minimize(sphere, vec(1, 1), bounds = Some(bounds)).toOption.get
    assertEquals(result.stopReason, CMAESStop.AllCoordinatesFixed)
    assertEquals(result.work.evaluations, 1)
    assertEquals(result.best.get.value, 2.0)
    assertEquals(result.distribution.covariance.rows, 0)
    assert(CMAES.start(vec(0, 1), bounds = Some(bounds)).isLeft)
    assert(CMAES.start(vec(Double.NaN)).isLeft)
    assert(CMAES.start(DVec.zeros(0)).isLeft)
    assert(CMAES.start(vec(1), CMAESConfig(initialStepSize = 0)).isLeft)
    assert(CMAES.start(vec(1), CMAESConfig(populationSize = 1)).isLeft)
    assert(CMAES.start(vec(1), CMAESConfig(maxEvaluations = 0)).isLeft)
    assert(CMAES.start(vec(1), CMAESConfig(targetValue = Some(Double.NaN))).isLeft)

  test("sampling is capped without objective calls; vanished representable steps are explicit"):
    var calls = 0
    val bounds = MatrixBoxBounds.uniform(2, -1e-100, 1e-100).toOption.get
    val result = CMAES
      .minimize(_ => { calls += 1; 0.0 }, vec(0, 0), CMAESConfig(maxSamplingAttempts = 5), Some(bounds))
      .toOption
      .get
    assertEquals(result.stopReason, CMAESStop.SamplingLimit)
    assertEquals(result.work.sampled, 5L)
    assertEquals(result.work.rejected, 5L)
    assertEquals(calls, 0)
    val stuck = CMAES.minimize(sphere, vec(1e100), CMAESConfig(initialStepSize = 1e-100)).toOption.get
    assertEquals(stuck.stopReason, CMAESStop.NoRepresentableStep)
    assertEquals(stuck.work.evaluations, 0)

  test("rank invariance, positive definite covariance, and snapshot ownership across repeated updates"):
    val a = start(CMAESConfig(seed = 91, maxGenerations = 15, stagnationGenerations = 0))
    val b = start(CMAESConfig(seed = 91, maxGenerations = 15, stagnationGenerations = 0))
    val before = a.snapshot
    var generation = 0
    while generation < 15 do
      val x = batch(a)
      val y = batch(b)
      a.tell(x, x.points.map(sphere)).toOption.get
      b.tell(y, y.points.map(v => math.log1p(sphere(v)))).toOption.get
      assertEquals(a.snapshot.distribution.mean.toSeq, b.snapshot.distribution.mean.toSeq)
      val c = a.snapshot.distribution.covariance
      val eig = Eigen.eigSymmetric(c, EigenSelection.All).toOption.get
      assert(eig.eigenvalues.toSeq.forall(_ > 0.0))
      assertEquals(c(0, 1), c(1, 0))
      generation += 1
    assertEquals(before.distribution.mean.toSeq, Seq(2.0, -3.0))
    assertEquals(before.distribution.covariance(0, 0), 1.0)
    assertEquals(a.stopReason, Some(CMAESStop.GenerationLimit))

  test("rotated ill-conditioned and nonsmooth analytic objectives reach independent targets"):
    val rotated: DVec => Double = x => {
      val a = (x(0) + x(1)) / math.sqrt(2)
      val b = (x(0) - x(1)) / math.sqrt(2)
      1e4 * a * a + b * b
    }
    for objective <- Vector(rotated, (x: DVec) => math.abs(x(0)) + 3.0 * math.abs(x(1))) do
      val result = CMAES
        .minimize(
          objective,
          vec(2, -3),
          CMAESConfig(seed = 99, targetValue = Some(1e-8), maxEvaluations = 20000, stagnationGenerations = 200)
        )
        .toOption
        .get
      assertEquals(result.stopReason, CMAESStop.TargetReached)
      assert(objective(result.best.get.coordinates) <= 1e-8)

  test("IPOP resets state with distinct seeds, grows populations, and keeps one global budget"):
    val result = CMAES
      .minimize(
        _ => 1.0,
        vec(0, 0),
        CMAESConfig(seed = 12, populationSize = 4, maxEvaluations = 29, stagnationGenerations = 1),
        restarts = CMAESRestarts(maxRestarts = 5)
      )
      .toOption
      .get
    assertEquals(result.work.evaluations, 29)
    assertEquals(result.runs.map(_.populationSize), Vector(4, 8, 16))
    assertEquals(result.runs.map(_.seed).distinct.size, result.runs.size)
    assertEquals(result.runs.map(_.work.evaluations).sum, result.work.evaluations)
    assertEquals(result.stopReason, CMAESStop.EvaluationLimit)
    assertEquals(result.best.get.value, 1.0)

  test("trace is bounded; control failures and zero generation budgets do not invent results"):
    val result = CMAES
      .minimize(sphere, vec(1, 2), CMAESConfig(maxGenerations = 8), control = CMAESControl(traceCapacity = 3))
      .toOption
      .get
    assertEquals(result.trace.size, 3)
    assertEquals(result.trace.last.evaluations, result.work.evaluations)
    val failed = CMAES
      .minimize(sphere, vec(1, 2), control = CMAESControl(cancelled = () => throw new IllegalStateException("control")))
      .toOption
      .get
    assertEquals(failed.stopReason, CMAESStop.ControlFailure("control"))
    assertEquals(failed.work.evaluations, 0)
    val zero = CMAES.minimize(sphere, vec(1, 2), CMAESConfig(maxGenerations = 0)).toOption.get
    assertEquals(zero.stopReason, CMAESStop.GenerationLimit)
    assertEquals(zero.best, None)

  test("a progress stop or failure vetoes a restart at a terminal generation"):
    for throws <- Vector(false, true) do
      var calls = 0
      val result = CMAES
        .minimize(
          _ => { calls += 1; 1.0 },
          vec(0, 0),
          CMAESConfig(populationSize = 4, stagnationGenerations = 1),
          control = CMAESControl(progress =
            Some(p =>
              if p.generation < 2 then true
              else if throws then throw new IllegalStateException("stop restart")
              else false
            )
          ),
          restarts = CMAESRestarts(maxRestarts = 3)
        )
        .toOption
        .get
      assertEquals(calls, 8)
      assertEquals(result.work.evaluations, 8)
      assertEquals(result.runs.size, 1)
      assertEquals(result.stopReason, if throws then CMAESStop.ControlFailure("stop restart") else CMAESStop.Cancelled)
