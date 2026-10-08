package gale.optim

import gale.linalg.DVec

/** cmaes 0.12.0 with positive weights and canonical initial generation index; see tools/optim-cmaes/compare.py. */
final class CMAESReferenceSuite extends munit.FunSuite:
  test("independent first-generation covariance and CSA reference, seed 7"):
    val session = CMAES
      .start(
        DVec.fromSeq(Seq(1.0, -2.0, 0.5)),
        CMAESConfig(seed = 7L, initialStepSize = 0.3, populationSize = 6, eigenUpdatePeriod = 1)
      )
      .toOption
      .get
    val batch = session.ask().toOption.get.get
    val fitness = batch.points.map(x => x(0) * x(0) + 2 * x(1) * x(1) + 4 * x(2) * x(2))
    session.tell(batch, fitness).toOption.get
    val d = session.snapshot.distribution
    val expectedMean = Vector(1.1813249865317117, -1.509671688252211, 0.24291686017388808)
    val expectedC = Vector(0.9581169244325932, 0.20390053144268167, -0.08687705899350641, 0.20390053144268167,
      1.407427364555064, -0.27085547919618363, -0.08687705899350641, -0.27085547919618363, 1.020647809945344)
    assertEqualsDouble(d.stepSize, 0.3352696019824501, 2e-12)
    for i <- 0 until 3 do
      assertEqualsDouble(d.mean(i), expectedMean(i), 2e-12)
      for j <- 0 until 3 do assertEqualsDouble(d.covariance(i, j), expectedC(i * 3 + j), 2e-12)

  test("independent first-generation covariance and CSA reference, seed 42"):
    val session = CMAES
      .start(
        DVec.fromSeq(Seq(1.0, -2.0, 0.5)),
        CMAESConfig(seed = 42L, initialStepSize = 0.3, populationSize = 6, eigenUpdatePeriod = 1)
      )
      .toOption
      .get
    val batch = session.ask().toOption.get.get
    val fitness = batch.points.map(x => x(0) * x(0) + 2 * x(1) * x(1) + 4 * x(2) * x(2))
    session.tell(batch, fitness).toOption.get
    val d = session.snapshot.distribution
    val expectedMean = Vector(1.2620393747466938, -1.6171382410246629, 0.07003850530724276)
    val expectedC = Vector(1.0347106807277517, 0.23761385330247786, -0.26598490045622064, 0.23761385330247786,
      1.2018852528863881, -0.3782008311171554, -0.26598490045622064, -0.3782008311171554, 1.2875654662836977)
    assertEqualsDouble(d.stepSize, 0.3469055038737974, 2e-12)
    for i <- 0 until 3 do
      assertEqualsDouble(d.mean(i), expectedMean(i), 2e-12)
      for j <- 0 until 3 do assertEqualsDouble(d.covariance(i, j), expectedC(i * 3 + j), 2e-12)
  test("canonical first-generation h-sigma gating suppresses the rank-one path"):
    // In 1D with mu=1, cmu=0. For 1.95 < y < 2.1 the canonical h-sigma gate
    // is false but an off-by-one generation denominator would make it true.
    // Hence pc=0 and C=1-c1+c1*cc*(2-cc), independent of the particular y.
    val session =
      CMAES.start(DVec.zeros(1), CMAESConfig(seed = 100L, initialStepSize = 1.0, populationSize = 2)).toOption.get
    val batch = session.ask().toOption.get.get
    val selected = batch.points.map(_(0)).max
    assert(selected > 1.95 && selected < 2.1)
    session.tell(batch, batch.points.map(x => -x(0))).toOption.get
    assertEqualsDouble(session.snapshot.distribution.covariance(0, 0), 0.9740436715226631, 2e-15)
