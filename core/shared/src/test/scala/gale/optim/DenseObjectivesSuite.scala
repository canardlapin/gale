package gale.optim

import gale.TestAccess
import gale.linalg.{DMat, DVec}

final class DenseObjectivesSuite extends munit.FunSuite:
  private def matrix(rows: Vector[Vector[Double]]): DMat =
    DMat.tabulate(rows.size, rows.head.size)((row, column) => rows(row)(column))
  private def column(values: Double*): DMat = DMat.fromArrayRowMajor(values.size, 1, values.toArray)

  test("least squares agrees with direct objective and gradient arithmetic"):
    val objective = DenseObjectives
      .leastSquares(matrix(Vector(Vector(1.0, 2.0), Vector(3.0, 4.0))), DVec.fromSeq(Seq(1.0, -1.0)))
      .toOption
      .get
    val evaluated = objective.evaluate(column(0.5, -0.25)).toOption.get
    assertEqualsDouble(evaluated.value, 0.8125, 1e-15)
    assertEqualsDouble(evaluated.gradient(0, 0), 1.75, 1e-15)
    assertEqualsDouble(evaluated.gradient(1, 0), 2.0, 1e-15)
    assert(GradientCheck.directional(objective, column(0.5, -0.25), column(0.3, -0.8)).toOption.get.passed)

  test("logistic is stable at saturated margins and matches independent derivative"):
    val design = matrix(Vector(Vector(1.0), Vector(-1.0)))
    val objective = DenseObjectives.logistic(design, DVec.fromSeq(Seq(1.0, -1.0)), l2 = 0.2).toOption.get
    val evaluated = objective.evaluate(column(0.7)).toOption.get
    val tail = math.exp(-0.7)
    val expectedLoss = math.log1p(tail) + 0.1 * 0.7 * 0.7
    val expectedGradient = -tail / (1.0 + tail) + 0.14
    assertEqualsDouble(evaluated.value, expectedLoss, 1e-14)
    assertEqualsDouble(evaluated.gradient(0, 0), expectedGradient, 1e-14)
    assert(GradientCheck.directional(objective, column(0.7), column(-0.4)).toOption.get.passed)
    val saturated = DenseObjectives
      .logistic(matrix(Vector(Vector(1000.0), Vector(-1000.0))), DVec.fromSeq(Seq(1.0, -1.0)))
      .toOption
      .get
    for weight <- Vector(column(1.0), column(-1.0)) do
      val result = saturated.evaluate(weight).toOption.get
      assert(result.value.isFinite)
      assert(result.gradient(0, 0).isFinite)

  test("constructors reject invalid data and evaluation rejects invalid parameter arithmetic"):
    assert(DenseObjectives.leastSquares(DMat.zeros(0, 1), DVec.zeros(0)).isLeft)
    assert(DenseObjectives.leastSquares(DMat.zeros(2, 1), DVec.zeros(1)).isLeft)
    assert(DenseObjectives.logistic(DMat.zeros(1, 1), DVec.fromSeq(Seq(0.0))).isLeft)
    assert(DenseObjectives.logistic(DMat.zeros(1, 1), DVec.fromSeq(Seq(1.0)), l2 = -1.0).isLeft)
    assert(DenseObjectives.leastSquares(matrix(Vector(Vector(Double.NaN))), DVec.fromSeq(Seq(0.0))).isLeft)
    val overflow =
      DenseObjectives.leastSquares(matrix(Vector(Vector(Double.MaxValue))), DVec.fromSeq(Seq(0.0))).toOption.get
    assert(overflow.evaluate(column(Double.MaxValue)).isLeft)
    val valid = DenseObjectives.leastSquares(matrix(Vector(Vector(1.0))), DVec.fromSeq(Seq(0.0))).toOption.get
    assert(valid.evaluate(DMat.zeros(1, 2)).isLeft)

  test("objectives own design and response snapshots"):
    val design = matrix(Vector(Vector(1.0), Vector(2.0)))
    val response = DVec.fromSeq(Seq(3.0, 4.0))
    val leastSquares = DenseObjectives.leastSquares(design, response).toOption.get
    val logistic = DenseObjectives.logistic(design, DVec.fromSeq(Seq(1.0, -1.0))).toOption.get
    TestAccess.dmatStorage(design)(0) = 99.0
    TestAccess.dvecStorage(response)(0) = -99.0
    val ls = leastSquares.evaluate(column(1.0)).toOption.get
    val log = logistic.evaluate(column(1.0)).toOption.get
    assertEqualsDouble(ls.value, 2.0, 1e-15)
    assertEqualsDouble(log.value, (math.log1p(math.exp(-1.0)) + math.log1p(math.exp(2.0))) / 2.0, 1e-15)
