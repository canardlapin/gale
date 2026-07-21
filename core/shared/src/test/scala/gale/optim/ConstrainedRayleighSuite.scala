package gale.optim

import gale.linalg.{DMat, Matrix, Vec}

class ConstrainedRayleighSuite extends munit.FunSuite:
  test("nonnegative solve matches an independent two-dimensional angle oracle"):
    val numerator = matrix(Vector(Vector(3.0, -2.0), Vector(-2.0, 2.0)))
    val denominator = matrix(Vector(Vector(2.0, 0.4), Vector(0.4, 1.0)))
    val result = ProjectedRayleigh.solveNonnegative(numerator, denominator).toOption.get
    val expected = angleOracle(numerator, denominator, 200000)

    assert(result.converged)
    assertEqualsDouble(result.root, expected, 2e-7)
    assert(result.direction(0) >= 0.0)
    assert(result.direction(1) >= 0.0)
    assertEqualsDouble(quadratic(result.direction, denominator), 1.0, 1e-9)
    assert(result.certificate.stationarityResidual <= 1e-7)
    assert(result.certificate.constraintViolation <= 1e-12)

  test("inactive nonnegative constraint recovers the ordinary positive eigenvector"):
    val numerator = matrix(Vector(Vector(4.0, 1.0), Vector(1.0, 2.0)))
    val denominator = DMat.eye(2)
    val result = ProjectedRayleigh.solveNonnegative(numerator, denominator).toOption.get
    val expected = 3.0 + Math.sqrt(2.0)

    assert(result.converged)
    assertEqualsDouble(result.root, expected, 1e-8)
    assert(result.direction(0) > 0.0)
    assert(result.direction(1) > 0.0)

  test("configuration and denominator failures are typed"):
    val identity = DMat.eye(2)
    val invalidConfig = ProjectedRayleigh.solveNonnegative(
      identity,
      identity,
      ProjectedRayleighConfig(tolerance = -1.0)
    )
    val indefinite = matrix(Vector(Vector(1.0, 0.0), Vector(0.0, -1.0)))
    val denominatorFailure = ProjectedRayleigh.solve(
      identity,
      indefinite,
      RayleighCone.NonnegativeOrthant,
      Vec(0.0, 1.0)
    )

    assert(invalidConfig.left.toOption.exists(_.isInstanceOf[ConstrainedRayleighError.InvalidConfiguration]))
    assertEquals(
      denominatorFailure.left.toOption,
      Some(ConstrainedRayleighError.NonPositiveDenominator(-1.0))
    )

  private def angleOracle(numerator: DMat, denominator: DMat, intervals: Int): Double =
    var best = Double.NegativeInfinity
    var index = 0
    while index <= intervals do
      val angle = 0.5 * Math.PI * index.toDouble / intervals.toDouble
      val direction = Vec(Math.cos(angle), Math.sin(angle))
      best = Math.max(best, quadratic(direction, numerator) / quadratic(direction, denominator))
      index += 1
    best

  private def quadratic(value: gale.linalg.DVec, matrix: DMat): Double =
    var result = 0.0
    var row = 0
    while row < matrix.rows do
      var projected = 0.0
      var column = 0
      while column < matrix.cols do
        projected += matrix(row, column) * value(column)
        column += 1
      result += value(row) * projected
      row += 1
    result

  private def matrix(rows: Vector[Vector[Double]]): DMat =
    val out = Matrix.newBuilder(rows.length, rows.head.length)
    var row = 0
    while row < rows.length do
      var column = 0
      while column < rows.head.length do
        out(row, column) = rows(row)(column)
        column += 1
      row += 1
    out.result()
