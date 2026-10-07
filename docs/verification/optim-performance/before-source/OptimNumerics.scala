package gale.optim

import gale.linalg.DMat

private[optim] object OptimNumerics:
  def validate(value: DMat, rows: Int): Either[FirstOrderError, Unit] =
    if rows <= 0 || value.cols <= 0 then
      Left(FirstOrderError.InvalidConfiguration("optimization variables must be nonempty"))
    else if value.rows != rows then Left(FirstOrderError.ShapeMismatch("objective", rows, value.rows))
    else finiteMatrix(value)

  def like(value: DMat, expected: DMat): Either[FirstOrderError, Unit] =
    if value.rows != expected.rows || value.cols != expected.cols then
      Left(
        FirstOrderError.InvalidConfiguration(
          s"expected ${expected.rows}x${expected.cols}, got ${value.rows}x${value.cols}"
        )
      )
    else finiteMatrix(value)

  def finiteMatrix(value: DMat): Either[FirstOrderError, Unit] =
    var r = 0
    while r < value.rows do
      var c = 0
      while c < value.cols do
        if !value(r, c).isFinite then
          return Left(FirstOrderError.NonFiniteValue("numerical value", r * value.cols + c, value(r, c)))
        c += 1
      r += 1
    Right(())

  def finiteScalar(value: Double, context: String): Either[FirstOrderError, Unit] =
    if value.isFinite then Right(()) else Left(FirstOrderError.NonFiniteValue(context, 0, value))

  def array(value: DMat): Array[Double] =
    val result = new Array[Double](value.rows * value.cols)
    value.copyRowMajorTo(result)
    result

  def matrix(values: Array[Double], rows: Int, cols: Int): DMat = DMat.fromArrayRowMajor(rows, cols, values)

  def affine(at: DMat, direction: DMat, step: Double): Either[FirstOrderError, DMat] =
    val result = DMat.tabulate(at.rows, at.cols)((r, c) => at(r, c) + step * direction(r, c))
    finiteMatrix(result).map(_ => result)

  def difference(left: DMat, right: DMat): Array[Double] =
    val result = array(left)
    var i = 0
    while i < result.length do
      result(i) -= right(i / left.cols, i % left.cols)
      i += 1
    result

  def dot(left: Array[Double], right: Array[Double]): Double =
    var sum = 0.0
    var i = 0
    while i < left.length do
      sum += left(i) * right(i)
      i += 1
    sum

  def dot(left: DMat, right: DMat): Double =
    var sum = 0.0
    var r = 0
    while r < left.rows do
      var c = 0
      while c < left.cols do
        sum += left(r, c) * right(r, c)
        c += 1
      r += 1
    sum

  def normInf(value: DMat): Double =
    var result = 0.0
    var r = 0
    while r < value.rows do
      var c = 0
      while c < value.cols do
        result = Math.max(result, Math.abs(value(r, c)))
        c += 1
      r += 1
    result

  def identical(left: DMat, right: DMat): Boolean =
    left.rows == right.rows && left.cols == right.cols &&
      (left.eq(right) || {
        var same = true
        var r = 0
        while r < left.rows && same do
          var c = 0
          while c < left.cols && same do
            same = left(r, c) == right(r, c)
            c += 1
          r += 1
        same
      })

  def mapping(at: DMat, next: DMat, step: Double): Either[FirstOrderError, (Double, Double)] =
    var residual = 0.0
    var resolution = 0.0
    var r = 0
    while r < at.rows do
      var c = 0
      while c < at.cols do
        residual = Math.max(residual, Math.abs(at(r, c) - next(r, c)) / step)
        resolution = Math.max(resolution, 2.0 * (Math.ulp(Math.max(Math.abs(at(r, c)), Math.abs(next(r, c)))) / step))
        c += 1
      r += 1
    if !residual.isFinite || !resolution.isFinite then
      Left(FirstOrderError.NumericalFailure("gradient mapping overflow"))
    else Right((residual, resolution))

  def result(
      point: DMat,
      objective: Double,
      residual: Double,
      iterations: Int,
      status: FirstOrderStoppingStatus,
      settings: FirstOrderSettings,
      execution: OptimizationExecution,
      step: Double,
      resolution: Double = 0.0,
      change: Double = 0.0
  ): FirstOrderSolution =
    val certificate = FirstOrderCertificate(
      ValueSummary.from(point),
      None,
      objective,
      residual,
      0.0,
      change,
      iterations,
      settings,
      resolution,
      0.0,
      Some(point),
      None
    )
    FirstOrderSolution(point, None, objective, status, certificate, execution.counts, step, 0.0, execution.trace)
