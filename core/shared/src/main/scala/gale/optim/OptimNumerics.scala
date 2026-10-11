package gale.optim

import gale.linalg.{DMat, value}
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Shared numerics for the optimization solvers.
  *
  * Every elementwise reduction visits logical elements in row-major order. Matrices whose storage is contiguous
  * row-major (every owned matrix, including the `n × 1` columns the solvers pass around) take a flat-index fast path
  * that visits the same elements in the same order, so both paths produce identical bits.
  */
private[optim] object OptimNumerics:
  /** Logical row-major element `i` of `value` is `value.data(value.offset.value + i)`. */
  private inline def packed(value: DMat): Boolean = value.isContiguousRowMajor

  /** Equal shapes, both [[packed]]: flat index `i` addresses the same logical element of each. */
  private inline def packedPair(left: DMat, right: DMat): Boolean =
    left.rows == right.rows && left.cols == right.cols && packed(left) && packed(right)

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
    if packed(value) then
      val data = value.data
      val base = value.offset.value
      val size = value.rows * value.cols
      var i = 0
      while i < size do
        val x = data(base + i)
        if !x.isFinite then return Left(FirstOrderError.NonFiniteValue("numerical value", i, x))
        i += 1
      return Right(())
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

  inline def tabulate(size: Int)(inline f: Int => Double): Array[Double] =
    val out = new Array[Double](size)
    var i = 0
    while i < size do
      out(i) = f(i)
      i += 1
    out

  def allFinite(values: Array[Double]): Boolean =
    var i = 0
    while i < values.length do
      if !values(i).isFinite then return false
      i += 1
    true

  /** `at + step * direction`, failing on the first non-finite element in row-major order. */
  def affine(at: DMat, direction: DMat, step: Double): Either[FirstOrderError, DMat] =
    if packedPair(at, direction) then
      val size = at.rows * at.cols
      val out = DoubleArray.alloc(size)
      val a = at.data
      val aBase = at.offset.value
      val d = direction.data
      val dBase = direction.offset.value
      var i = 0
      while i < size do
        val x = a(aBase + i) + step * d(dBase + i)
        if !x.isFinite then return Left(FirstOrderError.NonFiniteValue("numerical value", i, x))
        out(i) = x
        i += 1
      Right(DMat.fromDoubleArrayOwned(at.rows, at.cols, out))
    else
      val result = DMat.tabulate(at.rows, at.cols)((r, c) => at(r, c) + step * direction(r, c))
      finiteMatrix(result).map(_ => result)

  def difference(left: DMat, right: DMat): Array[Double] =
    val result = array(left)
    if packedPair(left, right) then
      val data = right.data
      val base = right.offset.value
      var i = 0
      while i < result.length do
        result(i) -= data(base + i)
        i += 1
    else
      var i = 0
      while i < result.length do
        result(i) -= right(i / left.cols, i % left.cols)
        i += 1
    result

  /** Streaming scaled sum of squares, avoiding both squared-value overflow and a sqrt per element. */
  final class ScaledNorm:
    private var scale = 0.0
    private var squares = 1.0
    def add(value: Double): Unit =
      val magnitude = Math.abs(value)
      if !magnitude.isFinite then scale += magnitude
      else if magnitude > 0.0 then
        if scale < magnitude then
          val ratio = scale / magnitude
          squares = 1.0 + squares * ratio * ratio
          scale = magnitude
        else
          val ratio = magnitude / scale
          squares += ratio * ratio
    def result: Double = scale * Math.sqrt(squares)

  def norm(values: Array[Double]): Double =
    val accumulated = new ScaledNorm
    var i = 0
    while i < values.length do
      accumulated.add(values(i))
      i += 1
    accumulated.result

  def dot(left: Array[Double], right: Array[Double]): Double =
    var sum = 0.0
    var i = 0
    while i < left.length do
      sum += left(i) * right(i)
      i += 1
    sum

  def dot(left: DMat, right: DMat): Double =
    var sum = 0.0
    if packedPair(left, right) then
      val l = left.data
      val lBase = left.offset.value
      val rData = right.data
      val rBase = right.offset.value
      val size = left.rows * left.cols
      var i = 0
      while i < size do
        sum += l(lBase + i) * rData(rBase + i)
        i += 1
      return sum
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
    if packed(value) then
      val data = value.data
      val base = value.offset.value
      val size = value.rows * value.cols
      var i = 0
      while i < size do
        result = Math.max(result, Math.abs(data(base + i)))
        i += 1
      return result
    var r = 0
    while r < value.rows do
      var c = 0
      while c < value.cols do
        result = Math.max(result, Math.abs(value(r, c)))
        c += 1
      r += 1
    result

  def identical(left: DMat, right: DMat): Boolean =
    if left.rows != right.rows || left.cols != right.cols then false
    else if left.eq(right) then true
    else if packedPair(left, right) then
      val l = left.data
      val lBase = left.offset.value
      val rData = right.data
      val rBase = right.offset.value
      val size = left.rows * left.cols
      var i = 0
      while i < size && l(lBase + i) == rData(rBase + i) do i += 1
      i == size
    else
      var same = true
      var r = 0
      while r < left.rows && same do
        var c = 0
        while c < left.cols && same do
          same = left(r, c) == right(r, c)
          c += 1
        r += 1
      same

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
