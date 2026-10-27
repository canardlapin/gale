package gale.bench

import gale.linalg.*
import gale.syntax.all.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import scala.compiletime.uninitialized

/** Compare identical subtract/divide arithmetic. Use -prof gc for allocation; no performance conclusion follows merely
  * from compiling this benchmark.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 3, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(1)
@State(Scope.Thread)
class MarginJmh:
  @Param(Array("256"))
  var n: Int = 0

  @Param(Array("rows", "cols"))
  var axisName: String = "cols"

  @Param(Array("rowMajor", "transposed"))
  var layout: String = "rowMajor"

  private var input: DMat = uninitialized
  private var center: Margin[?] = uninitialized
  private var scale: Margin[?] = uninitialized
  private var centerRows: Margin.ByRow = uninitialized
  private var scaleRows: Margin.ByRow = uninitialized
  private var centerCols: Margin.ByCol = uninitialized
  private var scaleCols: Margin.ByCol = uninitialized
  private var destination: DMatBuilder = uninitialized
  private var axis: Axis = Axis.Cols

  @Setup(Level.Trial)
  def setup(): Unit =
    axis = if axisName == "rows" then Axis.Rows else Axis.Cols
    val base = Matrix.tabulate(n, n)((i, j) => i.toDouble * 0.125 - j.toDouble * 0.375)
    input = if layout == "transposed" then base.t else base
    center = Vec.tabulate(n)(i => i.toDouble * 0.125).along(axis)
    scale = Vec.tabulate(n)(i => 1.0 + i.toDouble * 0.01).along(axis)
    centerRows = center.toDVec.along(Axis.Rows)
    scaleRows = scale.toDVec.along(Axis.Rows)
    centerCols = center.toDVec.along(Axis.Cols)
    scaleCols = scale.toDVec.along(Axis.Cols)
    destination = DMatBuilder.zeros(n, n)
    val expected = explicitLoop()
    val actual = fused()
    for i <- 0 until n; j <- 0 until n do
      require(java.lang.Double.doubleToLongBits(expected(i, j)) == java.lang.Double.doubleToLongBits(actual(i, j)))

  @Benchmark def chained(): DMat = (input - center).pointwise / scale

  @Benchmark def fused(): DMat =
    if axis == Axis.Rows then input.standardize(centerRows, scaleRows)
    else input.standardize(centerCols, scaleCols)

  @Benchmark def into(): DMatBuilder =
    if axis == Axis.Rows then input.standardizeInto(centerRows, scaleRows, destination)
    else input.standardizeInto(centerCols, scaleCols, destination)
    destination

  @Benchmark def explicitLoop(): DMat =
    val out = DMatBuilder.zeros(input.rows, input.cols)
    var i = 0
    while i < input.rows do
      var j = 0
      while j < input.cols do
        val k = if axis == Axis.Rows then i else j
        out(i, j) = (input(i, j) - center(k)) / scale(k)
        j += 1
      i += 1
    out.result()
