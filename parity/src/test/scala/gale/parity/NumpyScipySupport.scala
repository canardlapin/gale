package gale.parity

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Matrix
import gale.linalg.Vec
import gale.sparse.CSR
import gale.sparse.Sparse

/** Comparison helpers for NumPy / SciPy reference fixtures.
  *
  * Same mixed absolute/relative rule as [[ParitySupport]]:
  * `|x − y| ≤ tol · max(1, |x|, |y|)`, with NaN and signed-infinity agreement.
  */
object NumpyScipySupport:

  def galeMatrix(data: Array[Array[Double]]): DMat =
    Matrix.tabulate(data.length, if data.isEmpty then 0 else data(0).length)((i, j) => data(i)(j))

  def galeVector(data: Array[Double]): DVec =
    Vec(data.toIndexedSeq*)

  def galeCsr(data: Array[Array[Double]]): CSR =
    val builder = Sparse.coo(data.length, data(0).length)
    var i = 0
    while i < data.length do
      var j = 0
      while j < data(0).length do
        if data(i)(j) != 0.0 then builder.add(i, j, data(i)(j))
        j += 1
      i += 1
    builder.toCSR()

  def isClose(x: Double, y: Double, tol: Double): Boolean =
    if x.isNaN && y.isNaN then true
    else if x.isPosInfinity && y.isPosInfinity then true
    else if x.isNegInfinity && y.isNegInfinity then true
    else if !x.isFinite || !y.isFinite then false
    else
      val scale = math.max(1.0, math.max(math.abs(x), math.abs(y)))
      math.abs(x - y) <= tol * scale

  def assertScalarClose(g: Double, ref: Double, tol: Double, clue: => String): Unit =
    if !isClose(g, ref, tol) then
      throw new AssertionError(s"$clue: gale=$g numpy/scipy=$ref |Δ|=${math.abs(g - ref)} tol=$tol")

  def assertVecClose(g: DVec, ref: Array[Double], tol: Double, clue: => String): Unit =
    if g.length != ref.length then
      throw new AssertionError(s"$clue: length gale=${g.length} numpy/scipy=${ref.length}")
    var i = 0
    while i < g.length do
      if !isClose(g(i), ref(i), tol) then
        throw new AssertionError(
          s"$clue: [$i] gale=${g(i)} numpy/scipy=${ref(i)} |Δ|=${math.abs(g(i) - ref(i))} tol=$tol"
        )
      i += 1

  def assertArrayClose(g: Array[Double], ref: Array[Double], tol: Double, clue: => String): Unit =
    if g.length != ref.length then
      throw new AssertionError(s"$clue: length gale=${g.length} numpy/scipy=${ref.length}")
    var i = 0
    while i < g.length do
      if !isClose(g(i), ref(i), tol) then
        throw new AssertionError(
          s"$clue: [$i] gale=${g(i)} numpy/scipy=${ref(i)} |Δ|=${math.abs(g(i) - ref(i))} tol=$tol"
        )
      i += 1

  def assertMatClose(g: DMat, ref: Array[Array[Double]], tol: Double, clue: => String): Unit =
    val rows = ref.length
    val cols = if rows == 0 then 0 else ref(0).length
    if g.rows != rows || g.cols != cols then
      throw new AssertionError(s"$clue: shape gale=${g.rows}x${g.cols} numpy/scipy=${rows}x$cols")
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        if !isClose(g(i, j), ref(i)(j), tol) then
          throw new AssertionError(
            s"$clue: ($i,$j) gale=${g(i, j)} numpy/scipy=${ref(i)(j)} |Δ|=${math.abs(g(i, j) - ref(i)(j))} tol=$tol"
          )
        j += 1
      i += 1

  def frobenius(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  def residual2(a: DMat, x: DVec, b: DVec): Double =
    (b - (a * x)).norm2
