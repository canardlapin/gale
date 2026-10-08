package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Named, kernel-backed elementwise numerics and log-domain reductions.
  *
  * Elementwise operations are explicit in Gale: `exp(A)` is never spelled as an
  * operator on the matrix itself. [[gale.syntax.all]]'s `A.pointwise.map(f)`
  * remains the general mechanism for an arbitrary function; the functions here
  * are its fast paths for common ones. `Numerics.exp(A)` returns '''exactly'''
  * the same values as `A.pointwise.map(math.exp)` (the same scalar function is
  * applied to each entry), but runs as a primitive loop without boxing.
  *
  * Every result is a newly owned contiguous vector or row-major matrix; inputs
  * may be any view (slice, transpose, strided). Scalar semantics are those of
  * `java.lang.Math` (for example `log` of a negative number is NaN), except
  * [[sigmoid]], which is evaluated in an overflow-free form.
  *
  * The log-domain functions use a max shift and so never overflow for finite
  * inputs:
  *
  *   - `logSumExp` of an empty input is `-Inf`. Otherwise NaN anywhere gives
  *     NaN; else any `+Inf` gives `+Inf`; else all `-Inf` gives `-Inf`.
  *   - `softmax` and `logSoftmax` return all-NaN for an input (or line) that
  *     contains NaN or `+Inf`, or whose entries are all `-Inf`; an empty input
  *     gives an empty result.
  *
  * Breeze naming trap: Breeze's `softmax(v)` returns the '''scalar'''
  * log-sum-exp. Its Gale spelling is [[logSumExp]]; Gale's [[softmax]] is the
  * normalized-exponential vector.
  */
object Numerics:

  def exp(x: DVec): DVec = mapVec(x)(DoubleKernels.dexpInto)
  def exp(a: DMat): DMat = mapMat(a)(DoubleKernels.dexpInto)

  def log(x: DVec): DVec = mapVec(x)(DoubleKernels.dlogInto)
  def log(a: DMat): DMat = mapMat(a)(DoubleKernels.dlogInto)

  /** `log(1 + x)`, accurate for small `|x|`. */
  def log1p(x: DVec): DVec = mapVec(x)(DoubleKernels.dlog1pInto)
  def log1p(a: DMat): DMat = mapMat(a)(DoubleKernels.dlog1pInto)

  /** `exp(x) - 1`, accurate for small `|x|`. */
  def expm1(x: DVec): DVec = mapVec(x)(DoubleKernels.dexpm1Into)
  def expm1(a: DMat): DMat = mapMat(a)(DoubleKernels.dexpm1Into)

  /** Logistic sigmoid `1 / (1 + exp(-x))`, evaluated as `1/(1+t)` for `x >= 0`
    * and `t/(1+t)` otherwise with `t = exp(-|x|)`, so it never overflows:
    * `sigmoid(+Inf) = 1`, `sigmoid(-Inf) = 0`, `sigmoid(NaN) = NaN`.
    */
  def sigmoid(x: DVec): DVec = mapVec(x)(DoubleKernels.dsigmoidInto)
  def sigmoid(a: DMat): DMat = mapMat(a)(DoubleKernels.dsigmoidInto)

  /** `log(sum exp(x_i))` over the whole vector, computed with a max shift. */
  def logSumExp(x: DVec): Double =
    DoubleKernels.dlogSumExp(x.length, x.data, x.offset.value, x.stride.value)

  /** `log(sum exp(a_ij))` over every entry of the matrix. */
  def logSumExp(a: DMat): Double =
    if a.isContiguousRowMajor || a.isContiguousColMajor then
      DoubleKernels.dlogSumExp(a.rows * a.cols, a.data, a.offset.value, 1)
    else logSumExp(DVec.fromDoubleArrayOwned(a.toDoubleArrayCopyRowMajor))

  /** Per-axis log-sum-exp: one value per row (`Axis.Rows`) or per column
    * (`Axis.Cols`). An empty line gives `-Inf`.
    */
  def logSumExp(a: DMat, axis: Axis): DVec =
    val lines = a.axisLines(axis)
    val length = a.axisLength(axis)
    val lineStep = a.axisLineStep(axis)
    val elementStep = a.axisElementStep(axis)
    val out = DVec.newBuilder(lines)
    var line = 0
    while line < lines do
      out(line) = DoubleKernels.dlogSumExp(length, a.data, a.offset.value + line * lineStep, elementStep)
      line += 1
    out.result()

  /** Normalized exponential `exp(x_i) / sum exp(x_j)`, computed with a max
    * shift; the result sums to one up to rounding.
    */
  def softmax(x: DVec): DVec = mapVec(x)(DoubleKernels.dsoftmaxInto)

  /** Softmax over every entry of the matrix jointly (the entries of the result
    * sum to one). Use the [[Axis]] overload for one distribution per line.
    */
  def softmax(a: DMat): DMat = wholeMat(a)(DoubleKernels.dsoftmaxInto)

  /** Softmax of each row (`Axis.Rows`) or each column (`Axis.Cols`)
    * independently; the result has the shape of `a`.
    */
  def softmax(a: DMat, axis: Axis): DMat = perLine(a, axis)(DoubleKernels.dsoftmaxInto)

  /** `x_i - logSumExp(x)`, evaluated as `(x_i - m) - log(sum exp(x_j - m))`. */
  def logSoftmax(x: DVec): DVec = mapVec(x)(DoubleKernels.dlogSoftmaxInto)

  /** Log-softmax over every entry of the matrix jointly. */
  def logSoftmax(a: DMat): DMat = wholeMat(a)(DoubleKernels.dlogSoftmaxInto)

  /** Log-softmax of each row (`Axis.Rows`) or each column (`Axis.Cols`). */
  def logSoftmax(a: DMat, axis: Axis): DMat = perLine(a, axis)(DoubleKernels.dlogSoftmaxInto)

  private type LineKernel = (Int, DoubleArray, Int, Int, DoubleArray, Int, Int) => Unit

  private inline def mapVec(x: DVec)(inline kernel: LineKernel): DVec =
    val out = DoubleArray.alloc(x.length)
    kernel(x.length, x.data, x.offset.value, x.stride.value, out, 0, 1)
    DVec.fromDoubleArrayOwned(out)

  /** Elementwise map into an owned row-major result: one call when `a` is
    * row-major contiguous, otherwise one strided call per row.
    */
  private inline def mapMat(a: DMat)(inline kernel: LineKernel): DMat =
    val rows = a.rows
    val cols = a.cols
    val out = DoubleArray.alloc(rows * cols)
    if a.isContiguousRowMajor then kernel(rows * cols, a.data, a.offset.value, 1, out, 0, 1)
    else
      var row = 0
      while row < rows do
        kernel(cols, a.data, a.offset.value + row * a.rowStride.value, a.colStride.value, out, row * cols, 1)
        row += 1
    DMat.fromDoubleArrayOwned(rows, cols, out)

  /** A whole-matrix (non-elementwise) line kernel over the row-major copy. */
  private inline def wholeMat(a: DMat)(inline kernel: LineKernel): DMat =
    val rows = a.rows
    val cols = a.cols
    val out = DoubleArray.alloc(rows * cols)
    if a.isContiguousRowMajor then kernel(rows * cols, a.data, a.offset.value, 1, out, 0, 1)
    else
      val source = a.toDoubleArrayCopyRowMajor
      kernel(rows * cols, source, 0, 1, out, 0, 1)
    DMat.fromDoubleArrayOwned(rows, cols, out)

  private inline def perLine(a: DMat, axis: Axis)(inline kernel: LineKernel): DMat =
    val rows = a.rows
    val cols = a.cols
    val out = DoubleArray.alloc(rows * cols)
    val lines = a.axisLines(axis)
    val length = a.axisLength(axis)
    val lineStep = a.axisLineStep(axis)
    val elementStep = a.axisElementStep(axis)
    // Output is row-major: a row line starts at `row * cols` with unit stride, a
    // column line starts at `col` with stride `cols`.
    val outLineStep = if axis == Axis.Rows then cols else 1
    val outElementStep = if axis == Axis.Rows then 1 else cols
    var line = 0
    while line < lines do
      kernel(length, a.data, a.offset.value + line * lineStep, elementStep, out, line * outLineStep, outElementStep)
      line += 1
    DMat.fromDoubleArrayOwned(rows, cols, out)
