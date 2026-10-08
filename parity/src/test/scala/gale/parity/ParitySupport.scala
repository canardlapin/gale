package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Matrix
import gale.linalg.Vec

/** Shared fixtures and comparison helpers for the Breeze parity suites.
  *
  * Every input is generated once as a plain `Array[Array[Double]]` / `Array[Double]`
  * and then handed to '''both''' libraries via identical `tabulate` closures, so
  * gale and Breeze see bit-for-bit identical data. Any disagreement is therefore a
  * genuine numerical difference between the two implementations, never a
  * construction artifact. Comparisons use a mixed absolute/relative tolerance:
  * `|x − y| ≤ tol · max(1, |x|, |y|)`.
  */
object ParitySupport:

  // ---------------------------------------------------------------------------
  // Data generation (library-agnostic, deterministic)
  // ---------------------------------------------------------------------------

  /** An `rows × cols` array of entries drawn uniformly from `[-1, 1)`. */
  def matrixData(rows: Int, cols: Int, seed: Long): Array[Array[Double]] =
    val rng = new scala.util.Random(seed)
    Array.tabulate(rows, cols)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  /** A length-`n` array of entries drawn uniformly from `[-1, 1)`. */
  def vectorData(n: Int, seed: Long): Array[Double] =
    val rng = new scala.util.Random(seed)
    Array.tabulate(n)(_ => rng.nextDouble() * 2.0 - 1.0)

  /** An `n × n` strictly diagonally dominant (hence well-conditioned, nonsingular)
    * matrix: off-diagonals in `[-1, 1)`, each diagonal set to the row's absolute
    * off-diagonal sum plus one.
    */
  def diagonallyDominant(n: Int, seed: Long): Array[Array[Double]] =
    val a = matrixData(n, n, seed)
    var i = 0
    while i < n do
      var sum = 0.0
      var j = 0
      while j < n do
        if j != i then sum += math.abs(a(i)(j))
        j += 1
      a(i)(i) = sum + 1.0
      i += 1
    a

  /** An `n × n` symmetric positive-definite matrix `B Bᵀ + n·I`, built into an
    * exactly-symmetric array (lower triangle computed, then mirrored) so Breeze's
    * `requireSymmetricMatrix` accepts it.
    */
  def spd(n: Int, seed: Long): Array[Array[Double]] =
    val b = matrixData(n, n, seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        var s = 0.0
        var k = 0
        while k < n do
          s += b(i)(k) * b(j)(k)
          k += 1
        if i == j then s += n.toDouble
        a(i)(j) = s
        a(j)(i) = s
        j += 1
      i += 1
    a

  /** A symmetric matrix with entries in `[-1, 1)` (lower triangle mirrored). */
  def symmetric(n: Int, seed: Long): Array[Array[Double]] =
    val src = matrixData(n, n, seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        a(i)(j) = src(i)(j)
        a(j)(i) = src(i)(j)
        j += 1
      i += 1
    a

  /** A symmetric matrix `Q diag(spectrum) Qᵀ` with a prescribed spectrum and a
    * random orthonormal `Q` (from a gale QR of a random matrix). Built into an
    * exactly-symmetric array. Lets a suite compare recovered eigenvalues against a
    * known reference and probe repeated/clustered eigenvalues.
    */
  def withSpectrum(spectrum: Array[Double], seed: Long): Array[Array[Double]] =
    val n = spectrum.length
    val q = orthonormal(n, seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        var s = 0.0
        var k = 0
        while k < n do
          s += q(i)(k) * spectrum(k) * q(j)(k)
          k += 1
        a(i)(j) = s
        a(j)(i) = s
        j += 1
      i += 1
    a

  /** An `n × n` orthonormal matrix as a row-major array, via gale's QR `Q`. */
  private def orthonormal(n: Int, seed: Long): Array[Array[Double]] =
    val q = galeMatrix(matrixData(n, n, seed)).qr.q
    Array.tabulate(n, n)((i, j) => q(i, j))

  // ---------------------------------------------------------------------------
  // Hardening fixtures: ill-conditioned inputs and non-contiguous views
  // ---------------------------------------------------------------------------

  /** Unit roundoff spacing `ε = 2⁻⁵²` used to scale hardening tolerances. */
  val Eps: Double = Math.ulp(1.0)

  /** The `n × n` Hilbert matrix `1 / (i + j + 1)`: SPD, `κ₂ ≈ e^{3.5n}`. */
  def hilbert(n: Int): Array[Array[Double]] =
    Array.tabulate(n, n)((i, j) => 1.0 / (i + j + 1).toDouble)

  /** A geometric spectrum from `1` down to `1 / kappa` (descending), length `n`. */
  def geometricSpectrum(n: Int, kappa: Double): Array[Double] =
    if n == 1 then Array(1.0)
    else Array.tabulate(n)(i => math.pow(kappa, -i.toDouble / (n - 1).toDouble))

  /** [[diagonallyDominant]] with each row divided by its diagonal: the unit
    * diagonal keeps `det` representable at large `n` while `κ` stays small.
    */
  def unitDiagonallyDominant(n: Int, seed: Long): Array[Array[Double]] =
    val a = diagonallyDominant(n, seed)
    Array.tabulate(n)(i => a(i).map(_ / a(i)(i)))

  def transpose(data: Array[Array[Double]]): Array[Array[Double]] =
    if data.isEmpty then data
    else Array.tabulate(data(0).length, data.length)((i, j) => data(j)(i))

  /** Embed `data` at `(top, left)` in a larger `NaN`-padded buffer, so a slice of
    * it is a genuinely strided view whose neighbours would poison any kernel that
    * read outside the window.
    */
  def embedded(data: Array[Array[Double]], top: Int, left: Int, bottom: Int, right: Int): Array[Array[Double]] =
    val rows = data.length
    val cols = if rows == 0 then 0 else data(0).length
    Array.tabulate(top + rows + bottom, left + cols + right): (i, j) =>
      val r = i - top
      val c = j - left
      if r >= 0 && r < rows && c >= 0 && c < cols then data(r)(c) else Double.NaN

  /** Three non-contiguous gale views of the logical matrix `data`: a transposed
    * view (unit row stride), a strided interior slice (row stride > cols), and a
    * transposed strided slice (neither stride is 1 along a row).
    */
  def galeViews(data: Array[Array[Double]]): List[(String, DMat)] =
    val rows = data.length
    val cols = if rows == 0 then 0 else data(0).length
    val transposedView = galeMatrix(transpose(data)).t
    val slice = galeMatrix(embedded(data, 2, 3, 1, 4)).slice(2, 2 + rows, 3, 3 + cols)
    val transposedSlice =
      galeMatrix(embedded(transpose(data), 1, 2, 3, 1)).slice(1, 1 + cols, 2, 2 + rows).t
    List("transposed" -> transposedView, "strided slice" -> slice, "transposed slice" -> transposedSlice)

  /** A stride-3 gale vector view of `data` (a column of a `NaN`-padded matrix). */
  def stridedVector(data: Array[Double]): DVec =
    galeMatrix(Array.tabulate(data.length, 3)((i, j) => if j == 1 then data(i) else Double.NaN)).col(1)

  // ---------------------------------------------------------------------------
  // Norm-wise error measures for conditioning-scaled comparisons
  // ---------------------------------------------------------------------------

  def maxAbs(data: Array[Array[Double]]): Double =
    data.foldLeft(0.0)((m, row) => row.foldLeft(m)((mm, x) => math.max(mm, math.abs(x))))

  /** `max |g − b| / max |b|`, the entrywise-max relative distance. */
  def relDiff(g: DVec, b: BDV[Double]): Double =
    if g.length != b.length then throw new AssertionError(s"length gale=${g.length} breeze=${b.length}")
    var num = 0.0
    var den = 0.0
    var i = 0
    while i < g.length do
      num = math.max(num, math.abs(g(i) - b(i)))
      den = math.max(den, math.abs(b(i)))
      i += 1
    if den == 0.0 then num else num / den

  def relDiff(g: DMat, b: BDM[Double]): Double =
    if g.rows != b.rows || g.cols != b.cols then
      throw new AssertionError(s"shape gale=${g.rows}x${g.cols} breeze=${b.rows}x${b.cols}")
    var num = 0.0
    var den = 0.0
    var i = 0
    while i < g.rows do
      var j = 0
      while j < g.cols do
        num = math.max(num, math.abs(g(i, j) - b(i, j)))
        den = math.max(den, math.abs(b(i, j)))
        j += 1
      i += 1
    if den == 0.0 then num else num / den

  /** Norm-wise backward error `‖b − A x‖∞ / (‖A‖∞ ‖x‖∞ + ‖b‖∞)` (Rigal–Gaches),
    * meaningful even where the forward error is not (`κ ε ≳ 1`).
    */
  def backwardError(a: Array[Array[Double]], x: Int => Double, b: Array[Double]): Double =
    val m = a.length
    val n = if m == 0 then 0 else a(0).length
    var resid = 0.0
    var aNorm = 0.0
    var i = 0
    while i < m do
      var s = 0.0
      var rowSum = 0.0
      var j = 0
      while j < n do
        s += a(i)(j) * x(j)
        rowSum += math.abs(a(i)(j))
        j += 1
      resid = math.max(resid, math.abs(b(i) - s))
      aNorm = math.max(aNorm, rowSum)
      i += 1
    val xNorm = (0 until n).foldLeft(0.0)((acc, j) => math.max(acc, math.abs(x(j))))
    val bNorm = b.foldLeft(0.0)((acc, v) => math.max(acc, math.abs(v)))
    resid / (aNorm * xNorm + bNorm)

  def assertBelow(value: Double, bound: Double, clue: => String): Unit =
    if !(value <= bound) then throw new AssertionError(s"$clue: value=$value bound=$bound")

  /** `‖G_c − B_c (B_cᵀ G_c)‖_F` for orthonormal eigenvector columns `cols` of
    * gale (`G`) and Breeze (`B`): the sine of the principal angles between the two
    * invariant subspaces (Frobenius), insensitive to signs and to rotations within
    * a cluster. Unlike `√(1 − cos²)` it has no `√ε` floor.
    */
  def subspaceSine(g: DMat, b: BDM[Double], cols: Range): Double =
    val n = g.rows
    val k = cols.length
    val proj = Array.ofDim[Double](k, k) // proj(p)(q) = ⟨b_p, g_q⟩
    for p <- 0 until k; q <- 0 until k do
      var s = 0.0
      var i = 0
      while i < n do
        s += b(i, cols(p)) * g(i, cols(q))
        i += 1
      proj(p)(q) = s
    var sum = 0.0
    var i = 0
    while i < n do
      var q = 0
      while q < k do
        var r = g(i, cols(q))
        var p = 0
        while p < k do
          r -= b(i, cols(p)) * proj(p)(q)
          p += 1
        sum += r * r
        q += 1
      i += 1
    math.sqrt(sum)

  /** Group ascending `values` into clusters whose internal gaps are `≤ gap`; each
    * cluster comes with its separation from the rest of the spectrum (`∞` if none).
    */
  def spectralClusters(values: IndexedSeq[Double], gap: Double): List[(Range, Double)] =
    val ranges = List.newBuilder[Range]
    var start = 0
    var i = 1
    while i < values.length do
      if values(i) - values(i - 1) > gap then
        ranges += (start until i)
        start = i
      i += 1
    if values.nonEmpty then ranges += (start until values.length)
    ranges.result().map: r =>
      val below = if r.start > 0 then values(r.start) - values(r.start - 1) else Double.PositiveInfinity
      val above = if r.end < values.length then values(r.end) - values(r.end - 1) else Double.PositiveInfinity
      (r, math.min(below, above))

  // ---------------------------------------------------------------------------
  // Conversions to each library
  // ---------------------------------------------------------------------------

  def galeMatrix(data: Array[Array[Double]]): DMat =
    Matrix.tabulate(data.length, if data.isEmpty then 0 else data(0).length)((i, j) => data(i)(j))

  def breezeMatrix(data: Array[Array[Double]]): BDM[Double] =
    BDM.tabulate(data.length, if data.isEmpty then 0 else data(0).length)((i, j) => data(i)(j))

  def galeVector(data: Array[Double]): DVec =
    Vec(data.toIndexedSeq*)

  def breezeVector(data: Array[Double]): BDV[Double] =
    BDV(data.clone())

  // ---------------------------------------------------------------------------
  // Comparison assertions (throw on failure; suites call from within `test`)
  // ---------------------------------------------------------------------------

  private def isClose(x: Double, y: Double, tol: Double): Boolean =
    if x.isNaN && y.isNaN then true
    else if x.isPosInfinity && y.isPosInfinity then true
    else if x.isNegInfinity && y.isNegInfinity then true
    else if !x.isFinite || !y.isFinite then false
    else
      val scale = math.max(1.0, math.max(math.abs(x), math.abs(y)))
      math.abs(x - y) <= tol * scale

  def assertScalarClose(g: Double, b: Double, tol: Double, clue: => String): Unit =
    if !isClose(g, b, tol) then
      throw new AssertionError(s"$clue: gale=$g breeze=$b |Δ|=${math.abs(g - b)} tol=$tol")

  def assertVecClose(g: DVec, b: BDV[Double], tol: Double, clue: => String): Unit =
    if g.length != b.length then
      throw new AssertionError(s"$clue: length gale=${g.length} breeze=${b.length}")
    var i = 0
    while i < g.length do
      if !isClose(g(i), b(i), tol) then
        throw new AssertionError(s"$clue: [$i] gale=${g(i)} breeze=${b(i)} |Δ|=${math.abs(g(i) - b(i))} tol=$tol")
      i += 1

  def assertMatClose(g: DMat, b: BDM[Double], tol: Double, clue: => String): Unit =
    if g.rows != b.rows || g.cols != b.cols then
      throw new AssertionError(s"$clue: shape gale=${g.rows}x${g.cols} breeze=${b.rows}x${b.cols}")
    var i = 0
    while i < g.rows do
      var j = 0
      while j < g.cols do
        if !isClose(g(i, j), b(i, j), tol) then
          throw new AssertionError(s"$clue: ($i,$j) gale=${g(i, j)} breeze=${b(i, j)} |Δ|=${math.abs(g(i, j) - b(i, j))} tol=$tol")
        j += 1
      i += 1

  /** Compare two Breeze matrices (used for gale-side reconstructions expressed as
    * Breeze products, or Breeze-vs-Breeze invariants).
    */
  def assertBreezeMatClose(x: BDM[Double], y: BDM[Double], tol: Double, clue: => String): Unit =
    if x.rows != y.rows || x.cols != y.cols then
      throw new AssertionError(s"$clue: shape ${x.rows}x${x.cols} vs ${y.rows}x${y.cols}")
    var i = 0
    while i < x.rows do
      var j = 0
      while j < x.cols do
        if !isClose(x(i, j), y(i, j), tol) then
          throw new AssertionError(s"$clue: ($i,$j) ${x(i, j)} vs ${y(i, j)} |Δ|=${math.abs(x(i, j) - y(i, j))} tol=$tol")
        j += 1
      i += 1
