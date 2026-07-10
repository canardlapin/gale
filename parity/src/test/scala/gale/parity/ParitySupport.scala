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
