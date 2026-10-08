package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import breeze.linalg.`*` as Every
import breeze.linalg.argmax as bArgmax
import breeze.linalg.argmin as bArgmin
import breeze.linalg.max as bMax
import breeze.linalg.min as bMin
import breeze.linalg.norm as bNorm
import breeze.linalg.softmax as bLogSumExp
import breeze.linalg.sum as bSum
import breeze.numerics.abs as bAbs
import breeze.numerics.exp as bExp
import breeze.numerics.expm1 as bExpm1
import breeze.numerics.log as bLog
import breeze.numerics.log1p as bLog1p
import breeze.numerics.sigmoid as bSigmoid
import breeze.stats.mean as bMean
import gale.linalg.*
import gale.parity.ParitySupport.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAllNoShrink

/** Breeze differential tests for the everyday reductions, norms and elementwise
  * numerics: `sum`, `sumExact`, `mean`, `max`, `min`, `argmax`, `argmin` (whole
  * and per [[Axis]]), vector and matrix norms, and [[Numerics]].
  *
  * Every property runs over contiguous, strided, sliced and transposed Gale views
  * of the same logical data, so layout-dependent kernels are compared against
  * the same Breeze reference. Sizes are `1..64`; fixed `n = 1K` and `64K` cases
  * cover the unrolled kernels. Replay a failure with [[ParitySeed]].
  *
  * Tolerances (`ε = 2^-52`):
  *   - `max`, `min`, `argmax`, `argmin`, `normInf`, and the per-axis extrema:
  *     exact. They select an entry and do no arithmetic.
  *   - `exp`, `log`, `log1p`, `expm1`: exact. Gale and `breeze.numerics` both
  *     apply `java.lang.Math` to each entry.
  *   - `sigmoid`: 4 ulps. Breeze evaluates `1/(1+exp(-x))`; Gale uses
  *     `exp(x)/(1+exp(x))` for `x < 0` to avoid overflow.
  *   - `sum`, `mean`, per-axis sums: `|Δ| ≤ 2 n ε Σ|x|` (divided by `n` for a
  *     mean), the forward bound of two differently associated summations. As a
  *     relative error this is `2 n ε` times the summation condition number
  *     `Σ|x| / |Σx|`. `sumExact` is additionally compared bit for bit with the
  *     correctly rounded `BigDecimal` sum.
  *   - `norm1`, `norm2`, Frobenius, matrix 1/∞ norms: relative `4 n ε`, since
  *     every term is non-negative and the summation is well conditioned.
  *   - `logSumExp`, `logSoftmax`: absolute `4 (n + 2) ε max(1, |r|)`; softmax
  *     probabilities: absolute `4 (n + max|x|) ε`, because the Breeze reference
  *     `exp(x - softmax(x))` loses `ε |x|` in the subtraction.
  *
  * Divergences from Breeze 2.1.0 are asserted as Gale's documented semantics,
  * with Breeze's observed behaviour pinned beside it (tests named
  * `divergence: …`); they are listed in the Breeze migration guide.
  *   1. `argmax`/`argmin` with NaN: Gale returns the first NaN; Breeze skips NaN
  *      unless it is the first entry.
  *   1. `argmin` ties: Gale returns the first minimum, Breeze the last.
  *   1. Matrix `argmax`/`argmin`: Gale breaks ties and locates NaN in row-major
  *      order; Breeze breaks ties in column-major order and skips NaN.
  *   1. Signed zeros: Gale's `max`/`min` return the first of `0.0`/`-0.0`;
  *      Breeze's `max` returns `0.0` and `min` returns `-0.0`.
  *   1. Empty input: Gale's `max`, `mean`, per-axis `max` throw
  *      `LinAlgError.EmptyInput`; Breeze returns `-Inf`, `0.0`, `-Inf`. `min` and
  *      `argmax` throw in both (Breeze: `IllegalArgumentException`).
  *   1. `mean` overflow: Gale is `sum / n` (NumPy), so `mean(MaxValue, MaxValue)`
  *      is `+Inf`; Breeze's running mean returns `MaxValue`.
  *   1. `norm2`/Frobenius: Gale scales, so `1e300` entries give a finite norm and
  *      `1e-300` entries a nonzero one; Breeze overflows to `+Inf` and underflows
  *      to `0.0`.
  *   1. `logSumExp` with `+Inf`: Gale gives `+Inf`; Breeze's `softmax` gives
  *      `-Inf`.
  *   1. `sigmoid` below `-log(Double.MaxValue) ≈ -709.78`: Breeze's
  *      `1/(1+exp(-x))` overflows to exactly `0.0`; Gale returns the subnormal
  *      `≈ exp(x)` down to `x ≈ -745`.
  */
class ReductionsNumericsParitySuite extends ScalaCheckSuite:
  override def scalaCheckInitialSeed =
    ParitySeed.initial("QnkrVBoujb2NO7KGCM4XAAba3j4fCn0AeJaFKoZ1tcE=")

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(40).withWorkers(1)

  private val Eps = math.ulp(1.0)
  private val sizeGen = Gen.choose(1, 64)
  private val seedGen = Gen.choose(1L, 9_000_000L)
  private val vectorCaseGen: Gen[(Int, Long)] = Gen.zip(sizeGen, seedGen)
  private val matrixCaseGen: Gen[(Int, Int, Long)] = Gen.zip(sizeGen, sizeGen, seedGen)
  private val logScaleGen: Gen[Double] = Gen.oneOf(1.0, 30.0, 1000.0)

  // ---------------------------------------------------------------------------
  // Gale views of one logical vector / matrix
  // ---------------------------------------------------------------------------

  /** Contiguous, strided (a column of a row-major matrix), offset slice, and a
    * strided row of a transposed view: all hold exactly `data`.
    */
  private def vectorViews(data: Array[Double]): Seq[(String, DVec)] =
    val n = data.length
    Seq(
      "contiguous" -> galeVector(data),
      "strided" -> Matrix.tabulate(n, 3)((i, j) => if j == 1 then data(i) else 99.0).col(1),
      "sliced" -> Vec((Array(7.0, 7.0) ++ data ++ Array(7.0)).toIndexedSeq*).slice(2, 2 + n),
      "transposed-row" -> Matrix.tabulate(n, 2)((i, j) => if j == 0 then data(i) else -5.0).t.row(0)
    )

  /** Row-major, column-major (transposed view), interior slice of a padded
    * matrix, and a transposed slice (non-unit strides on both axes).
    */
  private def matrixViews(data: Array[Array[Double]]): Seq[(String, DMat)] =
    val rows = data.length
    val cols = data(0).length
    val transposed = Array.tabulate(cols, rows)((i, j) => data(j)(i))
    val padded = Matrix.tabulate(rows + 2, cols + 3)((i, j) =>
      if i >= 1 && i < rows + 1 && j >= 2 && j < cols + 2 then data(i - 1)(j - 2) else 42.0
    )
    val paddedTransposed = Matrix.tabulate(cols + 3, rows + 1)((i, j) =>
      if i >= 2 && i < cols + 2 && j < rows then transposed(i - 2)(j) else -42.0
    )
    Seq(
      "row-major" -> galeMatrix(data),
      "col-major" -> galeMatrix(transposed).t,
      "sliced" -> padded.slice(1, rows + 1, 2, cols + 2),
      "sliced-transposed" -> paddedTransposed.slice(2, cols + 2, 0, rows).t
    )

  // ---------------------------------------------------------------------------
  // Tolerance helpers
  // ---------------------------------------------------------------------------

  private def absSum(xs: Iterable[Double]): Double = xs.foldLeft(0.0)((a, x) => a + math.abs(x))

  private def assertSumClose(g: Double, b: Double, n: Int, sumAbs: Double, clue: => String): Unit =
    val bound = 2.0 * math.max(n, 1) * Eps * sumAbs
    if !(math.abs(g - b) <= bound) then
      fail(s"$clue: gale=$g breeze=$b |Δ|=${math.abs(g - b)} bound=2nεΣ|x|=$bound")

  private def assertRelClose(g: Double, b: Double, rel: Double, clue: => String): Unit =
    if !(g == b || math.abs(g - b) <= rel * math.max(math.abs(g), math.abs(b))) then
      fail(s"$clue: gale=$g breeze=$b relΔ=${math.abs(g - b) / math.max(math.abs(g), math.abs(b))} tol=$rel")

  private def assertUlps(g: Double, b: Double, ulps: Int, clue: => String): Unit =
    val same = g == b || (g.isNaN && b.isNaN)
    if !same && !(math.abs(g - b) <= ulps * math.ulp(math.max(math.abs(g), math.abs(b)))) then
      fail(s"$clue: gale=$g breeze=$b |Δ|=${math.abs(g - b)} > $ulps ulp")

  private def assertExact(g: Double, b: Double, clue: => String): Unit =
    if !(java.lang.Double.compare(g, b) == 0 || (g == b)) then fail(s"$clue: gale=$g breeze=$b (exact)")

  /** Breeze's `1/(1+exp(-x))` overflows to exactly `0.0` once `exp(-x)` is
    * infinite; below that point Gale is compared with `exp(x)`, which equals
    * `sigmoid(x)` to within an ulp there (divergence 9).
    */
  private val BreezeSigmoidFloor = -math.log(Double.MaxValue)

  private def assertSigmoid(g: Double, x: Double, clue: => String): Unit =
    if x < BreezeSigmoidFloor then
      assertExact(bSigmoid(x), 0.0, s"breeze $clue")
      assertUlps(g, math.exp(x), 4, s"$clue vs exp(x)")
    else assertUlps(g, bSigmoid(x), 4, clue)

  /** The exact sum of `xs`, rounded once to the nearest `Double`. */
  private def correctlyRoundedSum(xs: Iterable[Double]): Double =
    xs.foldLeft(java.math.BigDecimal.ZERO)((acc, x) => acc.add(new java.math.BigDecimal(x))).doubleValue

  private def bvec(v: DVec): BDV[Double] = BDV.tabulate(v.length)(v(_))

  // ---------------------------------------------------------------------------
  // Vector reductions and norms
  // ---------------------------------------------------------------------------

  private def checkVectorReductions(data: Array[Double], clue: String): Unit =
    val n = data.length
    val b = breezeVector(data)
    val sumAbs = absSum(data)
    val exact = correctlyRoundedSum(data)
    for (layout, g) <- vectorViews(data) do
      val c = s"$clue $layout"
      assertSumClose(g.sum, bSum(b), n, sumAbs, s"sum $c")
      assertExact(g.sumExact, exact, s"sumExact vs correctly rounded $c")
      assertSumClose(g.sumExact, bSum(b), n, sumAbs, s"sumExact vs breeze $c")
      assertSumClose(g.mean, bMean(b), n, sumAbs / n, s"mean $c")
      assertExact(g.max, bMax(b), s"max $c")
      assertExact(g.min, bMin(b), s"min $c")
      assertEquals(g.argmax, bArgmax(b), s"argmax $c")
      assertEquals(g.argmin, bArgmin(b), s"argmin $c")
      assertRelClose(g.norm1, bNorm(b, 1.0), 4.0 * n * Eps, s"norm1 $c")
      assertRelClose(g.norm2, bNorm(b), 4.0 * n * Eps, s"norm2 $c")
      assertExact(g.normInf, bNorm(b, Double.PositiveInfinity), s"normInf $c")

  property("vector sum, sumExact, mean, extrema, arg-extrema and norms match Breeze") {
    forAllNoShrink(vectorCaseGen) { (sample: (Int, Long)) =>
      val (n, seed) = sample
      checkVectorReductions(vectorData(n, seed), s"n=$n seed=$seed")
    }
  }

  test("vector reductions at n = 1K and 64K") {
    for n <- Seq(1024, 65536) do
      checkVectorReductions(vectorData(n, 17L * n), s"n=$n")
      // Cancellation-heavy: large alternating entries plus small noise.
      val noisy = vectorData(n, 31L * n).zipWithIndex.map((x, i) => (if i % 2 == 0 then 1e8 else -1e8) + x)
      checkVectorReductions(noisy, s"cancelling n=$n")
  }

  // ---------------------------------------------------------------------------
  // Matrix reductions, per-axis forms, and norms
  // ---------------------------------------------------------------------------

  private def checkMatrixReductions(data: Array[Array[Double]], clue: String): Unit =
    val rows = data.length
    val cols = data(0).length
    val n = rows * cols
    val b = breezeMatrix(data)
    val flat = data.flatten
    val sumAbs = absSum(flat)
    val exact = correctlyRoundedSum(flat)
    val absB: BDM[Double] = bAbs(b)

    val colSums: BDV[Double] = bSum(b(::, Every)).t
    val rowSums: BDV[Double] = bSum(b(Every, ::))
    val colMeans: BDV[Double] = bMean(b(::, Every)).t
    val rowMeans: BDV[Double] = bMean(b(Every, ::))
    val colMax: BDV[Double] = bMax(b(::, Every)).t
    val rowMax: BDV[Double] = bMax(b(Every, ::))
    val colMin: BDV[Double] = bMin(b(::, Every)).t
    val rowMin: BDV[Double] = bMin(b(Every, ::))
    val absColSums: BDV[Double] = bSum(absB(::, Every)).t
    val absRowSums: BDV[Double] = bSum(absB(Every, ::))
    val norm1Ref: Double = bMax(absColSums)
    val normInfRef: Double = bMax(absRowSums)
    val frobRef = bNorm(b.toDenseVector)

    for (layout, g) <- matrixViews(data) do
      val c = s"$clue $layout"
      assertSumClose(g.sum, bSum(b), n, sumAbs, s"sum $c")
      assertExact(g.sumExact, exact, s"sumExact vs correctly rounded $c")
      assertSumClose(g.mean, bMean(b), n, sumAbs / n, s"mean $c")
      assertExact(g.max, bMax(b), s"max $c")
      assertExact(g.min, bMin(b), s"min $c")
      assertEquals(g.argmax, bArgmax(b), s"argmax $c")
      assertEquals(g.argmin, bArgmin(b), s"argmin $c")

      val gColSums = g.sum(Axis.Cols)
      val gRowSums = g.sum(Axis.Rows)
      val gColMeans = g.mean(Axis.Cols)
      val gRowMeans = g.mean(Axis.Rows)
      assertEquals(gColSums.length, cols, s"sum(Cols) length $c")
      assertEquals(gRowSums.length, rows, s"sum(Rows) length $c")
      for j <- 0 until cols do
        val lineAbs = absSum((0 until rows).map(data(_)(j)))
        assertSumClose(gColSums(j), colSums(j), rows, lineAbs, s"sum(Cols)[$j] $c")
        assertSumClose(gColMeans(j), colMeans(j), rows, lineAbs / rows, s"mean(Cols)[$j] $c")
        assertExact(g.max(Axis.Cols)(j), colMax(j), s"max(Cols)[$j] $c")
        assertExact(g.min(Axis.Cols)(j), colMin(j), s"min(Cols)[$j] $c")
      for i <- 0 until rows do
        val lineAbs = absSum(data(i))
        assertSumClose(gRowSums(i), rowSums(i), cols, lineAbs, s"sum(Rows)[$i] $c")
        assertSumClose(gRowMeans(i), rowMeans(i), cols, lineAbs / cols, s"mean(Rows)[$i] $c")
        assertExact(g.max(Axis.Rows)(i), rowMax(i), s"max(Rows)[$i] $c")
        assertExact(g.min(Axis.Rows)(i), rowMin(i), s"min(Rows)[$i] $c")

      assertRelClose(g.norm1, norm1Ref, 4.0 * rows * Eps, s"norm1 (max abs column sum) $c")
      assertRelClose(g.normInf, normInfRef, 4.0 * cols * Eps, s"normInf (max abs row sum) $c")
      assertRelClose(g.normFrobenius, frobRef, 4.0 * n * Eps, s"normFrobenius $c")

  property("matrix sum, mean, extrema, per-axis reductions and norms match Breeze") {
    forAllNoShrink(matrixCaseGen) { (sample: (Int, Int, Long)) =>
      val (rows, cols, seed) = sample
      checkMatrixReductions(matrixData(rows, cols, seed), s"${rows}x$cols seed=$seed")
    }
  }

  test("matrix reductions at 32x32 (1K) and 256x256 (64K)") {
    for n <- Seq(32, 256) do checkMatrixReductions(matrixData(n, n, 23L * n), s"${n}x$n")
    checkMatrixReductions(matrixData(1, 1024, 5L), "1x1024")
    checkMatrixReductions(matrixData(1024, 1, 6L), "1024x1")
  }

  // ---------------------------------------------------------------------------
  // Elementwise numerics and log-domain functions
  // ---------------------------------------------------------------------------

  private def checkElementwise(data: Array[Double], clue: String): Unit =
    val b = breezeVector(data)
    val positive = data.map(x => math.abs(x) + 1e-3)
    val bPositive = breezeVector(positive)
    val aboveMinusOne = data.map(x => x * 0.999)
    val bAboveMinusOne = breezeVector(aboveMinusOne)
    val viewsWithPositive = vectorViews(positive).map(_._2)
    val viewsAboveMinusOne = vectorViews(aboveMinusOne).map(_._2)
    for (((layout, g), gPos), gAbove) <- vectorViews(data).zip(viewsWithPositive).zip(viewsAboveMinusOne) do
      val c = s"$clue $layout"
      val ge = Numerics.exp(g)
      val gs = Numerics.sigmoid(g)
      val gm = Numerics.expm1(g)
      val gl = Numerics.log(gPos)
      val gl1 = Numerics.log1p(gAbove)
      val be = bExp(b)
      val bm = bExpm1(b)
      val bl = bLog(bPositive)
      val bl1 = bLog1p(bAboveMinusOne)
      for i <- data.indices do
        assertExact(ge(i), be(i), s"exp[$i] $c")
        assertExact(gm(i), bm(i), s"expm1[$i] $c")
        assertExact(gl(i), bl(i), s"log[$i] $c")
        assertExact(gl1(i), bl1(i), s"log1p[$i] $c")
        assertSigmoid(gs(i), data(i), s"sigmoid[$i] $c")

  private def checkLogDomain(data: Array[Double], clue: String): Unit =
    val n = data.length
    val b = breezeVector(data)
    val lse = bLogSumExp(b)
    val maxAbs = data.foldLeft(0.0)((a, x) => math.max(a, math.abs(x)))
    val lseTol = 4.0 * (n + 2) * Eps
    val probTol = 4.0 * (n + maxAbs) * Eps
    for (layout, g) <- vectorViews(data) do
      val c = s"$clue $layout"
      assertScalarClose(Numerics.logSumExp(g), lse, lseTol, s"logSumExp vs Breeze softmax $c")
      val p = Numerics.softmax(g)
      val lp = Numerics.logSoftmax(g)
      for i <- 0 until n do
        val ref = math.exp(data(i) - lse)
        if !(math.abs(p(i) - ref) <= probTol) then
          fail(s"softmax[$i] $c: gale=${p(i)} breeze exp(x - softmax(x))=$ref tol=$probTol")
        assertScalarClose(lp(i), data(i) - lse, probTol, s"logSoftmax[$i] $c")
      assertSumClose(p.sum, 1.0, n, 1.0, s"softmax sums to one $c")

  property("exp, log, log1p, expm1 and sigmoid match breeze.numerics") {
    forAllNoShrink(vectorCaseGen) { (sample: (Int, Long)) =>
      val (n, seed) = sample
      checkElementwise(vectorData(n, seed).map(_ * 40.0), s"n=$n seed=$seed")
    }
  }

  property("logSumExp matches Breeze softmax(v) and softmax matches exp(v - softmax(v))") {
    forAllNoShrink(vectorCaseGen, logScaleGen) { (sample: (Int, Long), scale: Double) =>
      val (n, seed) = sample
      checkLogDomain(vectorData(n, seed).map(_ * scale), s"n=$n seed=$seed scale=$scale")
    }
  }

  test("elementwise and log-domain functions at n = 1K and 64K") {
    for n <- Seq(1024, 65536) do
      checkElementwise(vectorData(n, 3L * n).map(_ * 40.0), s"n=$n")
      checkLogDomain(vectorData(n, 7L * n).map(_ * 30.0), s"n=$n")
  }

  private def checkMatrixLogDomain(data: Array[Array[Double]], clue: String): Unit =
    val rows = data.length
    val cols = data(0).length
    val b = breezeMatrix(data)
    val maxAbs = data.flatten.foldLeft(0.0)((a, x) => math.max(a, math.abs(x)))
    val rowLse: BDV[Double] = bLogSumExp(b(Every, ::))
    val colLse: BDV[Double] = bLogSumExp(b(::, Every)).t
    val allLse = bLogSumExp(b.toDenseVector)
    for (layout, g) <- matrixViews(data) do
      val c = s"$clue $layout"
      assertScalarClose(Numerics.logSumExp(g), allLse, 4.0 * (rows * cols + 2) * Eps, s"logSumExp(all) $c")
      assertVecClose(Numerics.logSumExp(g, Axis.Rows), rowLse, 4.0 * (cols + 2) * Eps, s"logSumExp(Rows) $c")
      assertVecClose(Numerics.logSumExp(g, Axis.Cols), colLse, 4.0 * (rows + 2) * Eps, s"logSumExp(Cols) $c")
      val pRows = Numerics.softmax(g, Axis.Rows)
      val pCols = Numerics.softmax(g, Axis.Cols)
      val pAll = Numerics.softmax(g)
      val lpRows = Numerics.logSoftmax(g, Axis.Rows)
      val lpCols = Numerics.logSoftmax(g, Axis.Cols)
      val gExp = Numerics.exp(g)
      val gSig = Numerics.sigmoid(g)
      val rowTol = 4.0 * (cols + maxAbs) * Eps
      val colTol = 4.0 * (rows + maxAbs) * Eps
      val allTol = 4.0 * (rows * cols + maxAbs) * Eps
      for i <- 0 until rows; j <- 0 until cols do
        val x = data(i)(j)
        assertScalarClose(pRows(i, j), math.exp(x - rowLse(i)), rowTol, s"softmax(Rows)($i,$j) $c")
        assertScalarClose(pCols(i, j), math.exp(x - colLse(j)), colTol, s"softmax(Cols)($i,$j) $c")
        assertScalarClose(pAll(i, j), math.exp(x - allLse), allTol, s"softmax(all)($i,$j) $c")
        assertScalarClose(lpRows(i, j), x - rowLse(i), rowTol, s"logSoftmax(Rows)($i,$j) $c")
        assertScalarClose(lpCols(i, j), x - colLse(j), colTol, s"logSoftmax(Cols)($i,$j) $c")
        assertExact(gExp(i, j), math.exp(x), s"exp($i,$j) $c")
        assertSigmoid(gSig(i, j), x, s"sigmoid($i,$j) $c")

  property("matrix logSumExp, softmax and logSoftmax (whole and per axis) match Breeze") {
    forAllNoShrink(matrixCaseGen, logScaleGen) { (sample: (Int, Int, Long), scale: Double) =>
      val (rows, cols, seed) = sample
      checkMatrixLogDomain(
        matrixData(rows, cols, seed).map(_.map(_ * scale)),
        s"${rows}x$cols seed=$seed scale=$scale"
      )
    }
  }

  // ---------------------------------------------------------------------------
  // Explicit edge cases where Gale and Breeze agree
  // ---------------------------------------------------------------------------

  private val Inf = Double.PositiveInfinity
  private val NaN = Double.NaN

  test("length-1 inputs agree for every reduction and log-domain function") {
    for x <- Seq(-2.5, 0.0, 3.0, 1000.0, -1000.0) do
      val g = Vec(x)
      val b = BDV(x)
      assertExact(g.sum, bSum(b), s"sum $x")
      assertExact(g.mean, bMean(b), s"mean $x")
      assertExact(g.max, bMax(b), s"max $x")
      assertExact(g.min, bMin(b), s"min $x")
      assertEquals(g.argmax, 0)
      assertEquals(bArgmax(b), 0)
      assertExact(g.norm1, bNorm(b, 1.0), s"norm1 $x")
      assertExact(g.norm2, bNorm(b), s"norm2 $x")
      assertExact(Numerics.logSumExp(g), bLogSumExp(b), s"logSumExp $x")
      assertExact(Numerics.softmax(g)(0), 1.0, s"softmax $x")
      val m = Matrix(1, 1)(x)
      assertExact(m.normFrobenius, bNorm(BDV(x)), s"frobenius $x")
      assertEquals(m.argmax, (0, 0))
  }

  test("infinities follow IEEE arithmetic in both libraries") {
    val withInf = BDV(1.0, Inf, -3.0)
    val g = Vec(1.0, Inf, -3.0)
    assertExact(g.sum, bSum(withInf), "sum with +Inf")
    assertExact(g.max, bMax(withInf), "max with +Inf")
    assertEquals(g.argmax, bArgmax(withInf))
    assertExact(g.norm1, bNorm(withInf, 1.0), "norm1 with +Inf")
    assertExact(g.normInf, bNorm(withInf, Inf), "normInf with +Inf")
    assertExact(Vec(1.0, -Inf).normInf, bNorm(BDV(1.0, -Inf), Inf), "normInf with -Inf")
    assertExact(Vec(Inf, -Inf).norm2, bNorm(BDV(Inf, -Inf)), "norm2 of two infinities")
    assert(Vec(Inf, -Inf).sum.isNaN && bSum(BDV(Inf, -Inf)).isNaN, "+Inf + -Inf is NaN")
    assert(Vec(Inf, -Inf).sumExact.isNaN, "sumExact +Inf + -Inf is NaN")
    // Ties at +Inf: both return the first occurrence.
    assertEquals(Vec(1.0, Inf, Inf).argmax, 1)
    assertEquals(bArgmax(BDV(1.0, Inf, Inf)), 1)
  }

  test("NaN propagates through sums, max, min and norms in both libraries") {
    for at <- 0 until 3 do
      val data = Array(1.0, 2.0, 3.0)
      data(at) = NaN
      val g = galeVector(data)
      val b = breezeVector(data)
      val c = s"NaN at $at"
      assert(g.sum.isNaN && bSum(b).isNaN, s"sum $c")
      assert(g.mean.isNaN && bMean(b).isNaN, s"mean $c")
      assert(g.max.isNaN && bMax(b).isNaN, s"max $c")
      assert(g.min.isNaN && bMin(b).isNaN, s"min $c")
      assert(g.norm1.isNaN && bNorm(b, 1.0).isNaN, s"norm1 $c")
      assert(g.norm2.isNaN && bNorm(b).isNaN, s"norm2 $c")
      assert(g.normInf.isNaN && bNorm(b, Inf).isNaN, s"normInf $c")
      assert(Numerics.logSumExp(g).isNaN && bLogSumExp(b).isNaN, s"logSumExp $c")
      val matrix = Matrix(1, 3)(data*)
      assert(matrix.max.isNaN && bMax(breezeMatrix(Array(data))).isNaN, s"matrix max $c")
      assert(matrix.normFrobenius.isNaN, s"frobenius $c")
    assert(Numerics.sigmoid(Vec(NaN))(0).isNaN && bSigmoid(NaN).isNaN, "sigmoid NaN")
  }

  test("argmax ties, all -Inf, and log-sum-exp of all -Inf agree") {
    assertEquals(Vec(1.0, 3.0, 3.0, 2.0).argmax, 1)
    assertEquals(bArgmax(BDV(1.0, 3.0, 3.0, 2.0)), 1)
    val allNegInf = Vec(-Inf, -Inf, -Inf)
    val bAllNegInf = BDV(-Inf, -Inf, -Inf)
    assertExact(allNegInf.max, bMax(bAllNegInf), "max all -Inf")
    assertEquals(allNegInf.argmax, bArgmax(bAllNegInf))
    assertExact(Numerics.logSumExp(allNegInf), bLogSumExp(bAllNegInf), "logSumExp all -Inf")
    // softmax of all -Inf: Gale documents all-NaN; Breeze's exp(v - softmax(v))
    // is exp(-Inf + Inf) = NaN as well.
    val p = Numerics.softmax(allNegInf)
    val bp: BDV[Double] = bExp(bAllNegInf - bLogSumExp(bAllNegInf))
    for i <- 0 until 3 do assert(p(i).isNaN && bp(i).isNaN, s"softmax all -Inf [$i]")
    // A single -Inf among finite entries contributes probability exactly 0.
    assertExact(Numerics.logSumExp(Vec(-Inf, 1.0)), bLogSumExp(BDV(-Inf, 1.0)), "logSumExp(-Inf, 1)")
    assertExact(Numerics.softmax(Vec(-Inf, 1.0))(0), 0.0, "softmax(-Inf, 1)[0]")
  }

  test("logSumExp and softmax stay finite for inputs of magnitude 1000") {
    for data <- Seq(Array(1000.0, 1000.0, -1000.0), Array(-1000.0, -1000.0), Array(1000.0, 999.0, 0.0, -1000.0)) do
      checkLogDomain(data, data.mkString("[", ",", "]"))
    assertScalarClose(Numerics.logSumExp(Vec(1000.0, 1000.0)), 1000.0 + math.log(2.0), 4 * Eps, "lse(1000,1000)")
    assertExact(Numerics.softmax(Vec(1000.0, 1000.0, -1000.0))(2), 0.0, "softmax underflow to 0")
  }

  test("signed zeros compare equal and argmax returns the first") {
    assertEquals(Vec(-0.0, 0.0).argmax, 0)
    assertEquals(bArgmax(BDV(-0.0, 0.0)), 0)
    assertEquals(Vec(0.0, -0.0).argmin, 0)
    assertExact(Vec(-0.0, 0.0).sum, bSum(BDV(-0.0, 0.0)), "sum of signed zeros")
    assert(Vec(-0.0, 0.0).max == 0.0, "max of signed zeros is zero")
  }

  test("sigmoid saturates without overflow at extreme arguments") {
    for x <- Seq(-800.0, -745.0, -700.0, -40.0, 0.0, 40.0, 800.0, -Inf, Inf) do
      assertSigmoid(Numerics.sigmoid(Vec(x))(0), x, s"sigmoid($x)")
  }

  test("empty inputs: sums, norms and logSumExp agree") {
    val g = Vec.zeros(0)
    val b = BDV.zeros[Double](0)
    assertExact(g.sum, bSum(b), "sum empty")
    assertExact(g.sumExact, 0.0, "sumExact empty")
    assertExact(g.norm1, bNorm(b, 1.0), "norm1 empty")
    assertExact(g.norm2, bNorm(b), "norm2 empty")
    assertExact(g.normInf, bNorm(b, Inf), "normInf empty")
    assertExact(Numerics.logSumExp(g), bLogSumExp(b), "logSumExp empty is -Inf")
    val m03 = Matrix.zeros(0, 3)
    val bm03 = BDM.zeros[Double](0, 3)
    val colSums: BDV[Double] = bSum(bm03(::, Every)).t
    assertVecClose(m03.sum(Axis.Cols), colSums, 0.0, "per-column sum of 0x3")
    assertExact(m03.sum, 0.0, "sum 0x3")
    assertExact(m03.normFrobenius, 0.0, "frobenius 0x3")
  }

  // ---------------------------------------------------------------------------
  // Divergences: Gale's documented semantics, Breeze's observed behaviour pinned
  // ---------------------------------------------------------------------------

  test("divergence: argmax/argmin with NaN return the first NaN; Breeze skips non-leading NaN") {
    // (NaN position, Breeze argmax, Breeze argmin) observed for (1, 2, 3) with one NaN.
    val breezeObserved = Seq((0, 0, 1), (1, 2, 0), (2, 1, 0))
    for (at, breezeArgmax, breezeArgmin) <- breezeObserved do
      val data = Array(1.0, 2.0, 3.0)
      data(at) = NaN
      assertEquals(galeVector(data).argmax, at, s"gale argmax NaN at $at")
      assertEquals(galeVector(data).argmin, at, s"gale argmin NaN at $at")
      assertEquals(bArgmax(breezeVector(data)), breezeArgmax, s"breeze argmax NaN at $at")
      assertEquals(bArgmin(breezeVector(data)), breezeArgmin, s"breeze argmin NaN at $at")
  }

  test("divergence: argmin ties return the first minimum; Breeze returns the last") {
    assertEquals(Vec(1.0, 0.0, 0.0, 2.0).argmin, 1)
    assertEquals(bArgmin(BDV(1.0, 0.0, 0.0, 2.0)), 2)
  }

  test("divergence: matrix argmax ties and NaN use row-major order; Breeze column-major, skipping NaN") {
    val ties = Array(Array(1.0, 5.0), Array(5.0, 2.0))
    assertEquals(galeMatrix(ties).argmax, (0, 1))
    assertEquals(bArgmax(breezeMatrix(ties)), (1, 0))
    for (layout, g) <- matrixViews(ties) do assertEquals(g.argmax, (0, 1), s"ties $layout")
    val withNaN = Array(Array(1.0, 5.0), Array(NaN, 2.0))
    for (layout, g) <- matrixViews(withNaN) do assertEquals(g.argmax, (1, 0), s"NaN $layout")
    assertEquals(bArgmax(breezeMatrix(withNaN)), (0, 1))
  }

  test("divergence: max/min of signed zeros return the first; Breeze max is +0 and min is -0") {
    assertEquals(1.0 / Vec(-0.0, 0.0).max, -Inf)
    assertEquals(1.0 / Vec(0.0, -0.0).min, Inf)
    assertEquals(1.0 / bMax(BDV(-0.0, 0.0)), Inf)
    assertEquals(1.0 / bMin(BDV(0.0, -0.0)), -Inf)
    assertEquals(bArgmin(BDV(0.0, -0.0)), 1)
  }

  test("divergence: empty max/min/argmax/mean throw EmptyInput; Breeze returns -Inf/0.0 or throws IAE") {
    val g = Vec.zeros(0)
    val b = BDV.zeros[Double](0)
    intercept[LinAlgError.EmptyInput](g.max)
    intercept[LinAlgError.EmptyInput](g.min)
    intercept[LinAlgError.EmptyInput](g.argmax)
    intercept[LinAlgError.EmptyInput](g.argmin)
    intercept[LinAlgError.EmptyInput](g.mean)
    intercept[LinAlgError.EmptyInput](Matrix.zeros(0, 0).max)
    intercept[LinAlgError.EmptyInput](Matrix.zeros(0, 3).max(Axis.Cols))
    intercept[LinAlgError.EmptyInput](Matrix.zeros(0, 3).mean(Axis.Cols))
    assertEquals(bMax(b), -Inf)
    assertEquals(bMean(b), 0.0)
    intercept[IllegalArgumentException](bMin(b))
    intercept[IllegalArgumentException](bArgmax(b))
    val bm03 = BDM.zeros[Double](0, 3)
    val colMax: BDV[Double] = bMax(bm03(::, Every)).t
    assertEquals(colMax.toArray.toSeq, Seq(-Inf, -Inf, -Inf))
  }

  test("divergence: mean is sum / n and overflows; Breeze's running mean does not") {
    val big = Vec(Double.MaxValue, Double.MaxValue)
    assertEquals(big.mean, Inf)
    assertEquals(bMean(BDV(Double.MaxValue, Double.MaxValue)), Double.MaxValue)
  }

  test("divergence: norm2 and Frobenius are scaled; Breeze overflows and underflows") {
    val huge = Array(1e300, 1e300, 1e300, 1e300)
    assertRelClose(galeVector(huge).norm2, 2e300, 4 * Eps, "gale norm2 1e300")
    assertEquals(bNorm(breezeVector(huge)), Inf)
    val hugeMatrix = Array(Array(1e300, 1e300), Array(1e300, 1e300))
    for (layout, g) <- matrixViews(hugeMatrix) do
      assertRelClose(g.normFrobenius, 2e300, 4 * Eps, s"gale frobenius 1e300 $layout")
    assertEquals(bNorm(breezeMatrix(hugeMatrix).toDenseVector), Inf)
    val tiny = Array(1e-300, 1e-300)
    assertRelClose(galeVector(tiny).norm2, math.sqrt(2.0) * 1e-300, 4 * Eps, "gale norm2 1e-300")
    assertEquals(bNorm(breezeVector(tiny)), 0.0)
  }

  test("divergence: logSumExp with +Inf is +Inf; Breeze softmax returns -Inf") {
    assertEquals(Numerics.logSumExp(Vec(1.0, Inf)), Inf)
    assertEquals(Numerics.logSumExp(Vec(Inf, Inf)), Inf)
    assertEquals(bLogSumExp(BDV(1.0, Inf)), -Inf)
    assertEquals(bLogSumExp(BDV(Inf, Inf)), -Inf)
    // Gale's softmax of a line with a +Inf maximum is all NaN (documented).
    assert(Numerics.softmax(Vec(1.0, Inf)).toSeq.forall(_.isNaN))
  }

  test("divergence: sigmoid below -709.78 is subnormal; Breeze flushes to 0.0") {
    for x <- Seq(-710.0, -730.0, -740.0) do
      val g = Numerics.sigmoid(Vec(x))(0)
      assert(g > 0.0, s"gale sigmoid($x) = $g")
      assertExact(bSigmoid(x), 0.0, s"breeze sigmoid($x)")
  }
