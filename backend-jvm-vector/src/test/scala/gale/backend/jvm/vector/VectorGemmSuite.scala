package gale.backend.jvm.vector

import gale.backend.Backend
import gale.backend.DenseDoubleKernel
import gale.backend.PureBackend
import gale.backend.PureDenseDoubleKernel
import gale.linalg.* // DMat, Matrix, Offset/Rows/Cols/Stride, and the `.value` extensions
import gale.platform.DoubleArray

/** Correctness parity for the Vector-API (SIMD) `gemm`: it must agree with the pure
  * reference within a small tolerance (reassociation from SIMD lane order + FMA makes
  * the result law-equivalent, not bit-identical). Lives in a `gale.*` subpackage so it
  * can read the `private[gale]` `DMat.data` handle and drive the kernels directly.
  */
class VectorGemmSuite extends munit.FunSuite:
  private val AbsTol = 1e-12
  private val RelTol = 1e-10

  // Deterministic, asymmetric fills (no RNG) — small integers so exact products are
  // representable and the SIMD-vs-pure gap is pure reassociation rounding.
  private def mkA(rows: Int, cols: Int): DMat =
    Matrix.tabulate(rows, cols)((i, j) => ((i * 7 + j * 13) % 5 - 2).toDouble)

  private def mkB(rows: Int, cols: Int): DMat =
    Matrix.tabulate(rows, cols)((i, j) => ((i * 11 + j * 5) % 7 - 3).toDouble)

  /** Drive `kernel.gemm` for `alpha·a·b + beta·c` in place on `c`. */
  private def runGemm(kernel: DenseDoubleKernel, a: DMat, b: DMat, alpha: Double, beta: Double, c: DMat): Unit =
    kernel.gemm(
      a.rows,
      b.cols,
      a.cols,
      alpha,
      a.data,
      a.offset.value,
      a.rowStride.value,
      a.colStride.value,
      b.data,
      b.offset.value,
      b.rowStride.value,
      b.colStride.value,
      beta,
      c.data,
      c.offset.value,
      c.rowStride.value,
      c.colStride.value
    )

  private def assertClose(actual: DMat, expected: DMat)(using munit.Location): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var i = 0
    while i < expected.rows do
      var j = 0
      while j < expected.cols do
        val a = actual(i, j)
        val e = expected(i, j)
        val tolerance = AbsTol + RelTol * math.max(math.abs(a), math.abs(e))
        assert(a.isFinite == e.isFinite, s"finiteness mismatch at ($i,$j): vector=$a pure=$e")
        assert(math.abs(a - e) <= tolerance,
          s"mismatch at ($i,$j): vector=$a pure=$e (|Δ|=${math.abs(a - e)}, tol=$tolerance)")
        j += 1
      i += 1

  private def runGemv(
      kernel: DenseDoubleKernel,
      a: DMat,
      x: DVec,
      alpha: Double,
      beta: Double,
      y: Array[Double],
      yOffset: Int = 0,
      yStride: Int = 1
  ): Unit =
    kernel.gemv(
      a.rows,
      a.cols,
      alpha,
      a.data,
      a.offset.value,
      a.rowStride.value,
      a.colStride.value,
      x.data,
      x.offset.value,
      x.stride.value,
      beta,
      DoubleArray.adopt(y), // no copy: the kernel must write the caller's array
      yOffset,
      yStride
    )

  private def assertArrayClose(actual: Array[Double], expected: Array[Double])(using munit.Location): Unit =
    assertEquals(actual.length, expected.length)
    actual.indices.foreach { i =>
      val tolerance = AbsTol + RelTol * math.max(math.abs(actual(i)), math.abs(expected(i)))
      assert(math.abs(actual(i) - expected(i)) <= tolerance,
        s"mismatch at $i: vector=${actual(i)} pure=${expected(i)}")
    }

  /** Bitwise equality. The column-contiguous SIMD gemv applies each element's FMAs in
    * the pure column-major kernel's order, so on the JVM the two agree exactly.
    */
  private def assertArrayBits(actual: Array[Double], expected: Array[Double])(using munit.Location): Unit =
    assertEquals(actual.length, expected.length)
    actual.indices.foreach { i =>
      assertEquals(
        java.lang.Double.doubleToRawLongBits(actual(i)),
        java.lang.Double.doubleToRawLongBits(expected(i)),
        s"bits differ at $i: vector=${actual(i)} pure=${expected(i)}"
      )
    }

  /** Parity for `C := A·B` (alpha=1, beta=0 — exactly what the seam issues). */
  private def checkGemm(rows: Int, cols: Int, shared: Int)(using munit.Location): Unit =
    val a = mkA(rows, shared)
    val b = mkB(shared, cols)
    val cVec = DMat.zeros(rows, cols)
    val cPure = DMat.zeros(rows, cols)
    runGemm(VectorDenseDoubleKernel, a, b, 1.0, 0.0, cVec)
    runGemm(PureDenseDoubleKernel, a, b, 1.0, 0.0, cPure)
    assertClose(cVec, cPure)

  // Sizes chosen to stress lane handling: 1, 3, 5, 7, 11, 17 are non-multiples of the
  // typical 2/4/8 lane count, exercising the scalar tail; 64 exercises whole lane groups
  // and the pure kernel's blocked path.
  test("gemm parity: 1x1x1")(checkGemm(1, 1, 1))
  test("gemm parity: 3x3x3")(checkGemm(3, 3, 3))
  test("gemm parity: 5x7x11 (rows x cols x shared)")(checkGemm(5, 7, 11))
  test("gemm parity: 17x17x17")(checkGemm(17, 17, 17))
  test("gemm parity: 64x64x64")(checkGemm(64, 64, 64))
  test("gemm parity: rectangular 13x20x9")(checkGemm(13, 20, 9))
  test("gemm parity: rectangular 20x13x33")(checkGemm(20, 13, 33))

  test("gemm parity: non-unit alpha/beta over a seeded C (scale + accumulate path)"):
    val rows = 9
    val cols = 11
    val shared = 7
    val a = mkA(rows, shared)
    val b = mkB(shared, cols)
    // Same deterministic starting C for both kernels, so beta-scaling is compared too.
    def seedC(): DMat = Matrix.tabulate(rows, cols)((i, j) => ((i * 3 + j) % 4 - 1).toDouble)
    val cVec = seedC()
    val cPure = seedC()
    runGemm(VectorDenseDoubleKernel, a, b, 2.5, -1.5, cVec)
    runGemm(PureDenseDoubleKernel, a, b, 2.5, -1.5, cPure)
    assertClose(cVec, cPure)

  test("gemv parity: padded rows, scalar tail, and non-unit alpha/beta"):
    val rows = 7
    val cols = 11
    val a = submatrixView(rows, cols, cols + 5)((i, j) => ((i * 17 + j * 7) % 13 - 6).toDouble / 3.0)
    val x = Vec.tabulate(cols)(j => (j - 4).toDouble / 5.0)
    val vectorY = Array.tabulate(rows)(i => i.toDouble - 2.0)
    val pureY = vectorY.clone()
    runGemv(VectorDenseDoubleKernel, a, x, 2.5, -1.25, vectorY)
    runGemv(PureDenseDoubleKernel, a, x, 2.5, -1.25, pureY)
    assertArrayClose(vectorY, pureY)

  test("gemv parity: transposed A with strided x and y falls back to pure"):
    val rows = 9
    val cols = 6
    val a = mkA(cols, rows).t
    val xBacking = Array.fill(cols * 2)(Double.NaN)
    var j = 0
    while j < cols do
      xBacking(j * 2) = (j - 3).toDouble
      j += 1
    val x = new DVec(
      DoubleArray.fromArray(xBacking), Offset.unsafe(0), Length.unsafe(cols), Stride.unsafe(2)
    )
    val vectorY = Array.fill(rows * 3)(7.0)
    val pureY = vectorY.clone()
    runGemv(VectorDenseDoubleKernel, a, x, -0.75, 0.5, vectorY, yOffset = 1, yStride = 3)
    runGemv(PureDenseDoubleKernel, a, x, -0.75, 0.5, pureY, yOffset = 1, yStride = 3)
    assertArrayClose(vectorY, pureY)

  /** Row-major SIMD `gemv` on an interior submatrix of a NaN-padded parent (row and
    * column offsets), with `x` and `y` at non-zero offsets. Column counts are chosen
    * so the two-accumulator pair loop, the single-vector loop and the scalar tail
    * all run on 2-, 4- and 8-lane species.
    */
  private def checkGemvRowMajor(
      rows: Int,
      cols: Int,
      alpha: Double,
      beta: Double,
      rowOff: Int = 3,
      colOff: Int = 5
  )(using munit.Location): Unit =
    val parentCols = colOff + cols + 3
    val lead = 1
    val backing = Array.fill(lead + (rowOff + rows + 1) * parentCols)(Double.NaN)
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        backing(lead + (rowOff + i) * parentCols + colOff + j) = ((i * 13 + j * 7) % 11 - 5).toDouble / 4.0 + i * 1e-3
        j += 1
      i += 1
    val a = new DMat(
      DoubleArray.fromArray(backing),
      Offset.unsafe(lead + rowOff * parentCols + colOff),
      Rows.unsafe(rows),
      Cols.unsafe(cols),
      Stride.unsafe(parentCols),
      Stride.unsafe(1)
    )
    val xBacking = Array.fill(cols + 4)(Double.NaN)
    var j = 0
    while j < cols do
      xBacking(3 + j) = ((j * 5) % 9 - 4).toDouble / 3.0
      j += 1
    val x = new DVec(DoubleArray.fromArray(xBacking), Offset.unsafe(3), Length.unsafe(cols), Stride.unsafe(1))
    val vectorY = Array.tabulate(rows + 4)(k => if k < 3 || k >= rows + 3 then 99.0 else (k % 7).toDouble - 3.0)
    val pureY = vectorY.clone()
    runGemv(VectorDenseDoubleKernel, a, x, alpha, beta, vectorY, yOffset = 3)
    runGemv(PureDenseDoubleKernel, a, x, alpha, beta, pureY, yOffset = 3)
    assertArrayClose(vectorY, pureY)
    assertEquals(Seq(vectorY(0), vectorY(1), vectorY(2), vectorY(rows + 3)), Seq(99.0, 99.0, 99.0, 99.0))

  test("gemv parity: 4x16 offset view (one pair step on 8 lanes)")(checkGemvRowMajor(4, 16, 1.0, 0.0))
  test("gemv parity: 9x75 offset view, alpha/beta")(checkGemvRowMajor(9, 75, -0.5, 2.0))
  test("gemv parity: 33x67 offset view")(checkGemvRowMajor(33, 67, 1.0, 0.0))
  test("gemv parity: 64x64 offset view, beta = 1")(checkGemvRowMajor(64, 64, 1.5, 1.0))
  test("gemv parity: 131x129 offset view, alpha/beta")(checkGemvRowMajor(131, 129, 0.75, -1.25))
  test("gemv parity: 6x7 offset view, columns below one pair")(checkGemvRowMajor(6, 7, 1.0, 0.0, rowOff = 0, colOff = 1))

  /** `gemv` on the transpose view of an offset, padded row-major parent: the
    * column-contiguous SIMD path. Padding is NaN, so a read past a column fails.
    */
  private def checkGemvTransposed(
      rows: Int,
      cols: Int,
      alpha: Double,
      beta: Double,
      xStride: Int = 1
  )(using munit.Location): Unit =
    val parentCols = rows + 3
    val lead = 5
    val backing = Array.fill(lead + cols * parentCols)(Double.NaN)
    var j = 0
    while j < cols do
      var i = 0
      while i < rows do
        backing(lead + j * parentCols + i) = ((i * 13 + j * 7) % 11 - 5).toDouble / 4.0
        i += 1
      j += 1
    val parent = new DMat(
      DoubleArray.fromArray(backing),
      Offset.unsafe(lead),
      Rows.unsafe(cols),
      Cols.unsafe(rows),
      Stride.unsafe(parentCols),
      Stride.unsafe(1)
    )
    val a = parent.t
    assertEquals(a.rowStride.value, 1)
    assertEquals(a.colStride.value, parentCols)
    val xBacking = Array.fill(cols * xStride)(Double.NaN)
    j = 0
    while j < cols do
      xBacking(j * xStride) = ((j * 5) % 9 - 4).toDouble / 3.0
      j += 1
    val x = new DVec(
      DoubleArray.fromArray(xBacking), Offset.unsafe(0), Length.unsafe(cols), Stride.unsafe(xStride)
    )
    // y sits at offset 2 between sentinels that must survive.
    val vectorY = Array.tabulate(rows + 3)(i => if i < 2 || i >= rows + 2 then 99.0 else (i % 7).toDouble - 3.0)
    val pureY = vectorY.clone()
    runGemv(VectorDenseDoubleKernel, a, x, alpha, beta, vectorY, yOffset = 2)
    runGemv(PureDenseDoubleKernel, a, x, alpha, beta, pureY, yOffset = 2)
    assertArrayBits(vectorY, pureY)
    assertEquals(Seq(vectorY(0), vectorY(1), vectorY(rows + 2)), Seq(99.0, 99.0, 99.0))

  // Row counts straddle every lane count (tails of 1..7); column counts straddle
  // the four-column block (tails of 0..3).
  test("gemvT parity: 1x1")(checkGemvTransposed(1, 1, 1.0, 0.0))
  test("gemvT parity: 7x5")(checkGemvTransposed(7, 5, 1.0, 0.0))
  test("gemvT parity: 17x18, alpha/beta")(checkGemvTransposed(17, 18, -0.75, 0.5))
  test("gemvT parity: 33x19, beta = 1")(checkGemvTransposed(33, 19, 2.0, 1.0))
  test("gemvT parity: 64x67, strided x")(checkGemvTransposed(64, 67, 1.0, 0.0, xStride = 3))
  test("gemvT parity: 131x129, strided x, alpha/beta")(checkGemvTransposed(131, 129, 1.5, -2.0, xStride = 2))

  // Eight-column passes followed by a four-column pass and single-column tails.
  test("gemvT parity: 9x12 (one eight-block, one four-block)")(checkGemvTransposed(9, 12, 1.0, 0.0))
  test("gemvT parity: 21x15, alpha/beta (eight + four + three)")(checkGemvTransposed(21, 15, 0.5, -1.0))

  test("gemvT: beta = 0 ignores a non-finite prior y"):
    val n = 9
    val a = mkA(n, n).t
    val x = Vec.tabulate(n)(j => (j - 4).toDouble)
    val vectorY = Array.fill(n)(Double.NaN)
    val pureY = Array.fill(n)(0.0)
    runGemv(VectorDenseDoubleKernel, a, x, 1.0, 0.0, vectorY)
    runGemv(PureDenseDoubleKernel, a, x, 1.0, 0.0, pureY)
    assertArrayBits(vectorY, pureY)

  test("seam: threshold-selected transposed gemv under VectorBackend matches pure"):
    val n = 160
    val a = mkA(n, n + 3)
    val x = Vec.tabulate(n)(i => ((i * 7) % 11 - 5).toDouble)
    val viaVector = a.t.*(x)(using VectorBackend)
    val viaPure = a.t.*(x)(using PureBackend)
    assertEquals(viaVector.length, n + 3)
    assertArrayBits(viaVector.toSeq.toArray, viaPure.toSeq.toArray)

  /** `syrk` on an interior row-major view of a NaN-padded parent into a padded `C`
    * prefilled with `prior`, whose surround must survive: the routed kernel and the
    * forced SIMD path (whatever the shape guard says) vs pure within tolerance,
    * plus an exactly symmetric result. A NaN `prior` checks that the assign-only
    * `syrk` never reads the old `C`.
    */
  private def checkSyrk(m: Int, k: Int, random: Boolean = false, prior: Double = 99.0)(using munit.Location): Unit =
    val rng = new scala.util.Random(31L * m + k)
    val parentCols = k + 3
    val lead = 2
    val backing = Array.fill(lead + m * parentCols)(Double.NaN)
    var r = 0
    while r < m do
      var col = 0
      while col < k do
        backing(lead + r * parentCols + col) =
          if random then rng.nextDouble() * 2.0 - 1.0 else ((r * 7 + col * 13) % 5 - 2).toDouble
        col += 1
      r += 1
    val cRowStride = k + 2
    val cOffset = 1
    def run(syrk: Array[Double] => Unit): Array[Double] =
      val c = Array.fill(cOffset + k * cRowStride + 1)(prior)
      syrk(c)
      c
    val pure = run(c =>
      PureDenseDoubleKernel.syrk(m, k, DoubleArray.fromArray(backing), lead, parentCols, DoubleArray.adopt(c), cOffset, cRowStride)
    )
    val routed = run(c =>
      VectorDenseDoubleKernel.syrk(m, k, DoubleArray.fromArray(backing), lead, parentCols, DoubleArray.adopt(c), cOffset, cRowStride)
    )
    val simd = run(c => VectorDenseDoubleKernel.syrkPackedSimd(m, k, backing, lead, parentCols, c, cOffset, cRowStride))
    for (vector, label) <- Seq(routed -> "routed", simd -> "simd") do
      var i = 0
      while i < k do
        var j = 0
        while j < k do
          val v = vector(cOffset + i * cRowStride + j)
          val e = pure(cOffset + i * cRowStride + j)
          val tolerance = AbsTol + RelTol * math.max(math.abs(v), math.abs(e))
          assert(math.abs(v - e) <= tolerance, s"$label syrk mismatch at ($i,$j): vector=$v pure=$e")
          assertEquals(v, vector(cOffset + j * cRowStride + i), s"$label syrk asymmetric at ($i,$j)")
          j += 1
        i += 1
      // Everything outside the k×k block keeps its sentinel.
      var idx = 0
      while idx < vector.length do
        val inside = idx >= cOffset && (idx - cOffset) / cRowStride < k && (idx - cOffset) % cRowStride < k
        if !inside then
          assertEquals(
            java.lang.Double.doubleToRawLongBits(vector(idx)),
            java.lang.Double.doubleToRawLongBits(prior),
            s"$label syrk wrote outside C at $idx"
          )
        idx += 1

  test("syrk parity: 7x4 (edge-only, below one tile column)")(checkSyrk(7, 4))
  test("syrk parity: 1x1")(checkSyrk(1, 1))
  test("syrk parity: 5x2")(checkSyrk(5, 2))
  test("syrk parity: 40x9 (whole tiles)")(checkSyrk(40, 9))
  test("syrk parity: 33x17 (tiles plus edges)")(checkSyrk(33, 17))
  test("syrk parity: 600x13 random (three shared-dimension slices)")(checkSyrk(600, 13, random = true))
  test("syrk parity: 513x31 random (slice boundary plus one row)")(checkSyrk(513, 31, random = true))
  test("syrk parity: 512x7 random (two whole slices)")(checkSyrk(512, 7, random = true))
  test("syrk parity: 600x1 (single column, edge-only)")(checkSyrk(600, 1, random = true))
  test("syrk parity: 300x2 (k below one tile, two slices)")(checkSyrk(300, 2, random = true))
  test("syrk parity: 257x3 (one tile, slice boundary plus one row)")(checkSyrk(257, 3, random = true))
  test("syrk parity: 300x7 over a NaN-filled C (assign-only)")(checkSyrk(300, 7, random = true, prior = Double.NaN))
  test("syrk parity: 4x40 short-wide")(checkSyrk(4, 40, random = true))
  test("syrk parity: 8x40 short-wide")(checkSyrk(8, 40, random = true))

  test("syrk: m = 0 assigns the zero Gram matrix over a non-finite C"):
    val k = 5
    for simd <- Seq(false, true) do
      val c = Array.fill(k * k)(Double.NaN)
      if simd then VectorDenseDoubleKernel.syrkPackedSimd(0, k, new Array[Double](0), 0, k, c, 0, k)
      else VectorDenseDoubleKernel.syrk(0, k, DoubleArray.alloc(0), 0, k, DoubleArray.adopt(c), 0, k)
      assert(c.forall(_ == 0.0), s"simd=$simd: " + c.mkString(","))

  test("syrk: k = 0 neither throws nor writes"):
    val a = Array.fill(5)(1.0)
    for simd <- Seq(false, true) do
      val c = Array.fill(3)(7.0)
      if simd then VectorDenseDoubleKernel.syrkPackedSimd(5, 0, a, 0, 1, c, 1, 1)
      else VectorDenseDoubleKernel.syrk(5, 0, DoubleArray.fromArray(a), 0, 1, DoubleArray.adopt(c), 1, 1)
      assertEquals(c.toSeq, Seq(7.0, 7.0, 7.0), s"simd=$simd")

  test("syrk shape guard: small and short-wide AᵀA stay pure, tall goes SIMD"):
    if VectorDenseDoubleKernel.preferredLaneCount >= 2 then
      val lanes = VectorDenseDoubleKernel.preferredLaneCount
      assert(!VectorDenseDoubleKernel.syrkUsesSimd(64, 64))
      assert(VectorDenseDoubleKernel.syrkUsesSimd(math.max(128, 4 * lanes), 64))
      assert(!VectorDenseDoubleKernel.syrkUsesSimd(4 * lanes - 1, 4096))
      assert(VectorDenseDoubleKernel.syrkUsesSimd(4 * lanes, 4096))
    else assert(!VectorDenseDoubleKernel.syrkUsesSimd(1024, 256))

  test("seam: threshold-selected AᵀA under VectorBackend matches pure and is symmetric"):
    val a = Matrix.tabulate(300, 130)((i, j) => ((i * 7 + j * 13) % 11 - 5).toDouble / 3.0)
    val viaVector = a.t.*(a)(using VectorBackend)
    val viaPure = a.t.*(a)(using PureBackend)
    assertClose(viaVector, viaPure)
    for i <- 0 until viaVector.rows; j <- 0 until i do
      assertEquals(viaVector(i, j), viaVector(j, i), s"asymmetric at ($i,$j)")

  test("gemm parity: transposed (column-strided) B operand falls back to pure"):
    // A row-major matrix's transpose view has non-unit column stride, so the SIMD gemm
    // must route it to the pure kernel internally — result still correct.
    val rows = 6
    val cols = 8
    val shared = 5
    val a = mkA(rows, shared)
    val bt = mkB(cols, shared).t // (shared x cols) view, colStride = shared != 1
    assert(bt.colStride.value != 1, "test setup: B operand must be column-strided")
    val cVec = DMat.zeros(rows, cols)
    val cPure = DMat.zeros(rows, cols)
    runGemm(VectorDenseDoubleKernel, a, bt, 1.0, 0.0, cVec)
    runGemm(PureDenseDoubleKernel, a, bt, 1.0, 0.0, cPure)
    assertClose(cVec, cPure)

  test("gemm parity: transposed (column-strided) A operand falls back to pure"):
    val rows = 7
    val cols = 6
    val shared = 9
    val at = mkA(shared, rows).t // (rows x shared) view, colStride = shared != 1
    assert(at.colStride.value != 1, "test setup: A operand must be column-strided")
    val b = mkB(shared, cols)
    val cVec = DMat.zeros(rows, cols)
    val cPure = DMat.zeros(rows, cols)
    runGemm(VectorDenseDoubleKernel, at, b, 1.0, 0.0, cVec)
    runGemm(PureDenseDoubleKernel, at, b, 1.0, 0.0, cPure)
    assertClose(cVec, cPure)

  // --- end-to-end through the coarse gemm seam (Matrix.*) --------------------------

  test("seam: threshold-selected product under VectorBackend matches pure"):
    // Two-or-more-lane runtimes route 128^3 (above the 64^3 crossover) to the packed SIMD kernel.
    val n = 128
    val a = mkA(n, n)
    val b = mkB(n, n)
    val viaVector = a.*(b)(using VectorBackend)
    val viaPure = a.*(b)(using PureBackend)
    assertClose(viaVector, viaPure)

  test("seam: below-threshold product stays on the pure path under VectorBackend"):
    val n = 8
    val a = mkA(n, n)
    val b = mkB(n, n)
    val viaVector = a.*(b)(using VectorBackend)
    val viaPure = a.*(b)(using PureBackend)
    assertClose(viaVector, viaPure)

  test("seam: threshold-selected gemv under VectorBackend matches pure"):
    val n = 128
    val a = mkA(n, n)
    val x = Vec.tabulate(n)(i => ((i * 7) % 11 - 5).toDouble)
    val viaVector = a.*(x)(using VectorBackend)
    val viaPure = a.*(x)(using PureBackend)
    assertEquals(viaVector.length, viaPure.length)
    assertArrayClose(viaVector.toSeq.toArray, viaPure.toSeq.toArray)

  test("seam: transposed operand above threshold is correct under VectorBackend"):
    val n = 40
    val a = mkA(n, n)
    val bt = mkB(n, n).t // column-strided -> kernel falls back to pure internally
    val viaVector = a.*(bt)(using VectorBackend)
    val viaPure = a.*(bt)(using PureBackend)
    assertClose(viaVector, viaPure)

  test("given import resolves Backend to jvm-vector"):
    import gale.backend.jvm.vector.given
    assertEquals(summon[Backend].name, "jvm-vector")
    assert(summon[Backend].acceleratesGemm)

  test("threshold requires at least two SIMD lanes"):
    if VectorDenseDoubleKernel.preferredLaneCount >= 2 then
      assertEquals(VectorThresholds.nativeGemmMinFlops, 64L * 64L * 64L)
      assertEquals(VectorThresholds.nativeGemvMinWork, 128L * 128L)
    else
      assertEquals(VectorThresholds.nativeGemmMinFlops, Long.MaxValue)
      assertEquals(VectorThresholds.nativeGemvMinWork, Long.MaxValue)

  test("independent BigDecimal oracle across scale and long inner dimensions"):
    val rng = new scala.util.Random(8675309L)
    for scale <- Seq(1e-100, 1.0, 1e100) do
      val rows = 5
      val shared = 257
      val cols = 7
      val a = Matrix.tabulate(rows, shared)((_, _) => (rng.nextDouble() * 2.0 - 1.0) * scale)
      val b = Matrix.tabulate(shared, cols)((_, _) => (rng.nextDouble() * 2.0 - 1.0) / scale)
      // Drive the SIMD kernel directly: this shape is far below the seam threshold, so
      // going through `a * b` would only ever validate the pure branch against the oracle.
      val actual = DMat.zeros(rows, cols)
      runGemm(VectorDenseDoubleKernel, a, b, 1.0, 0.0, actual)
      var i = 0
      while i < rows do
        var j = 0
        while j < cols do
          var k = 0
          var expected = BigDecimal(0)
          while k < shared do
            expected += BigDecimal(a(i, k)) * BigDecimal(b(k, j))
            k += 1
          val e = expected.toDouble
          val tolerance = 5e-12 * math.max(1.0, math.abs(e))
          assert(actual(i, j).isFinite, s"non-finite result at scale=$scale ($i,$j)")
          assert(math.abs(actual(i, j) - e) <= tolerance,
            s"BigDecimal mismatch at scale=$scale ($i,$j): ${actual(i, j)} vs $e")
          j += 1
        i += 1

  // --- review-hardening: cases the exact-integer fixtures above never exercise --------

  /** A `rows × cols` view (unit column stride) over a wider `rows × parentCols` row-major
    * backing, so `rowStride = parentCols > cols`. Cells beyond `cols` in each row are NaN
    * poison: a kernel that mistook `rowStride == cols` and read into the stride padding
    * would produce NaN and fail parity.
    */
  private def submatrixView(rows: Int, cols: Int, parentCols: Int)(f: (Int, Int) => Double): DMat =
    val backing = new Array[Double](rows * parentCols)
    var i = 0
    while i < rows do
      var j = 0
      while j < parentCols do
        backing(i * parentCols + j) = if j < cols then f(i, j) else Double.NaN
        j += 1
      i += 1
    new DMat(
      DoubleArray.fromArray(backing),
      Offset.unsafe(0),
      Rows.unsafe(rows),
      Cols.unsafe(cols),
      Stride.unsafe(parentCols),
      Stride.unsafe(1)
    )

  test("gemm parity: submatrix views (colStride==1, rowStride>cols) read only real data"):
    val rows = 5
    val cols = 6
    val shared = 7
    val a = submatrixView(rows, shared, shared + 3)((i, j) => ((i * 7 + j * 13) % 5 - 2).toDouble)
    val b = submatrixView(shared, cols, cols + 4)((i, j) => ((i * 11 + j * 5) % 7 - 3).toDouble)
    val cVec = DMat.zeros(rows, cols)
    val cPure = DMat.zeros(rows, cols)
    runGemm(VectorDenseDoubleKernel, a, b, 1.0, 0.0, cVec)
    runGemm(PureDenseDoubleKernel, a, b, 1.0, 0.0, cPure)
    assertClose(cVec, cPure)

  test("gemm parity: random inputs at larger shared bound the reassociation error under Tol"):
    // Exact-integer fixtures never round; irrational inputs at shared=200 make the SIMD
    // lane-order vs pure summation gap real — this asserts it stays within 1e-9.
    val rng = new scala.util.Random(1234567L)
    val rows = 12
    val cols = 15
    val shared = 200
    val a = Matrix.tabulate(rows, shared)((_, _) => rng.nextDouble() * 2.0 - 1.0)
    val b = Matrix.tabulate(shared, cols)((_, _) => rng.nextDouble() * 2.0 - 1.0)
    val cVec = DMat.zeros(rows, cols)
    val cPure = DMat.zeros(rows, cols)
    runGemm(VectorDenseDoubleKernel, a, b, 1.0, 0.0, cVec)
    runGemm(PureDenseDoubleKernel, a, b, 1.0, 0.0, cPure)
    assertClose(cVec, cPure)

  test("gemm parity: column-strided C output falls back to pure"):
    // A transposed view as the destination has colStride != 1, exercising the third leg of
    // the SIMD dispatch guard (cColStride != 1 -> pure fallback).
    val rows = 6
    val cols = 8
    val shared = 5
    val a = mkA(rows, shared)
    val b = mkB(shared, cols)
    val cVec = DMat.zeros(cols, rows).t // rows×cols view, colStride = rows != 1
    val cPure = DMat.zeros(cols, rows).t
    assert(cVec.colStride.value != 1, "test setup: C must be column-strided")
    runGemm(VectorDenseDoubleKernel, a, b, 1.0, 0.0, cVec)
    runGemm(PureDenseDoubleKernel, a, b, 1.0, 0.0, cPure)
    assertClose(cVec, cPure)
