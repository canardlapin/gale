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
      DoubleArray.fromArray(y),
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

  test("gemv parity: strided x/y and transposed A fall back to pure"):
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
    // Two-or-more-lane runtimes route 128^3 to the measured packed SIMD kernel.
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
      assertEquals(VectorThresholds.nativeGemmMinFlops, 128L * 128L * 128L)
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
