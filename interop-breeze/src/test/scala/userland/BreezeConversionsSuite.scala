package userland

import breeze.linalg.CSCMatrix
import breeze.linalg.DenseMatrix
import breeze.linalg.DenseVector
import gale.interop.breeze.*
import gale.linalg.*
import gale.sparse.*

/** Round-trip fidelity and the copy/view contract for the Breeze interop bridge.
  *
  * Lives in a non-`gale.interop` package so it exercises the module exactly as a
  * downstream user would (`import gale.interop.breeze.*`) and so the library name
  * `breeze` is not shadowed by the module's own package.
  */
class BreezeConversionsSuite extends munit.FunSuite:

  private def bits(x: Double): Long = java.lang.Double.doubleToRawLongBits(x)

  /** Bit-exact scalar equality (conversions copy raw doubles, never recompute). */
  private def assertBits(a: Double, b: Double, clue: => String): Unit =
    assert(bits(a) == bits(b), s"$clue: $a != $b (bits ${bits(a)} vs ${bits(b)})")

  private def elem(i: Int, j: Int): Double =
    (i * 7 + j * 3) * 0.1 + i * 0.001 - j * 0.017 + 0.5

  private def galeDense(m: Int, n: Int): DMat =
    Matrix.tabulate(m, n)((i, j) => elem(i, j))

  private def breezeDense(m: Int, n: Int): DenseMatrix[Double] =
    DenseMatrix.tabulate(m, n)((i, j) => elem(i, j))

  // ---------------------------------------------------------------------------
  // Dense matrix round trips
  // ---------------------------------------------------------------------------

  test("DMat -> Breeze copy -> DMat: shape + bit-exact, contiguous") {
    val a  = galeDense(4, 6)
    val b  = toBreezeCopy(a)
    assertEquals(b.rows, 4)
    assertEquals(b.cols, 6)
    for i <- 0 until 4; j <- 0 until 6 do assertBits(b(i, j), a(i, j), s"->breeze ($i,$j)")
    val a2 = fromBreezeCopy(b)
    assertEquals(a2.rows, 4)
    assertEquals(a2.cols, 6)
    for i <- 0 until 4; j <- 0 until 6 do assertBits(a2(i, j), a(i, j), s"round-trip ($i,$j)")
  }

  test("DMat -> Breeze copy: non-contiguous (transposed) gale input") {
    val a  = galeDense(5, 3)
    val at = a.t // transposed view: 3x5, column-contiguous, not row-major
    val b  = toBreezeCopy(at)
    assertEquals((b.rows, b.cols), (3, 5))
    for i <- 0 until 3; j <- 0 until 5 do assertBits(b(i, j), a(j, i), s"transposed ($i,$j)")
  }

  test("Breeze -> DMat copy: transposed and sliced Breeze inputs") {
    val b = breezeDense(5, 4)

    val gt = fromBreezeCopy(b.t) // 4x5 transposed view
    assertEquals((gt.rows, gt.cols), (4, 5))
    for i <- 0 until 4; j <- 0 until 5 do assertBits(gt(i, j), b(j, i), s"t ($i,$j)")

    val slice = b(1 until 4, 1 until 3) // 3x2 view with offset + majorStride
    val gs    = fromBreezeCopy(slice)
    assertEquals((gs.rows, gs.cols), (3, 2))
    for i <- 0 until 3; j <- 0 until 2 do assertBits(gs(i, j), b(1 + i, 1 + j), s"slice ($i,$j)")
  }

  // ---------------------------------------------------------------------------
  // Dense matrix zero-copy view (Breeze -> gale) + aliasing contract
  // ---------------------------------------------------------------------------

  test("unsafeFromBreezeView: matches elementwise for column-major, transpose, and slice") {
    val b = breezeDense(4, 5)

    val v1 = unsafeFromBreezeView(b)
    for i <- 0 until 4; j <- 0 until 5 do assertBits(v1(i, j), b(i, j), s"cm ($i,$j)")

    val v2 = unsafeFromBreezeView(b.t)
    assertEquals((v2.rows, v2.cols), (5, 4))
    for i <- 0 until 5; j <- 0 until 4 do assertBits(v2(i, j), b(j, i), s"t ($i,$j)")

    val slice = b(1 until 4, 2 until 5)
    val v3    = unsafeFromBreezeView(slice)
    assertEquals((v3.rows, v3.cols), (3, 3))
    for i <- 0 until 3; j <- 0 until 3 do assertBits(v3(i, j), b(1 + i, 2 + j), s"slice ($i,$j)")
  }

  test("copy/view aliasing contract (dense matrix)") {
    val b    = breezeDense(3, 3)
    val view = unsafeFromBreezeView(b)
    val copy = fromBreezeCopy(b)
    val before = b(1, 2)

    b(1, 2) = before + 123.0 // mutate the Breeze storage
    assertBits(view(1, 2), before + 123.0, "view aliases the Breeze write")
    assertBits(copy(1, 2), before, "copy is independent of the Breeze write")
  }

  test("gale -> Breeze is copy-only: mutating the Breeze copy does not touch gale") {
    val a = galeDense(3, 3)
    val b = toBreezeCopy(a)
    val original = a(0, 0)
    b(0, 0) = original + 999.0
    assertBits(a(0, 0), original, "gale is unaffected by mutating its Breeze copy")
  }

  test("dense copies preserve signed zero, infinities, and a NaN payload bit-for-bit") {
    val payloadNaN = java.lang.Double.longBitsToDouble(0x7ff8000000000042L)
    val b = new DenseMatrix[Double](
      2,
      3,
      Array(0.0, -0.0, Double.PositiveInfinity, Double.NegativeInfinity, payloadNaN, 1.25)
    )
    val g = fromBreezeCopy(b)
    val roundTrip = toBreezeCopy(g)
    for i <- 0 until b.rows; j <- 0 until b.cols do
      assertBits(roundTrip(i, j), b(i, j), s"special-value round trip ($i,$j)")
  }

  // ---------------------------------------------------------------------------
  // Dense vector round trips + view
  // ---------------------------------------------------------------------------

  test("DVec <-> Breeze DenseVector: round trip incl. strided gale input") {
    val v = Vec.tabulate(7)(i => i * 0.3 - 1.1)
    val b = toBreezeCopy(v)
    assertEquals(b.length, 7)
    for i <- 0 until 7 do assertBits(b(i), v(i), s"->breeze [$i]")
    val v2 = fromBreezeCopy(b)
    for i <- 0 until 7 do assertBits(v2(i), v(i), s"round-trip [$i]")

    // A gale column of a row-major matrix is a strided (stride>1) view.
    val col = galeDense(4, 3).col(1)
    val bc  = toBreezeCopy(col)
    for i <- 0 until 4 do assertBits(bc(i), col(i), s"strided col [$i]")
  }

  test("unsafeFromBreezeView (vector) aliases; fromBreezeCopy does not") {
    val b    = DenseVector.tabulate(5)(i => i * 1.25)
    val view = unsafeFromBreezeView(b)
    val copy = fromBreezeCopy(b)
    val before = b(3)
    b(3) = before + 42.0
    assertBits(view(3), before + 42.0, "vector view aliases")
    assertBits(copy(3), before, "vector copy independent")
  }

  test("negative-stride Breeze vectors are copyable but cannot become Gale views") {
    val reversed = new DenseVector[Double](Array(1.0, 2.0, 3.0, 4.0, 5.0), 4, -1, 5)
    val copy = fromBreezeCopy(reversed)
    assertEquals(copy.toSeq, Seq(5.0, 4.0, 3.0, 2.0, 1.0))
    intercept[IllegalArgumentException]:
      unsafeFromBreezeView(reversed)
  }

  // ---------------------------------------------------------------------------
  // Sparse round trips (copy-only)
  // ---------------------------------------------------------------------------

  private val sparseEntries = Seq(
    (0, 0, 1.5), (0, 3, -2.25), (1, 1, 3.0), (2, 0, 0.75), (2, 2, -4.5), (3, 3, 6.125)
  )

  private def galeCsr(rows: Int, cols: Int): CSR =
    val builder = Sparse.coo(rows, cols)
    sparseEntries.foreach { case (r, c, v) => builder.add(r, c, v) }
    builder.toCSR()

  test("CSR -> Breeze CSCMatrix -> CSR/CSC: structure + bit-exact values") {
    val rows = 4
    val cols = 4
    val csr  = galeCsr(rows, cols)
    val bcsc = toBreezeCopy(csr)
    assertEquals((bcsc.rows, bcsc.cols), (rows, cols))
    for i <- 0 until rows; j <- 0 until cols do
      assertBits(bcsc(i, j), csr(i, j), s"csr->breeze ($i,$j)")

    val backCsc = fromBreezeToCsc(bcsc)
    val backCsr = fromBreezeToCsr(bcsc)
    for i <- 0 until rows; j <- 0 until cols do
      assertBits(backCsc(i, j), csr(i, j), s"breeze->csc ($i,$j)")
      assertBits(backCsr(i, j), csr(i, j), s"breeze->csr ($i,$j)")
  }

  test("CSC -> Breeze CSCMatrix preserves the same matrix") {
    val csc  = galeCsr(4, 4).toCSC
    val bcsc = toBreezeCopy(csc)
    for i <- 0 until 4; j <- 0 until 4 do assertBits(bcsc(i, j), csc(i, j), s"csc->breeze ($i,$j)")
  }

  test("Breeze CSCMatrix -> gale round trip") {
    val builder = new CSCMatrix.Builder[Double](3, 3)
    builder.add(0, 0, 2.5)
    builder.add(2, 1, -1.75)
    builder.add(1, 2, 3.25)
    val bcsc = builder.result()

    val csc = fromBreezeToCsc(bcsc)
    for i <- 0 until 3; j <- 0 until 3 do assertBits(csc(i, j), bcsc(i, j), s"round ($i,$j)")
  }
