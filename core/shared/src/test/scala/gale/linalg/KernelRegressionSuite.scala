package gale.linalg

import gale.TestAccess
import gale.platform.DoubleArray
import gale.sparse.Sparse

class KernelRegressionSuite extends munit.FunSuite:
  // Item 1: with beta == 0 the destination must be assigned, not read-and-scaled.
  // A NaN already sitting in the output buffer would poison the result via
  // 0.0 * NaN == NaN if the kernel reads y before overwriting it.
  test("mulInto with beta==0 ignores NaN already in the destination (all layouts)") {
    // A == [[1, 2], [3, 4]] realised three different storage ways; A * [1, 1] == [3, 7].
    val expected = Seq(3.0, 7.0)
    val x = Vec(1.0, 1.0)

    // (a) contiguous row-major -> dgemvRowMajor
    val rowMajor = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)

    // (b) column-major view -> dgemvColMajor (A.t of a row-major matrix)
    val colMajor = Matrix.dense(2, 2)(1.0, 3.0, 2.0, 4.0).t

    // (c) strided submatrix view -> dgemv (rowStride != 1 and colStride != 1)
    val backing = TestAccess.doubleArray(1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0)
    val strided = TestAccess.mat(backing, offset = 0, rows = 2, cols = 2, rowStride = 4, colStride = 2)

    for (label, a) <- Seq("row-major" -> rowMajor, "col-major" -> colMajor, "strided" -> strided) do
      val y = MutableDVec.zeros(2)
      y(0) = Double.NaN
      y(1) = Double.NaN
      a.mulInto(x, y)
      assert(y(0).isFinite && y(1).isFinite, s"$label produced non-finite output: ${y.asVec.toSeq}")
      assertEquals(y.asVec.toSeq, expected, s"$label mismatch")
  }

  // Item 1 continued: dgemm (matrix-matrix) has the same beta==0 read-not-assign hazard.
  test("dgemm with beta==0 ignores NaN already in the destination") {
    val a = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
    val b = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val c = DMat.zeros(2, 2)
    TestAccess.poisonWithNaN(c)
    TestAccess.gemm(a, b, c)
    assertEquals(c.valuesRowMajor, Seq(1.0, 2.0, 3.0, 4.0))
  }

  // Item 3: writing into a destination that shares storage with the source
  // vector corrupts the computation mid-flight. The kernels must reject it
  // rather than silently produce garbage.
  test("DMat.mulInto rejects a destination aliasing the source vector") {
    val a = Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0,
      7.0, 8.0, 9.0
    )
    val y = MutableDVec.zeros(3)
    y(0) = 1.0
    y(1) = 2.0
    y(2) = 3.0
    val aliased = y.asVec // shares y's backing storage
    intercept[LinAlgError.UnsupportedOperation] {
      a.mulInto(aliased, y)
    }
  }

  test("CSR.mulInto and tMulInto reject a destination aliasing the source vector") {
    val a =
      Sparse
        .coo(3, 3)
        .add(0, 0, 1.0)
        .add(1, 1, 2.0)
        .add(2, 0, 3.0)
        .add(2, 2, 4.0)
        .toCSR()
    val y = MutableDVec.zeros(3)
    y(0) = 1.0
    y(1) = 2.0
    y(2) = 3.0
    intercept[LinAlgError.UnsupportedOperation] {
      a.mulInto(y.asVec, y)
    }
    val z = MutableDVec.zeros(3)
    z(0) = 1.0
    z(1) = 2.0
    z(2) = 3.0
    intercept[LinAlgError.UnsupportedOperation] {
      a.tMulInto(z.asVec, z)
    }
  }

  // The row-major kernel tiles four rows and the column-major kernel sweeps four
  // columns; every shape here straddles both tiles and their scalar tails. Integer
  // data keeps every partial sum exact, so any summation order must match exactly.
  test("dgemvRowMajor/dgemvColMajor tiles and tails match a naive product exactly") {
    import gale.kernel.DoubleKernels
    def entry(i: Int, j: Int): Double = ((i * 7 + j * 3) % 9 - 4).toDouble
    for
      rows <- 1 to 11
      cols <- 1 to 11
      (alpha, beta) <- Seq((1.0, 0.0), (2.0, -1.0), (-0.5, 1.0))
      yStride <- Seq(1, 2)
    do
      val xs = Array.tabulate(cols)(j => ((j * 5) % 7 - 3).toDouble)
      val prior = Array.tabulate(rows)(i => (i % 5 - 2).toDouble)
      val expected = Array.tabulate(rows) { i =>
        var s = 0.0
        var j = 0
        while j < cols do
          s += entry(i, j) * xs(j)
          j += 1
        if beta == 0.0 then alpha * s else alpha * s + beta * prior(i)
      }
      def run(label: String, kernel: (DoubleArray, DoubleArray) => Unit, x: DoubleArray): Unit =
        val y = TestAccess.filled(1 + rows * yStride, Double.NaN)
        var i = 0
        while i < rows do
          if beta != 0.0 then y(1 + i * yStride) = prior(i)
          i += 1
        kernel(x, y)
        i = 0
        while i < rows do
          assertEquals(y(1 + i * yStride), expected(i), s"$label ${rows}x$cols alpha=$alpha beta=$beta yStride=$yStride row $i")
          i += 1
        // Untouched gap cells keep their sentinel.
        if yStride > 1 then assert(y(2).isNaN, s"$label wrote a stride gap")

      // Row-major parent padded by 2 columns, offset 3.
      val rowMajor = TestAccess.filled(3 + rows * (cols + 2), Double.NaN)
      // Column-major parent padded by 2 rows, offset 3.
      val colMajor = TestAccess.filled(3 + cols * (rows + 2), Double.NaN)
      for i <- 0 until rows; j <- 0 until cols do
        rowMajor(3 + i * (cols + 2) + j) = entry(i, j)
        colMajor(3 + j * (rows + 2) + i) = entry(i, j)
      val xUnit = TestAccess.doubleArray(xs*)
      val xStrided = TestAccess.filled(1 + cols * 3, Double.NaN)
      for j <- 0 until cols do xStrided(1 + j * 3) = xs(j)

      run(
        "row-major",
        (x, y) => DoubleKernels.dgemvRowMajor(rows, cols, alpha, rowMajor, 3, cols + 2, x, 0, beta, y, 1, yStride),
        xUnit
      )
      run(
        "col-major",
        (x, y) => DoubleKernels.dgemvColMajor(rows, cols, alpha, colMajor, 3, rows + 2, x, 1, 3, beta, y, 1, yStride),
        xStrided
      )
  }
