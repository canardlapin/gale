package gale.linalg

import gale.TestAccess

class MutableMatrixSuite extends munit.FunSuite:
  test("DMatBuilder fills row-major storage and transfers it without a public mutable alias") {
    val builder = DMat.newBuilder(2, 3)
    var i = 0
    while i < builder.size do
      builder.writeLinear(i, i.toDouble + 1.0)
      i += 1
    builder(1, 2) = 9.0

    val matrix = builder.result()

    assertEquals(matrix.valuesRowMajor, Seq(1.0, 2.0, 3.0, 4.0, 5.0, 9.0))
    intercept[LinAlgError.UnsupportedOperation](builder.update(0, 0, -1.0))
    intercept[LinAlgError.UnsupportedOperation](builder.result())
    assertEquals(matrix(0, 0), 1.0)
  }

  test("DMatBuilder supports empty shapes and validates logical indices") {
    val empty = Matrix.newBuilder(0, 3).result()
    assertEquals((empty.rows, empty.cols), (0, 3))

    val builder = DMatBuilder.zeros(2, 2)
    intercept[LinAlgError.IndexOutOfBounds](builder.update(-1, 0, 1.0))
    intercept[LinAlgError.IndexOutOfBounds](builder.update(0, 2, 1.0))
    intercept[LinAlgError.IndexOutOfBounds](builder.writeLinear(4, 1.0))
  }

  test("consumeQR transfers builder storage, matches result then QR, and closes every operation") {
    val values = Seq(
      1.0, 0.0, 2.0,
      1.0, 1.0, -1.0,
      1.0, 2.0, 0.5,
      1.0, 3.0, 1.5,
      1.0, 4.0, -0.5
    )
    def filledBuilder(): DMatBuilder =
      val builder = DMatBuilder.zeros(5, 3)
      var i = 0
      while i < values.length do
        builder.writeLinear(i, values(i))
        i += 1
      builder

    val options = QROptions(pivoting = QRPivoting.Column)
    val expected = filledBuilder().result().qr(options)
    val builder = filledBuilder()
    val transferred = TestAccess.dmatBuilderStorage(builder)
    val actual = builder.consumeQR(options, DenseWorkspace.empty)

    assert(TestAccess.sameStorage(transferred, TestAccess.dmatStorage(actual.r)))
    assertEquals(actual.r.valuesRowMajor, expected.r.valuesRowMajor)
    assertEquals(actual.columnPermutation.toIndexSeq, expected.columnPermutation.toIndexSeq)
    assertEquals(actual.diagnostics, expected.diagnostics)
    assertEquals(actual.q.valuesRowMajor, expected.q.valuesRowMajor)
    val rhs = Vec(1.0, 2.0, -1.0, 0.5, 3.0)
    assertEquals(actual.solveLeastSquares(rhs).orThrow.toSeq, expected.solveLeastSquares(rhs).orThrow.toSeq)
    assertEquals(actual.normalizedCovariance.orThrow.valuesRowMajor, expected.normalizedCovariance.orThrow.valuesRowMajor)

    intercept[LinAlgError.UnsupportedOperation](builder(0, 0))
    intercept[LinAlgError.UnsupportedOperation](builder.update(0, 0, -1.0))
    intercept[LinAlgError.UnsupportedOperation](builder.writeLinear(0, -1.0))
    intercept[LinAlgError.UnsupportedOperation](builder.fill(0.0))
    intercept[LinAlgError.UnsupportedOperation](builder.result())
    intercept[LinAlgError.UnsupportedOperation](builder.consumeQR)
    assertEquals(actual.r.valuesRowMajor, expected.r.valuesRowMajor)
  }

  test("every consumeQR overload transfers ownership and preserves QR semantics") {
    val rows = 11
    val cols = 5
    val values = Seq.tabulate(rows * cols): index =>
      val row = index / cols
      val col = index % cols
      math.sin((row + 1.0) * (col + 2.0) * 0.0625) + (if row == col then 2.0 else 0.0)
    val matrix = Matrix.dense(rows, cols, values)
    val pivoted = QROptions(pivoting = QRPivoting.Column)

    def builder(): DMatBuilder =
      val out = DMatBuilder.zeros(rows, cols)
      var index = 0
      while index < values.length do
        out.writeLinear(index, values(index))
        index += 1
      out

    val cases = Seq[(DMatBuilder, DMatBuilder => QR, QROptions)](
      (builder(), _.consumeQR, QROptions.Default),
      (builder(), _.consumeQR(DenseWorkspace.empty), QROptions.Default),
      (builder(), _.consumeQR(pivoted), pivoted),
      (builder(), _.consumeQR(pivoted, DenseWorkspace.empty), pivoted),
    )

    for (source, consume, options) <- cases do
      val transferred = TestAccess.dmatBuilderStorage(source)
      val expected = matrix.qr(options)
      val actual = consume(source)

      assert(TestAccess.sameStorage(transferred, TestAccess.dmatStorage(actual.r)))
      assertEquals(actual.r.valuesRowMajor, expected.r.valuesRowMajor)
      assertEquals(actual.columnPermutation.toIndexSeq, expected.columnPermutation.toIndexSeq)
      assertEquals(actual.diagnostics, expected.diagnostics)
      assertEquals(actual.q.valuesRowMajor, expected.q.valuesRowMajor)
      intercept[LinAlgError.UnsupportedOperation](source.result())
      intercept[LinAlgError.UnsupportedOperation](source.consumeQR)
      assertEquals(actual.r.valuesRowMajor, expected.r.valuesRowMajor)
  }

  test("consumeQR handles empty dimensions while closing the builder") {
    for (rows, cols) <- Seq((0, 0), (0, 3), (3, 0)) do
      val builder = DMatBuilder.zeros(rows, cols)
      val transferred = TestAccess.dmatBuilderStorage(builder)
      val actual = builder.consumeQR(QROptions(pivoting = QRPivoting.Column), DenseWorkspace.empty)

      assertEquals(actual.r.shape, Shape(Rows(rows), Cols(cols)))
      assertEquals(actual.diagnostics.rank, Some(0))
      assertEquals(actual.columnPermutation.toIndexSeq, (0 until cols).toIndexedSeq)
      assert(TestAccess.sameStorage(transferred, TestAccess.dmatStorage(actual.r)))
      intercept[LinAlgError.UnsupportedOperation](builder.result())
      intercept[LinAlgError.UnsupportedOperation](builder.consumeQR)
  }
