package gale.linalg

class MutableMatrixSuite extends munit.FunSuite:
  test("DMatBuilder fills row-major storage and transfers it without a public mutable alias") {
    val builder = DMat.newBuilder(2, 3)
    var i = 0
    while i < builder.size do
      builder.updateRowMajor(i, i.toDouble + 1.0)
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
    intercept[LinAlgError.IndexOutOfBounds](builder.updateRowMajor(4, 1.0))
  }
