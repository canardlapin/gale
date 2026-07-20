package gale.linalg

import gale.TestAccess

class DenseTransformSuite extends munit.FunSuite:

  private def assertMatrixEquals(actual: DMat, expected: DMat): Unit =
    assertEquals((actual.rows, actual.cols), (expected.rows, expected.cols))
    assertEquals(actual.valuesRowMajor, expected.valuesRowMajor)

  test("DMatBuilder.from makes one owned copy and result adopts it") {
    val source = DMat.tabulate(2, 3)((row, col) => row * 10.0 + col)
    val builder = DMatBuilder.from(source)
    val builderStorage = TestAccess.dmatBuilderStorage(builder)

    assertEquals(builder(1, 2), 12.0)
    builder(0, 0) = 99.0
    val result = builder.result()

    assertEquals(source(0, 0), 0.0)
    assertEquals(result(0, 0), 99.0)
    assert(!TestAccess.sameStorage(TestAccess.dmatStorage(source), TestAccess.dmatStorage(result)))
    assert(TestAccess.sameStorage(builderStorage, TestAccess.dmatStorage(result)))
    intercept[LinAlgError.UnsupportedOperation](builder.update(0, 0, -1.0))
  }

  test("toBuilder copies transposed and strided logical order exactly once") {
    val source = DMat.tabulate(4, 5)((row, col) => row * 10.0 + col)
    val view = source.t.slice(1, 5, 1, 4)
    val builder = view.toBuilder
    val builderStorage = TestAccess.dmatBuilderStorage(builder)
    val result = builder.result()
    val expected = DMat.tabulate(4, 3)((row, col) => source(col + 1, row + 1))

    assertMatrixEquals(result, expected)
    assert(result.isContiguousRowMajor)
    assert(!TestAccess.sameStorage(TestAccess.dmatStorage(source), TestAccess.dmatStorage(result)))
    assert(TestAccess.sameStorage(builderStorage, TestAccess.dmatStorage(result)))
  }

  test("matrix slice is a strided view and validates half-open ranges") {
    val source = DMat.tabulate(4, 5)((row, col) => row * 10.0 + col)
    val view = source.slice(1, 4, 2, 5)

    assertEquals((view.rows, view.cols), (3, 3))
    assertEquals(view.valuesRowMajor, Seq(12.0, 13.0, 14.0, 22.0, 23.0, 24.0, 32.0, 33.0, 34.0))
    assert(TestAccess.sameStorage(TestAccess.dmatStorage(source), TestAccess.dmatStorage(view)))
    assertEquals(view.rowStride, source.rowStride)
    assertEquals(view.colStride, source.colStride)

    intercept[LinAlgError.InvalidArgument](source.slice(-1, 2, 0, 2))
    intercept[LinAlgError.InvalidArgument](source.slice(2, 1, 0, 2))
    intercept[LinAlgError.InvalidArgument](source.slice(0, 5, 0, 2))
    intercept[LinAlgError.InvalidArgument](source.slice(0, 2, 4, 3))
    intercept[LinAlgError.InvalidArgument](source.slice(0, 2, 0, 6))

    val noRows = source.slice(2, 2, 1, 4)
    assertEquals((noRows.rows, noRows.cols), (0, 3))
    val noCols = DMat.zeros(3, 0).slice(1, 3, 0, 0)
    assertEquals((noCols.rows, noCols.cols), (2, 0))
  }

  test("row and column gathers own storage and preserve repeats on strided input") {
    val source = DMat.tabulate(4, 5)((row, col) => row * 10.0 + col)
    val view = source.t.slice(1, 5, 1, 4)

    val rows = view.gatherRows(IndexedSeq(3, 1, 3))
    assertEquals(rows.valuesRowMajor, Seq(14.0, 24.0, 34.0, 12.0, 22.0, 32.0, 14.0, 24.0, 34.0))
    assert(rows.isContiguousRowMajor)
    assert(!TestAccess.sameStorage(TestAccess.dmatStorage(view), TestAccess.dmatStorage(rows)))

    val columns = view.gatherColumns(IndexedSeq(2, 0, 2))
    assertEquals(columns.valuesRowMajor, Seq(31.0, 11.0, 31.0, 32.0, 12.0, 32.0, 33.0, 13.0, 33.0, 34.0, 14.0, 34.0))
    assert(columns.isContiguousRowMajor)
    assert(!TestAccess.sameStorage(TestAccess.dmatStorage(view), TestAccess.dmatStorage(columns)))

    intercept[LinAlgError.IndexOutOfBounds](view.gatherRows(IndexedSeq(0, 4)))
    intercept[LinAlgError.IndexOutOfBounds](view.gatherRows(IndexedSeq(-1)))
    intercept[LinAlgError.IndexOutOfBounds](view.gatherColumns(IndexedSeq(3)))

    assertEquals(view.gatherRows(IndexedSeq.empty).shape, Shape(Rows(0), Cols(view.cols)))
    assertEquals(view.gatherColumns(IndexedSeq.empty).shape, Shape(Rows(view.rows), Cols(0)))
    assertEquals(DMat.zeros(0, 3).gatherRows(IndexedSeq.empty).shape, Shape(Rows(0), Cols(3)))
  }

  test("DVec gather and builder copy work for strided vectors") {
    val source = TestAccess.stridedCopy(Vec(2.0, 4.0, 6.0, 8.0), stride = 3)
    val gathered = source.gather(IndexedSeq(3, 1, 3, 0))
    assertEquals(gathered.toSeq, Seq(8.0, 4.0, 8.0, 2.0))
    assert(!TestAccess.sameStorage(TestAccess.dvecStorage(source), TestAccess.dvecStorage(gathered)))

    val builder = source.toBuilder
    val builderStorage = TestAccess.dvecBuilderStorage(builder)
    builder(1) = 40.0
    val copied = builder.result()
    assertEquals(copied.toSeq, Seq(2.0, 40.0, 6.0, 8.0))
    assertEquals(source.toSeq, Seq(2.0, 4.0, 6.0, 8.0))
    assert(TestAccess.sameStorage(builderStorage, TestAccess.dvecStorage(copied)))
    assert(!TestAccess.sameStorage(TestAccess.dvecStorage(source), TestAccess.dvecStorage(copied)))

    intercept[LinAlgError.IndexOutOfBounds](source.gather(IndexedSeq(4)))
    intercept[LinAlgError.IndexOutOfBounds](source.gather(IndexedSeq(-1)))
    assertEquals(source.gather(IndexedSeq.empty).length, 0)
    assertEquals(DVec.zeros(0).gather(IndexedSeq.empty).length, 0)
  }

  test("diagonal addition and average symmetrization return owned matrices") {
    val rectangular = DMat.dense(2, 3, Seq(1.0, 2.0, 3.0, 4.0, 5.0, 6.0))
    val regularized = rectangular.t.addToDiagonal(10.0)
    assertEquals(regularized.valuesRowMajor, Seq(11.0, 4.0, 2.0, 15.0, 3.0, 6.0))
    assertEquals(rectangular.valuesRowMajor, Seq(1.0, 2.0, 3.0, 4.0, 5.0, 6.0))
    assert(!TestAccess.sameStorage(TestAccess.dmatStorage(rectangular), TestAccess.dmatStorage(regularized)))

    val asymmetric = DMat.dense(3, 3, Seq(1.0, 2.0, 8.0, 4.0, 5.0, 10.0, 6.0, 12.0, 9.0))
    val symmetric = asymmetric.symmetrizedAverage
    assertEquals(symmetric.valuesRowMajor, Seq(1.0, 3.0, 7.0, 3.0, 5.0, 11.0, 7.0, 11.0, 9.0))
    assertEquals(symmetric.valuesRowMajor, symmetric.t.valuesRowMajor)
    assert(!TestAccess.sameStorage(TestAccess.dmatStorage(asymmetric), TestAccess.dmatStorage(symmetric)))
    intercept[LinAlgError.NonSquareMatrix](rectangular.symmetrizedAverage)
    assertEquals(DMat.zeros(0, 0).symmetrizedAverage.shape, Shape(Rows(0), Cols(0)))
  }

  test("row-major traversal and primitive copy honor transposed strided order") {
    val source = DMat.tabulate(4, 5)((row, col) => row * 10.0 + col)
    val view = source.t.slice(1, 5, 1, 4)
    val expected = Array(11.0, 21.0, 31.0, 12.0, 22.0, 32.0, 13.0, 23.0, 33.0, 14.0, 24.0, 34.0)

    val visited = new Array[Double](expected.length)
    var count = 0
    view.foreachRowMajor: value =>
      visited(count) = value
      count += 1
    assertEquals(count, expected.length)
    assertEquals(visited.toSeq, expected.toSeq)

    val destination = Array.fill(expected.length + 2)(-1.0)
    view.copyRowMajorTo(destination, 1)
    assertEquals(destination.head, -1.0)
    assertEquals(destination.last, -1.0)
    assertEquals(destination.slice(1, destination.length - 1).toSeq, expected.toSeq)
    intercept[LinAlgError.InvalidArgument](view.copyRowMajorTo(new Array[Double](expected.length - 1)))
    intercept[LinAlgError.InvalidArgument](view.copyRowMajorTo(new Array[Double](expected.length), -1))
    DMat.zeros(0, 3).copyRowMajorTo(Array.emptyDoubleArray)
  }

  test("vector traversal and primitive copy honor stride and destination bounds") {
    val vector = TestAccess.stridedCopy(Vec(2.0, 4.0, 6.0, 8.0), stride = 2)
    val visited = new Array[Double](vector.length)
    var count = 0
    vector.foreachValue: value =>
      visited(count) = value
      count += 1
    assertEquals(visited.toSeq, Seq(2.0, 4.0, 6.0, 8.0))

    val destination = Array.fill(6)(-1.0)
    vector.copyTo(destination, 1)
    assertEquals(destination.toSeq, Seq(-1.0, 2.0, 4.0, 6.0, 8.0, -1.0))
    intercept[LinAlgError.InvalidArgument](vector.copyTo(new Array[Double](3)))
    intercept[LinAlgError.InvalidArgument](vector.copyTo(new Array[Double](4), -1))
    DVec.zeros(0).copyTo(Array.emptyDoubleArray)
  }
