package gale.interop.ravel

import gale.linalg.*
import munit.FunSuite
import ravel.*
import ravel.DType.given

final class RavelConversionsSuite extends FunSuite:
  test("composed offset and negative-stride time blocks produce owned row-major copies") {
    val input = NDArray.tabulate[Double](3, 6)((row, time) => row * 100.0 + time)
    val view = input.narrow(0, 1, 2).slice(1, Slice(1, 6, 2)).reverse(1).transpose
    val matrix = fromRavelCopy(view)
    assertEquals((matrix.rows, matrix.cols), (3, 2))
    assertEquals(matrix.valuesRowMajor, Seq(105.0, 205.0, 103.0, 203.0, 101.0, 201.0))
    assertEquals(matrix.rowStride.value, 2)
    assertEquals(matrix.colStride.value, 1)
    assert(toRavelCopy(matrix).sameElements(view))
  }

  test("zero-sized and broadcast shapes retain their logical dimensions") {
    for array <- Seq(NDArray.zeros[Double](0, 3), NDArray.zeros[Double](2, 0)) do
      val matrix = fromRavelCopy(array)
      assertEquals((matrix.rows, matrix.cols), (array.shape(0), array.shape(1)))
      assert(toRavelCopy(matrix).sameElements(array))
    val repeated = NDArray.fromSeq(ravel.Shape(1, 2), Seq(2.0, 3.0)).broadcastTo(ravel.Shape(3, 2))
    assertEquals(fromRavelCopy(repeated).valuesRowMajor, Seq(2.0, 3.0, 2.0, 3.0, 2.0, 3.0))
  }

  test("rank-one conversion copies reversed logical values") {
    val array =
      NDArray.fromSeq(ravel.Shape(4), Seq(1.0, 2.0, 3.0, 4.0)).reverse(0)
    val vector = fromRavelCopy(array)
    assertEquals(vector.toSeq, Seq(4.0, 3.0, 2.0, 1.0))
    assert(toRavelCopy(vector).sameElements(array))
  }

  test("rank-two conversion copies transposed logical values") {
    val array =
      NDArray.tabulate[Double](2, 3)((row, column) => row * 10.0 + column).transpose
    val matrix = fromRavelCopy(array)
    assertEquals((matrix.rows, matrix.cols), (3, 2))
    assertEquals(matrix(2, 1), 12.0)
    assert(toRavelCopy(matrix).sameElements(array))
  }

  test("Gale matrix views materialize in Ravel logical order") {
    val matrix = DMat.tabulate(2, 3)((row, column) => row * 10.0 + column).t
    val array = toRavelCopy(matrix)
    assertEquals(array.shape.toString, "(3, 2)")
    assertEquals(array(2, 1), 12.0)
  }
