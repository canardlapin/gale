package gale.interop.ravel

import gale.linalg.*
import munit.FunSuite
import ravel.*
import ravel.DType.given

final class RavelConversionsSuite extends FunSuite:
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
