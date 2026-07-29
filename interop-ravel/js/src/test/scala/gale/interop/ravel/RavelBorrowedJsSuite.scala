package gale.interop.ravel

import munit.FunSuite
import ravel.*
import ravel.DType.given
import ravel.js.JsInterop
import scala.scalajs.js.typedarray.Float64Array

final class RavelBorrowedJsSuite extends FunSuite:
  test("borrowed rank-one overload copies into a vector") {
    val values = new Float64Array(4)
    values(0) = 1.0
    values(1) = 2.0
    values(2) = 3.0
    values(3) = 4.0
    val borrowed = JsInterop.unsafeBorrow(values, Shape(4))
    val vector = fromRavelCopy(borrowed)
    values(0) = 99.0
    assertEquals(vector.toSeq, Seq(1.0, 2.0, 3.0, 4.0))
  }

  test("copying a borrowed typed array into Gale removes the external alias") {
    val values = new Float64Array(4)
    values(0) = 1.0
    values(1) = 2.0
    values(2) = 3.0
    values(3) = 4.0
    val borrowed = JsInterop.unsafeBorrow(values, Shape(2, 2))
    val matrix = fromRavelCopy(borrowed)
    values(0) = 99.0
    assertEquals(matrix(0, 0), 1.0)
  }
