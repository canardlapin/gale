package gale.interop.ravel

import munit.FunSuite
import ravel.*
import ravel.DType.given
import ravel.jvm.JvmInterop

final class RavelBorrowedJvmSuite extends FunSuite:
  test("borrowed rank-one overload copies into a vector") {
    val values = Array(1.0, 2.0, 3.0, 4.0)
    val borrowed = JvmInterop.unsafeBorrow(values, Shape(4))
    val vector = fromRavelCopy(borrowed)
    values(0) = 99.0
    assertEquals(vector.toSeq, Seq(1.0, 2.0, 3.0, 4.0))
  }

  test("copying a borrowed array into Gale removes the external alias") {
    val values = Array(1.0, 2.0, 3.0, 4.0)
    val borrowed = JvmInterop.unsafeBorrow(values, Shape(2, 2))
    val matrix = fromRavelCopy(borrowed)
    values(0) = 99.0
    assertEquals(matrix(0, 0), 1.0)
  }
