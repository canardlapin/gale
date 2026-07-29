package gale.interop.ravel

import gale.linalg.{DMat, LinAlgError}
import munit.FunSuite
import ravel.{CanonicalArray, NDArray}
import ravel.DType.given

final class AccessConventionSuite extends FunSuite:
  test("canonical refinement adds linear access without hiding coordinates") {
    val array =
      NDArray.tabulate[Double](2, 3)((row, column) => row * 10.0 + column)
    val canonical = CanonicalArray.require(array)

    assertEquals(canonical(1, -1), 12.0)
    assertEquals(canonical.readLinear(5), 12.0)
  }

  test("Ravel and Gale retain their intentional negative-index boundary") {
    val array =
      NDArray.tabulate[Double](2, 3)((row, column) => row * 10.0 + column)
    assertEquals(array(-1, -1), 12.0)

    val matrix = fromRavelCopy(array)
    assertEquals(matrix(1, 2), 12.0)
    intercept[LinAlgError.IndexOutOfBounds](matrix(-1, 2))
  }

  test("Gale builder separates coordinate and linear writes") {
    val builder = DMat.newBuilder(2, 2)
    builder(0, 1) = 2.0
    builder.writeLinear(3, 4.0)

    val matrix = builder.result()
    assertEquals(matrix(0, 1), 2.0)
    assertEquals(matrix(1, 1), 4.0)
  }
