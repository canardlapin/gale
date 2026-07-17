package gale.backend.jvm.`native`

import gale.backend.jvm.`native`.NativeDMat.toNative
import gale.linalg.Matrix

import java.lang.foreign.Arena

class NativeDMatSuite extends munit.FunSuite:
  private val matrix = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

  private def assertMatrixEquals(actual: gale.linalg.DMat, expected: gale.linalg.DMat): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        assertEquals(actual(i, j), expected(i, j))
        j += 1
      i += 1

  test("owned row-major and column-major matrices round-trip"):
    for layout <- Layout.values do
      val native = NativeDMat.allocate(2, 3, layout)
      try
        var i = 0
        while i < 2 do
          var j = 0
          while j < 3 do
            native(i, j) = matrix(i, j)
            j += 1
          i += 1
        assertMatrixEquals(native.toHeap, matrix)
      finally native.close()
      assert(!native.isAlive)
      intercept[IllegalStateException](native(0, 0))

  test("borrowed toNative lifetime follows caller arena"):
    val arena = Arena.ofConfined()
    given Arena = arena
    val native = matrix.toNative(Layout.ColMajor)
    assertMatrixEquals(native.toHeap, matrix)
    native.close() // borrowed: must not close the caller's arena
    assert(native.isAlive)
    arena.close()
    assert(!native.isAlive)
    intercept[IllegalStateException](native.toHeap)

  test("shape, leading dimension, and indices are checked"):
    intercept[IllegalArgumentException](NativeDMat.allocate(-1, 2))
    val native = NativeDMat.allocate(2, 2)
    try
      intercept[IndexOutOfBoundsException](native(2, 0))
      intercept[IndexOutOfBoundsException](native(0, -1))
    finally native.close()
