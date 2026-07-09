package gale.linalg

class ShapeOverflowSuite extends munit.FunSuite:
  // Item 4: rows * cols must not silently wrap around Int. A 65536x65536 matrix
  // has 2^32 elements, which overflows Int to 0 and would allocate a zero-length
  // buffer while claiming the full shape.
  test("Shape.size does not overflow Int") {
    val shape = Shape(Rows(65536), Cols(65536))
    assert(shape.size != 0L, s"Shape.size overflowed to ${shape.size}")
    assert(shape.size == 65536L * 65536L, s"Shape.size = ${shape.size}")
  }

  test("DMat.zeros rejects an element count that overflows Int") {
    intercept[LinAlgError.InvalidArgument] {
      DMat.zeros(65536, 65536)
    }
  }

  test("Matrix.dense rejects an element count that overflows Int") {
    intercept[LinAlgError.InvalidArgument] {
      Matrix.dense(50000, 50000, Seq.empty)
    }
  }

  test("DMat.fromArrayRowMajor rejects an element count that overflows Int") {
    intercept[LinAlgError.InvalidArgument] {
      DMat.fromArrayRowMajor(65536, 65536, Array.empty)
    }
  }

  test("valid shapes remain constructible") {
    val a = DMat.zeros(3, 4)
    assert(a.shape.size == 12L, s"unexpected size ${a.shape.size}")
  }
