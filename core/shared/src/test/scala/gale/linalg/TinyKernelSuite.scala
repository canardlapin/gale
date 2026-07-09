package gale.linalg

class TinyKernelSuite extends munit.FunSuite:
  test("Mat2 uses unrolled matrix-vector multiply and determinant") {
    val A = Mat2(1.0, 2.0, 3.0, 4.0)
    val y = A * Vec2(5.0, 6.0)

    assertEquals(y, Vec2(17.0, 39.0))
    assertEquals(A.det, -2.0)
    assertEquals(A.toDMat.valuesRowMajor, Seq(1.0, 2.0, 3.0, 4.0))
  }

  test("Mat3 uses unrolled matrix-vector multiply and determinant") {
    val A = Mat3(
      1.0, 2.0, 3.0,
      0.0, 1.0, 4.0,
      5.0, 6.0, 0.0
    )
    val y = A * Vec3(1.0, 2.0, 3.0)

    assertEquals(y, Vec3(14.0, 14.0, 17.0))
    assertEquals(A.det, 1.0)
  }

  test("Mat4 uses unrolled matrix-vector multiply and determinant expansion") {
    val A = Mat4(
      1.0, 0.0, 0.0, 0.0,
      0.0, 2.0, 0.0, 0.0,
      0.0, 0.0, 3.0, 0.0,
      0.0, 0.0, 0.0, 4.0
    )
    val y = A * Vec4(1.0, 2.0, 3.0, 4.0)

    assertEquals(y, Vec4(1.0, 4.0, 9.0, 16.0))
    assertEquals(A.det, 24.0)
    assertEquals(A.toDMat.valuesRowMajor, Seq(
      1.0, 0.0, 0.0, 0.0,
      0.0, 2.0, 0.0, 0.0,
      0.0, 0.0, 3.0, 0.0,
      0.0, 0.0, 0.0, 4.0
    ))
  }
