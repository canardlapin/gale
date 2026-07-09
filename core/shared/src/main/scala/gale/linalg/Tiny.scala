package gale.linalg

final case class Vec2(x0: Double, x1: Double):
  def toDVec: DVec =
    Vec(x0, x1)

final case class Vec3(x0: Double, x1: Double, x2: Double):
  def toDVec: DVec =
    Vec(x0, x1, x2)

final case class Vec4(x0: Double, x1: Double, x2: Double, x3: Double):
  def toDVec: DVec =
    Vec(x0, x1, x2, x3)

final case class Mat2(
    a00: Double,
    a01: Double,
    a10: Double,
    a11: Double
):
  def *(x: Vec2): Vec2 =
    Vec2(
      a00 * x.x0 + a01 * x.x1,
      a10 * x.x0 + a11 * x.x1
    )

  def det: Double =
    a00 * a11 - a01 * a10

  def toDMat: DMat =
    Matrix.dense(2, 2)(
      a00, a01,
      a10, a11
    )

final case class Mat3(
    a00: Double,
    a01: Double,
    a02: Double,
    a10: Double,
    a11: Double,
    a12: Double,
    a20: Double,
    a21: Double,
    a22: Double
):
  def *(x: Vec3): Vec3 =
    Vec3(
      a00 * x.x0 + a01 * x.x1 + a02 * x.x2,
      a10 * x.x0 + a11 * x.x1 + a12 * x.x2,
      a20 * x.x0 + a21 * x.x1 + a22 * x.x2
    )

  def det: Double =
    a00 * (a11 * a22 - a12 * a21) -
      a01 * (a10 * a22 - a12 * a20) +
      a02 * (a10 * a21 - a11 * a20)

  def toDMat: DMat =
    Matrix.dense(3, 3)(
      a00, a01, a02,
      a10, a11, a12,
      a20, a21, a22
    )

final case class Mat4(
    a00: Double,
    a01: Double,
    a02: Double,
    a03: Double,
    a10: Double,
    a11: Double,
    a12: Double,
    a13: Double,
    a20: Double,
    a21: Double,
    a22: Double,
    a23: Double,
    a30: Double,
    a31: Double,
    a32: Double,
    a33: Double
):
  def *(x: Vec4): Vec4 =
    Vec4(
      a00 * x.x0 + a01 * x.x1 + a02 * x.x2 + a03 * x.x3,
      a10 * x.x0 + a11 * x.x1 + a12 * x.x2 + a13 * x.x3,
      a20 * x.x0 + a21 * x.x1 + a22 * x.x2 + a23 * x.x3,
      a30 * x.x0 + a31 * x.x1 + a32 * x.x2 + a33 * x.x3
    )

  def det: Double =
    a00 * det3(a11, a12, a13, a21, a22, a23, a31, a32, a33) -
      a01 * det3(a10, a12, a13, a20, a22, a23, a30, a32, a33) +
      a02 * det3(a10, a11, a13, a20, a21, a23, a30, a31, a33) -
      a03 * det3(a10, a11, a12, a20, a21, a22, a30, a31, a32)

  def toDMat: DMat =
    Matrix.dense(4, 4)(
      a00, a01, a02, a03,
      a10, a11, a12, a13,
      a20, a21, a22, a23,
      a30, a31, a32, a33
    )

  private inline def det3(
      b00: Double,
      b01: Double,
      b02: Double,
      b10: Double,
      b11: Double,
      b12: Double,
      b20: Double,
      b21: Double,
      b22: Double
  ): Double =
    b00 * (b11 * b22 - b12 * b21) -
      b01 * (b10 * b22 - b12 * b20) +
      b02 * (b10 * b21 - b11 * b20)
