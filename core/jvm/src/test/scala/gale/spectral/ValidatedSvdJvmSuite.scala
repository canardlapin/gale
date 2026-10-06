package gale.spectral

import gale.linalg.{DMat, Matrix}
import munit.FunSuite

class ValidatedSvdJvmSuite extends FunSuite:
  private def rejectsConstructor(kind: Class[?], args: AnyRef*): Unit =
    val constructor = kind.getDeclaredConstructors.find(_.getParameterCount == args.size).get
    constructor.setAccessible(true)
    val error = intercept[java.lang.reflect.InvocationTargetException] {
      constructor.newInstance(args*)
    }
    assert(error.getCause.isInstanceOf[IllegalArgumentException])

  test("scalar constructor privacy bypass cannot forge invalid endpoints"):
    rejectsConstructor(classOf[RealInterval], Double.box(2.0), Double.box(1.0))
    rejectsConstructor(classOf[RealInterval], Double.box(Double.NaN), Double.box(1.0))
    rejectsConstructor(classOf[RealInterval], Double.box(0.0), Double.box(Double.PositiveInfinity))

  test("matrix constructor privacy bypass cannot forge invalid enclosure geometry"):
    val one = Matrix.dense(1, 1)(1.0)
    rejectsConstructor(classOf[MatrixEnclosure], one, DMat.zeros(2, 1))
    rejectsConstructor(classOf[MatrixEnclosure], one, DMat.zeros(1, 1))
    rejectsConstructor(classOf[MatrixEnclosure], Matrix.dense(1, 1)(Double.NaN), one)
