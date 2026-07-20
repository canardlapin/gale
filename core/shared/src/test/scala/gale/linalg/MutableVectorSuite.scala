package gale.linalg

class MutableVectorSuite extends munit.FunSuite:
  test("DVecBuilder transfers primitive-loop storage and closes permanently") {
    val builder = DVec.newBuilder(3)
    var i = 0
    while i < builder.length do
      builder(i) = i.toDouble + 1.0
      i += 1
    builder(2) = 9.0

    val vector = builder.result()

    assertEquals(vector.toSeq, Seq(1.0, 2.0, 9.0))
    intercept[LinAlgError.UnsupportedOperation](builder.update(0, -1.0))
    intercept[LinAlgError.UnsupportedOperation](builder.result())
    assertEquals(vector(0), 1.0)
  }

  test("DVecBuilder supports empty vectors and validates logical indices") {
    val empty = Vec.newBuilder(0).result()
    assertEquals(empty.length, 0)

    val builder = DVecBuilder.zeros(2)
    builder.fill(4.0)
    assertEquals(builder(1), 4.0)
    intercept[LinAlgError.IndexOutOfBounds](builder.update(-1, 1.0))
    intercept[LinAlgError.IndexOutOfBounds](builder.update(2, 1.0))
  }
