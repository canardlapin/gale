package gale.sparse

import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.IndexArray

class SparseVectorSuite extends munit.FunSuite:
  private def sv(length: Int, entries: (Int, Double)*): SparseVector =
    SparseVector.fromEntries(length, entries)

  test("entries are sorted and duplicates resolve under each policy") {
    val entries = Seq(3 -> 1.0, 0 -> 2.0, 3 -> 4.0, 1 -> -1.0, 3 -> 0.5)
    val summed = SparseVector.fromEntries(5, entries)
    assertEquals(summed.activeIndices, IndexedSeq(0, 1, 3))
    assertEquals(summed.activeValues, IndexedSeq(2.0, -1.0, 5.5))
    val last = SparseVector.fromEntries(5, entries, duplicates = DuplicatePolicy.Last)
    assertEquals(last.activeValues, IndexedSeq(2.0, -1.0, 0.5))
    assertEquals(
      SparseVector.tryFromEntries(5, entries, duplicates = DuplicatePolicy.Error),
      Left(LinAlgError.InvalidArgument("duplicate sparse entry at index 3"))
    )
    intercept[LinAlgError.InvalidArgument] {
      SparseVector.fromEntries(5, entries, duplicates = DuplicatePolicy.Error)
    }
  }

  test("Last keeps input order among equal indices even when unsorted") {
    val v = SparseVector.fromEntries(4, Seq(2 -> 9.0, 1 -> 1.0, 2 -> 7.0, 0 -> 3.0, 2 -> 8.0), DuplicatePolicy.Last)
    assertEquals(v(2), 8.0)
  }

  test("explicit zeros are stored until compact, which also drops negative zero") {
    val v = sv(6, 1 -> 0.0, 2 -> 3.0, 4 -> -0.0, 5 -> 1.0, 5 -> -1.0)
    assertEquals(v.activeSize, 4)
    assertEquals(v.activeIndices, IndexedSeq(1, 2, 4, 5))
    val compacted = v.compact
    assertEquals(compacted.activeIndices, IndexedSeq(2))
    assertEquals(compacted.toDense.toSeq, v.toDense.toSeq)
    assert(compacted.compact eq compacted)
    assertEquals(SparseVector.fromDense(Vec(0.0, -0.0, 2.0)).activeIndices, IndexedSeq(2))
  }

  test("checked construction reports range, length, and array-shape errors as values") {
    assertEquals(SparseVector.tryFromEntries(3, Seq(3 -> 1.0)), Left(LinAlgError.IndexOutOfBounds(3, 3)))
    assertEquals(SparseVector.tryFromEntries(3, Seq(-1 -> 1.0)), Left(LinAlgError.IndexOutOfBounds(-1, 3)))
    assertEquals(
      SparseVector.tryFromEntries(-1, Nil),
      Left(LinAlgError.InvalidArgument("sparse vector length must be non-negative, got -1"))
    )
    assert(SparseVector.tryFromArrays(3, Array(0, 1), Array(1.0)).isLeft)
    assertEquals(
      SparseVector.tryFromArrays(3, Array(2, 0), Array(5.0, 6.0)).map(_.toDense.toSeq),
      Right(Seq(6.0, 0.0, 5.0))
    )
    intercept[IllegalArgumentException](SparseVector.zeros(-1))
  }

  test("tryFromArrays does not retain caller arrays") {
    val idx = Array(0, 2)
    val vals = Array(1.0, 2.0)
    val v = SparseVector.tryFromArrays(3, idx, vals).toOption.get
    idx(0) = 1
    vals(0) = 99.0
    assertEquals(v.toDense.toSeq, Seq(1.0, 0.0, 2.0))
  }

  test("non-finite values follow SparseValuePolicy") {
    val entries = Seq(0 -> Double.NaN, 1 -> Double.PositiveInfinity)
    val allowed = SparseVector.fromEntries(2, entries)
    assert(allowed(0).isNaN)
    assertEquals(allowed(1), Double.PositiveInfinity)
    assertEquals(
      SparseVector.tryFromEntries(2, entries, valuePolicy = SparseValuePolicy.RequireFinite),
      Left(LinAlgError.InvalidArgument("non-finite sparse value NaN at index 0"))
    )
    assert(
      SparseVector
        .tryFromEntries(2, Seq(1 -> Double.NegativeInfinity), valuePolicy = SparseValuePolicy.RequireFinite)
        .isLeft
    )
  }

  test("length-0 vector: empty products and norms, typed failure for max and min") {
    val empty = SparseVector.zeros(0)
    assertEquals(empty.activeSize, 0)
    assertEquals(empty.toDense.length, 0)
    assertEquals(empty.dot(empty), 0.0)
    assertEquals(empty.sum, 0.0)
    assertEquals(empty.norm1, 0.0)
    assertEquals(empty.norm2, 0.0)
    assertEquals(empty.normInf, 0.0)
    val maxError = intercept[LinAlgError.EmptyInput](empty.max)
    assertEquals(maxError.getMessage, "SparseVector.max requires a non-empty input")
    assertEquals(intercept[LinAlgError.EmptyInput](empty.min), LinAlgError.EmptyInput("SparseVector.min"))
    assertEquals((empty + empty).activeSize, 0)
  }

  test("length-1 vectors with and without the active entry") {
    val one = sv(1, 0 -> -2.0)
    val none = SparseVector.zeros(1)
    assertEquals(one.max, -2.0)
    assertEquals(one.min, -2.0)
    assertEquals(none.max, 0.0)
    assertEquals(none.min, 0.0)
    assertEquals(one.dot(none), 0.0)
    assertEquals((one - none).toDense.toSeq, Seq(-2.0))
    assertEquals((none - one).toDense.toSeq, Seq(2.0))
  }

  test("max and min include implicit zeros only when one exists") {
    val negatives = sv(3, 0 -> -1.0, 2 -> -3.0)
    assertEquals(negatives.max, 0.0)
    assertEquals(negatives.min, -3.0)
    val full = sv(2, 0 -> -1.0, 1 -> -3.0)
    assertEquals(full.max, -1.0)
    val positives = sv(3, 1 -> 4.0)
    assertEquals(positives.min, 0.0)
    assertEquals(positives.max, 4.0)
  }

  test("NaN propagates through max, min and norms") {
    val v = sv(4, 0 -> 1.0, 2 -> Double.NaN)
    assert(v.max.isNaN)
    assert(v.min.isNaN)
    assert(v.norm1.isNaN)
    assert(v.norm2.isNaN)
    assert(v.normInf.isNaN)
    assert(v.sum.isNaN)
  }

  test("norm2 neither overflows nor underflows and maps infinities to +Infinity") {
    val huge = sv(4, 0 -> 1e200, 3 -> -1e200)
    assertEqualsDouble(huge.norm2, 1e200 * math.sqrt(2.0), 1e186)
    val tiny = sv(4, 1 -> 3e-200, 2 -> 4e-200)
    assertEqualsDouble(tiny.norm2, 5e-200, 1e-213)
    val infinite = sv(3, 0 -> Double.NegativeInfinity, 2 -> Double.PositiveInfinity)
    assertEquals(infinite.norm2, Double.PositiveInfinity)
    assertEquals(infinite.normInf, Double.PositiveInfinity)
  }

  test("products skip implicit zeros, so a non-finite value facing one contributes nothing") {
    val a = sv(3, 0 -> Double.NaN, 1 -> 2.0)
    val b = sv(3, 1 -> 3.0)
    assertEquals(a.dot(b), 6.0)
    assertEquals(b.dot(Vec(Double.NaN, 1.0, Double.PositiveInfinity)), 3.0)
  }

  test("merge covers disjoint, overlapping and identical patterns; cancellation stays stored") {
    val a = sv(6, 0 -> 1.0, 2 -> 2.0)
    val disjoint = sv(6, 1 -> 5.0, 5 -> 6.0)
    assertEquals((a + disjoint).activeIndices, IndexedSeq(0, 1, 2, 5))
    assertEquals((a - disjoint).toDense.toSeq, Seq(1.0, -5.0, 2.0, 0.0, 0.0, -6.0))
    val overlapping = sv(6, 2 -> 2.0, 3 -> 1.0)
    val difference = a - overlapping
    assertEquals(difference.activeIndices, IndexedSeq(0, 2, 3))
    assertEquals(difference(2), 0.0)
    assertEquals(difference.compact.activeIndices, IndexedSeq(0, 3))
    val identical = a - a
    assertEquals(identical.activeIndices, a.activeIndices)
    assertEquals(identical.compact.activeSize, 0)
    assertEquals((a + a).activeValues, IndexedSeq(2.0, 4.0))
    assertEquals(a.dot(disjoint), 0.0)
    assertEquals(a.dot(overlapping), 4.0)
  }

  test("scaling by zero keeps the stored pattern") {
    val v = sv(4, 1 -> 3.0, 3 -> -1.0) * 0.0
    assertEquals(v.activeIndices, IndexedSeq(1, 3))
    assertEquals(v.compact.activeSize, 0)
  }

  test("dot accepts strided dense views") {
    val matrix = Matrix.dense(3, 2)(
      1.0, 10.0,
      2.0, 20.0,
      3.0, 30.0
    )
    val column = matrix.col(1)
    assertEquals(column.stride.value, 2)
    assertEquals(sv(3, 0 -> 1.0, 2 -> 2.0).dot(column), 70.0)
  }

  test("axpyInto updates only active positions of the destination") {
    val y = MutableDVec.from(Vec(1.0, 1.0, 1.0, 1.0))
    sv(4, 1 -> 2.0, 3 -> -1.0).axpyInto(3.0, y)
    assertEquals(y.toVec.toSeq, Seq(1.0, 7.0, 1.0, -2.0))
  }

  test("length mismatches and out-of-range access are typed errors") {
    val a = SparseVector.zeros(3)
    val b = SparseVector.zeros(4)
    assertEquals(intercept[LinAlgError.VectorLengthMismatch](a.dot(b)), LinAlgError.VectorLengthMismatch(3, 4))
    intercept[LinAlgError.VectorLengthMismatch](a + b)
    intercept[LinAlgError.VectorLengthMismatch](a - b)
    intercept[LinAlgError.VectorLengthMismatch](a.dot(Vec(1.0)))
    intercept[LinAlgError.VectorLengthMismatch](a.axpyInto(1.0, MutableDVec.zeros(2)))
    assertEquals(intercept[LinAlgError.IndexOutOfBounds](a(3)), LinAlgError.IndexOutOfBounds(3, 3))
    intercept[LinAlgError.IndexOutOfBounds](a(-1))
  }

  test("foreachActive visits entries in index order") {
    val seen = scala.collection.mutable.ArrayBuffer.empty[(Int, Double)]
    sv(5, 4 -> 1.0, 0 -> 2.0, 2 -> 0.0).foreachActive((i, v) => seen += i -> v)
    assertEquals(seen.toSeq, Seq(0 -> 2.0, 2 -> 0.0, 4 -> 1.0))
  }

  test("toString lists active entries") {
    // Non-integral values: Scala.js prints integral doubles without ".0".
    assertEquals(sv(5, 3 -> -0.25, 1 -> 1.5).toString, "SparseVector(5)(1 -> 1.5, 3 -> -0.25)")
  }

  test("CSR and CSC sparse slices keep stored zeros and validate their index") {
    val csr = Sparse.coo(2, 3).add(0, 0, 1.0).add(0, 2, 2.0).add(1, 1, 3.0).toCSR()
    val zeroed = csr.rebind(Array(1.0, 0.0, 3.0)).toOption.get
    assertEquals(zeroed.rowSparse(0).activeIndices, IndexedSeq(0, 2))
    assertEquals(zeroed.rowSparse(0).activeValues, IndexedSeq(1.0, 0.0))
    val csc = zeroed.toCSC
    assertEquals(csc.colSparse(2).activeIndices, IndexedSeq(0))
    assertEquals(csc.colSparse(2).activeValues, IndexedSeq(0.0))
    intercept[LinAlgError.IndexOutOfBounds](csr.rowSparse(2))
    intercept[LinAlgError.IndexOutOfBounds](csc.colSparse(-1))
  }

  test("rowSparse sorts and sums the duplicate columns of a non-canonical CSR") {
    val raw = new CSR(
      1,
      4,
      IndexArray.fromArray(Array(0, 4)),
      IndexArray.fromArray(Array(3, 1, 3, 0)),
      DoubleArray.fromArray(Array(1.0, 2.0, 4.0, 5.0))
    )
    val row = raw.rowSparse(0)
    assertEquals(row.activeIndices, IndexedSeq(0, 1, 3))
    assertEquals(row.activeValues, IndexedSeq(5.0, 2.0, 5.0))
    assertEquals((raw * sv(4, 3 -> 1.0)).toSeq, Seq(5.0))
  }

  test("CSR * SparseVector equals the dense product and checks its shape") {
    val csr = Sparse.coo(2, 3).add(0, 0, 1.0).add(0, 2, 2.0).add(1, 1, 3.0).toCSR()
    val x = sv(3, 2 -> 4.0, 1 -> -1.0)
    assertEquals((csr * x).toSeq, Seq(8.0, -3.0))
    intercept[LinAlgError.DimensionMismatch](csr * SparseVector.zeros(2))
  }
