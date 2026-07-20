package gale.sparse

import gale.TestAccess
import gale.linalg.LinAlgError
import gale.linalg.Vec
import scala.collection.mutable.ArrayBuffer

class CompressedPatternSuite extends munit.FunSuite:
  private def samplePattern: CSRPattern =
    CSRPattern
      .checked(
        rows = 3,
        cols = 4,
        rowOffsets = Array(0, 2, 3, 5),
        columnIndices = Array(0, 2, 1, 0, 3)
      )
      .toOption
      .get

  test("checked CSR patterns reject invalid shapes, offsets, and indices") {
    assert(CSRPattern.checked(-1, 2, Array(0), Array.emptyIntArray).isLeft)
    assert(CSRPattern.checked(2, 2, Array(0, 0), Array.emptyIntArray).isLeft)
    assert(CSRPattern.checked(2, 2, Array(1, 1, 1), Array.emptyIntArray).isLeft)
    assert(CSRPattern.checked(2, 2, Array(0, 2, 1), Array(0, 1)).isLeft)
    assert(CSRPattern.checked(2, 2, Array(0, 1, 1), Array(0, 1)).isLeft)
    assert(CSRPattern.checked(1, 2, Array(0, 1), Array(2)).isLeft)
    assert(CSRPattern.checked(1, 2, Array(0, 1), Array(-1)).isLeft)
  }

  test("checked CSC patterns use the same compressed invariant model") {
    val pattern = CSCPattern.checked(3, 2, Array(0, 2, 3), Array(0, 2, 1)).toOption.get
    assertEquals(pattern.rows, 3)
    assertEquals(pattern.cols, 2)
    assertEquals(pattern.nnz, 3)
    assertEquals(pattern.columnNnz(0), 2)
    assertEquals(pattern.columnNnz(1), 1)
    assert(pattern.hasCanonicalFormat)
    assert(CSCPattern.checked(3, 2, Array(0, 2, 3), Array(0, 3, 1)).isLeft)
  }

  test("pattern traversal and lane queries preserve deterministic storage order") {
    val pattern = samplePattern
    val positions = ArrayBuffer.empty[(Int, Int)]
    pattern.foreachStoredPosition((row, col) => positions += ((row, col)))
    assertEquals(positions.toSeq, Seq((0, 0), (0, 2), (1, 1), (2, 0), (2, 3)))
    assertEquals(pattern.rowNnz(0), 2)
    assertEquals(pattern.rowNnz(1), 1)
    assertEquals(pattern.rowNnz(2), 2)
    assertEquals(pattern.columnIndexAt(3), 0)

    val rowTwo = ArrayBuffer.empty[(Int, Int)]
    pattern.foreachRow(2)((col, storedIndex) => rowTwo += ((col, storedIndex)))
    assertEquals(rowTwo.toSeq, Seq((0, 3), (3, 4)))
    intercept[LinAlgError.IndexOutOfBounds](pattern.rowNnz(3))
    intercept[LinAlgError.IndexOutOfBounds](pattern.columnIndexAt(5))
  }

  test("transpose shares immutable structure and preserves structural equality") {
    val csr = samplePattern
    val expected =
      CSCPattern.checked(4, 3, Array(0, 2, 3, 5), Array(0, 2, 1, 0, 3)).toOption.get
    assertEquals(csr.t, expected)
    assertEquals(csr.t.hashCode(), expected.hashCode())
    assertEquals(csr.t.t, csr)
    assert(csr.t.storage eq csr.storage)

    val nonCanonical = CSRPattern.checked(1, 4, Array(0, 3), Array(2, 2, 1)).toOption.get
    assert(!nonCanonical.hasCanonicalFormat)
  }

  test("empty rectangular patterns bind and transpose without special cases") {
    val zeroRows = CSRPattern.checked(0, 3, Array(0), Array.emptyIntArray).toOption.get
    val matrix = zeroRows.bind(Array.emptyDoubleArray).toOption.get
    assertEquals(matrix.rows, 0)
    assertEquals(matrix.cols, 3)
    assertEquals(matrix.nnz, 0)
    assertEquals(zeroRows.t.rows, 3)
    assertEquals(zeroRows.t.cols, 0)

    val zeroCols = CSRPattern.checked(2, 0, Array(0, 0, 0), Array.emptyIntArray).toOption.get
    assertEquals(zeroCols.bind(Vec()).toOption.get.toDense().valuesRowMajor, Seq.empty)
  }

  test("pattern and value factories own their inputs") {
    val offsets = Array(0, 1, 2)
    val indices = Array(0, 1)
    val pattern = CSRPattern.checked(2, 2, offsets, indices).toOption.get
    offsets(1) = 0
    indices(0) = 1

    val values = Array(2.0, 3.0)
    val matrix = pattern.bind(values).toOption.get
    values(0) = 99.0
    assertEquals(matrix.toDense().valuesRowMajor, Seq(2.0, 0.0, 0.0, 3.0))

    val strided = TestAccess.stridedCopy(Vec(5.0, 7.0), stride = 3)
    val fromStrided = pattern.bind(strided).toOption.get
    assertEquals(fromStrided.toDense().valuesRowMajor, Seq(5.0, 0.0, 0.0, 7.0))
  }

  test("one pattern rebinds independent values without copying structure") {
    val pattern = samplePattern
    val first = pattern.bind(Array(1.0, 2.0, 3.0, 4.0, 5.0)).toOption.get
    val second = pattern.bind(Array(5.0, 4.0, 3.0, 2.0, 1.0)).toOption.get
    assert(first.rowPtr.asInstanceOf[AnyRef] eq pattern.storage.majorOffsets.asInstanceOf[AnyRef])
    assert(first.colIdx.asInstanceOf[AnyRef] eq pattern.storage.minorIndices.asInstanceOf[AnyRef])
    assert(second.rowPtr.asInstanceOf[AnyRef] eq pattern.storage.majorOffsets.asInstanceOf[AnyRef])
    assert(second.colIdx.asInstanceOf[AnyRef] eq pattern.storage.minorIndices.asInstanceOf[AnyRef])
    assertEquals(first.pattern, second.pattern)
    assertEquals(first(0, 0), 1.0)
    assertEquals(second(0, 0), 5.0)
    assert(pattern.bind(Array(1.0)).isLeft)
  }

  test("single-owner values builders transfer storage and close") {
    val pattern = samplePattern
    val builder = pattern.valuesBuilder()
    var i = 0
    while i < builder.length do
      builder(i) = (i + 1).toDouble
      i += 1
    val matrix = builder.result()
    assertEquals(matrix.toDense().valuesRowMajor, Seq(1.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 4.0, 0.0, 0.0, 5.0))
    intercept[LinAlgError.UnsupportedOperation](builder(0))
    intercept[LinAlgError.UnsupportedOperation](builder.result())

    val cscBuilder = pattern.t.valuesBuilder(4.0)
    val csc = cscBuilder.result()
    assertEquals(csc.nnz, pattern.nnz)
    assertEquals(csc(0, 0), 4.0)
  }

  test("numeric updates retain explicit zeros until an explicit prune") {
    val matrix = samplePattern.bind(Array(1.0, 2.0, 3.0, 4.0, 5.0)).toOption.get
    val zeroed = matrix.mapValues(_ => 0.0)
    assertEquals(zeroed.nnz, matrix.nnz)
    assertEquals(zeroed.pattern, matrix.pattern)
    assert(!zeroed.hasCanonicalFormat, "explicit numeric zeros are intentionally non-canonical")
    assertEquals(zeroed.toDense().valuesRowMajor, Seq.fill(12)(0.0))

    val pruned = zeroed.pruneZeros
    assertEquals(pruned.nnz, 0)
    assertNotEquals(pruned.pattern, matrix.pattern)

    val rebound = matrix.rebind(Array.fill(matrix.nnz)(0.0)).toOption.get
    assertEquals(rebound.nnz, matrix.nnz)
    assertEquals(rebound.pattern, matrix.pattern)
    assertEquals((matrix * 0.0).nnz, matrix.nnz)
  }

  test("four-cycle adjacency reuses one graph pattern across edge weights") {
    val cycle =
      CSRPattern
        .checked(
          4,
          4,
          Array(0, 2, 4, 6, 8),
          Array(1, 3, 0, 2, 1, 3, 0, 2)
        )
        .toOption
        .get
    val adjacency = cycle.bind(Array.fill(8)(1.0)).toOption.get
    val weighted = cycle.bind(Array(2.0, 4.0, 2.0, 3.0, 3.0, 5.0, 4.0, 5.0)).toOption.get

    assertEquals((adjacency * Vec(1.0, 2.0, 3.0, 4.0)).toSeq, Seq(6.0, 4.0, 6.0, 4.0))
    assertEquals(adjacency.pattern, weighted.pattern)
    assertEquals(adjacency.pattern.t, adjacency.t.pattern)
  }
