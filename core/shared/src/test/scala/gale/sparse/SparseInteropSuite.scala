package gale.sparse

import gale.TestAccess
import gale.linalg.LinAlgError

class SparseInteropSuite extends munit.FunSuite:

  private def storedEntries(matrix: COO): Vector[(Int, Int, Double)] =
    val out = Vector.newBuilder[(Int, Int, Double)]
    matrix.foreachStoredEntry((row, col, value) => out += ((row, col, value)))
    out.result()

  private def storedEntries(matrix: CSR): Vector[(Int, Int, Double)] =
    val out = Vector.newBuilder[(Int, Int, Double)]
    matrix.foreachStoredEntry((row, col, value) => out += ((row, col, value)))
    out.result()

  test("COO primitive traversal preserves physical order, duplicates, zeros, and non-finite values") {
    val coo = TestAccess.coo(
      rows = 3,
      cols = 4,
      entries = Seq(
        (2, 1, 0.0),
        (0, 3, Double.NaN),
        (2, 1, Double.PositiveInfinity),
        (1, 0, -2.0)
      ),
      canonical = false
    )
    val rows = new Array[Int](coo.nnz)
    val cols = new Array[Int](coo.nnz)
    val values = new Array[Double](coo.nnz)
    var count = 0
    coo.foreachStoredEntry: (row, col, value) =>
      rows(count) = row
      cols(count) = col
      values(count) = value
      count += 1

    assertEquals(count, 4)
    assertEquals(rows.toSeq, Seq(2, 0, 2, 1))
    assertEquals(cols.toSeq, Seq(1, 3, 1, 0))
    assertEquals(values(0), 0.0)
    assert(values(1).isNaN)
    assertEquals(values(2), Double.PositiveInfinity)
    assertEquals(values(3), -2.0)
  }

  test("CSR primitive traversal is deterministic row-major storage order") {
    val builder = Sparse.cooChecked(3, 3).toOption.get
    assert(builder.tryAdd(2, 1, 4.0).isRight)
    assert(builder.tryAdd(0, 2, 1.0).isRight)
    assert(builder.tryAdd(0, 2, 3.0).isRight)
    assert(builder.tryAdd(1, 0, 0.0).isRight)
    assert(builder.tryAdd(2, 0, -2.0).isRight)

    val csr = builder.tryToCSR(DuplicatePolicy.Sum).toOption.get
    assert(csr.hasCanonicalFormat)
    assertEquals(
      storedEntries(csr),
      Vector((0, 2, 4.0), (2, 0, -2.0), (2, 1, 4.0))
    )
  }

  test("checked builder factory and tryAdd make shape and index failures total") {
    Sparse.cooChecked(-1, 3) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected invalid shape Left, got $other")
    COOBuilder.checked(3, -1) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected invalid shape Left, got $other")

    val builder = Sparse.cooChecked(2, 3).toOption.get
    assertEquals(builder.nnz, 0)
    assert(builder.tryAdd(-1, 0, 1.0).left.exists(_.isInstanceOf[LinAlgError.IndexOutOfBounds]))
    assert(builder.tryAdd(2, 0, 1.0).left.exists(_.isInstanceOf[LinAlgError.IndexOutOfBounds]))
    assert(builder.tryAdd(0, -1, 1.0).left.exists(_.isInstanceOf[LinAlgError.IndexOutOfBounds]))
    assert(builder.tryAdd(0, 3, 1.0).left.exists(_.isInstanceOf[LinAlgError.IndexOutOfBounds]))
    assertEquals(builder.nnz, 0)

    val zero = Sparse.cooChecked(0, 3).toOption.get
    assert(zero.tryAdd(0, 0, 1.0).isLeft)
    val zeroCsr = zero.tryToCSR().toOption.get
    assertEquals((zeroCsr.rows, zeroCsr.cols, zeroCsr.nnz), (0, 3, 0))
    assert(zeroCsr.hasCanonicalFormat)
    val noColumns = Sparse.cooChecked(3, 0).toOption.get.tryToCSR().toOption.get
    assertEquals((noColumns.rows, noColumns.cols, noColumns.nnz), (3, 0, 0))
    assertEquals(Sparse.cooChecked(0, 0).toOption.get.tryToCOO().toOption.get.nnz, 0)

    val maximumRows = Sparse.cooChecked(Int.MaxValue, 0).toOption.get
    assert(maximumRows.tryToCSR().left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    val maximumRowsCsc = maximumRows.tryToCSC().toOption.get
    assertEquals((maximumRowsCsc.rows, maximumRowsCsc.cols, maximumRowsCsc.nnz), (Int.MaxValue, 0, 0))

    val maximumColumns = Sparse.cooChecked(0, Int.MaxValue).toOption.get
    assert(maximumColumns.tryToCSC().left.exists(_.isInstanceOf[LinAlgError.InvalidArgument]))
    val maximumColumnsCsr = maximumColumns.tryToCSR().toOption.get
    assertEquals((maximumColumnsCsr.rows, maximumColumnsCsr.cols, maximumColumnsCsr.nnz), (0, Int.MaxValue, 0))
  }

  test("non-finite policy is explicit, allowed by default, and strict on request") {
    val permissive = Sparse.cooChecked(1, 3).toOption.get
    assertEquals(permissive.valuePolicy, SparseValuePolicy.AllowNonFinite)
    assert(permissive.tryAdd(0, 0, Double.NaN).isRight)
    assert(permissive.tryAdd(0, 1, Double.PositiveInfinity).isRight)
    assert(permissive.tryAdd(0, 2, Double.NegativeInfinity).isRight)
    val values = storedEntries(permissive.tryToCOO().toOption.get).map(_._3)
    assert(values(0).isNaN)
    assertEquals(values(1), Double.PositiveInfinity)
    assertEquals(values(2), Double.NegativeInfinity)

    val strict = Sparse
      .cooChecked(1, 3, SparseValuePolicy.RequireFinite)
      .toOption
      .get
    assert(strict.tryAdd(0, 0, Double.NaN).isLeft)
    assert(strict.tryAdd(0, 1, Double.PositiveInfinity).isLeft)
    assert(strict.tryAdd(0, 2, Double.NegativeInfinity).isLeft)
    assertEquals(strict.nnz, 0)
    assert(strict.tryAdd(0, 0, 2.0).isRight)
    intercept[LinAlgError.InvalidArgument](strict.add(0, 1, Double.NaN))
    assertEquals(strict.nnz, 1)

    val legacy = Sparse.coo(1, 1).add(0, 0, Double.NaN).toCOO()
    assert(storedEntries(legacy).head._3.isNaN)
  }

  test("checked duplicate policies sum, retain the stable last value, or return Left") {
    val builder = Sparse.cooChecked(2, 2).toOption.get
    assert(builder.tryAdd(1, 1, 1.0).isRight)
    assert(builder.tryAdd(0, 0, 9.0).isRight)
    assert(builder.tryAdd(1, 1, 7.0).isRight)
    assert(builder.tryAdd(0, 1, 4.0).isRight)
    assert(builder.tryAdd(1, 1, 3.0).isRight)

    assertEquals(builder.tryToCOO(DuplicatePolicy.Sum).toOption.get(1, 1), 11.0)
    assertEquals(builder.tryToCOO(DuplicatePolicy.Last).toOption.get(1, 1), 3.0)
    builder.tryToCOO(DuplicatePolicy.Error) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected duplicate Left, got $other")
    assertEquals(builder.nnz, 5)

    assert(builder.tryCanonicalize(DuplicatePolicy.Error).isLeft)
    assertEquals(builder.nnz, 5, "failed canonicalization must not mutate the builder")
    assertEquals(builder.tryToCOO(DuplicatePolicy.Last).toOption.get(1, 1), 3.0)
  }

  test("checked CSR and CSC finalization is canonical and round-trips by streaming") {
    val result = for
      builder <- Sparse.cooChecked(3, 4, SparseValuePolicy.RequireFinite)
      _ <- builder.tryAdd(2, 3, 5.0)
      _ <- builder.tryAdd(0, 1, 2.0)
      _ <- builder.tryAdd(2, 0, 4.0)
      _ <- builder.tryAdd(0, 1, -2.0)
      _ <- builder.tryAdd(1, 2, 7.0)
      csr <- builder.tryToCSR(DuplicatePolicy.Sum)
    yield csr
    val csr = result.toOption.get
    assert(csr.hasCanonicalFormat)
    assertEquals(storedEntries(csr), Vector((1, 2, 7.0), (2, 0, 4.0), (2, 3, 5.0)))

    val reconstructed = Sparse.cooChecked(csr.rows, csr.cols).toOption.get
    csr.foreachStoredEntry: (row, col, value) =>
      assert(reconstructed.tryAdd(row, col, value).isRight)
    val roundTrip = reconstructed.tryToCSR(DuplicatePolicy.Error).toOption.get
    assert(roundTrip.hasCanonicalFormat)
    assertEquals(storedEntries(roundTrip), storedEntries(csr))

    val csc = reconstructed.tryToCSC(DuplicatePolicy.Error).toOption.get
    assert(csc.hasCanonicalFormat)
    assertEquals(csc.toCSR.toDense().valuesRowMajor, csr.toDense().valuesRowMajor)
  }

  test("checked Error finalization composes without exceptions") {
    val result = for
      builder <- Sparse.cooChecked(2, 2)
      _ <- builder.tryAdd(0, 0, 1.0)
      _ <- builder.tryAdd(0, 0, 2.0)
      csr <- builder.tryToCSR(DuplicatePolicy.Error)
    yield csr
    result match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected duplicate Left, got $other")
  }
