package gale.sparse

import gale.linalg.*

/** Plain regression suite for Matrix Market I/O that stays in core.
  *
  * The property-based sparse matvec law moved to the gale-laws module
  * (`gale.laws.MatrixMarketLawSuite`), driven through [[gale.laws.SparseLaws]].
  * The write/read round-trip below is a concrete regression, so it stays here.
  */
class MatrixMarketLawSuite extends munit.FunSuite:
  test("Matrix Market coordinate round-trips CSR matrices") {
    val A =
      Sparse
        .coo(3, 3)
        .add(0, 0, 1.0)
        .add(1, 2, -2.5)
        .add(2, 1, 4.0)
        .toCSR()

    val roundTrip = MatrixMarket.readCoordinate(MatrixMarket.writeCoordinate(A))

    assertEquals(roundTrip.toDense().valuesRowMajor, A.toDense().valuesRowMajor)
  }
