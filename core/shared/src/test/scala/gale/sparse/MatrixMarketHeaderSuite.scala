package gale.sparse

import gale.linalg.*

class MatrixMarketHeaderSuite extends munit.FunSuite:
  // Item 7: only `matrix coordinate real|integer general` is supported. Every
  // other banner (symmetric storage, pattern/complex fields, dense array format)
  // must be rejected rather than silently mis-parsed as coordinate/real/general.
  test("symmetric symmetry is rejected") {
    val text =
      """%%MatrixMarket matrix coordinate real symmetric
        |2 2 1
        |1 1 3.0
        |""".stripMargin
    intercept[LinAlgError.UnsupportedRepresentation](MatrixMarket.readCoordinate(text))
  }

  test("pattern field is rejected") {
    val text =
      """%%MatrixMarket matrix coordinate pattern general
        |2 2 1
        |1 1
        |""".stripMargin
    intercept[LinAlgError.UnsupportedRepresentation](MatrixMarket.readCoordinate(text))
  }

  test("complex field is rejected") {
    val text =
      """%%MatrixMarket matrix coordinate complex general
        |2 2 1
        |1 1 3.0 0.0
        |""".stripMargin
    intercept[LinAlgError.UnsupportedRepresentation](MatrixMarket.readCoordinate(text))
  }

  test("dense array format is rejected") {
    val text =
      """%%MatrixMarket matrix array real general
        |2 2
        |1.0
        |2.0
        |3.0
        |4.0
        |""".stripMargin
    intercept[LinAlgError.UnsupportedRepresentation](MatrixMarket.readCoordinate(text))
  }

  test("coordinate real general parses, skips comments, and sums duplicates") {
    val text =
      """%%MatrixMarket matrix coordinate real general
        |% a comment line
        |% another comment line
        |2 2 3
        |1 1 1.5
        |1 1 2.5
        |2 2 4.0
        |""".stripMargin
    val csr = MatrixMarket.readCoordinate(text)
    assertEquals(csr(0, 0), 4.0) // 1.5 + 2.5 summed
    assertEquals(csr(1, 1), 4.0)
    assertEquals(csr(0, 1), 0.0)
    assertEquals(csr.nnz, 2)
  }
