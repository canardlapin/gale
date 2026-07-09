package gale.linalg

class FactorizationTypeSuite extends munit.FunSuite:
  test("PivotVector copies input and reports bounds errors") {
    val raw = Array(2, 0, 1)
    val pivots = PivotVector.fromArray(raw)
    raw(0) = 99

    assertEquals(pivots.length, 3)
    assertEquals(pivots(0), 2)
    assertEquals(pivots.toArray.toSeq, Seq(2, 0, 1))
    intercept[LinAlgError.IndexOutOfBounds] {
      pivots(3)
    }
  }

  test("diagnostics expose success and optional rank metadata") {
    val ok = FactorizationDiagnostics(rank = Some(3))
    val failed = FactorizationDiagnostics(info = 2)

    assert(ok.isSuccess)
    assertEquals(ok.rank, Some(3))
    assert(!failed.isSuccess)
  }

  test("dense matrices expose typed decomposition entry points") {
    val A = Matrix.eye(2)

    assert(A.lu.isRight)
    assert(A.cholesky.isRight)
    assertEquals(A.qr.diagnostics.rank, Some(2))
  }
