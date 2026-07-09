package gale.linalg

class QRSuite extends munit.FunSuite:
  test("QR reconstructs a tall dense matrix") {
    val A = Matrix.dense(3, 2)(
      1.0, 0.0,
      1.0, 1.0,
      1.0, 2.0
    )

    val qr = A.qr
    val reconstructed = qr.q * qr.r
    val qtq = qr.q.t * qr.q

    assertMatrixClose(reconstructed, A, 1e-10)
    assertMatrixClose(qtq, Matrix.eye(3), 1e-10)
    assertEquals(qr.diagnostics.rank, Some(2))
  }

  test("QR least-squares solve recovers full-rank coefficients") {
    val A = Matrix.dense(3, 2)(
      1.0, 0.0,
      1.0, 1.0,
      1.0, 2.0
    )
    val b = Vec(1.0, 3.0, 5.0)

    val x = A.leastSquares(b).orThrow
    val r = A * x - b

    assert(norm(r) < 1e-10)
    assert(math.abs(x(0) - 1.0) < 1e-10)
    assert(math.abs(x(1) - 2.0) < 1e-10)
  }

  test("QR least-squares reports rank deficiency") {
    val A = Matrix.dense(3, 2)(
      1.0, 2.0,
      2.0, 4.0,
      3.0, 6.0
    )
    val b = Vec(1.0, 2.0, 3.0)

    assertEquals(A.leastSquares(b), Left(LinAlgError.RankDeficient(1, 2)))
  }

  // K5: QR keeps compact reflectors + tau (shape m x min(m, n)); the dense Q is
  // rebuilt on demand and still satisfies the reconstruction and orthogonality
  // laws q * r == A and qᵀq == I.
  test("QR stores compact reflectors and materialises an orthogonal Q on demand") {
    val A = Matrix.dense(4, 2)(
      1.0, 0.0,
      1.0, 1.0,
      1.0, 2.0,
      1.0, 3.0
    )

    val qr = A.qr
    assertEquals(qr.reflectors.rows, 4)
    assertEquals(qr.reflectors.cols, 2) // min(m, n)
    assertEquals(qr.tau.length, 2)

    assertMatrixClose(qr.q * qr.r, A, 1e-10)
    assertMatrixClose(qr.q.t * qr.q, Matrix.eye(4), 1e-10)
    assertEquals(qr.diagnostics.rank, Some(2))
  }

  // K5: solveLeastSquares applies Qᵀ implicitly through the reflectors (never
  // forming Q). The answer must match the normal-equations solution on a random
  // full-rank tall system.
  test("QR least-squares via implicit Qᵀ matches the normal-equations solution") {
    val rng = new scala.util.Random(20260709L)
    val m = 14
    val n = 5
    val A = Matrix.dense(m, n, Seq.fill(m * n)(rng.nextDouble() * 2.0 - 1.0))
    val b = DVec.fromSeq(Seq.fill(m)(rng.nextDouble() * 2.0 - 1.0))

    val x = A.leastSquares(b).orThrow

    // Reference: solve the normal equations AᵀA x = Aᵀb directly.
    val ref = (A.t * A).solve(A.t * b).orThrow
    var i = 0
    while i < n do
      assert(math.abs(x(i) - ref(i)) < 1e-8, s"coefficient $i: ${x(i)} != ${ref(i)}")
      i += 1

    // The residual is orthogonal to the column space at the optimum.
    val residual = A * x - b
    assert((A.t * residual).norm2 < 1e-8)
  }

  private def assertMatrixClose(actual: DMat, expected: DMat, tolerance: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        assert(math.abs(actual(i, j) - expected(i, j)) < tolerance)
        j += 1
      i += 1
