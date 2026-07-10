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

    val qr = A.qr
    assertEquals(
      qr.solveLeastSquares(b),
      Left(LinAlgError.RankDeficient(1, 2)),
      clues(qr.r.valuesRowMajor, DenseDecompositions.rankToleranceFromMatrix(qr.r), qr.diagnostics)
    )
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

  test("blocked QR reconstructs and stays orthogonal across a partial final panel") {
    val rng = new scala.util.Random(2026071001L)
    val m = 129
    val n = 97
    val A = Matrix.dense(m, n, Seq.fill(m * n)(rng.nextDouble() * 2.0 - 1.0))
    val workspace = DenseWorkspace.empty

    val qr = A.qrWith(workspace)
    val reconstructed = qr.q * qr.r
    val qtq = qr.q.t * qr.q

    assert(workspace.workCapacity > m, "large QR did not acquire reusable block workspace")
    assertMatrixRelative(reconstructed, A, rel = 2e-12, abs = 2e-13)
    assertMatrixRelative(qtq, Matrix.eye(m), rel = 2e-12, abs = 2e-12)
    assertEquals(qr.diagnostics.rank, Some(n))
  }

  test("blocked least squares recovers a known solution") {
    val rng = new scala.util.Random(2026071002L)
    val m = 192
    val n = 96
    val A = Matrix.dense(m, n, Seq.fill(m * n)(rng.nextDouble() * 2.0 - 1.0))
    val expected = DVec.fromSeq(Seq.fill(n)(rng.nextDouble() * 2.0 - 1.0))
    val b = A * expected

    val actual = A.leastSquares(b).orThrow
    var i = 0
    while i < n do
      assert(math.abs(actual(i) - expected(i)) <= 2e-11, s"coefficient $i: ${actual(i)} != ${expected(i)}")
      i += 1
    assert((A * actual - b).norm2 <= 2e-10)
  }

  test("blocked QR reconstructs a wide matrix") {
    val rng = new scala.util.Random(2026071004L)
    val m = 97
    val n = 129
    val A = Matrix.dense(m, n, Seq.fill(m * n)(rng.nextDouble() * 2.0 - 1.0))
    val qr = A.qr

    assertMatrixRelative(qr.q * qr.r, A, rel = 3e-12, abs = 3e-13)
    assertMatrixRelative(qr.q.t * qr.q, Matrix.eye(m), rel = 3e-12, abs = 3e-12)
    assertEquals(qr.diagnostics.rank, Some(m))
  }

  test("blocked QR detects an exactly dependent column in the final panel") {
    val rng = new scala.util.Random(2026071005L)
    val m = 128
    val n = 96
    val independent = Array.fill(m * (n - 1))(rng.nextDouble() * 2.0 - 1.0)
    val A = Matrix.tabulate(m, n): (i, j) =>
      if j == n - 1 then 2.0 * independent(i * (n - 1))
      else independent(i * (n - 1) + j)

    assertEquals(A.qr.diagnostics.rank, Some(n - 1))
  }

  test("Householder construction remains finite at extreme representable scales") {
    val base = Matrix.dense(5, 3)(
      1.0, -2.0, 0.5,
      -3.0, 4.0, 1.5,
      2.0, 1.0, -1.0,
      0.25, -0.75, 2.0,
      -1.5, 0.5, 3.0
    )

    for scale <- Seq(1e155, 1e-170) do
      val A = Matrix.tabulate(base.rows, base.cols)((i, j) => scale * base(i, j))
      val qr = A.qr
      val reconstructed = qr.q * qr.r

      assert(qr.q.valuesRowMajor.forall(_.isFinite), s"non-finite Q at scale=$scale")
      assert(qr.r.valuesRowMajor.forall(_.isFinite), s"non-finite R at scale=$scale")
      assertMatrixRelative(reconstructed, A, rel = 3e-14, abs = 3e-14 * scale)
  }

  test("Householder normalization does not overflow when x0 minus beta would") {
    val A = Matrix.dense(3, 2)(
      1.0e308, 1.0,
      1.0e307, 2.0,
      -1.0e307, 3.0
    )
    val qr = A.qr

    assert(qr.q.valuesRowMajor.forall(_.isFinite))
    assert(qr.r.valuesRowMajor.forall(_.isFinite))
    assertMatrixRelative(qr.q * qr.r, A, rel = 5e-15, abs = 5e-14)
    assertMatrixRelative(qr.q.t * qr.q, Matrix.eye(3), rel = 5e-14, abs = 5e-14)
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

  private def assertMatrixRelative(actual: DMat, expected: DMat, rel: Double, abs: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        val a = actual(i, j)
        val e = expected(i, j)
        val tolerance = abs + rel * math.abs(e)
        assert(a.isFinite, s"non-finite [$i,$j]: $a")
        assert(math.abs(a - e) <= tolerance, s"[$i,$j] $a != $e (tol=$tolerance)")
        j += 1
      i += 1
