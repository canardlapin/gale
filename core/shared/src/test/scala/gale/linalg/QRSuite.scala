package gale.linalg

import gale.TestAccess
import gale.platform.PlatformMath.fma

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

  test("column-pivoted QR reconstructs the permuted input and unpermutes coefficients") {
    val A = Matrix.dense(4, 3)(
      1.0e-3, 0.0, 0.0,
      0.0, 10.0, 0.0,
      0.0, 0.0, 1.0,
      0.0, 0.0, 0.0
    )
    val expected = Vec(2.0, -3.0, 4.0)
    val options = QROptions(pivoting = QRPivoting.Column)

    val qr = A.qr(options)
    val permutedInput = Matrix.tabulate(A.rows, A.cols): (row, col) =>
      A(row, qr.columnPermutation(col))
    val actual = qr.solveLeastSquares(A * expected).orThrow

    assert(!qr.columnPermutation.isIdentity)
    assertMatrixClose(qr.q * qr.r, permutedInput, 1e-11)
    assertVectorClose(actual, expected, 1e-10)
    assertEquals(qr.diagnostics.rank, Some(3))
  }

  test("row-first pivot norms exactly match a column-wise reference") {
    val rng = new scala.util.Random(2026080301L)
    val random = Matrix.dense(31, 8, Seq.fill(31 * 8)(rng.nextDouble() * 2.0 - 1.0))
    val tiedAndDeficient = Matrix.dense(7, 6)(
      1.0, 0.0, 0.0, 2.0, 0.0, 1.0 + 4.0e-16,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    val overflowSensitive = Matrix.dense(6, 3)(
      1.0e308, 0.0, 0.0,
      1.0e308, 0.0, 0.0,
      0.0, 1.0e150, 0.0,
      0.0, 0.0, 1.0e-150,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0
    )

    for fixture <- Seq(random, tiedAndDeficient, overflowSensitive) do
      val reference = columnWisePivotedReference(fixture)
      val actual = fixture.qr(QROptions(pivoting = QRPivoting.Column))

      assertEquals(actual.columnPermutation.toIndexSeq, reference.permutation.toIndexedSeq)
      assertEquals(actual.diagnostics.rank, Some(reference.rank))
      assertEquals(actual.r.valuesRowMajor, reference.r.toIndexedSeq)

      val repeated = fixture.qr(QROptions(pivoting = QRPivoting.Column))
      assertEquals(repeated.columnPermutation.toIndexSeq, actual.columnPermutation.toIndexSeq)
      assertEquals(repeated.r.valuesRowMajor, actual.r.valuesRowMajor)
      assertEquals(repeated.diagnostics, actual.diagnostics)
  }

  test("pivoted QR permutation and rank are invariant under safe uniform scaling") {
    val base = Matrix.tabulate(17, 5): (row, col) =>
      math.sin((row + 1).toDouble * (col + 2).toDouble * 0.125) +
        (if row == col then 0.5 else 0.0)
    val expected = base.qr(QROptions(pivoting = QRPivoting.Column))

    for scale <- Seq(1.0e-150, 1.0e150) do
      val scaled = Matrix.tabulate(base.rows, base.cols)((row, col) => scale * base(row, col))
      val actual = scaled.qr(QROptions(pivoting = QRPivoting.Column))
      val reference = columnWisePivotedReference(scaled)

      assertEquals(actual.columnPermutation.toIndexSeq, expected.columnPermutation.toIndexSeq)
      assertEquals(actual.columnPermutation.toIndexSeq, reference.permutation.toIndexedSeq)
      assertEquals(actual.diagnostics.rank, expected.diagnostics.rank)
      assertEquals(actual.diagnostics.rank, Some(reference.rank))
      assertMatrixRelative(actual.q * actual.r, permuteColumns(scaled, actual.columnPermutation), 4e-13, 1e-300)
  }

  test("pivoted QR preserves column-wise non-finite pivot selection") {
    val fixture = Matrix.dense(4, 4)(
      Double.NaN, 3.0, Double.PositiveInfinity, 0.0,
      0.0, 4.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0
    )
    val reference = columnWisePivotedReference(fixture)

    val actual = fixture.qr(QROptions(pivoting = QRPivoting.Column))
    assertEquals(actual.columnPermutation.toIndexSeq, reference.permutation.toIndexedSeq)
    assertEquals(actual.diagnostics.rank, Some(reference.rank))
  }

  test("caller rank tolerance controls the explicit rank decision") {
    val A = Matrix.dense(3, 2)(
      1.0, 0.0,
      0.0, 1.0e-8,
      0.0, 0.0
    )

    val automatic = A.qr(QROptions(pivoting = QRPivoting.Column))
    val thresholded = A.qr(
      QROptions(
        pivoting = QRPivoting.Column,
        rankTolerance = Some(1.0e-6)
      )
    )

    assertEquals(automatic.diagnostics.rank, Some(2))
    assertEquals(thresholded.diagnostics.rank, Some(1))
    assertEquals(thresholded.diagnostics.rankTolerance, Some(1.0e-6))
    assert(thresholded.solveLeastSquares(Vec(1.0, 1.0e-8, 0.0)).left.exists(_.isInstanceOf[LinAlgError.RankDeficient]))
  }

  test("pivoted QR solves matrix right-hand sides and exposes Q applications") {
    val A = Matrix.dense(5, 3)(
      1.0, 0.0, 2.0,
      1.0, 1.0, -1.0,
      1.0, 2.0, 0.5,
      1.0, 3.0, 1.5,
      1.0, 4.0, -0.5
    )
    val expected = Matrix.dense(3, 2)(
      1.0, -2.0,
      0.5, 3.0,
      -1.0, 4.0
    )
    val observations = A * expected
    val qr = A.qr(QROptions(pivoting = QRPivoting.Column))

    val actual = qr.solveLeastSquares(observations).orThrow
    val transformed = qr.applyQT(observations).orThrow

    assertMatrixClose(actual, expected, 1e-10)
    assertMatrixClose(transformed, qr.q.t * observations, 1e-10)
    assertMatrixClose(qr.applyQ(transformed).orThrow, observations, 1e-10)
  }

  test("matrix RHS Q applications and least squares reduce to independent columns") {
    val rng = new scala.util.Random(2026080201L)
    val m = 37
    val p = 6
    val design = Matrix.dense(m, p, Seq.fill(m * p)(rng.nextDouble() * 2.0 - 1.0))
    val qr = design.qr(QROptions(pivoting = QRPivoting.Column))

    for q <- Seq(1, 8, 17) do
      val rhs = Matrix.dense(m, q, Seq.fill(m * q)(rng.nextDouble() * 2.0 - 1.0))
      val appliedQt = qr.applyQT(rhs).orThrow
      val reappliedQ = qr.applyQ(appliedQt).orThrow
      val solved = qr.solveLeastSquares(rhs).orThrow
      val appliedQtColumns = Seq.tabulate(q): col =>
        qr.applyQT(Matrix.tabulate(m, 1)((row, _) => rhs(row, col))).orThrow
      val appliedQColumns = Seq.tabulate(q): col =>
        qr.applyQ(Matrix.tabulate(m, 1)((row, _) => appliedQt(row, col))).orThrow
      val appliedQtByColumn = Matrix.tabulate(m, q)((row, col) => appliedQtColumns(col)(row, 0))
      val appliedQByColumn = Matrix.tabulate(m, q)((row, col) => appliedQColumns(col)(row, 0))
      val solvedByColumn = Matrix.tabulate(p, q): (row, col) =>
        qr.solveLeastSquares(rhs.col(col)).orThrow(row)

      assertMatrixRelative(appliedQt, qr.q.t * rhs, rel = 8e-13, abs = 8e-14)
      assertMatrixRelative(appliedQt, appliedQtByColumn, rel = 8e-13, abs = 8e-14)
      assertMatrixRelative(reappliedQ, rhs, rel = 2e-12, abs = 2e-13)
      assertMatrixRelative(reappliedQ, appliedQByColumn, rel = 8e-13, abs = 8e-14)
      assertMatrixRelative(solved, solvedByColumn, rel = 8e-13, abs = 8e-14)

      val appliedAgain = qr.applyQT(rhs).orThrow
      val solvedAgain = qr.solveLeastSquares(rhs).orThrow
      assertEquals(appliedAgain.valuesRowMajor, appliedQt.valuesRowMajor)
      assertEquals(solvedAgain.valuesRowMajor, solved.valuesRowMajor)
      assert(!TestAccess.sameStorage(TestAccess.dmatStorage(rhs), TestAccess.dmatStorage(appliedQt)))
      assert(!TestAccess.sameStorage(TestAccess.dmatStorage(rhs), TestAccess.dmatStorage(solved)))
  }

  test("matrix RHS paths accept a non-contiguous transpose without changing its storage") {
    val m = 41
    val p = 6
    val q = 9
    val design = Matrix.tabulate(m, p)((row, col) =>
      math.sin((row + 1).toDouble * (col + 2).toDouble * 0.03125) + (if row == col then 2.0 else 0.0)
    )
    val backing = Matrix.tabulate(q, m)((row, col) => math.cos((row + 3).toDouble * (col + 1).toDouble * 0.015625))
    val strided = backing.t
    val owned = Matrix.tabulate(m, q)(strided.apply)
    val before = backing.valuesRowMajor
    val qr = design.qr(QROptions(pivoting = QRPivoting.Column))

    assert(!strided.isContiguousRowMajor)
    assertMatrixRelative(qr.applyQT(strided).orThrow, qr.applyQT(owned).orThrow, rel = 0.0, abs = 0.0)
    assertMatrixRelative(
      qr.solveLeastSquares(strided).orThrow,
      qr.solveLeastSquares(owned).orThrow,
      rel = 0.0,
      abs = 0.0
    )
    assertEquals(backing.valuesRowMajor, before)
  }

  test("matrix RHS paths preserve rank decisions at exact and near deficiency") {
    val exact = Matrix.dense(5, 3)(
      1.0, 2.0, 3.0,
      2.0, 4.0, 6.0,
      3.0, 6.0, 9.0,
      4.0, 8.0, 12.0,
      5.0, 10.0, 15.0
    ).qr(QROptions(pivoting = QRPivoting.Column))
    val rhs = Matrix.tabulate(5, 4)((row, col) => row.toDouble - 0.5 * col)
    assertEquals(exact.diagnostics.rank, Some(1))
    assert(exact.solveLeastSquares(rhs).left.exists(_.isInstanceOf[LinAlgError.RankDeficient]))

    val near = Matrix.dense(5, 3)(
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0e-10,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0
    ).qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(1.0e-8)))
    assertEquals(near.diagnostics.rank, Some(2))
    assert(near.solveLeastSquares(rhs).left.exists(_.isInstanceOf[LinAlgError.RankDeficient]))
    assertMatrixRelative(near.applyQ(near.applyQT(rhs).orThrow).orThrow, rhs, rel = 1e-15, abs = 1e-15)
  }

  test("matrix RHS least squares remains finite and scale-equivariant at extreme magnitudes") {
    val base = Matrix.tabulate(19, 5)((row, col) =>
      math.sin((row + 1).toDouble * (col + 1).toDouble) + (if row == col then 1.0 else 0.0)
    )
    val expected = Matrix.tabulate(5, 7)((row, col) => (row - 2).toDouble * 0.125 + col.toDouble * 0.03125)

    for scale <- Seq(1.0e150, 1.0e-150) do
      val design = Matrix.tabulate(base.rows, base.cols)((row, col) => scale * base(row, col))
      val observations = design * expected
      val actual = design.qr(QROptions(pivoting = QRPivoting.Column)).solveLeastSquares(observations).orThrow

      assert(actual.valuesRowMajor.forall(_.isFinite), s"non-finite matrix solution at scale=$scale")
      assertMatrixRelative(actual, expected, rel = 4e-13, abs = 4e-14)
  }

  test("QR residualization and normalized covariance satisfy independent identities") {
    val A = Matrix.dense(5, 2)(
      1.0, 0.0,
      1.0, 1.0,
      1.0, 2.0,
      1.0, 3.0,
      1.0, 4.0
    )
    val data = Matrix.dense(5, 2)(
      2.0, 1.0,
      -1.0, 3.0,
      0.5, -2.0,
      4.0, 0.0,
      1.0, 5.0
    )
    val qr = A.qr(QROptions(pivoting = QRPivoting.Column))

    val residual = qr.residualize(data).orThrow
    val covariance = qr.normalizedCovariance.orThrow

    assertMatrixClose(A.t * residual, Matrix.zeros(A.cols, data.cols), 1e-10)
    assertMatrixClose((A.t * A) * covariance, Matrix.eye(A.cols), 1e-10)
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

  private final case class ReferencePivotedQR(
      permutation: Array[Int],
      r: Array[Double],
      rank: Int
  )

  /** Deliberately slow test-only reference for the pre-optimization contract:
    * recompute every candidate norm in its own strided column scan, then rescan
    * the selected column while constructing its reflector.
    */
  private def columnWisePivotedReference(A: DMat): ReferencePivotedQR =
    val m = A.rows
    val n = A.cols
    val limit = math.min(m, n)
    val r = A.valuesRowMajor.toArray
    val reflectors = Array.fill(m * limit)(0.0)
    val tau = Array.fill(limit)(0.0)
    val permutation = Array.tabulate(n)(identity)
    val dots = Array.fill(n)(0.0)

    var k = 0
    while k < limit do
      var pivot = k
      var bestNorm = -1.0
      var col = k
      while col < n do
        val norm = referenceNorm(r, m, n, k, col)
        if norm > bestNorm then
          bestNorm = norm
          pivot = col
        col += 1

      if pivot != k then
        var row = 0
        while row < m do
          val left = row * n + k
          val right = row * n + pivot
          val value = r(left)
          r(left) = r(right)
          r(right) = value
          row += 1
        val original = permutation(k)
        permutation(k) = permutation(pivot)
        permutation(pivot) = original

      val diagIndex = k * n + k
      val x0 = r(diagIndex)
      val norm = referenceNorm(r, m, n, k, k)
      reflectors(k * limit + k) = 1.0
      if norm > 0.0 then
        val beta = if x0 >= 0.0 then -norm else norm
        val x0OverBeta = x0 / beta
        val denominatorOverBeta = x0OverBeta - 1.0
        tau(k) = 1.0 - x0OverBeta
        r(diagIndex) = beta
        var row = k + 1
        while row < m do
          val index = row * n + k
          reflectors(row * limit + k) = (r(index) / beta) / denominatorOverBeta
          r(index) = 0.0
          row += 1
      else
        tau(k) = 0.0
        var row = k + 1
        while row < m do
          r(row * n + k) = 0.0
          row += 1

      col = k + 1
      while col < n do
        dots(col) = 0.0
        col += 1
      var row = k
      while row < m do
        val vi = reflectors(row * limit + k)
        col = k + 1
        while col < n do
          dots(col) = fma(vi, r(row * n + col), dots(col))
          col += 1
        row += 1
      col = k + 1
      while col < n do
        dots(col) *= tau(k)
        col += 1
      row = k
      while row < m do
        val vi = reflectors(row * limit + k)
        col = k + 1
        while col < n do
          val index = row * n + col
          r(index) = fma(-vi, dots(col), r(index))
          col += 1
        row += 1
      k += 1

    var maxDiag = 0.0
    k = 0
    while k < limit do
      maxDiag = math.max(maxDiag, math.abs(r(k * n + k)))
      k += 1
    val tolerance = 2.0 * math.max(m, n).toDouble * 2.220446049250313e-16 * maxDiag
    var rank = 0
    k = 0
    while k < limit do
      if math.abs(r(k * n + k)) > tolerance then rank += 1
      k += 1
    ReferencePivotedQR(permutation, r, rank)

  private def referenceNorm(values: Array[Double], rows: Int, cols: Int, fromRow: Int, col: Int): Double =
    var scale = 0.0
    var ssq = 1.0
    var row = fromRow
    while row < rows do
      val value = values(row * cols + col)
      if value != 0.0 then
        val abs = math.abs(value)
        if scale < abs then
          val ratio = scale / abs
          ssq = 1.0 + ssq * ratio * ratio
          scale = abs
        else
          val ratio = abs / scale
          ssq += ratio * ratio
      row += 1
    scale * math.sqrt(ssq)

  private def permuteColumns(A: DMat, permutation: ColumnPermutation): DMat =
    Matrix.tabulate(A.rows, A.cols)((row, col) => A(row, permutation(col)))

  private def assertVectorClose(actual: DVec, expected: DVec, tolerance: Double): Unit =
    assertEquals(actual.length, expected.length)
    var i = 0
    while i < actual.length do
      assert(math.abs(actual(i) - expected(i)) <= tolerance, s"index $i: ${actual(i)} != ${expected(i)}")
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
