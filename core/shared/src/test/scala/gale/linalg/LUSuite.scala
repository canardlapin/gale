package gale.linalg

class LUSuite extends munit.FunSuite:
  test("LU solves a dense square system with small residual") {
    val A = Matrix.dense(3, 3)(
      1.0, 2.0, 3.0,
      4.0, 5.0, 6.0,
      7.0, 8.0, 10.0
    )
    val b = Vec(3.0, 3.0, 4.0)

    val x = A.solve(b).orThrow
    val r = A * x - b

    assert(norm(r) < 1e-10)
  }

  test("LU uses partial pivoting") {
    val A = Matrix.dense(2, 2)(
      0.0, 1.0,
      2.0, 3.0
    )
    val b = Vec(1.0, 5.0)

    val lu = A.lu.orThrow
    val x = lu.solve(b).orThrow

    assertEquals(lu.pivots.toArray.toSeq, Seq(1, 0))
    assert(math.abs(x(0) - 1.0) < 1e-12)
    assert(math.abs(x(1) - 1.0) < 1e-12)
  }

  test("LU determinant comes from packed factors and pivot parity") {
    val A = Matrix.dense(2, 2)(
      1.0, 2.0,
      3.0, 4.0
    )

    assert(math.abs(A.lu.orThrow.det.orThrow + 2.0) < 1e-12)
  }

  test("singular and non-square matrices return diagnostics") {
    val singular = Matrix.dense(2, 2)(
      1.0, 2.0,
      2.0, 4.0
    )
    val rectangular = Matrix.zeros(2, 3)

    assertEquals(singular.lu, Left(LinAlgError.SingularMatrix(1)))
    assert(rectangular.lu.left.exists(_.isInstanceOf[LinAlgError.NonSquareMatrix]))
  }

  test("blocked matrix solve agrees with per-column vector solves and accepts a transposed RHS") {
    // n = 100 spans many triangular blocks, and 7 RHS columns cover the
    // four-column groups plus a remainder.
    val (n, k) = (100, 7)
    val rng = new scala.util.Random(41)
    val a = Matrix.tabulate(n, n)((i, j) => if i == j then n.toDouble else rng.nextDouble() * 2.0 - 1.0)
    val b = Matrix.tabulate(n, k)((_, _) => rng.nextDouble() * 2.0 - 1.0)
    val transposedB = Matrix.tabulate(k, n)((c, i) => b(i, c)).t
    val lu = a.lu.orThrow
    val x = lu.solve(b).orThrow
    for c <- 0 until k do
      val column = lu.solve(Vec.tabulate(n)(i => b(i, c))).orThrow
      for i <- 0 until n do
        assert(math.abs(x(i, c) - column(i)) <= 1e-13 * (1.0 + math.abs(column(i))), s"x($i, $c)")
    val fromTransposed = lu.solve(transposedB).orThrow
    for i <- 0 until n; c <- 0 until k do
      assertEquals(fromTransposed(i, c), x(i, c))
      assertEquals(transposedB(i, c), b(i, c))
  }

  test("inverse matches a hand-computed 2x2 and agrees with solve(I) under row pivoting") {
    val small = Matrix.dense(2, 2)(4.0, 7.0, 2.0, 6.0).inverse.orThrow
    val expected = Seq(0.6, -0.7, -0.2, 0.4)
    for i <- 0 until 2; j <- 0 until 2 do assertEqualsDouble(small(i, j), expected(i * 2 + j), 1e-15)

    // Random entries force genuine row swaps; n = 100 spans several column blocks.
    for n <- Seq(1, 3, 17, 100) do
      val rng = new scala.util.Random(53 + n)
      val a = Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)
      val before = a.valuesRowMajor
      val lu = a.lu.orThrow
      if n >= 3 then assert((0 until n).exists(i => lu.pivots(i) != i), s"n=$n should pivot")
      val x = a.inverse.orThrow
      val viaSolve = a.solve(Matrix.eye(n)).orThrow
      // Backward-stable bounds: residual ‖AX − I‖ ≲ c·n·ε·‖A‖·‖X‖, and the two
      // computed inverses agree to that relative level times the condition number.
      val bound = 32.0 * n * 2.220446049250313e-16 * a.normInf * x.normInf
      val residual = a * x - Matrix.eye(n)
      assert(residual.normInf <= bound, s"n=$n residual ${residual.normInf} > $bound")
      val diff = (x - viaSolve).normInf
      assert(diff <= bound * x.normInf, s"n=$n inverse vs solve(I): $diff")
      assertEquals(a.valuesRowMajor, before)
      assertEquals(lu.inverse.orThrow.valuesRowMajor, x.valuesRowMajor)
  }

  test("inverse reports singular, non-square, and empty inputs") {
    val singular = Matrix.dense(2, 2)(1.0, 2.0, 2.0, 4.0)
    assertEquals(singular.inverse, Left(LinAlgError.SingularMatrix(1)))
    assert(Matrix.zeros(2, 3).inverse.left.exists(_.isInstanceOf[LinAlgError.NonSquareMatrix]))
    val empty = Matrix.zeros(0, 0).inverse.orThrow
    assertEquals((empty.rows, empty.cols), (0, 0))
  }

  /** The one-column-at-a-time right-looking elimination the paired-column LU
    * must reproduce exactly: `(packed row-major, pivots, parity)` or the
    * singular column.
    */
  private def sequentialLu(a: DMat): Either[Int, (Seq[Double], Seq[Int], Int)] =
    val n = a.rows
    val p = Array.tabulate(n * n)(i => a(i / n, i % n))
    val pivots = Array.tabulate(n)(identity)
    var parity = 1
    var singular = -1
    var k = 0
    while singular < 0 && k < n do
      var pivot = k
      var maxAbs = math.abs(p(k * n + k))
      for i <- k + 1 until n do
        if math.abs(p(i * n + k)) > maxAbs then
          maxAbs = math.abs(p(i * n + k))
          pivot = i
      if maxAbs == 0.0 || maxAbs.isNaN then singular = k
      else
        if pivot != k then
          for j <- 0 until n do
            val t = p(k * n + j); p(k * n + j) = p(pivot * n + j); p(pivot * n + j) = t
          val t = pivots(k); pivots(k) = pivots(pivot); pivots(pivot) = t
          parity = -parity
        for i <- k + 1 until n do
          p(i * n + k) = p(i * n + k) / p(k * n + k)
          for j <- k + 1 until n do p(i * n + j) = p(i * n + j) - p(i * n + k) * p(k * n + j)
        k += 1
    if singular >= 0 then Left(singular) else Right((p.toSeq, pivots.toSeq, parity))

  test("paired-column LU reproduces sequential elimination bit for bit, including singular columns") {
    def bits(xs: Seq[Double]) = xs.map(java.lang.Double.doubleToRawLongBits)
    for n <- Seq(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 37, 64); seed <- Seq(1, 2) do
      val rng = new scala.util.Random(1000 * n + seed)
      val a = Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)
      val lu = a.lu.orThrow
      val (packed, pivots, parity) = sequentialLu(a).toOption.get
      assertEquals(bits(lu.packed.valuesRowMajor), bits(packed), s"n=$n seed=$seed")
      assertEquals(lu.pivots.toIndexSeq, pivots.toIndexedSeq, s"n=$n seed=$seed")
      assertEquals(lu.parity, parity, s"n=$n seed=$seed")
      // A zero column stays exactly zero under elimination, so the factorization
      // stops there: an even and an odd paired column...
      for c <- Seq(0, 1, 2, n - 1).distinct if c < n do
        val singular = Matrix.tabulate(n, n)((i, j) => if j == c then 0.0 else a(i, j))
        assertEquals(singular.lu, Left(LinAlgError.SingularMatrix(c)), s"n=$n zero column $c")
        assertEquals(sequentialLu(singular).left.toOption, Some(c))
      // ...and a NaN diagonal is reported at its column.
      if n >= 3 then
        val poisoned = Matrix.tabulate(n, n)((i, j) => if i >= 1 && j == 1 then Double.NaN else a(i, j))
        assertEquals(poisoned.lu.left.toOption, sequentialLu(poisoned).left.toOption.map(LinAlgError.SingularMatrix(_)))
  }

  test("paired-column LU keeps the first-row tie-break on tied pivot magnitudes") {
    // ±1 entries tie every pivot candidate in the first column and many later
    // ones; the first row at or below k with the largest magnitude must win.
    def bits(xs: Seq[Double]) = xs.map(java.lang.Double.doubleToRawLongBits)
    for n <- Seq(2, 3, 4, 5, 8, 9, 16, 33); seed <- 1 to 4 do
      val rng = new scala.util.Random(77 * n + seed)
      val a = Matrix.tabulate(n, n)((_, _) => if rng.nextBoolean() then 1.0 else -1.0)
      (a.lu, sequentialLu(a)) match
        case (Right(lu), Right((packed, pivots, parity))) =>
          assertEquals(bits(lu.packed.valuesRowMajor), bits(packed), s"n=$n seed=$seed")
          assertEquals(lu.pivots.toIndexSeq, pivots.toIndexedSeq, s"n=$n seed=$seed")
          assertEquals(lu.parity, parity, s"n=$n seed=$seed")
        case (Left(error), Left(k)) => assertEquals(error, LinAlgError.SingularMatrix(k), s"n=$n seed=$seed")
        case (got, want)            => fail(s"n=$n seed=$seed: $got vs $want")
  }

  test("inverse of an ill-conditioned Hilbert matrix meets the backward bound and agrees with solve(I)") {
    val eps = 2.220446049250313e-16
    for n <- 8 to 10 do
      val a = Matrix.tabulate(n, n)((i, j) => 1.0 / (i + j + 1).toDouble)
      val x = a.inverse.orThrow
      val viaSolve = a.solve(Matrix.eye(n)).orThrow
      // ‖A‖·‖X‖ is the condition estimate (≈ 1e10 to 1e13 here).
      val kappaScale = a.normInf * x.normInf
      assert(kappaScale > 1e9, s"n=$n should be ill-conditioned: $kappaScale")
      val residual = (a * x - Matrix.eye(n)).normInf
      assert(residual <= 32.0 * n * eps * kappaScale, s"n=$n residual $residual vs κ-scale $kappaScale")
      val diff = (x - viaSolve).normInf
      assert(diff <= 32.0 * n * eps * kappaScale * x.normInf, s"n=$n inverse vs solve(I): $diff")
  }

  test("matrix solve columns equal vector solves bit for bit within one triangular block") {
    for n <- 1 to 8; k <- Seq(1, 3, 4, 5, 9) do
      val rng = new scala.util.Random(13 * n + k)
      val a = Matrix.tabulate(n, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)
      val b = Matrix.tabulate(n, k)((_, _) => rng.nextDouble() * 2.0 - 1.0)
      val lu = a.lu.orThrow
      val x = lu.solve(b).orThrow
      for c <- 0 until k do
        val column = lu.solve(Vec.tabulate(n)(i => b(i, c))).orThrow
        for i <- 0 until n do assertEquals(x(i, c), column(i), s"n=$n k=$k ($i, $c)")
  }
