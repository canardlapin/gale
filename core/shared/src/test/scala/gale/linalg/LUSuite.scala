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
