package gale.linalg

/** Tests for the `DMat.pinv` facade (Moore–Penrose pseudo-inverse via the
  * economy SVD, NumPy-default cutoff) and the total `DMat.kron` Kronecker
  * product. `pinv` is checked through the four Moore–Penrose conditions —
  * `A A⁺ A = A`, `A⁺ A A⁺ = A⁺`, `(A A⁺)ᵀ = A A⁺`, `(A⁺ A)ᵀ = A⁺ A` — on
  * full-rank and rank-deficient fixtures; `kron` through known values, shape
  * algebra, the mixed-product identity, and strided (transpose-view) operands.
  */
class PinvKronSuite extends munit.FunSuite:

  // --- helpers ---------------------------------------------------------------

  private def randomMat(m: Int, n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(m, n)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  private def frob(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  private def asymmetry(a: DMat): Double =
    frob(a - a.t)

  /** Assert the four Moore–Penrose conditions at `tol` (scaled by ‖A‖). */
  private def assertMoorePenrose(a: DMat, p: DMat, tol: Double, clue: String): Unit =
    assertEquals(p.rows, a.cols, s"$clue: pinv rows")
    assertEquals(p.cols, a.rows, s"$clue: pinv cols")
    val scale = math.max(1.0, frob(a))
    assert(frob(a * p * a - a) < tol * scale, s"$clue: A A+ A != A (${frob(a * p * a - a)})")
    assert(frob(p * a * p - p) < tol * math.max(1.0, frob(p)), s"$clue: A+ A A+ != A+ (${frob(p * a * p - p)})")
    assert(asymmetry(a * p) < tol * scale, s"$clue: A A+ not symmetric (${asymmetry(a * p)})")
    assert(asymmetry(p * a) < tol * scale, s"$clue: A+ A not symmetric (${asymmetry(p * a)})")

  // --- pinv: full rank -------------------------------------------------------

  test("pinv full-rank tall, wide, and square satisfy the Moore–Penrose conditions") {
    for (m, n, seed) <- Seq((8, 4, 11L), (4, 9, 12L), (7, 7, 13L)) do
      val a = randomMat(m, n, seed)
      val p = a.pinv.orThrow
      assertMoorePenrose(a, p, 1e-10, s"${m}x$n seed=$seed")
  }

  test("pinv of a full-rank tall matrix is the least-squares inverse: A+ A = I") {
    val a = randomMat(10, 4, 21L)
    val p = a.pinv.orThrow
    assert(frob(p * a - Matrix.eye(4)) < 1e-10, s"A+ A != I: ${frob(p * a - Matrix.eye(4))}")
  }

  test("pinv of a well-conditioned square matrix is its inverse") {
    // Diagonally dominant, hence comfortably invertible.
    val n = 6
    val base = randomMat(n, n, 31L)
    val a = Matrix.tabulate(n, n)((i, j) => if i == j then base(i, j) + 4.0 else base(i, j))
    val p = a.pinv.orThrow
    assert(frob(a * p - Matrix.eye(n)) < 1e-10)
    assert(frob(p * a - Matrix.eye(n)) < 1e-10)
  }

  // --- pinv: rank deficiency -------------------------------------------------

  test("pinv rank-deficient satisfies the Moore–Penrose conditions") {
    // rank 2 by construction (6x5 from two outer products).
    val x1 = (0 until 6).map(i => 1.0 + 0.3 * i)
    val y1 = (0 until 5).map(j => 2.0 - 0.4 * j)
    val x2 = (0 until 6).map(i => math.sin(i + 1.0))
    val y2 = (0 until 5).map(j => math.cos(j + 0.5))
    val a = Matrix.tabulate(6, 5)((i, j) => 3.0 * x1(i) * y1(j) + x2(i) * y2(j))
    val p = a.pinv.orThrow
    assertMoorePenrose(a, p, 1e-9, "rank-2 6x5")
  }

  test("pinv of a rank-1 outer product is the closed form v uT / (|u|^2 |v|^2)") {
    val u = Vec(1.0, 2.0, 2.0) // |u|^2 = 9
    val v = Vec(3.0, 4.0)      // |v|^2 = 25
    val a = Matrix.tabulate(3, 2)((i, j) => u(i) * v(j))
    val p = a.pinv.orThrow
    val expected = Matrix.tabulate(2, 3)((i, j) => v(i) * u(j) / (9.0 * 25.0))
    assert(frob(p - expected) < 1e-12, s"closed form mismatch: ${frob(p - expected)}")
  }

  test("pinv of the zero matrix is the zero matrix, and 1x1 edges hold") {
    val z = DMat.zeros(3, 5).pinv.orThrow
    assertEquals(z.rows, 5)
    assertEquals(z.cols, 3)
    assert(frob(z) == 0.0)
    val inv = Matrix.dense(1, 1)(2.0).pinv.orThrow
    assert(math.abs(inv(0, 0) - 0.5) < 1e-15)
    val zero1 = Matrix.dense(1, 1)(0.0).pinv.orThrow
    assertEquals(zero1(0, 0), 0.0)
  }

  test("pinv structural violation: empty dimension is Left") {
    assert(DMat.zeros(0, 4).pinv.isLeft)
  }

  // --- kron ------------------------------------------------------------------

  test("kron known values: 2x2 blocks") {
    val a = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
    val b = Matrix.dense(2, 2)(0.0, 5.0, 6.0, 7.0)
    val k = a.kron(b)
    assertEquals(k.rows, 4)
    assertEquals(k.cols, 4)
    val expected = Matrix.dense(4, 4)(
      0.0, 5.0, 0.0, 10.0,
      6.0, 7.0, 12.0, 14.0,
      0.0, 15.0, 0.0, 20.0,
      18.0, 21.0, 24.0, 28.0
    )
    assert(frob(k - expected) == 0.0, s"kron mismatch")
  }

  test("kron shape algebra and identity blocks") {
    val a = randomMat(2, 3, 41L)
    val b = randomMat(4, 5, 42L)
    val k = a.kron(b)
    assertEquals(k.rows, 8)
    assertEquals(k.cols, 15)
    // I2 (x) B is block-diagonal with B blocks.
    val ib = Matrix.eye(2).kron(b)
    assertEquals(ib.rows, 8)
    assertEquals(ib.cols, 10)
    var i = 0
    while i < 4 do
      var j = 0
      while j < 5 do
        assertEquals(ib(i, j), b(i, j))
        assertEquals(ib(4 + i, 5 + j), b(i, j))
        assertEquals(ib(i, 5 + j), 0.0)
        assertEquals(ib(4 + i, j), 0.0)
        j += 1
      i += 1
  }

  test("kron mixed-product identity: (A kron B)(C kron D) = (AC) kron (BD)") {
    val a = randomMat(3, 2, 51L)
    val b = randomMat(2, 4, 52L)
    val c = randomMat(2, 3, 53L)
    val d = randomMat(4, 2, 54L)
    val lhs = a.kron(b) * c.kron(d)
    val rhs = (a * c).kron(b * d)
    assert(frob(lhs - rhs) < 1e-12, s"mixed product: ${frob(lhs - rhs)}")
  }

  test("kron reads strided transpose views correctly") {
    val a = randomMat(3, 2, 61L)
    val b = randomMat(2, 3, 62L)
    val viaViews = a.t.kron(b.t)
    val viaCopies = Matrix.tabulate(2, 3)((i, j) => a(j, i)).kron(Matrix.tabulate(3, 2)((i, j) => b(j, i)))
    assert(frob(viaViews - viaCopies) == 0.0)
  }

  test("kron with empty operands is total") {
    val a = randomMat(2, 3, 71L)
    val e = a.kron(DMat.zeros(0, 4))
    assertEquals(e.rows, 0)
    assertEquals(e.cols, 12)
  }
