package gale.linalg

class DctBasisSuite extends munit.FunSuite:
  private val a = math.sqrt(2.0 + math.sqrt(2.0)) / math.sqrt(8.0)
  private val b = math.sqrt(2.0 - math.sqrt(2.0)) / math.sqrt(8.0)
  // Radical values: independent of the cosine formula used by the implementation.
  private val four = Vector(
    Vector(0.5, a, 0.5, b),
    Vector(0.5, b, -0.5, -a),
    Vector(0.5, -b, -0.5, a),
    Vector(0.5, -a, 0.5, -b)
  )

  test("four-sample basis matches the closed-form radical oracle") {
    val basis = DctBasis.columns(4, 0, 4).orThrow
    assertEquals(basis.rows, 4)
    assertEquals(basis.cols, 4)
    for row <- 0 until 4; col <- 0 until 4 do assertEqualsDouble(basis(row, col), four(row)(col), 8e-16)
  }

  test("selected columns preserve component identity and optional unscaled normalization") {
    for first <- 0 to 4; count <- 0 to (4 - first) do
      for norm <- DctNormalization.values do
        val basis = DctBasis.columns(4, first, count, norm).orThrow
        assertEquals(basis.rows, 4)
        assertEquals(basis.cols, count)
        for row <- 0 until 4; col <- 0 until count do
          val component = first + col
          val scale = norm match
            case DctNormalization.Orthonormal => 1.0
            case DctNormalization.Unscaled    => if component == 0 then 2.0 else math.sqrt(2.0)
          assertEqualsDouble(basis(row, col), four(row)(component) * scale, 1e-15)
  }

  test("orthonormal columns satisfy independent scalar inner products") {
    for n <- Vector(1, 2, 3, 7, 31) do
      val basis = DctBasis.columns(n, 0, n).orThrow
      for left <- 0 until n; right <- 0 until n do
        var dot = 0.0
        var row = 0
        while row < n do
          dot += basis(row, left) * basis(row, right)
          row += 1
        assertEqualsDouble(dot, if left == right then 1.0 else 0.0, 1e-14)
  }

  test("basis diagonalizes the path Laplacian including its boundary rows") {
    // The path Laplacian has endpoint degree one and interior degree two.
    // This identity checks half-sample phase, ordering, and both boundaries.
    for n <- Vector(2, 3, 8, 29) do
      val basis = DctBasis.columns(n, 0, n).orThrow
      for col <- 0 until n do
        val eigenvalue = 4.0 * math.pow(math.sin(math.Pi * col.toDouble / (2.0 * n)), 2)
        for row <- 0 until n do
          val left = if row == 0 then 0.0 else basis(row, col) - basis(row - 1, col)
          val right = if row == n - 1 then 0.0 else basis(row, col) - basis(row + 1, col)
          // Phase evaluation error grows with the highest component, hence with n.
          assertEqualsDouble(left + right, eigenvalue * basis(row, col), 8.0 * math.ulp(1.0) * n)
  }

  test("one sample and empty ranges retain dimensions, including an enormous empty matrix") {
    assertEquals(DctBasis.columns(1, 0, 1).orThrow(0, 0), 1.0)
    for first <- Vector(0, Int.MaxValue) do
      val empty = DctBasis.columns(Int.MaxValue, first, 0).orThrow
      assertEquals(empty.rows, Int.MaxValue)
      assertEquals(empty.cols, 0)
      assertEquals(empty.valuesRowMajor.size, 0)
  }

  test("invalid and overflowing requests fail before allocation") {
    val requests = Vector(
      (0, 0, 0),
      (-1, 0, 0),
      (4, -1, 1),
      (4, 0, -1),
      (4, 5, 0),
      (4, 3, 2),
      (Int.MaxValue, Int.MaxValue, 1),
      (Int.MaxValue, Int.MaxValue, Int.MaxValue),
      (Int.MaxValue, 0, 2),
      (50000, 0, 50000)
    )
    requests.foreach { case (samples, first, count) =>
      DctBasis.columns(samples, first, count) match
        case Left(_: LinAlgError.InvalidArgument) => ()
        case other => fail(s"expected InvalidArgument for ($samples,$first,$count), got $other")
    }
  }
