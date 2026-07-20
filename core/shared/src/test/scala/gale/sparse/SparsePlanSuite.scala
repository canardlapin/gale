package gale.sparse

import gale.linalg.*

class SparsePlanSuite extends munit.FunSuite:
  private def csr(rows: Int, cols: Int, entries: (Int, Int, Double)*): CSR =
    val builder = Sparse.coo(rows, cols)
    entries.foreach { case (row, col, value) => builder.add(row, col, value) }
    builder.toCSR()

  private def assertDenseApprox(actual: DMat, expected: DMat, tolerance: Double = 1e-12): Unit =
    assertEquals(actual.shape, expected.shape)
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assert(
          math.abs(actual(row, col) - expected(row, col)) <= tolerance,
          s"($row,$col): ${actual(row, col)} != ${expected(row, col)}"
        )
        col += 1
      row += 1

  private def denseProduct(left: DMat, right: DMat): DMat =
    Matrix.tabulate(left.rows, right.cols): (row, col) =>
      var sum = 0.0
      var inner = 0
      while inner < left.cols do
        sum += left(row, inner) * right(inner, col)
        inner += 1
      sum

  private def positions(pattern: CSRPattern): Set[(Int, Int)] =
    val out = Set.newBuilder[(Int, Int)]
    pattern.foreachStoredPosition((row, col) => out += ((row, col)))
    out.result()

  test("union analysis replays changing values into one destination and preserves explicit cancellation") {
    val left0 = csr(2, 4, (0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0))
    val right0 = csr(2, 4, (0, 1, 4.0), (0, 2, -2.0), (1, 3, 5.0))
    val plan = CSRUnionPlan.analyze(left0.pattern, right0.pattern).toOption.get
    assert(plan.resultPattern.hasCanonicalFormat)
    assertEquals(positions(plan.resultPattern), Set((0, 0), (0, 1), (0, 2), (1, 1), (1, 3)))

    val destination = plan.newDestination()
    assertEquals(plan.evaluateInto(left0, right0, destination), Right(()))
    val first = destination.snapshot()
    assertEquals(first.nnz, 5)
    assertEquals(first(0, 2), 0.0)
    assert(!first.hasCanonicalFormat, "an explicit numeric zero must remain stored")
    assert(first.pattern.hasCanonicalFormat)
    assertEquals(first.pruneZeros.nnz, 4)
    assertDenseApprox(first.toDense(), left0.toDense() + right0.toDense())
    val firstDense = first.toDense()

    val left1 = plan.leftPattern.bind(Array(10.0, 20.0, 30.0)).toOption.get
    val right1 = plan.rightPattern.bind(Array(-1.0, 2.0, -3.0)).toOption.get
    assertEquals(plan.evaluateInto(left1, right1, destination, leftScale = 0.5, rightScale = 2.0), Right(()))
    val second = destination.snapshot()
    val expected = Matrix.tabulate(2, 4): (row, col) =>
      0.5 * left1(row, col) + 2.0 * right1(row, col)
    assertDenseApprox(second.toDense(), expected)
    assertDenseApprox(first.toDense(), firstDense)
  }

  test("union allocating replay shares result structure and zero scales suppress source reads") {
    val left0 = csr(1, 2, (0, 0, Double.NaN))
    val right0 = csr(1, 2, (0, 1, 3.0))
    val plan = CSRUnionPlan.analyze(left0.pattern, right0.pattern).toOption.get
    val result = plan.evaluate(left0, right0, leftScale = 0.0, rightScale = 1.0).toOption.get
    assertEquals(result(0, 0), 0.0)
    assertEquals(result(0, 1), 3.0)
    assertEquals(result.pattern, plan.resultPattern)
    assertEquals(result.nnz, 2)
  }

  test("union rejects noncanonical analysis and any input or destination pattern changed after analysis") {
    val canonical = csr(1, 3, (0, 0, 1.0), (0, 2, 2.0))
    val other = csr(1, 3, (0, 1, 3.0))
    val plan = CSRUnionPlan.analyze(canonical.pattern, other.pattern).toOption.get

    val changed = csr(1, 3, (0, 0, 1.0), (0, 1, 2.0))
    assert(plan.evaluateInto(changed, other, plan.newDestination()).isLeft)

    // Even structurally equal but independently owned storage is a changed
    // symbolic input. Rebind through the analyzed pattern for accepted updates.
    val equalButIndependent = csr(1, 3, (0, 0, 8.0), (0, 2, 9.0))
    assert(plan.evaluateInto(equalButIndependent, other, plan.newDestination()).isLeft)

    val anotherPlan = CSRUnionPlan.analyze(canonical.pattern, canonical.pattern).toOption.get
    assert(plan.evaluateInto(canonical, other, anotherPlan.newDestination()).isLeft)

    val noncanonical = CSRPattern.checked(1, 3, Array(0, 2), Array(2, 0)).toOption.get
    assert(CSRUnionPlan.analyze(noncanonical, canonical.pattern).isLeft)
    assert(CSRUnionPlan.analyze(canonical.pattern, csr(2, 3).pattern).isLeft)
  }

  test("product analysis precomputes contributions and replay matches a dense oracle") {
    val left0 = csr(2, 3, (0, 0, 2.0), (0, 2, 1.0), (1, 1, 4.0))
    val right0 = csr(3, 2, (0, 0, 5.0), (1, 1, 6.0), (2, 0, -10.0), (2, 1, 7.0))
    val plan = CSRProductPlan.analyze(left0.pattern, right0.pattern).toOption.get
    assertEquals(plan.contributionCount, 4)
    assert(plan.resultPattern.hasCanonicalFormat)
    assertEquals(positions(plan.resultPattern), Set((0, 0), (0, 1), (1, 1)))

    val destination = plan.newDestination()
    assertEquals(plan.evaluateInto(left0, right0, destination), Right(()))
    val first = destination.snapshot()
    assertEquals(first.nnz, 3)
    assertEquals(first(0, 0), 0.0)
    assert(!first.hasCanonicalFormat)
    assertDenseApprox(first.toDense(), denseProduct(left0.toDense(), right0.toDense()))
    val firstDense = first.toDense()

    val left1 = plan.leftPattern.bind(Array(1.0, 2.0, 3.0)).toOption.get
    val right1 = plan.rightPattern.bind(Array(4.0, 5.0, 6.0, 7.0)).toOption.get
    assertEquals(plan.evaluateInto(left1, right1, destination, scale = -0.5), Right(()))
    val second = destination.snapshot()
    val unscaled = denseProduct(left1.toDense(), right1.toDense())
    val expected = Matrix.tabulate(2, 2)((r, c) => -0.5 * unscaled(r, c))
    assertDenseApprox(second.toDense(), expected)
    assertDenseApprox(first.toDense(), firstDense)
  }

  test("product rejects incompatible, noncanonical, changed input, and changed destination patterns") {
    val left = csr(2, 3, (0, 0, 1.0), (1, 2, 2.0))
    val right = csr(3, 2, (0, 0, 3.0), (2, 1, 4.0))
    val plan = CSRProductPlan.analyze(left.pattern, right.pattern).toOption.get

    assert(CSRProductPlan.analyze(left.pattern, csr(4, 2).pattern).isLeft)
    val noncanonical = CSRPattern.checked(2, 3, Array(0, 2, 2), Array(2, 0)).toOption.get
    assert(CSRProductPlan.analyze(noncanonical, right.pattern).isLeft)

    val changedLeft = csr(2, 3, (0, 1, 1.0), (1, 2, 2.0))
    assert(plan.evaluateInto(changedLeft, right, plan.newDestination()).isLeft)
    val changedRight = csr(3, 2, (0, 1, 3.0), (2, 1, 4.0))
    assert(plan.evaluateInto(left, changedRight, plan.newDestination()).isLeft)

    val other = CSRProductPlan.analyze(left.pattern, csr(3, 3, (0, 0, 1.0)).pattern).toOption.get
    assert(plan.evaluateInto(left, right, other.newDestination()).isLeft)
  }

  test("zero-sized symbolic union and product are total") {
    val empty00 = CSRPattern.checked(0, 0, Array(0), Array.empty[Int]).toOption.get
    val union = CSRUnionPlan.analyze(empty00, empty00).toOption.get
    val emptyMatrix = empty00.bind(Array.empty[Double]).toOption.get
    val unionDestination = union.newDestination()
    assertEquals(union.evaluateInto(emptyMatrix, emptyMatrix, unionDestination), Right(()))
    assertEquals(unionDestination.snapshot().nnz, 0)

    val left = CSRPattern.checked(2, 0, Array(0, 0, 0), Array.empty[Int]).toOption.get
    val right = CSRPattern.checked(0, 4, Array(0), Array.empty[Int]).toOption.get
    val product = CSRProductPlan.analyze(left, right).toOption.get
    assertEquals(product.resultPattern.rows, 2)
    assertEquals(product.resultPattern.cols, 4)
    assertEquals(product.resultPattern.nnz, 0)
    assertEquals(product.contributionCount, 0)
    val leftMatrix = left.bind(Array.empty[Double]).toOption.get
    val rightMatrix = right.bind(Array.empty[Double]).toOption.get
    assertEquals(product.evaluateInto(leftMatrix, rightMatrix, product.newDestination()), Right(()))
  }

  test("randomized union and Boolean-product pattern laws plus numeric replay hold") {
    var seed = 0
    while seed < 60 do
      val rng = new scala.util.Random(2026072000L + seed)
      val rows = 1 + rng.nextInt(6)
      val inner = 1 + rng.nextInt(6)
      val cols = 1 + rng.nextInt(6)

      def randomMatrix(r: Int, c: Int): CSR =
        val builder = Sparse.coo(r, c)
        var i = 0
        while i < r do
          var j = 0
          while j < c do
            if rng.nextDouble() < 0.35 then builder.add(i, j, rng.nextDouble() * 4.0 - 2.0)
            j += 1
          i += 1
        builder.toCSR()

      val unionLeft = randomMatrix(rows, cols)
      val unionRight = randomMatrix(rows, cols)
      val unionPlan = CSRUnionPlan.analyze(unionLeft.pattern, unionRight.pattern).toOption.get
      assertEquals(
        positions(unionPlan.resultPattern),
        positions(unionLeft.pattern) union positions(unionRight.pattern),
        s"union seed=$seed"
      )
      val unionResult = unionPlan.evaluate(unionLeft, unionRight).toOption.get
      assertDenseApprox(unionResult.toDense(), unionLeft.toDense() + unionRight.toDense())

      val productLeft = randomMatrix(rows, inner)
      val productRight = randomMatrix(inner, cols)
      val productPlan = CSRProductPlan.analyze(productLeft.pattern, productRight.pattern).toOption.get
      val leftPositions = positions(productLeft.pattern)
      val rightPositions = positions(productRight.pattern)
      val expectedPattern =
        (for
          row <- 0 until rows
          col <- 0 until cols
          if (0 until inner).exists(k => leftPositions.contains((row, k)) && rightPositions.contains((k, col)))
        yield (row, col)).toSet
      assertEquals(positions(productPlan.resultPattern), expectedPattern, s"product pattern seed=$seed")
      val productResult = productPlan.evaluate(productLeft, productRight).toOption.get
      assertDenseApprox(productResult.toDense(), denseProduct(productLeft.toDense(), productRight.toDense()), 1e-10)
      seed += 1
  }
