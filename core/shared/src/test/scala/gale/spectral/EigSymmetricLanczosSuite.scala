package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.LinearOperator
import gale.linalg.Matrix

/** Tests for the iterative partial symmetric eigensolver
  * `Eigen.eigSymmetric(op, n, Count, options)` — block Krylov with soft locking,
  * thick restarts, and full reorthogonalization. Convergence is checked against
  * dense or analytic oracles; repeated eigenspaces are checked as projectors.
  */
class EigSymmetricLanczosSuite extends munit.FunSuite:

  private def randomSymmetric(n: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        val v = rng.nextDouble() * 2.0 - 1.0
        a(i)(j) = v
        a(j)(i) = v
        j += 1
      i += 1
    Matrix.tabulate(n, n)((r, c) => a(r)(c))

  private def laplacian(n: Int): DMat =
    Matrix.tabulate(n, n): (i, j) =>
      if i == j then 2.0
      else if math.abs(i - j) == 1 then -1.0
      else 0.0

  /** Matrix-free diagonal operator `diag(d)`. */
  private def diagonalOperator(d: Array[Double]): DoubleLinearOperator =
    val n = d.length
    LinearOperator.fromFunction(n, n): (x, into) =>
      var i = 0
      while i < n do
        into(i) = d(i) * x(i)
        i += 1

  /** Reproduces BlockSymmetricEigen's portable deterministic probe stream. */
  private def deterministicProbe(n: Int, stream: Int): DVec =
    var state = 123456789 ^ (stream * 0x9e3779b9)
    DVec.tabulate(n): _ =>
      state = state * 1103515245 + 12345
      ((state >>> 9) & 0x7fffff).toDouble / 0x800000.toDouble * 2.0 - 1.0

  private def values(d: EigenDecomposition): IndexedSeq[Double] =
    (0 until d.eigenvalues.length).map(d.eigenvalues(_))

  private def denseSpectrum(a: DMat): IndexedSeq[Double] =
    val d = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
    (0 until d.size).map(d.eigenvalues(_))

  private def assertAscending(v: IndexedSeq[Double]): Unit =
    var i = 1
    while i < v.length do
      assert(v(i) >= v(i - 1) - 1e-10, s"not ascending at $i: $v")
      i += 1

  private def assertClose(a: IndexedSeq[Double], b: IndexedSeq[Double], tol: Double): Unit =
    assertEquals(a.length, b.length, s"length $a vs $b")
    a.zip(b).foreach { case (x, y) => assert(math.abs(x - y) < tol, s"$x != $y in $a vs $b") }

  private def projector(vectors: DMat, columns: IndexedSeq[Int]): DMat =
    val basis = DMat.tabulate(vectors.rows, columns.length)((row, col) => vectors(row, columns(col)))
    basis * basis.t

  private def projectorDistance(actual: DMat, expected: DMat): Double =
    var sum = 0.0
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        val delta = actual(row, col) - expected(row, col)
        sum += delta * delta
        col += 1
      row += 1
    math.sqrt(sum)

  private def assertProjectorClose(actual: DMat, expected: DMat, tolerance: Double): Unit =
    val distance = projectorDistance(actual, expected)
    assert(distance < tolerance, s"projector distance $distance >= $tolerance")

  test("block krylov: largest/smallest algebraic match the dense extremes") {
    val n = 60
    val a = randomSymmetric(n, 3131L)
    val full = denseSpectrum(a)
    val k = 6

    val largest = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestAlgebraic)).toOption.get
    assert(largest.diagnostics.allConverged, s"largest not converged: ${largest.diagnostics}")
    assertAscending(values(largest))
    assertClose(values(largest), full.takeRight(k), 1e-6)
    assert(largest.diagnostics.worstResidual < 1e-6)

    val smallest = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic)).toOption.get
    assert(smallest.diagnostics.allConverged, s"smallest not converged: ${smallest.diagnostics}")
    assertClose(values(smallest), full.take(k), 1e-6)
  }

  test("block krylov: BothEnds returns floor/ceil split of both extremes") {
    val n = 50
    val a = randomSymmetric(n, 5150L)
    val full = denseSpectrum(a)
    val k = 6
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.BothEnds)).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    val expected = (full.take(k / 2) ++ full.takeRight((k + 1) / 2)).sorted
    assertClose(values(d), expected, 1e-6)
  }

  test("block krylov: largest magnitude on a sign-straddling spectrum") {
    val n = 50
    val a = randomSymmetric(n, 2026L)
    val full = denseSpectrum(a)
    val k = 5
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestMagnitude)).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    val expected = full.sortBy(v => -math.abs(v)).take(k).sorted
    assertClose(values(d), expected, 1e-6)
  }

  test("block krylov: matrix-free diagonal operator resolves the smallest eigenvalues") {
    val n = 40
    val d = Array.tabulate(n)(i => (i + 1).toDouble) // eigenvalues 1..n
    val op = diagonalOperator(d)
    val k = 4
    val result = Eigen.eigSymmetric(op, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic)).toOption.get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(1.0, 2.0, 3.0, 4.0), 1e-7)
    // Eigenvectors satisfy op y = θ y.
    var c = 0
    while c < k do
      val y = result.eigenvectors.col(c)
      val r = op * y - y * result.eigenvalues(c)
      assert(r.norm2 < 1e-6)
      c += 1
  }

  test("block krylov: partial convergence returns Right with converged < requested") {
    val n = 40
    val a = randomSymmetric(n, 909L)
    val k = 6
    // Minimal subspace + a single restart + a punishing tolerance: one build
    // cannot resolve all k pairs, so the result is partial (never a Left).
    val options = SpectralOptions(tolerance = 1e-13, maxIterations = 1, subspaceDimension = Some(k + 1))
    val result = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestAlgebraic), options)
    val d = result.toOption.get
    assertEquals(d.diagnostics.requested, k)
    assert(d.diagnostics.converged < k, s"expected partial, got ${d.diagnostics.converged}")
    assert(!d.diagnostics.allConverged)
    assertEquals(d.size, d.diagnostics.converged)
    // requireConverged is the caller's opt-in to fail-fast.
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  test("block krylov: values-only omits eigenvectors but still converges") {
    val n = 50
    val a = randomSymmetric(n, 6161L)
    val k = 4
    val options = SpectralOptions(returnVectors = EigenVectors.ValuesOnly)
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.LargestAlgebraic), options).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assertEquals(d.eigenvectors.cols, 0)
    assertClose(values(d), denseSpectrum(a).takeRight(k), 1e-6)
  }

  test("block krylov: explicit start vector seeds the first block column") {
    val n = 40
    val a = randomSymmetric(n, 7171L)
    val k = 4
    val start = DVec.tabulate(n)(i => if i % 2 == 0 then 1.0 else 0.5)
    val options = SpectralOptions(startVector = Some(start))
    val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic), options).toOption.get
    assert(d.diagnostics.allConverged, s"not converged: ${d.diagnostics}")
    assertClose(values(d), denseSpectrum(a).take(k), 1e-6)
  }

  test("block krylov: structural violations return Left") {
    val n = 20
    val a = laplacian(n)
    val opts = SpectralOptions()

    // k out of range.
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(0, EigenOrder.LargestAlgebraic), opts).isLeft)
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(n, EigenOrder.LargestAlgebraic), opts).isLeft)
    // illegal order for symmetric.
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestRealPart), opts).isLeft)
    // square operator whose shape disagrees with n -> DimensionMismatch.
    Eigen.eigSymmetric(a, n + 1, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), opts) match
      case Left(_: LinAlgError.DimensionMismatch) => ()
      case other                                  => fail(s"expected DimensionMismatch, got $other")
    // genuinely non-square operator -> NonSquareMatrix (matches the dense path).
    val rect = LinearOperator.fromFunction(n, n + 1)((_, _) => ())
    Eigen.eigSymmetric(rect, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), opts) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected NonSquareMatrix, got $other")
    // illegal vector flag.
    val leftVecs = SpectralOptions(returnVectors = EigenVectors.Left)
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), leftVecs).isLeft)
    // wrong-length start vector.
    val badStart = SpectralOptions(startVector = Some(DVec.zeros(n - 1)))
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), badStart).isLeft)
    // zero start vector.
    val zeroStart = SpectralOptions(startVector = Some(DVec.zeros(n)))
    assert(Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), zeroStart).isLeft)
  }

  test("block krylov: repeated diagonal roots count algebraic multiplicity") {
    val d = Array(5.0, 5.0, 3.0, 2.0, 1.0)
    val op = diagonalOperator(d)
    // Seed exactly one direction in the repeated top eigenspace. The remaining
    // block columns must recover its independent companion.
    val start = DVec.tabulate(5)(i => if i == 0 then 1.0 else 0.0)
    val options = SpectralOptions(startVector = Some(start))
    val result = Eigen.eigSymmetric(op, 5, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), options).toOption.get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(3.0, 5.0, 5.0), 1e-9)
    assertEquals(result.diagnostics.residuals.length, 3)
    assert(result.diagnostics.worstResidual < 1e-9)
    assert(result.diagnostics.orthogonalityError < 1e-10)

    val repeatedColumns = (0 until result.size).filter(i => math.abs(result.eigenvalues(i) - 5.0) < 1e-8)
    assertEquals(repeatedColumns.length, 2)
    val expected = DMat.tabulate(5, 5)((i, j) => if i == j && i < 2 then 1.0 else 0.0)
    assertProjectorClose(projector(result.eigenvectors, repeatedColumns), expected, 1e-8)
  }

  test("block krylov: residual convergence does not certify a hidden spectral extreme") {
    val n = 8
    val start = DVec.tabulate(n)(_ => 1.0)
    val q0 = start * (1.0 / start.norm2)
    val probe = deterministicProbe(n, stream = 1)
    val q1Raw = probe - q0 * q0.dot(probe)
    val q1 = q1Raw * (1.0 / q1Raw.norm2)

    // Analytic spectrum: 0.8 on q0, 0.9 on q1, and 2.0 on their
    // six-dimensional orthogonal complement. The supplied start and Gale's
    // deterministic second block probe therefore span an invariant plane that
    // hides the requested largest algebraic eigenvalue.
    val op = LinearOperator.fromFunction(n, n): (x, into) =>
      val c0 = q0.dot(x)
      val c1 = q1.dot(x)
      var i = 0
      while i < n do
        into(i) = 2.0 * x(i) - 1.2 * c0 * q0(i) - 1.1 * c1 * q1(i)
        i += 1

    val trapped = Eigen
      .eigSymmetric(
        op,
        n,
        EigenSelection.Count(1, EigenOrder.LargestAlgebraic),
        SpectralOptions(startVector = Some(start), subspaceDimension = Some(2))
      )
      .toOption
      .get
    assertClose(values(trapped), IndexedSeq(0.9), 1e-12)
    assert(trapped.diagnostics.allConverged)
    assertEquals(
      trapped.diagnostics.convergenceStatus,
      SpectralConvergenceStatus.ResidualConverged
    )

    val fullSpace = Eigen
      .eigSymmetric(
        op,
        n,
        EigenSelection.Count(1, EigenOrder.LargestAlgebraic),
        SpectralOptions(startVector = Some(start), subspaceDimension = Some(n))
      )
      .toOption
      .get
    assertClose(values(fullSpace), IndexedSeq(2.0), 1e-12)
    assertEquals(
      fullSpace.diagnostics.convergenceStatus,
      SpectralConvergenceStatus.ExtremeCertified
    )
  }

  test("block krylov: repeated eigenspace converges without an n-dimensional projection") {
    val n = 120
    val diagonal = Array.tabulate(n): i =>
      if i < 8 then 100.0 else ((i - 8) % 3 + 1).toDouble
    var applications = 0
    val op = LinearOperator.fromFunction(n, n): (x, into) =>
      applications += 1
      var i = 0
      while i < n do
        into(i) = diagonal(i) * x(i)
        i += 1

    val options = SpectralOptions(
      tolerance = 1e-10,
      maxIterations = 4,
      subspaceDimension = Some(32)
    )
    val result = Eigen
      .eigSymmetric(op, n, EigenSelection.Count(8, EigenOrder.LargestAlgebraic), options)
      .toOption
      .get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq.fill(8)(100.0), 1e-8)
    val expected = DMat.tabulate(n, n)((i, j) => if i == j && i < 8 then 1.0 else 0.0)
    assertProjectorClose(projector(result.eigenvectors, 0 until 8), expected, 1e-8)
    assert(applications < n, s"solver used $applications operator applications for n=$n")
  }

  test("block krylov: four-cycle multiplicity and projector are matrix-free correct") {
    val n = 4
    val op = LinearOperator.fromFunction(n, n): (x, into) =>
      var i = 0
      while i < n do
        val left = (i + n - 1) % n
        val right = (i + 1) % n
        into(i) = 2.0 * x(i) - x(left) - x(right)
        i += 1

    val result = Eigen.eigSymmetric(op, n, EigenSelection.Count(3, EigenOrder.SmallestAlgebraic)).toOption.get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(0.0, 2.0, 2.0), 1e-9)

    val repeatedColumns = (0 until result.size).filter(i => math.abs(result.eigenvalues(i) - 2.0) < 1e-8)
    val scale = 1.0 / math.sqrt(2.0)
    val expectedBasis = DMat.tabulate(4, 2): (row, col) =>
      if col == 0 then
        if row == 0 then scale else if row == 2 then -scale else 0.0
      else if row == 1 then scale else if row == 3 then -scale else 0.0
    assertProjectorClose(
      projector(result.eigenvectors, repeatedColumns),
      expectedBasis * expectedBasis.t,
      1e-8
    )
  }

  test("block krylov: disconnected graph returns the complete repeated nullspace") {
    val n = 6
    val op = LinearOperator.fromFunction(n, n): (x, into) =>
      var i = 0
      while i < n do
        val local = i % 3
        var value = 0.0
        if local > 0 then
          value += x(i) - x(i - 1)
        if local < 2 then
          value += x(i) - x(i + 1)
        into(i) = value
        i += 1

    val result = Eigen.eigSymmetric(op, n, EigenSelection.Count(2, EigenOrder.SmallestAlgebraic)).toOption.get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(0.0, 0.0), 1e-9)

    val scale = 1.0 / math.sqrt(3.0)
    val expectedBasis = DMat.tabulate(n, 2): (row, col) =>
      if (col == 0 && row < 3) || (col == 1 && row >= 3) then scale else 0.0
    assertProjectorClose(
      projector(result.eigenvectors, IndexedSeq(0, 1)),
      expectedBasis * expectedBasis.t,
      1e-8
    )
  }

  test("block krylov: projector is invariant under graph reindexing") {
    val n = 4
    val permutation = Array(2, 0, 3, 1)
    def cycleEntry(i: Int, j: Int): Double =
      if i == j then 2.0
      else if (i + 1) % n == j || (j + 1) % n == i then -1.0
      else 0.0

    val permuted = LinearOperator.fromFunction(n, n): (x, into) =>
      var i = 0
      while i < n do
        var value = 0.0
        var j = 0
        while j < n do
          value += cycleEntry(permutation(i), permutation(j)) * x(j)
          j += 1
        into(i) = value
        i += 1

    val result = Eigen.eigSymmetric(permuted, n, EigenSelection.Count(3, EigenOrder.SmallestAlgebraic)).toOption.get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(0.0, 2.0, 2.0), 1e-9)

    val scale = 1.0 / math.sqrt(2.0)
    val originalBasis = DMat.tabulate(4, 2): (row, col) =>
      if col == 0 then
        if row == 0 then scale else if row == 2 then -scale else 0.0
      else if row == 1 then scale else if row == 3 then -scale else 0.0
    val originalProjector = originalBasis * originalBasis.t
    val expectedPermuted = DMat.tabulate(n, n)((i, j) => originalProjector(permutation(i), permutation(j)))
    assertProjectorClose(
      projector(result.eigenvectors, IndexedSeq(1, 2)),
      expectedPermuted,
      1e-8
    )
  }

  test("block krylov: BothEnds preserves repeated roots at both spectral ends") {
    val diagonal = Array(-5.0, -5.0, -1.0, 0.0, 2.0, 7.0, 7.0)
    val result = Eigen
      .eigSymmetric(diagonalOperator(diagonal), diagonal.length, EigenSelection.Count(4, EigenOrder.BothEnds))
      .toOption
      .get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(-5.0, -5.0, 7.0, 7.0), 1e-9)
  }

  test("block krylov: smallest magnitude preserves a repeated interior root") {
    val diagonal = Array(-4.0, -1.0, -1.0, 2.0, 6.0)
    val result = Eigen
      .eigSymmetric(diagonalOperator(diagonal), diagonal.length, EigenSelection.Count(2, EigenOrder.SmallestMagnitude))
      .toOption
      .get
    assert(result.diagnostics.allConverged, s"not converged: ${result.diagnostics}")
    assertClose(values(result), IndexedSeq(-1.0, -1.0), 1e-9)
  }

  test("block krylov: clustered spectrum converges via thick-restart growth") {
    // A 1e-6-wide cluster at the bottom is unresolvable in the initial minimal
    // subspace; the growth strategy (not a lossy restart) must separate it.
    val n = 40
    val d = Array.tabulate(n)(i => if i == 0 then 1.0 else if i == 1 then 1.0 + 1e-6 else (i + 1).toDouble)
    val op = diagonalOperator(d)
    val k = 2
    val options = SpectralOptions(tolerance = 1e-9, maxIterations = 50, subspaceDimension = Some(k + 2))
    val result = Eigen.eigSymmetric(op, n, EigenSelection.Count(k, EigenOrder.SmallestAlgebraic), options).toOption.get
    assert(result.diagnostics.allConverged, s"cluster not resolved: ${result.diagnostics}")
    assert(result.diagnostics.iterations > 1, s"expected growth, converged in ${result.diagnostics.iterations} build(s)")
    assertClose(values(result), IndexedSeq(1.0, 1.0 + 1e-6), 1e-8)
  }

  test("block krylov: zero convergence returns Right(empty), never a Left") {
    // An unsatisfiable tolerance: no pair can converge, so the result is an empty
    // Right with converged == 0 — the parity doc's zero-convergence boundary.
    val n = 30
    val a = randomSymmetric(n, 8888L)
    val options = SpectralOptions(tolerance = 0.0, maxIterations = 1, subspaceDimension = Some(5))
    val result = Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), options)
    val d = result.toOption.get
    assertEquals(d.size, 0)
    assertEquals(d.diagnostics.converged, 0)
    assert(!d.diagnostics.allConverged)
    assertEquals(d.diagnostics.requested, 3)
    d.requireConverged match
      case Left(_: LinAlgError.DidNotConverge) => ()
      case other                               => fail(s"expected DidNotConverge, got $other")
  }

  test("block krylov: shift-invert / target is deferred") {
    val n = 20
    val a = laplacian(n)
    val target = Some(SpectralTarget.ShiftInvert(0.5, LinearSolvePlan.Direct))
    val result = Eigen.eigSymmetric(a, n, EigenSelection.Count(3, EigenOrder.SmallestAlgebraic), SpectralOptions(), target)
    result match
      case Left(_: LinAlgError.UnsupportedOperation) => ()
      case other                                     => fail(s"expected UnsupportedOperation, got $other")
  }
