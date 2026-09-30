package gale.linalg

class BandedLogDetJetSuite extends munit.FunSuite:
  private def factor(bands: DMat, options: CholeskyOptions = CholeskyOptions.Default): BandedCholesky =
    BandedCholesky.factorLower(bands, options).fold(throw _, identity)

  private def derivatives(first: Vector[DMat], secondUpper: Vector[DMat]): BandedDerivatives =
    BandedDerivatives(first, secondUpper).fold(throw _, identity)

  private def jet(f: BandedCholesky, first: Vector[DMat], secondUpper: Vector[DMat]): BandedLogDetJet =
    f.logDetJet(derivatives(first, secondUpper)).fold(throw _, identity)

  private def pair(p: Int, q: Int): Int = q * (q + 1) / 2 + p

  /** Lower band of a dense symmetric matrix, with NaN in the ignored padding. */
  private def bandsOf(dense: DMat, b: Int): DMat =
    Matrix.tabulate(dense.rows, b + 1)((i, d) => if d <= i then dense(i, i - d) else Double.NaN)

  /** Symmetric matrix of bandwidth `b` from a deterministic seed. */
  private def symmetricBand(n: Int, b: Int, seed: Double, diagonal: Double): DMat =
    Matrix.tabulate(n, n): (i, j) =>
      val distance = math.abs(i - j)
      if distance > b then 0.0
      else if i == j then diagonal + 0.3 * math.sin(seed + i)
      else 0.4 * math.cos(seed * 1.7 + i + j) / (1.0 + distance)

  private def trace(m: DMat): Double = (0 until m.rows).map(i => m(i, i)).sum

  /** Independent dense oracle: ∂p = tr(A⁻¹Ap), ∂pq = tr(A⁻¹Apq) - tr(A⁻¹Ap A⁻¹Aq). */
  private def denseJet(a: DMat, first: Vector[DMat], second: (Int, Int) => DMat): (Vector[Double], DMat) =
    val chol = a.cholesky.fold(throw _, identity)
    val solved = first.map(ap => chol.solve(ap).fold(throw _, identity))
    val d = first.length
    val gradient = solved.map(trace)
    val hessian = Matrix.tabulate(d, d): (p, q) =>
      trace(chol.solve(second(p, q)).fold(throw _, identity)) - trace(solved(p) * solved(q))
    (gradient, hessian)

  private def assertJet(actual: BandedLogDetJet, gradient: Vector[Double], hessian: DMat, tolerance: Double)(using
      munit.Location
  ): Unit =
    assertEquals(actual.dimension, gradient.length)
    for p <- gradient.indices do
      assertEqualsDouble(actual.gradient(p), gradient(p), tolerance * math.max(1.0, math.abs(gradient(p))))
    for p <- gradient.indices; q <- gradient.indices do
      assertEqualsDouble(actual.hessian(p, q), hessian(p, q), tolerance * math.max(1.0, math.abs(hessian(p, q))))

  test("dense trace identities for every mixed second derivative and supported bandwidth"):
    for n <- List(1, 2, 9, 17); b <- List(0, 1, 2, n - 1).distinct.filter(_ < n); d <- 1 to 3 do
      val a = symmetricBand(n, b, 0.2, 4.0 + b)
      val first = Vector.tabulate(d)(p => symmetricBand(n, b, 1.1 + p, 0.5 * (p + 1)))
      val second = Vector.tabulate(d, d)((p, q) => symmetricBand(n, b, 2.3 + math.min(p, q) + 3.0 * math.max(p, q), 0.7))
      val secondUpper = Vector.tabulate(d * (d + 1) / 2): index =>
        val q = (0 until d).find(q => q * (q + 1) / 2 <= index && index < (q + 1) * (q + 2) / 2).get
        second(index - q * (q + 1) / 2)(q)
      val f = factor(bandsOf(a, b))
      val result = jet(f, first.map(bandsOf(_, b)), secondUpper.map(bandsOf(_, b)))
      val (gradient, hessian) = denseJet(a, first, (p, q) => second(math.min(p, q))(math.max(p, q)))
      assertJet(result, gradient, hessian, 1e-11)
      for p <- 0 until d; q <- 0 until d do assertEquals(result.hessian(p, q), result.hessian(q, p))

  test("central finite differences converge to the jet of a nonlinear banded family"):
    val n = 11
    val b = 2
    val d = 3
    val base = symmetricBand(n, b, 0.4, 5.0)
    val linear = Vector.tabulate(d)(p => symmetricBand(n, b, 1.3 + p, 0.9))
    val quadratic = Vector.tabulate(d, d)((p, q) => symmetricBand(n, b, 3.1 + p + q, 0.4))
    // A(θ) = A0 + Σ θp Bp + ½ Σ θp θq Cpq + Σ sin(θp) Ep, with Cpq = Cqp.
    val periodic = Vector.tabulate(d)(p => symmetricBand(n, b, 5.7 + p, 0.3))
    def at(theta: Vector[Double]): DMat =
      var out = base
      for p <- 0 until d do out = out + linear(p) * theta(p) + periodic(p) * math.sin(theta(p))
      for p <- 0 until d; q <- 0 until d do out = out + quadratic(p)(q) * (0.5 * theta(p) * theta(q))
      out
    def logDet(theta: Vector[Double]): Double = factor(bandsOf(at(theta), b)).logDet
    val theta0 = Vector(0.21, -0.35, 0.44)
    val first = Vector.tabulate(d): p =>
      var out = linear(p) + periodic(p) * math.cos(theta0(p))
      for q <- 0 until d do out = out + quadratic(p)(q) * theta0(q)
      out
    val second = Vector.tabulate(d, d): (p, q) =>
      val periodicTerm = if p == q then periodic(p) * -math.sin(theta0(p)) else Matrix.zeros(n, n)
      quadratic(p)(q) + periodicTerm
    val secondUpper = for q <- (0 until d).toVector; p <- 0 to q yield bandsOf(second(p)(q), b)
    val result = jet(factor(bandsOf(at(theta0), b)), first.map(bandsOf(_, b)), secondUpper)
    def shifted(p: Int, h: Double): Vector[Double] = theta0.updated(p, theta0(p) + h)
    def shifted2(p: Int, hp: Double, q: Int, hq: Double): Vector[Double] =
      val once = shifted(p, hp)
      once.updated(q, once(q) + hq)
    val errors = List(1e-2, 5e-3).map: h =>
      var worst = 0.0
      for p <- 0 until d do
        val fd = (logDet(shifted(p, h)) - logDet(shifted(p, -h))) / (2.0 * h)
        worst = math.max(worst, math.abs(fd - result.gradient(p)))
        for q <- 0 until d do
          val fd2 =
            if p == q then (logDet(shifted(p, h)) - 2.0 * result.value + logDet(shifted(p, -h))) / (h * h)
            else
              (logDet(shifted2(p, h, q, h)) - logDet(shifted2(p, h, q, -h)) -
                logDet(shifted2(p, -h, q, h)) + logDet(shifted2(p, -h, q, -h))) / (4.0 * h * h)
          worst = math.max(worst, math.abs(fd2 - result.hessian(p, q)))
      worst
    assert(errors(0) < 1e-3, s"finite-difference errors $errors")
    // Central differences are O(h²): halving h cuts the error by about four.
    assert(errors(0) > 3.0 * errors(1), s"finite-difference errors do not converge: $errors")

  test("congruence D(θ) A0 D(θ) has analytic linear log determinant and exactly cancelling Hessian"):
    for (n, b) <- List((1, 0), (6, 0), (12, 1), (12, 3), (8, 7)) do
      val d = 3
      val a0 = symmetricBand(n, b, 0.9, 3.0 + b)
      val w = Matrix.tabulate(n, d)((i, p) => 0.1 * math.sin(i + 2.0 * p + 0.5))
      val theta = Vector(0.3, -0.2, 0.15)
      val scale = Vector.tabulate(n)(i => math.exp((0 until d).map(p => theta(p) * w(i, p)).sum))
      val a = Matrix.tabulate(n, n)((i, j) => scale(i) * scale(j) * a0(i, j))
      def rate(i: Int, j: Int, p: Int): Double = w(i, p) + w(j, p)
      val first = Vector.tabulate(d)(p => Matrix.tabulate(n, n)((i, j) => rate(i, j, p) * a(i, j)))
      val secondUpper = for q <- (0 until d).toVector; p <- 0 to q yield
        Matrix.tabulate(n, n)((i, j) => rate(i, j, p) * rate(i, j, q) * a(i, j))
      val result = jet(factor(bandsOf(a, b)), first.map(bandsOf(_, b)), secondUpper.map(bandsOf(_, b)))
      val logDetA0 = factor(bandsOf(a0, b)).logDet
      assertEqualsDouble(result.value, logDetA0 + 2.0 * (0 until n).map(i => math.log(scale(i))).sum, 1e-11)
      for p <- 0 until d do
        val expected = 2.0 * (0 until n).map(i => w(i, p)).sum
        assert(first(p).toDoubleArrayCopyRowMajor.toArray.exists(_ != 0.0))
        assertEqualsDouble(result.gradient(p), expected, 1e-12)
        for q <- 0 until d do assertEqualsDouble(result.hessian(p, q), 0.0, 1e-12)

  test("empty and singleton systems"):
    val empty = factor(DMat.zeros(0, 1))
    val emptyJet = jet(empty, Vector(DMat.zeros(0, 1)), Vector(DMat.zeros(0, 1)))
    assertEquals(emptyJet.value, 0.0)
    assertEquals(emptyJet.gradient(0), 0.0)
    assertEquals(emptyJet.hessian(0, 0), 0.0)
    val single = factor(Matrix.tabulate(1, 1)((_, _) => 4.0))
    val singleJet = jet(single, Vector(Matrix.tabulate(1, 1)((_, _) => 2.0)), Vector(Matrix.tabulate(1, 1)((_, _) => 3.0)))
    // log a: first a'/a, second a''/a - (a'/a)².
    assertEqualsDouble(singleJet.gradient(0), 0.5, 1e-15)
    assertEqualsDouble(singleJet.hessian(0, 0), 0.75 - 0.25, 1e-15)

  test("strided views and NaN padding are read logically; factor identity and inputs are preserved"):
    val n = 10
    val b = 2
    val a = symmetricBand(n, b, 0.6, 5.0)
    val ap = symmetricBand(n, b, 1.9, 1.0)
    val app = symmetricBand(n, b, 2.8, 0.6)
    val f = factor(bandsOf(a, b))
    val factorSnapshot = DMatBuilder.from(f.lowerBands).result()
    val contiguous = jet(f, Vector(bandsOf(ap, b)), Vector(bandsOf(app, b)))
    // Transposed storage and a sliced view of a larger buffer.
    val transposed = DMatBuilder.from(bandsOf(ap, b).t).result().t
    val padded = Matrix.tabulate(n + 3, b + 4)((i, d) =>
      if i >= 2 && i < n + 2 && d >= 1 && d <= b + 1 then bandsOf(app, b)(i - 2, d - 1) else Double.NaN)
    val slice = padded.slice(2, n + 2, 1, b + 2)
    val sliceSnapshot = DMatBuilder.from(slice).result()
    val strided = jet(f, Vector(transposed), Vector(slice))
    assertEquals(strided.gradient(0), contiguous.gradient(0))
    assertEquals(strided.hessian(0, 0), contiguous.hessian(0, 0))
    assert(strided.factor eq f)
    assertEquals(strided.value, f.logDet)
    for i <- 0 until n; d <- 0 to b do
      assertEquals(f.lowerBands(i, d), factorSnapshot(i, d))
      if d <= i then assertEquals(slice(i, d), sliceSnapshot(i, d))
    // The factor still solves the original system after the jet.
    val x = Matrix.tabulate(n, 3)((i, j) => math.cos(i + 2.0 * j))
    val solved = f.solve(a * x).fold(throw _, identity)
    for i <- 0 until n; j <- 0 until 3 do assertEqualsDouble(solved(i, j), x(i, j), 1e-12)
    val v = Vec.tabulate(n)(i => math.sin(i + 0.5))
    val sv = f.solve(a * v).fold(throw _, identity)
    for i <- 0 until n do assertEqualsDouble(sv(i), v(i), 1e-12)

  test("malformed derivative counts, shapes and nonfinite active entries are typed refusals"):
    val band = Matrix.tabulate(4, 2)((i, d) => if d == 0 then 3.0 else if d <= i then -1.0 else Double.NaN)
    val f = factor(band)
    def invalid(result: Either[LinAlgError, ?]): Boolean = result match
      case Left(LinAlgError.InvalidArgument(_)) => true
      case _ => false
    assert(invalid(BandedDerivatives(Vector.empty, Vector.empty)))
    assert(invalid(BandedDerivatives(Vector.fill(4)(band), Vector.fill(10)(band))))
    assert(invalid(BandedDerivatives(Vector(band, band), Vector(band, band))))
    assert(BandedDerivatives(Vector(band), Vector(DMat.zeros(4, 3))).left.toOption.exists:
      case LinAlgError.DimensionMismatch(_, _) => true
      case _ => false)
    val wider = derivatives(Vector(DMat.zeros(4, 3)), Vector(DMat.zeros(4, 3)))
    assertEquals(f.logDetJet(wider).left.toOption, Some(LinAlgError.DimensionMismatch(band.shape, DMat.zeros(4, 3).shape)))
    val badFirst = DMatBuilder.from(band).result()
    val withNaN = Matrix.tabulate(4, 2)((i, d) => if i == 2 && d == 1 then Double.NaN else band(i, d))
    assert(invalid(f.logDetJet(derivatives(Vector(withNaN), Vector(band)))))
    assert(invalid(f.logDetJet(derivatives(Vector(badFirst), Vector(withNaN)))))
    val infinite = Matrix.tabulate(4, 2)((i, d) => if i == 3 && d == 0 then Double.PositiveInfinity else band(i, d))
    assert(invalid(f.logDetJet(derivatives(Vector(infinite), Vector(band)))))
    // Padding (d > i) may hold NaN; only row 0, column 1 is padding here.
    assert(band(0, 1).isNaN)
    assert(f.logDetJet(derivatives(Vector(band), Vector(band))).isRight)

  test("near-singular and extreme ridge scales match dense traces; overflow is a typed refusal"):
    val n = 12
    val b = 1
    // Discrete Laplacian (eigenvalues down to about 0.06) plus ridge lambda; dA/dλ = I, d²A/dλ² = 0.
    for lambda <- List(1e-10, 1e-6, 1.0, 1e6) do
      val a = Matrix.tabulate(n, n)((i, j) => if i == j then 2.0 + lambda else if math.abs(i - j) == 1 then -1.0 else 0.0)
      val eye = Matrix.eye(n)
      val result = jet(factor(bandsOf(a, b)), Vector(bandsOf(eye, b)), Vector(bandsOf(DMat.zeros(n, n), b)))
      val (gradient, hessian) = denseJet(a, Vector(eye), (_, _) => DMat.zeros(n, n))
      assertJet(result, gradient, hessian, 1e-10)
    // A barely accepted singular Laplacian: tiny positive last pivot, huge but finite derivatives.
    val singular = Matrix.tabulate(n, n)((i, j) =>
      if i == j then (if i == 0 || i == n - 1 then 1.0 else 2.0) + 1e-13 else if math.abs(i - j) == 1 then -1.0 else 0.0)
    val tiny = factor(bandsOf(singular, b))
    val eye = Matrix.eye(n)
    val nearJet = jet(tiny, Vector(bandsOf(eye, b)), Vector(bandsOf(DMat.zeros(n, n), b)))
    assert(nearJet.gradient(0) > 1e11 && nearJet.gradient(0).isFinite)
    assert(nearJet.hessian(0, 0) < -1e22 && nearJet.hessian(0, 0).isFinite)
    // Overflowing derivatives are refused, not clipped, and the factor is unchanged.
    val snapshot = DMatBuilder.from(tiny.lowerBands).result()
    val huge = Matrix.tabulate(n, b + 1)((i, d) => if d <= i then 1e300 else Double.NaN)
    assert(tiny.logDetJet(derivatives(Vector(huge), Vector(huge))).left.toOption.exists:
      case LinAlgError.InvalidArgument(message) => message.contains("nonfinite")
      case _ => false)
    for i <- 0 until n; d <- 0 to b do assertEquals(tiny.lowerBands(i, d), snapshot(i, d))
    // The accepted absolute pivot policy is inherited: a refused factor yields no jet, with no repair.
    val strict = CholeskyOptions.Default.copy(pivotTolerance = 1e-10)
    assertEquals(BandedCholesky.factorLower(bandsOf(singular, b), strict).left.toOption, Some(LinAlgError.NotPositiveDefinite(n - 1)))

  test("large narrow band where a dense N-by-N oracle is infeasible"):
    val n = 100000
    val b = 2
    val d = 2
    // Congruence of a diagonally dominant pentadiagonal A0 by D(θ) = diag(exp(θ·w_i)).
    def a0(i: Int, delta: Int): Double = if delta == 0 then 6.0 + math.sin(i) else -0.8 / delta
    def w(i: Int, p: Int): Double = 1e-3 * math.cos(i * 0.01 + p)
    val theta = Vector(0.4, -0.7)
    def s(i: Int): Double = math.exp((0 until d).map(p => theta(p) * w(i, p)).sum)
    def entry(i: Int, delta: Int, factorOf: (Int, Int) => Double): Double =
      if delta > i then Double.NaN else s(i) * s(i - delta) * a0(i, delta) * factorOf(i, i - delta)
    val bands = Matrix.tabulate(n, b + 1)((i, delta) => entry(i, delta, (_, _) => 1.0))
    val first = Vector.tabulate(d)(p => Matrix.tabulate(n, b + 1)((i, delta) => entry(i, delta, (x, y) => w(x, p) + w(y, p))))
    val secondUpper = for q <- (0 until d).toVector; p <- 0 to q yield
      Matrix.tabulate(n, b + 1)((i, delta) => entry(i, delta, (x, y) => (w(x, p) + w(y, p)) * (w(x, q) + w(y, q))))
    val result = jet(factor(bands), first, secondUpper)
    for p <- 0 until d do
      val expected = 2.0 * (0 until n).map(i => w(i, p)).sum
      assertEqualsDouble(result.gradient(p), expected, 1e-9)
      for q <- 0 until d do assertEqualsDouble(result.hessian(p, q), 0.0, 1e-9)
