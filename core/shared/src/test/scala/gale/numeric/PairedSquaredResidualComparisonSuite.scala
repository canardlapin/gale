package gale.numeric

class PairedSquaredResidualComparisonSuite extends munit.FunSuite:
  private def ok[A](value: Either[PairedResidualError, A]): A = value.fold(e => fail(e.message), identity)

  // This oracle evaluates two complete squared residual sums as exact dyadic
  // rationals. It does not use the implementation's paired difference formula.
  private final case class Dyadic(value: BigInt, exponent: Int):
    def +(that: Dyadic): Dyadic =
      val power = math.min(exponent, that.exponent)
      Dyadic((value << (exponent - power)) + (that.value << (that.exponent - power)), power)
    def unary_- : Dyadic = Dyadic(-value, exponent)
    def -(that: Dyadic): Dyadic = this + (-that)
    def *(that: Dyadic): Dyadic = Dyadic(value * that.value, exponent + that.exponent)
    def compare(that: Dyadic): Int =
      val power = math.min(exponent, that.exponent)
      (value << (exponent - power)).compare(that.value << (that.exponent - power))

  private val zero = Dyadic(BigInt(0), 0)

  private def exact(value: Double): Dyadic =
    assert(value.isFinite)
    val bits = java.lang.Double.doubleToRawLongBits(value)
    val exponent = ((bits >>> 52) & 0x7ffL).toInt
    val fraction = bits & 0x000fffffffffffffL
    val significand = if exponent == 0 then fraction else fraction | (1L << 52)
    Dyadic(BigInt(if bits < 0 then -significand else significand), if exponent == 0 then -1074 else exponent - 1075)

  private def residualEnergy(f: PairedResidualFixtures.Fixture, design: Array[Double], beta: Array[Double]): Dyadic =
    var sum = zero
    var row = 0
    while row < f.rows do
      var prediction = zero
      var col = 0
      while col < f.columns do
        prediction = prediction + exact(design(row * f.columns + col)) * exact(beta(col))
        col += 1
      val residual = exact(f.response(row)) - prediction
      sum = sum + residual * residual
      row += 1
    sum

  private def difference(f: PairedResidualFixtures.Fixture): Dyadic =
    residualEnergy(f, f.candidateDesign, f.candidateCoefficients) - residualEnergy(f, f.previousDesign, f.previousCoefficients)

  private def compare(f: PairedResidualFixtures.Fixture, scope: PairedResidualScope = PairedResidualScope.ProfileMinimum): PairedResidualComparison =
    val workspace = ok(PairedResidualWorkspace(f.rows, f.columns))
    ok(workspace.compare(f.response, f.previousDesign, f.candidateDesign, f.previousCoefficients, f.candidateCoefficients, scope))

  private def enclosed(result: PairedResidualComparison, oracle: Dyadic): Unit =
    assert(exact(result.modelDifference.lower).compare(oracle) <= 0, "exact RSS difference below lower endpoint")
    assert(exact(result.modelDifference.upper).compare(oracle) >= 0, "exact RSS difference above upper endpoint")
    result.direction match
      case ResidualDifferenceDirection.Decrease => assert(oracle.compare(zero) < 0)
      case ResidualDifferenceDirection.Increase => assert(oracle.compare(zero) > 0)
      case ResidualDifferenceDirection.Unresolved => ()

  PairedResidualFixtures.all.foreach: fixture =>
    test(s"bit-exact independent witness: ${fixture.name}"):
      val result = compare(fixture)
      enclosed(result, difference(fixture))
      assertEquals(result.certifiesProfileDecrease, fixture.profileDecrease)
      if fixture.name == "darwin-node2421-voxel-15.json" || fixture.name == "linux-node2421-voxel-2.json" then
        assertEquals(result.direction, ResidualDifferenceDirection.Increase)
        assert(!result.certifiesProfileDecrease)

  test("poor previous coefficients cannot falsely certify a profile decrease"):
    val f = PairedResidualFixtures.all.find(_.name == "synthetic-old-poor-beta.json").get
    val result = compare(f)
    assertEquals(result.direction, ResidualDifferenceDirection.Decrease)
    // Previous best RSS is exactly 0; candidate best RSS is exactly 1, although
    // these supplied previous coefficients have RSS 2.
    assert(result.modelDifference.contains(-1.0))
    result.profileMinimum match
      case ProfileMinimumComparison.Unresolved(ProfileMinimumUnresolved.DecreaseNotCertified(bound)) =>
        assert(bound.suboptimalityUpper >= 2.0)
        assert(bound.differenceUpper >= 1.0)
      case other => fail(s"expected the old-gap guard, got $other")

  test("a positive definite Gram need not have a positive Gershgorin lower bound"):
    val f = PairedResidualFixtures.all.find(_.name == "synthetic-insufficient-gram-bound.json").get
    // D=[[1,2],[0,1]], so Gram=[[1,2],[2,5]] has determinant 1 > 0.
    val result = compare(f)
    assertEquals(result.direction, ResidualDifferenceDirection.Decrease)
    result.profileMinimum match
      case ProfileMinimumComparison.Unresolved(ProfileMinimumUnresolved.NonPositiveGramLowerBound(lower)) => assert(lower <= 0.0)
      case other => fail(s"expected an unresolved sufficient bound, got $other")

  test("returned-model scope never implies an unrequested profile certificate"):
    val f = PairedResidualFixtures.all.find(_.name == "synthetic-profile-descent.json").get
    val result = compare(f, PairedResidualScope.ReturnedModels)
    assertEquals(result.direction, ResidualDifferenceDirection.Decrease)
    assertEquals(result.profileMinimum, ProfileMinimumComparison.NotRequested)
    assert(!result.certifiesProfileDecrease)

  test("common scalar energy offsets cannot erase a proven tiny uphill difference"):
    val f = PairedResidualFixtures.all.find(_.name == "darwin-node2421-voxel-15.json").get
    val result = compare(f)
    enclosed(result, difference(f))
    assertEquals(result.direction, ResidualDifferenceDirection.Increase)
    def roundedEnergy(design: Array[Double], beta: Array[Double]): Double =
      var energy = 0.0
      var row = 0
      while row < f.rows do
        var residual = f.response(row)
        var col = 0
        while col < f.columns do
          residual -= design(row * f.columns + col) * beta(col)
          col += 1
        energy += residual * residual
        row += 1
      energy
    for common <- Vector(1e8, 1e16) do
      val previous = common + roundedEnergy(f.previousDesign, f.previousCoefficients)
      val candidate = common + roundedEnergy(f.candidateDesign, f.candidateCoefficients)
      assertEqualsDouble(candidate - previous, 0.0, 0.0)
    // The generic API has no common energy constant or acceptance tolerance.
    assert(!result.certifiesProfileDecrease)

  test("uniform dyadic column scaling preserves the physical model comparison"):
    val f = PairedResidualFixtures.all.find(_.name == "synthetic-coupled-profile-descent.json").get
    val reference = difference(f)
    for power <- Vector(-40, -10, 10, 40) do
      val scale = math.pow(2.0, power.toDouble)
      val scaled = f.copy(previousDesign = f.previousDesign.map(_ * scale), candidateDesign = f.candidateDesign.map(_ * scale),
        previousCoefficients = f.previousCoefficients.map(_ / scale), candidateCoefficients = f.candidateCoefficients.map(_ / scale))
      assertEquals(difference(scaled).compare(reference), 0)
      val result = compare(scaled)
      enclosed(result, reference)
      assert(result.certifiesProfileDecrease)

  test("zero and subnormal changes remain conservative rather than manufacturing descent"):
    for value <- Vector(0.0, Double.MinPositiveValue, -Double.MinPositiveValue) do
      val f = PairedResidualFixtures.Fixture("tiny", 1, 1, Array(0.0), Array(1.0), Array(1.0), Array(0.0), Array(value), false)
      val result = compare(f)
      enclosed(result, difference(f))
      assertEquals(result.direction, ResidualDifferenceDirection.Unresolved)
      assert(!result.certifiesProfileDecrease)

  test("finite input overflow or an infinite outward endpoint refuses certification"):
    val workspace = ok(PairedResidualWorkspace(1, 1))
    assert(workspace.compare(Array(0.0), Array(Double.MaxValue), Array(Double.MaxValue), Array(2.0), Array(2.0)).isLeft)
    assert(workspace.compare(Array(Double.MaxValue), Array(0.0), Array(0.0), Array(0.0), Array(0.0)).isLeft)

  test("shape admission and every nonfinite input channel are typed failures"):
    assert(PairedResidualWorkspace(0, 1).isLeft)
    assert(PairedResidualWorkspace.requiredArrayCells(46341, 46341).isLeft)
    assertEquals(ok(PairedResidualWorkspace.requiredArrayCells(72, 3)), 150L)
    val workspace = ok(PairedResidualWorkspace(2, 1))
    assert(workspace.compare(Array(1.0), Array(1.0, 1.0), Array(1.0, 0.0), Array(0.0), Array(1.0)).isLeft)
    assertEquals(workspace.compare(null, Array(1.0, 1.0), Array(1.0, 0.0), Array(0.0), Array(1.0)),
      Left(PairedResidualError.MissingInput(PairedResidualInput.Response)))
    for bad <- Vector(Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity); channel <- 0 until 5 do
      val arrays = Array(Array(1.0, 1.0), Array(1.0, 1.0), Array(1.0, 0.0), Array(0.0), Array(1.0))
      arrays(channel)(0) = bad
      workspace.compare(arrays(0), arrays(1), arrays(2), arrays(3), arrays(4)) match
        case Left(PairedResidualError.NonFiniteInput(_, 0)) => ()
        case other => fail(s"expected nonfinite input failure, got $other")

  test("workspace reuse after failure cannot mutate an earlier immutable result"):
    val f = PairedResidualFixtures.all.find(_.name == "synthetic-profile-descent.json").get
    val workspace = ok(PairedResidualWorkspace(f.rows, f.columns))
    val first = ok(workspace.compare(f.response, f.previousDesign, f.candidateDesign, f.previousCoefficients, f.candidateCoefficients))
    val lower = java.lang.Double.doubleToRawLongBits(first.modelDifference.lower)
    val upper = java.lang.Double.doubleToRawLongBits(first.modelDifference.upper)
    val invalid = f.response.clone()
    invalid(0) = Double.NaN
    assert(workspace.compare(invalid, f.previousDesign, f.candidateDesign, f.previousCoefficients, f.candidateCoefficients).isLeft)
    val again = ok(workspace.compare(f.response, f.previousDesign, f.candidateDesign, f.previousCoefficients, f.candidateCoefficients))
    assertEquals(java.lang.Double.doubleToRawLongBits(first.modelDifference.lower), lower)
    assertEquals(java.lang.Double.doubleToRawLongBits(first.modelDifference.upper), upper)
    assertEquals(java.lang.Double.doubleToRawLongBits(again.modelDifference.lower), lower)
    assertEquals(java.lang.Double.doubleToRawLongBits(again.modelDifference.upper), upper)

  test("deterministic finite paired matrices enclose an independent exact RSS subtraction"):
    var sample = 0
    while sample < 128 do
      val old = Array.tabulate(12)(i => ((i * 37 + sample * 11) % 127 - 63).toDouble / 32.0)
      val next = old.zipWithIndex.map((x, i) => x + ((sample + i) % 7 - 3).toDouble / 1048576.0)
      val b0 = Array(0.25 + sample.toDouble / 1024.0, -0.75, 1.5)
      val b1 = Array(b0(0) - 1.0 / 65536.0, b0(1) + 1.0 / 32768.0, b0(2))
      val y = Array.tabulate(4)(i => (i * 17 - sample % 13).toDouble / 16.0)
      val f = PairedResidualFixtures.Fixture("generated", 4, 3, y, old, next, b0, b1, false)
      enclosed(compare(f, PairedResidualScope.ReturnedModels), difference(f))
      sample += 1
