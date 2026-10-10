package gale.numeric

class FDistributionSuite extends munit.FunSuite:

  private def dist(df1: Double, df2: Double): FDistribution =
    FDistribution(df1, df2).fold(e => fail(e.message), identity)

  private def assertRelative(actual: Double, expected: Double, tol: Double, clue: String): Unit =
    val scale = math.abs(expected)
    assertEqualsDouble(actual, expected, tol * scale, clue)

  // (f, df1, df2, pf(f, df1, df2, lower.tail = FALSE)) generated with R 4.3.3:
  //   Rscript -e 'cat(sprintf("%.17g\n", pf(f, df1, df2, lower.tail = FALSE)))'
  // for each row below (inputs printed with format(..., digits = 17)).
  private val moderate = List(
    (0.5, 1.0, 1.0, 0.60817344796939288),
    (2.5, 3.0, 10.0, 0.11903956265827811),
    (3.7, 2.5, 17.3, 0.037971339396467324),
    (1.2, 4.81, 23.67, 0.33906535671931687),
    (0.1, 0.3, 0.7, 0.53612697515608732),
    (0.9, 0.05, 0.02, 0.7087023656787903),
    (4.2, 12.5, 30.25, 0.00060936673002620493),
    (1.7, 15.8, 9.4, 0.20642465544544961),
    (20.0, 3.0, 12.4, 4.8723990430060612e-05),
    (1e5, 1.0, 2.0, 9.9998500024999572e-06),
    (1e-8, 3.0, 7.0, 0.99999999999847589),
    (1e-300, 1.0, 5.0, 1.0),
    (1e-6, 0.5, 2.5, 0.97718816371436956)
  )

  // Extreme upper tails: computed directly, never as 1 - cdf.
  private val extremeTails = List(
    (100.0, 2.2, 31.7, 7.5069447480771803e-15),
    (1000.0, 5.0, 40.0, 7.1609972573169714e-41),
    (50.0, 200.0, 300.0, 5.8222772739641108e-161),
    (30.0, 2.5, 1e6, 1.4242197791693736e-16)
  )

  // Degrees of freedom up to 1e6 (including one very large and one small, non-integer).
  private val largeDf = List(
    (5.0, 1.0, 1e6, 0.02534753835246906),
    (1.3, 7.5, 1e6, 0.24179424801322355),
    (3.0, 1e6, 7.5, 0.053163377358750243),
    (1.001, 1e6, 1e6, 0.30862554957108718),
    (1.01, 1e6, 1e6, 3.2597907372696257e-07),
    (0.99, 1e6, 1e6, 0.99999974848838913),
    (1.05, 1000.0, 1e6, 0.1326006449552202),
    (1.2, 1e6, 1000.5, 4.0612213314532931e-05),
    // Worst case of a 7147-point random R grid (one huge, one small df): about 3e-11 relative.
    (1.2294571910869099, 15.308207047395328, 864263.0223509145, 0.23867126963540616)
  )

  // Chi-square limits, R's pf with df = Inf.
  private val infiniteDf = List(
    (2.0, Double.PositiveInfinity, 7.3, 0.15778729176571696),
    (2.0, 3.5, Double.PositiveInfinity, 0.10110253028564252),
    (0.3, 40.5, Double.PositiveInfinity, 0.99999545976065851),
    (0.6, Double.PositiveInfinity, 12.25, 0.93456429584868206),
    (1.002, 1e6, Double.PositiveInfinity, 0.078718661386129637)
  )

  private def check(cases: List[(Double, Double, Double, Double)], tol: Double): Unit =
    cases.foreach { case (f, df1, df2, expected) =>
      assertRelative(dist(df1, df2).survival(f), expected, tol, s"pf($f, $df1, $df2, lower.tail = FALSE)")
    }

  test("survival matches R pf upper tail for moderate, non-integer degrees of freedom (relative 1e-12)") {
    check(moderate, 1e-12)
  }

  test("extreme upper tails keep relative accuracy (relative 1e-11)") {
    check(extremeTails, 1e-11)
  }

  test("degrees of freedom up to 1e6 (relative 1e-10)") {
    check(largeDf, 1e-10)
  }

  test("infinite degrees of freedom reduce to chi-square tails (relative 1e-10)") {
    check(infiniteDf, 1e-10)
  }

  test("cdf matches R pf lower tail, including small lower tails computed directly (relative 1e-12)") {
    // Rscript -e 'cat(sprintf("%.17g\n", c(pf(2.5, 3, 10), pf(1e-8, 3, 7), pf(0.02, 40, 60), pf(0.3, 12.5, 30.25))))'
    assertRelative(dist(3.0, 10.0).cdf(2.5), 0.88096043734172191, 1e-12, "pf(2.5, 3, 10)")
    assertRelative(dist(3.0, 7.0).cdf(1e-8), 1.5241708451375527e-12, 1e-12, "pf(1e-8, 3, 7)")
    assertRelative(dist(40.0, 60.0).cdf(0.02), 4.7469121437814091e-25, 1e-12, "pf(0.02, 40, 60)")
    assertRelative(dist(12.5, 30.25).cdf(0.3), 0.0139102465176243, 1e-12, "pf(0.3, 12.5, 30.25)")
  }

  test("cdf and survival are complementary") {
    val d = dist(4.81, 23.67)
    List(0.01, 0.5, 1.0, 1.2, 3.0, 9.0).foreach { f =>
      assertEqualsDouble(d.cdf(f) + d.survival(f), 1.0, 1e-15, s"f = $f")
    }
  }

  test("survival is monotone non-increasing across the support") {
    val d = dist(7.5, 1e6)
    val values = (0 to 400).map(i => d.survival(i * 0.01))
    values.sliding(2).foreach { pair => assert(pair(1) <= pair(0), s"$pair") }
  }

  test("boundary and degenerate arguments follow R pf") {
    val d = dist(3.0, 4.0)
    assertEquals(d.survival(-1.0), 1.0)
    assertEquals(d.survival(0.0), 1.0)
    assertEquals(d.survival(Double.PositiveInfinity), 0.0)
    assertEquals(d.cdf(0.0), 0.0)
    assert(d.survival(Double.NaN).isNaN)
    // Rscript -e 'pf(c(0.5, 1, 2), Inf, Inf, lower.tail = FALSE)' gives 1.0 0.5 0.0.
    val point = dist(Double.PositiveInfinity, Double.PositiveInfinity)
    assertEquals(point.survival(0.5), 1.0)
    assertEquals(point.survival(1.0), 0.5)
    assertEquals(point.survival(2.0), 0.0)
    // pf(2, 1e6, 1e6, lower.tail = FALSE) underflows to 0 in R as well.
    assertEquals(dist(1e6, 1e6).survival(2.0), 0.0)
    // df1 * f overflowing Double still yields a tiny, finite tail.
    val far = dist(1e6, 3.0).survival(1e305)
    assert(far >= 0.0 && far < 1e-300, s"$far")
  }

  test("invalid parameters and arguments are rejected") {
    assertEquals(FDistribution(0.0, 1.0), Left(DistributionError.InvalidParameter("df1", 0.0)))
    assertEquals(FDistribution(1.0, -2.0), Left(DistributionError.InvalidParameter("df2", -2.0)))
    assert(FDistribution(Double.NaN, 1.0).isLeft)
    assert(FDistribution.upperTail(Double.NaN, 1.0, 1.0).isLeft)
    assertEquals(
      FDistribution.upperTail(2.5, 3.0, 10.0).map(p => math.abs(p - 0.11903956265827811) < 1e-14),
      Right(true)
    )
  }

  test("incomplete beta tails match R pbeta (relative 1e-12)") {
    // Rscript -e 'cat(sprintf("%.17g\n", c(pbeta(0.3, 2.5, 7.25), pbeta(0.3, 2.5, 7.25, lower.tail = FALSE),
    //   pbeta(0.999, 0.5, 3.5, lower.tail = FALSE))))'
    val tails = IncompleteBeta.regularized(2.5, 7.25, 0.3).fold(e => fail(e.message), identity)
    assertRelative(tails.lower, 0.66048227358190736, 1e-12, "lower")
    assertRelative(tails.upper, 0.33951772641809269, 1e-12, "upper")
    val far = IncompleteBeta.regularized(0.5, 3.5, 0.999).fold(e => fail(e.message), identity)
    assertRelative(far.upper, 9.2066370916053418e-12, 1e-11, "far upper")
    assert(IncompleteBeta.regularized(1.0, 1.0, 1.5).isLeft)
    assert(IncompleteBeta.regularized(0.0, 1.0, 0.5).isLeft)
  }
