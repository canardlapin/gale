package gale.optim

import gale.linalg.DMat

class StandardTermsSuite extends munit.FunSuite:
  test("L1 and squared L2 have their analytic proximal maps"):
    val l1 = ProximalTerms.l1(2, 2.0).toOption.get
    val l2 = ProximalTerms.squaredL2(2, 3.0).toOption.get
    val x = DMat.tabulate(2, 1)((r, _) => if r == 0 then -3.0 else 0.5)
    val l1Out = l1.proximal(x, 0.5).toOption.get
    val l2Out = l2.proximal(x, 1.0).toOption.get
    assertEqualsDouble(l1Out(0, 0), -2.0, 0.0)
    assertEqualsDouble(l1Out(1, 0), 0.0, 0.0)
    assertEqualsDouble(l2Out(0, 0), -0.75, 0.0)
    assertEqualsDouble(l2Out(1, 0), 0.125, 0.0)

  test("weighted L1 requires exact matrix shape and does not mutate input"):
    val weights = DMat.tabulate(2, 2)((r, c) => (r + c + 1).toDouble)
    val term = ProximalTerms.weightedL1(weights).toOption.get
    val input = DMat.tabulate(2, 2)((r, c) => if r == 0 then -2.0 else c + 0.5)
    val out = term.proximal(input, 0.5).toOption.get
    assertEqualsDouble(out(0, 0), -1.5, 0.0)
    assertEqualsDouble(out(1, 1), 0.0, 0.0)
    assertEqualsDouble(input(0, 0), -2.0, 0.0)
    assert(term.proximal(DMat.zeros(2, 1), 1.0).isLeft)

  test("box indicator is finite exactly on feasible points and projection clamps"):
    val box = ProjectionSets.box(2, -1.0, 2.0).toOption.get
    val input = DMat.tabulate(2, 1)((r, _) => if r == 0 then -3.0 else 4.0)
    val projected = box.project(input).toOption.get
    assertEqualsDouble(projected(0, 0), -1.0, 0.0)
    assertEqualsDouble(projected(1, 0), 2.0, 0.0)
    assert(box.indicator.value(projected).contains(0.0))
    assert(box.indicator.value(input).isLeft)

  test("simplex projection is independently feasible in each column"):
    val simplex = ProjectionSets.simplex(3, 1.0).toOption.get
    val input = DMat.tabulate(3, 2)((r, c) => if c == 0 then Vector(2.0, -1.0, 0.0)(r) else Vector(-2.0, 3.0, 1.0)(r))
    val out = simplex.project(input).toOption.get
    var c = 0
    while c < 2 do
      val sum = out(0, c) + out(1, c) + out(2, c)
      assertEqualsDouble(sum, 1.0, 1e-12)
      assert(out(0, c) >= 0.0 && out(1, c) >= 0.0 && out(2, c) >= 0.0)
      c += 1
    assertEqualsDouble(out(0, 0), 1.0, 1e-12)
    assertEqualsDouble(out(1, 1), 1.0, 1e-12)

  test("Moreau L1 conjugate is clipping onto the dual box"):
    val l1 = ProximalTerms.l1(1, 2.0).toOption.get
    val conjugate = ProximalTerms.conjugate(l1)
    val input = DMat.tabulate(1, 3)((_, c) => Vector(-3.0, 1.0, 4.0)(c))
    val out = conjugate.proximalConjugate(input, 0.5).toOption.get
    assertEqualsDouble(out(0, 0), -2.0, 1e-12)
    assertEqualsDouble(out(0, 1), 1.0, 1e-12)
    assertEqualsDouble(out(0, 2), 2.0, 1e-12)

  test("catalog rejects invalid parameters and non-finite values"):
    assert(ProximalTerms.zero(0).isLeft)
    assert(ProximalTerms.l1(1, -1.0).isLeft)
    assert(ProjectionSets.box(1, 2.0, 1.0).isLeft)
    assert(ProjectionSets.simplex(2, 0.0).isLeft)
    assert(ProjectionSets.nonnegative(1).toOption.get.project(DMat.tabulate(1, 1)((_, _) => Double.NaN)).isLeft)

  test("simplex projection is translation stable and catalog conjugates avoid cancellation"):
    val projected = ProjectionSets.simplex(3).toOption.get.project(DMat.tabulate(3, 1)((_, _) => 1e15)).toOption.get
    for r <- 0 until 3 do assertEqualsDouble(projected(r, 0), 1.0 / 3.0, 1e-15)
    val dual = ProximalTerms
      .conjugate(ProximalTerms.l1(1, 1.0).toOption.get)
      .proximalConjugate(DMat.tabulate(1, 1)((_, _) => 1e16), 1.0)
      .toOption
      .get
    assertEqualsDouble(dual(0, 0), 1.0, 0.0)
    val l2 =
      ProximalTerms.squaredL2(1, 1e100).toOption.get.proximal(DMat.tabulate(1, 1)((_, _) => 1e50), 1e250).toOption.get
    assertEqualsDouble(l2(0, 0) / 1e-300, 1.0, 1e-14)
