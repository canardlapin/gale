package gale.numeric

class CenteredMomentsSuite extends munit.FunSuite:
  private def ok[E, A](v: Either[E, A]): A = v.fold(e => fail(e.toString), a => a)
  private def state(rows: Vector[Vector[Double]], weights: Vector[Long] = Vector.empty): CenteredMoments =
    val s = ok(CenteredMoments.make(rows.head.size))
    rows.zipWithIndex.foreach((row, i) => ok(s.add(row.toArray, if weights.isEmpty then 1L else weights(i))))
    s
  private def snapshot(s: CenteredMoments, ddof: Int = 1): (Vector[Double], Vector[Double]) =
    val x = ok(s.finish(ddof)); (x.mean.toSeq.toVector, x.covariance.valuesRowMajor.toVector)
  test("simple centered sample/population moments match independent values"):
    val s = state(Vector(Vector(1, 2), Vector(2, 4), Vector(3, 6)))
    val x = ok(s.finish()); assertEquals(x.mean.toSeq.toVector, Vector(2.0, 4.0))
    assertEquals(x.covariance.valuesRowMajor.toVector, Vector(1.0, 2.0, 2.0, 4.0))
    assertEquals(ok(s.finish(0)).covariance(0, 0), 2.0 / 3)
    assertEquals(x.frequency, 3L); assertEquals(x.positiveOccurrences, 3L); assertEquals(x.divisor, 2L)
  test("large offsets preserve a small exact centered covariance"):
    val rows = Vector(Vector(1e16), Vector(1e16 + 2), Vector(1e16 + 4))
    val x = ok(state(rows).finish())
    assertEquals(x.mean(0), 1e16 + 2); assertEquals(x.covariance(0, 0), 4.0)
  test("constant maximum finite inputs do not overflow unnormalized products or means"):
    val s = state(Vector.fill(3)(Vector(Double.MaxValue, -Double.MaxValue)))
    val x = ok(s.finish()); assertEquals(x.mean(0), Double.MaxValue); assertEquals(x.mean(1), -Double.MaxValue)
    assertEquals(x.covariance.valuesRowMajor.toVector, Vector.fill(4)(0.0))
  test("zero and rank deficient directions are preserved without fabricated noise"):
    val x = ok(state(Vector(Vector(0, 1, 2), Vector(0, 2, 4), Vector(0, 3, 6))).finish())
    assertEquals(x.covariance(0, 0), 0.0); assertEquals(x.covariance(1, 1), 1.0)
    assertEquals(x.covariance(1, 2), 2.0); assertEquals(x.covariance(2, 2), 4.0)
  test("integer frequency weights equal actual replication"):
    val rows: Vector[Vector[Double]] = Vector(Vector(1, 3), Vector(2, 4), Vector(5, 1));
    val weights = Vector(2L, 3L, 1L)
    val expanded = rows.zip(weights).flatMap((r, w) => Vector.fill(w.toInt)(r))
    assertEquals(snapshot(state(rows, weights)), snapshot(state(expanded)))
    assertEquals(ok(state(rows, weights).finish()).frequency, 6L)
  test("zero frequency contributes no moments but nonfinite inputs still refuse atomically"):
    val s = state(Vector(Vector(1), Vector(3))); val before = snapshot(s)
    ok(s.add(Array(999.0), 0)); assertEquals(snapshot(s), before)
    assert(s.add(Array(Double.NaN), 0).isLeft); assertEquals(snapshot(s), before)
  test("subnormal means use ties-to-even and covariance underflow is an explicit rounded readout"):
    val tiny = java.lang.Double.longBitsToDouble(1L)
    val x = ok(state(Vector(Vector(tiny), Vector(3 * tiny))).finish())
    assertEquals(x.mean(0), 2 * tiny); assertEquals(x.covariance(0, 0), 0.0)
    assertEquals(ok(state(Vector(Vector(0.0), Vector(tiny))).finish()).mean(0), 0.0)
  test("normalized covariance overflowing Double range refuses without changing state"):
    val s = state(Vector(Vector(Double.MaxValue), Vector(-Double.MaxValue)))
    assert(s.finish().isLeft); assertEquals(s.frequency, 2L)
    ok(s.add(Array(0.0))); assertEquals(s.frequency, 3L)
  test("arbitrary partitions and merge trees yield bit-identical exact state readout"):
    val rows = Vector.tabulate(24)(i => Vector(1e16 + 2 * (i % 7), -1e15 + 0.25 * (i % 11), i.toDouble - 12))
    val expected = snapshot(state(rows))
    val left = state(rows.take(8)); val middle = state(rows.slice(8, 17)); val right = state(rows.drop(17))
    val source = snapshot(middle); ok(left.merge(middle)); ok(left.merge(right));
    assertEquals(snapshot(left), expected); assertEquals(snapshot(middle), source)
    val a = state(rows.take(8)); val b = state(rows.slice(8, 17)); val c = state(rows.drop(17)); ok(b.merge(c));
    ok(a.merge(b)); assertEquals(snapshot(a), expected)
    assertEquals(snapshot(state(rows.reverse)), expected)
  test("copy and owned summaries do not alias later accumulation or source mutation"):
    val input = Array(1.0, 2.0); val s = ok(CenteredMoments.make(2)); ok(s.add(input)); input(0) = 999;
    ok(s.add(Array(3.0, 4.0)))
    val owned = ok(s.finish()); val copied = s.copy(); ok(copied.add(Array(100.0, 100.0)))
    assertEquals(owned.mean.toSeq.toVector, Vector(2.0, 3.0)); assertEquals(s.frequency, 2L)
    ok(s.add(Array(5.0, 6.0))); assertEquals(owned.mean.toSeq.toVector, Vector(2.0, 3.0))
  test("shape, nonfinite, negative frequency and self merge refusals leave the prefix unchanged"):
    val s = state(Vector(Vector(1, 2), Vector(3, 4))); val before = snapshot(s)
    assert(s.add(null).isLeft); assert(s.add(Array(1.0)).isLeft); assert(s.add(Array(1.0, 2.0), -1).isLeft)
    assert(s.add(Array(1.0, Double.PositiveInfinity)).isLeft); assert(s.merge(s).isLeft)
    assert(s.merge(state(Vector(Vector(1), Vector(2)))).isLeft); assertEquals(snapshot(s), before)
  test("maximal frequency has bounded state and overflow refusals are atomic"):
    val s = ok(CenteredMoments.make(1)); ok(s.add(Array(Double.MaxValue), CenteredMoments.MaxFrequency))
    val before = snapshot(s); assert(s.add(Array(0.0)).isLeft); assertEquals(snapshot(s), before)
    val another = ok(CenteredMoments.make(1)); ok(another.add(Array(0.0)))
    assert(s.merge(another).isLeft); assertEquals(snapshot(s), before); assertEquals(another.frequency, 1L)
  test("population/sample insufficiency and unsupported divisors are explicit"):
    val s = ok(CenteredMoments.make(1)); assert(s.finish(0).isLeft); assert(s.finish(1).isLeft)
    ok(s.add(Array(2.0))); assert(s.finish(1).isLeft); assertEquals(ok(s.finish(0)).covariance(0, 0), 0.0)
    assert(s.finish(2).isLeft); assert(s.finish(-1).isLeft)
  test("resource inspection is pure, conservative and shape checked before allocation"):
    assert(CenteredMoments.make(0).isLeft); assert(CenteredMoments.resources(50000).isLeft)
    val r = ok(CenteredMoments.resources(376)); assert(r.stateBytes > 0); assert(r.addScratchBytes > r.stateBytes)
    assert(r.readoutScratchBytes > r.ownedResultBytes); assertEquals(r.copyAdditionalBytes, r.stateBytes)
    val empty = ok(CenteredMoments.make(2)); val before = empty.resources; ok(empty.add(Array(1.0, 2.0)))
    assertEquals(empty.resources, before); assert(before.maximumProductSumBits > before.maximumSumBits)
  test("symmetric rounded covariance has nonnegative independent quadratic forms within declared tolerance"):
    val rows = Vector.tabulate(15)(i => Vector(i.toDouble - 5, (i % 4).toDouble, i.toDouble * 2))
    val c = ok(state(rows).finish()).covariance
    for i <- 0 until 3; j <- 0 until 3 do assertEquals(c(i, j), c(j, i))
    for v <- Vector(Vector(1.0, 0.0, 0.0), Vector(1.0, -2.0, 3.0), Vector(2.0, 0.0, -1.0)) do
      val q = (0 until 3).map(i => (0 until 3).map(j => v(i) * c(i, j) * v(j)).sum).sum
      assert(q >= -1e-12)
