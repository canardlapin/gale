package gale.golden

import gale.golden.GoldenTolerance.*

/** Mutations of the evidence itself: signs, repeated subspaces and rectangular
  * null completions must not become either false failures or false passes.
  */
class SpectralGoldenEvidenceSuite extends munit.FunSuite:
  private val u = Array(Array(1.0, 0.0), Array(0.0, 1.0), Array(0.0, 0.0))
  private val vt = Array(Array(1.0, 0.0), Array(0.0, 1.0))
  private def ratios(s: Array[Double], left: Array[Array[Double]], right: Array[Array[Double]]): Vector[Double] =
    svdFactorRatios((i, j) => if i == j then s(j) else 0.0, 3, 2, s, left(_)(_), right(_)(_)).map(_._2)
  private def projector(s: Array[Double], left: Array[Array[Double]]): Double =
    svdProjectorRatio(u(_)(_), left(_)(_), s, 3, 3)

  test("paired signs and common repeated-cluster rotations preserve evidence") {
    val flipped = u.map(_.clone())
    flipped(0)(0) = -1.0
    val signed = vt.map(_.clone())
    signed(0)(0) = -1.0
    val s = Array(3.0, 1.0)
    assert(ratios(s, flipped, signed).forall(_ <= 8.0))
    assert(projector(s, flipped) <= 8.0)
    val rotated = Array(Array(0.6, -0.8), Array(0.8, 0.6), Array(0.0, 0.0))
    val inverse = Array(Array(0.6, 0.8), Array(-0.8, 0.6))
    val repeated = Array(3.0, 3.0)
    assert(ratios(repeated, rotated, inverse).forall(_ <= 8.0))
    assert(projector(repeated, rotated) <= 8.0)
    assert(projector(s, rotated) > 8.0, "cross-cluster rotation must fail")
    assert(ratios(s, rotated, inverse).exists(_ > 8.0))
  }

  test("one-sided signs and corrupted trailing rectangular rows fail") {
    val flipped = u.map(_.clone())
    flipped(0)(0) = -1.0
    assert(ratios(Array(3.0, 1.0), flipped, vt).exists(_ > 8.0))
    val corrupt = u.map(_.clone())
    corrupt(2)(0) = 1e-3
    assert(projector(Array(3.0, 1.0), corrupt) > 8.0)
    assert(ratios(Array(3.0, 1.0), corrupt, vt).exists(_ > 8.0))
  }

  test("economy null completions may differ, and the ambient zero affects positive gaps") {
    val other = Array(Array(1.0, 0.0), Array(0.0, 0.0), Array(0.0, 1.0))
    val rankOne = Array(3.0, 0.0)
    assert(projector(rankOne, other) <= 8.0)
    assert(ratios(rankOne, other, vt).forall(_ <= 8.0))
    val angle = 1e-5
    val tilted = Array(Array(1.0, 0.0), Array(0.0, math.sqrt(1.0 - angle * angle)), Array(0.0, angle))
    val s = Array(1.0, 1e-10)
    val expected = angle * math.sqrt(1.0 - angle * angle) / (3.0 * Eps / s(1))
    assert(math.abs(projector(s, tilted) / expected - 1.0) < 1e-12)
  }
