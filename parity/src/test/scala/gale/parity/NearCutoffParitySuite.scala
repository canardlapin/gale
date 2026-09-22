package gale.parity

import gale.linalg.*
import gale.parity.NumpyScipyFixtures.*
import gale.parity.NumpyScipySupport.*

/** Near-cutoff rank, `pinv`, and 1-norm `cond` against NumPy / SciPy.
  *
  * Breeze and Gale use different default cutoffs, so the Breeze harness only
  * compares clear full-rank / exactly-deficient cases. This suite pins the
  * boundary against SciPy's MATLAB-compatible conventions:
  *
  *   - `scipy.linalg.pinv` default `rtol = max(m, n) · ε` — the same cutoff
  *     Gale documents for [[gale.linalg.DMat.pinv]];
  *   - `numpy.linalg.matrix_rank` with that same `rtol` (SVD rank);
  *   - `numpy.linalg.cond(A, 1)` as the exact 1-norm condition number.
  *
  * Gale `rankEstimate` is a QR rank with tolerance `2 · max(m, n) · ε · max|R_ii|`,
  * so it can disagree with NumPy's SVD rank on values that sit between the two
  * cutoffs. Clear cases must agree; near-cutoff cases must stay inside the
  * constructed spectral band. `conditionEstimate` is a Hager/Higham '''lower
  * bound''' on the 1-norm condition number, not bit equality with NumPy.
  */
class NearCutoffParitySuite extends munit.FunSuite:

  private val pinvTol = 1e-8
  private val condFactor = 10.0

  test("pinv matches SciPy pinv with the MATLAB/SciPy rcond away from ε-scale σ") {
    // A prescribed SVD with a kept σ ≈ 2·max(m,n)·ε makes A⁺ entries ~1/ε.
    // Independent SVDs then disagree in the 1e12 range even when they keep
    // the same values. The diagonal near-cutoff fixture isolates the policy
    // because that SVD is exact.
    for ref <- nearCutoff if ref.kind != "near_cutoff" do
      val g = galeMatrix(ref.a).pinv.orThrow
      assertMatClose(g, ref.pinv, pinvTol, s"${ref.name} pinv")
  }

  test("rankEstimate agrees with NumPy SVD rank on clear full-rank and exactly-deficient cases") {
    for ref <- nearCutoff if ref.kind == "clear" || ref.kind == "clear_deficient" do
      val galeRank = galeMatrix(ref.a).rankEstimate
      assertEquals(galeRank, ref.numpySvdRank, s"${ref.name}: Gale QR rank vs numpy.linalg.matrix_rank")
      assertEquals(galeRank, ref.definiteRank, s"${ref.name}: definite rank")
  }

  test("near-cutoff rank stays inside the constructed SVD band") {
    for ref <- nearCutoff if ref.kind == "near_cutoff" || ref.kind == "near_cutoff_diag" do
      val galeRank = galeMatrix(ref.a).rankEstimate
      val lo = ref.definiteRank
      val hi = ref.definiteRank + ref.nearCutoffCount
      assert(
        galeRank >= lo && galeRank <= hi,
        s"${ref.name}: Gale rank $galeRank outside [$lo, $hi]"
      )
      assert(
        ref.numpySvdRank >= lo && ref.numpySvdRank <= hi,
        s"${ref.name}: NumPy SVD rank ${ref.numpySvdRank} outside [$lo, $hi]"
      )
  }

  test("conditionEstimate is a 1-norm lower bound on NumPy cond(A, 1)") {
    for ref <- nearCutoff if ref.a.length == ref.a(0).length && !ref.kind.contains("near_cutoff") do
      val galeCond = galeMatrix(ref.a).conditionEstimate.orThrow
      if ref.cond1.isPosInfinity then
        assert(galeCond.isPosInfinity, s"${ref.name}: expected ∞, got $galeCond")
      else if ref.cond1.isFinite then
        assert(
          galeCond <= ref.cond1 * 1.05 + 1e-12,
          s"${ref.name}: Gale 1-norm estimate $galeCond exceeds NumPy cond_1 ${ref.cond1}"
        )
        assert(
          galeCond >= ref.cond1 / condFactor,
          s"${ref.name}: Gale 1-norm estimate $galeCond is more than ${condFactor}× below NumPy ${ref.cond1}"
        )
  }
