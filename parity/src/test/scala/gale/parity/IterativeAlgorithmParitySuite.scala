package gale.parity

import gale.linalg.*
import gale.parity.NumpyScipyFixtures.*
import gale.parity.NumpyScipySupport.*
import gale.solvers.*

/** Krylov '''algorithm''' parity versus SciPy `sparse.linalg`, not just
  * replaceability of the computed `x` by a dense direct solve.
  *
  * Breeze 2.1 has no first-class `cg` / `bicgstab` / `gmres` / `lsqr`. The
  * existing [[IterativeSolveParitySuite]] only checks that Gale's iterate
  * matches Breeze `\\`. This suite compares Gale to SciPy's solvers on the
  * same operator, right-hand side, relative residual, and zero start:
  *
  *   - both sides converge;
  *   - solutions agree;
  *   - residuals sit in the same band (the residual each algorithm reports);
  *   - iteration counts stay within a small factor — path-for-path Krylov
  *     vectors are '''not''' required, except unpreconditioned CG, which is
  *     unique and must finish within two steps of SciPy and at most `n + 2`
  *     iterations.
  *
  * Stopping conventions: Gale `ToleranceMode.RelativeToRhs` and SciPy
  * `rtol` / `atol=0` (LSQR `atol=btol=rtol`). That is the R/NumPy family
  * check (`solve` / `lsqr` in R's `Matrix` / `pracma` are the same methods).
  */
class IterativeAlgorithmParitySuite extends munit.FunSuite:

  private val solveTol = 1e-7
  private def configFor(ref: IterativeRef): SolverConfig =
    SolverConfig(tolerance = ref.rtol, maxIterations = 500, restart = 40)

  private def runGale(ref: IterativeRef): SolverResult =
    val b = galeVector(ref.b)
    val cfg = configFor(ref)
    val mode = ToleranceMode.RelativeToRhs
    ref.algorithm match
      case "cg"       => cg(galeCsr(ref.a), b, cfg, toleranceMode = mode)
      case "bicgstab" => bicgstab(galeCsr(ref.a), b, cfg, toleranceMode = mode)
      case "gmres"    => gmres(galeMatrix(ref.a), b, cfg, toleranceMode = mode)
      case "lsqr"     => lsqr(galeMatrix(ref.a), b, cfg, toleranceMode = mode)
      case other      => fail(s"unknown algorithm $other")

  private def assertIterationBand(galeIters: Int, ref: IterativeRef): Unit =
    val n = ref.a.length
    val slack = math.max(5, n / 2)
    val upper = math.max(3 * math.max(ref.iterations, 1), n + slack)
    assert(
      galeIters <= upper,
      s"${ref.name}: Gale iters=$galeIters exceeds band vs SciPy ${ref.iterations} (cap $upper)"
    )
    if ref.algorithm == "cg" then
      assert(galeIters <= n + 2, s"${ref.name}: CG should terminate in ≤ n+2, got $galeIters n=$n")
      assert(
        math.abs(galeIters - ref.iterations) <= 3,
        s"${ref.name}: unpreconditioned CG iters Gale=$galeIters SciPy=${ref.iterations}"
      )

  test("CG / BiCGSTAB / GMRES / LSQR match SciPy sparse.linalg on solution, residual, and iteration band") {
    for ref <- iterative do
      assert(ref.converged, s"${ref.name}: SciPy fixture did not converge")
      val result = runGale(ref)
      assert(
        result.converged,
        s"${ref.name}: Gale did not converge (iters=${result.iterations}, residual=${result.residual})"
      )
      assertVecClose(result.x, ref.x, solveTol, s"${ref.name} solution")
      val a = galeMatrix(ref.a)
      val b = galeVector(ref.b)
      if ref.algorithm == "lsqr" then
        val ar = a.t * (b - (a * result.x))
        assert(
          ar.norm2 <= ref.rtol * (a.t * b).norm2 * 10.0 + 1e-12,
          s"${ref.name}: Gale LSQR ||Aᵀr||=${ar.norm2}"
        )
      else
        val galeR = residual2(a, result.x, b)
        val scipyR = residual2(a, galeVector(ref.x), b)
        val band = math.max(ref.rtol * b.norm2 * 20.0, 1e-8)
        assert(galeR <= band, s"${ref.name}: Gale ||r||=$galeR exceeds $band")
        assert(scipyR <= band, s"${ref.name}: SciPy ||r||=$scipyR exceeds $band")
      assertIterationBand(result.iterations, ref)
  }
