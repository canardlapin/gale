package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.LU as BreezeLU
import breeze.linalg.cholesky
import breeze.linalg.cond
import breeze.linalg.det
import breeze.linalg.eigSym
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.spectral.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAllNoShrink

/** Generated-size companion to the fixed-seed [[FactorizationParitySuite]],
  * [[SvdQrParitySuite]] and [[SpectralParitySuite]] cases (which are kept).
  *
  * Each property draws a size in `1..64` and a data seed; the data come from the
  * same deterministic [[ParitySupport]] builders, so a failure is replayed from
  * the `n=… seed=…` clue alone. Shrinking is disabled (`forAllNoShrink`) so the
  * reported seed is exactly the failing one, and the suite's initial ScalaCheck
  * seed is pinned for reproducible runs. Tolerances follow
  * [[FactorizationHardeningParitySuite]]: `c · n · ε · κ` with `c = 32`
  * (`κ²` for least squares) and the Weyl / Davis–Kahan bounds for `eigSym`.
  */
class FactorizationPropertyParitySuite extends ScalaCheckSuite:
  override def scalaCheckInitialSeed =
    "Hk3Lr6mXq0dQ2s4V8yTz1bWc7nPf5gJ9aE-uR_xKoMN="

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(30).withWorkers(1)

  private val c = 32.0
  private val sizeGen = Gen.choose(1, 64)
  private val seedGen = Gen.choose(1L, 9_000_000L)
  private val caseGen: Gen[(Int, Long)] = Gen.zip(sizeGen, seedGen)

  private def forwardTol(n: Int, kappa: Double): Double = c * math.max(n, 1) * Eps * kappa

  property("LU: det, vector and multi-RHS solve, and P L U = A match breeze") {
    forAllNoShrink(caseGen, Gen.choose(1, 8)) { (sample: (Int, Long), rhsCols: Int) =>
      val (n, seed) = sample
      val data = unitDiagonallyDominant(n, seed)
      val clue = s"n=$n seed=$seed"
      val tol = forwardTol(n, cond(breezeMatrix(data)))
      val lu = galeMatrix(data).lu.orThrow
      val bd = det(breezeMatrix(data))
      assertBelow(math.abs(lu.det.orThrow - bd) / math.abs(bd), tol, s"det $clue")
      val bData = vectorData(n, seed + 1)
      assertBelow(relDiff(lu.solve(galeVector(bData)).orThrow, breezeMatrix(data) \ breezeVector(bData)), tol, s"solve $clue")
      val rhs = matrixData(n, rhsCols, seed + 2)
      assertBelow(
        relDiff(lu.solve(galeMatrix(rhs)).orThrow, breezeMatrix(data) \ breezeMatrix(rhs)),
        tol,
        s"multi-RHS($rhsCols) $clue"
      )
      // gale's (L·U)(i, ·) = A(pivots(i), ·); compare against breeze's P·L·U = A.
      val lower = Matrix.tabulate(n, n)((i, j) => if i == j then 1.0 else if j < i then lu.packed(i, j) else 0.0)
      val upper = Matrix.tabulate(n, n)((i, j) => if j >= i then lu.packed(i, j) else 0.0)
      val permutedA = BDM.tabulate(n, n)((i, j) => data(lu.pivots(i))(j))
      assertBelow(relDiff(lower * upper, permutedA), c * n * Eps, s"gale recon $clue")
      val luB = BreezeLU(breezeMatrix(data))
      assertBelow(relDiff(galeMatrix(data), luB.P * luB.L * luB.U), c * n * Eps, s"breeze recon $clue")
    }
  }

  property("Cholesky: lower factor and solve match breeze") {
    forAllNoShrink(caseGen) { (sample: (Int, Long)) =>
      val (n, seed) = sample
      val data = spd(n, seed)
      val clue = s"n=$n seed=$seed"
      val tol = forwardTol(n, cond(breezeMatrix(data)))
      val factor = galeMatrix(data).cholesky.orThrow
      val bL = cholesky(breezeMatrix(data))
      val gLower = Matrix.tabulate(n, n)((i, j) => if j <= i then factor.lower(i, j) else 0.0)
      assertBelow(relDiff(gLower, BDM.tabulate(n, n)((i, j) => if j <= i then bL(i, j) else 0.0)), tol, s"L $clue")
      val bData = vectorData(n, seed + 1)
      assertBelow(
        relDiff(factor.solve(galeVector(bData)).orThrow, breezeMatrix(data) \ breezeVector(bData)),
        tol,
        s"solve $clue"
      )
    }
  }

  property("QR and least squares on tall m×n: RᵀR = AᵀA, Q R = A, and A \\ b match breeze") {
    forAllNoShrink(caseGen, Gen.choose(1, 64)) { (sample: (Int, Long), extra: Int) =>
      val (n, seed) = sample
      val m = n + extra
      val data = matrixData(m, n, seed)
      val clue = s"${m}x$n seed=$seed"
      val ba = breezeMatrix(data)
      val gqr = galeMatrix(data).qr
      assertBelow(relDiff(gqr.q * gqr.r, ba), c * m * Eps, s"QR recon $clue")
      assertBelow(relDiff(gqr.r.t * gqr.r, ba.t * ba), c * m * Eps, s"RᵀR $clue")
      val k = cond(ba)
      val bData = vectorData(m, seed + 1)
      assertBelow(
        relDiff(galeMatrix(data).leastSquares(galeVector(bData)).orThrow, ba \ breezeVector(bData)),
        c * m * Eps * k * k,
        s"lstsq $clue κ=$k"
      )
    }
  }

  property("eigSym: eigenvalues within the Weyl bound, cluster subspaces within Davis–Kahan") {
    forAllNoShrink(caseGen) { (sample: (Int, Long)) =>
      val (n, seed) = sample
      val data = symmetric(n, seed)
      val clue = s"n=$n seed=$seed"
      val gd = Eigen.eigSymmetric(galeMatrix(data), EigenSelection.All, EigenVectors.Right).orThrow
      val es = eigSym(breezeMatrix(data))
      val bValues = (0 until n).map(es.eigenvalues(_))
      val norm2 = math.max(bValues.map(math.abs).max, Double.MinPositiveValue)
      val valueTol = c * n * Eps * norm2
      for i <- 0 until n do assertBelow(math.abs(gd.eigenvalues(i) - bValues(i)), valueTol, s"λ[$i] $clue")
      for (cluster, sep) <- spectralClusters(bValues, math.sqrt(Eps) * norm2) do
        val bound = math.sqrt(cluster.length.toDouble) * c * n * Eps * (if sep.isInfinite then 1.0 else norm2 / sep)
        assertBelow(subspaceSine(gd.eigenvectors, es.eigenvectors, cluster), bound, s"subspace $cluster $clue")
    }
  }
