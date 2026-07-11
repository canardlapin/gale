package gale.laws

import gale.linalg.*
import gale.spectral.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

/** Conformance suite driving the [[SpectralLaws]] bundle over seeded generators of
  * symmetric / general / rectangular / SPD-pencil inputs. Sizes are kept JS-sane
  * (`n ≤ ~9`) and the ScalaCheck iteration count small, so the suite stays fast on
  * both platforms while covering the residual, orthogonality, ordering,
  * membership, rank-deficiency, and generalized-problem law families.
  */
class SpectralLawSuite extends ScalaCheckSuite:

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(10).withWorkers(1)

  // --- seeded input builders -------------------------------------------------

  private def rng(seed: Long): scala.util.Random =
    new scala.util.Random(seed)

  private def general(n: Int, seed: Long): DMat =
    val r = rng(seed)
    Matrix.tabulate(n, n)((_, _) => r.nextDouble() * 2.0 - 1.0)

  private def symmetric(n: Int, seed: Long): DMat =
    val g = general(n, seed)
    Matrix.tabulate(n, n)((i, j) => 0.5 * (g(i, j) + g(j, i)))

  /** Well-conditioned SPD `M Mᵀ + n·I`. */
  private def spd(n: Int, seed: Long): DMat =
    val m = general(n, seed)
    val mmt = m * m.t
    Matrix.tabulate(n, n)((i, j) => if i == j then mmt(i, j) + n.toDouble else mmt(i, j))

  private def rect(rows: Int, cols: Int, seed: Long): DMat =
    val r = rng(seed)
    Matrix.tabulate(rows, cols)((_, _) => r.nextDouble() * 2.0 - 1.0)

  private def lowRank(m: Int, n: Int, r: Int, seed: Long): DMat =
    val rand = rng(seed)
    val us = Array.fill(r)(DVec.tabulate(m)(_ => rand.nextDouble() * 2.0 - 1.0))
    val vs = Array.fill(r)(DVec.tabulate(n)(_ => rand.nextDouble() * 2.0 - 1.0))
    val scales = Array.tabulate(r)(t => (r - t + 1).toDouble)
    Matrix.tabulate(m, n): (i, j) =>
      var s = 0.0
      var t = 0
      while t < r do
        s += scales(t) * us(t)(i) * vs(t)(j)
        t += 1
      s

  // --- generators ------------------------------------------------------------

  private val dimGen: Gen[Int] = Gen.choose(3, 9)
  private val seedGen: Gen[Long] = Gen.choose(1L, 9_000_000L)

  private val symGen: Gen[DMat] =
    for n <- dimGen; s <- seedGen yield symmetric(n, s)

  private val generalGen: Gen[DMat] =
    for n <- dimGen; s <- seedGen yield general(n, s)

  private val rectGen: Gen[DMat] =
    for
      m <- Gen.choose(3, 8)
      n <- Gen.choose(3, 8)
      s <- seedGen
    yield rect(m, n, s)

  private val pencilGen: Gen[(DMat, DMat)] =
    for
      n <- dimGen
      s1 <- seedGen
      s2 <- seedGen
    yield (symmetric(n, s1), spd(n, s2))

  private val gsvdGen: Gen[(DMat, DMat)] =
    for
      n <- Gen.choose(3, 6)
      mExtra <- Gen.choose(0, 3)
      pExtra <- Gen.choose(0, 3)
      s1 <- seedGen
      s2 <- seedGen
    yield (rect(n + mExtra, n, s1), rect(n + pExtra, n, s2))

  // --- helpers (test-side) ---------------------------------------------------

  private def valuesOf(d: EigenDecomposition): Seq[Double] =
    (0 until d.size).map(d.eigenvalues(_))

  private def referenceSingularValues(a: DMat): DVec =
    val gram = if a.rows >= a.cols then a.t * a else a * a.t
    val eigs = Eigen.eigSymmetric(gram, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get.eigenvalues
    val p = eigs.length
    DVec.tabulate(p)(i => math.sqrt(math.max(0.0, eigs(p - 1 - i))))

  /** Eigenvalues of the naively-reduced `L⁻¹ A L⁻ᵀ` (explicit L-inverse). */
  private def naiveReduction(a: DMat, b: DMat): Seq[Double] =
    val n = a.rows
    val l = b.cholesky.toOption.get.lower
    val aSym = Matrix.tabulate(n, n)((i, j) => if i >= j then a(i, j) else a(j, i))
    val linvCols = (0 until n).map: j =>
      TriangularSolve.lower(l, DVec.tabulate(n)(i => if i == j then 1.0 else 0.0)).toOption.get
    val linv = Matrix.tabulate(n, n)((i, j) => linvCols(j)(i))
    val c = linv * aSym * linv.t
    valuesOf(Eigen.eigSymmetric(c, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get)

  // ===========================================================================
  // Symmetric eigen
  // ===========================================================================

  property("eigSymmetric dense: residual, orthonormal V, ascending values") {
    forAll(symGen) { (a: DMat) =>
      val d = Eigen.eigSymmetric(a, EigenSelection.All).toOption.get
      SpectralLaws.symmetricResidual(a, d, 1e-8)
      SpectralLaws.orthonormalColumns(d.eigenvectors, 1e-8)
      SpectralLaws.ascending(d.eigenvalues, 1e-9)
    }
  }

  property("eigSymmetric dense: Count(k) selects the k algebraic extremes") {
    forAll(symGen) { (a: DMat) =>
      val n = a.rows
      val full = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get.eigenvalues
      val k = math.max(1, n / 2)
      for order <- Seq(EigenOrder.SmallestAlgebraic, EigenOrder.LargestAlgebraic) do
        val sub = Eigen.eigSymmetric(a, EigenSelection.Count(k, order), EigenVectors.ValuesOnly).toOption.get
        SpectralLaws.symmetricMembership(full, valuesOf(sub), k, order, 1e-8)
        SpectralLaws.ascending(sub.eigenvalues, 1e-9)
    }
  }

  property("eigSymmetric Lanczos: converged residual, ascending, subset of spectrum") {
    forAll(symGen) { (a: DMat) =>
      val n = a.rows
      val full = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get.eigenvalues
      val d = Eigen.eigSymmetric(a, n, EigenSelection.Count(2, EigenOrder.LargestAlgebraic)).toOption.get
      SpectralLaws.symmetricResidual(a, d, 1e-6)
      SpectralLaws.ascending(d.eigenvalues, 1e-9)
      SpectralLaws.subsetOfSpectrum(valuesOf(d), full, 1e-6)
    }
  }

  // ===========================================================================
  // Nonsymmetric eigen
  // ===========================================================================

  property("eigNonsymmetric dense All: residual, pairs adjacent, descending magnitude") {
    forAll(generalGen) { (a: DMat) =>
      val d = Eigen.eigNonsymmetric(a, EigenSelection.All).toOption.get
      SpectralLaws.nonsymmetricResidual(a, d, 1e-7)
      SpectralLaws.nonsymmetricOrdering(d, EigenOrder.LargestMagnitude, 1e-9)
    }
  }

  property("eigNonsymmetric dense Count: never-split membership and criterion order") {
    forAll(generalGen) { (a: DMat) =>
      val full = Eigen.eigNonsymmetric(a, EigenSelection.All).toOption.get
      val k = math.max(1, a.rows / 2)
      val orders = Seq(
        EigenOrder.LargestMagnitude,
        EigenOrder.SmallestMagnitude,
        EigenOrder.LargestRealPart,
        EigenOrder.SmallestRealPart
      )
      for order <- orders do
        val sub = Eigen.eigNonsymmetric(a, EigenSelection.Count(k, order)).toOption.get
        SpectralLaws.nonsymmetricMembership(full, sub, k, order, 1e-7)
        SpectralLaws.nonsymmetricOrdering(sub, order, 1e-9)
    }
  }

  property("eigNonsymmetric Arnoldi: converged residual, pairs, subset of spectrum") {
    forAll(generalGen) { (a: DMat) =>
      val n = a.rows
      val full = Eigen.eigNonsymmetric(a, EigenSelection.All).toOption.get
      val d = Eigen.eigNonsymmetric(a, n, EigenSelection.Count(1, EigenOrder.LargestMagnitude)).toOption.get
      SpectralLaws.nonsymmetricResidual(a, d, 1e-6)
      SpectralLaws.nonsymmetricOrdering(d, EigenOrder.LargestMagnitude, 1e-9)
      val fullMags = DVec.tabulate(full.size)(i => full.eigenvalue(i).magnitude)
      SpectralLaws.subsetOfSpectrum((0 until d.size).map(i => d.eigenvalue(i).magnitude), fullMags, 1e-6)
    }
  }

  // ===========================================================================
  // SVD
  // ===========================================================================

  property("svd: two-sided residual, orthonormal U/V, descending, largest membership") {
    forAll(rectGen) { (a: DMat) =>
      val p = math.min(a.rows, a.cols)
      val refSvals = referenceSingularValues(a)
      val scale = math.max(1.0, SpectralLaws.frobenius(a))
      val k = p - 1
      val svd = Svds.svd(a, SingularSelection.Count(k, SingularOrder.Largest)).toOption.get
      SpectralLaws.svdResidual(a, svd, 1e-6)
      SpectralLaws.orthonormalColumns(svd.u, 1e-7)
      SpectralLaws.orthonormalColumns(svd.vt.t, 1e-7)
      SpectralLaws.descending(svd.singularValues, 1e-9)
      SpectralLaws.singularMembership(refSvals, (0 until svd.size).map(svd.singularValues(_)), k, SingularOrder.Largest, 1e-6 * scale)

      // Smallest end: residual + descending, and the single smallest σ is small.
      val small = Svds.svd(a, SingularSelection.Count(1, SingularOrder.Smallest)).toOption.get
      SpectralLaws.svdResidual(a, small, 1e-6)
      SpectralLaws.descending(small.singularValues, 1e-9)
      assert(
        math.abs(small.singularValues(0) - refSvals(p - 1)) <= 1e-4 * scale + 1e-8,
        s"smallest σ ${small.singularValues(0)} vs ref ${refSvals(p - 1)}"
      )
    }
  }

  test("svd: exact-rank-deficient matrix converges and reports rank") {
    val a = lowRank(6, 4, 2, 123457L)
    // k = rank: exactly the two nonzero singular values, all converged.
    val svd2 = Svds.svd(a, SingularSelection.Count(2, SingularOrder.Largest)).toOption.get
    assert(svd2.diagnostics.allConverged, s"rank-2 solve not converged: ${svd2.diagnostics}")
    assertEquals(svd2.rank, 2)
    SpectralLaws.singularRankConsistent(svd2, 1e-8)
    SpectralLaws.svdResidual(a, svd2, 1e-6)
    // k beyond the rank drives the GKL pending-right path: it must still be a Right
    // (a near-zero singular value is legal, never a Left).
    Svds.svd(a, SingularSelection.Count(3, SingularOrder.Largest)) match
      case Right(svd3) => SpectralLaws.singularRankConsistent(svd3, 1e-8)
      case Left(err)   => fail(s"rank-deficient GKL should converge to a Right, got Left($err)")
  }

  // ===========================================================================
  // Generalized symmetric-definite eigen
  // ===========================================================================

  property("eigSymmetricGeneralized: residual, B-orthonormal, ascending, matches reduction") {
    forAll(pencilGen) { (pencil: (DMat, DMat)) =>
      val (a, b) = pencil
      val d = Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All).toOption.get
      SpectralLaws.generalizedResidual(a, b, d, 1e-7)
      SpectralLaws.bOrthonormal(d.eigenvectors, b, 1e-8)
      SpectralLaws.ascending(d.eigenvalues, 1e-9)
      SpectralLaws.sameSortedValues(valuesOf(d), naiveReduction(a, b), 1e-7)
    }
  }

  test("eigSymmetricGeneralized: indefinite B → Left(NotPositiveDefinite)") {
    val a = symmetric(3, 17L)
    val b = Matrix.tabulate(3, 3)((i, j) => if i == j then Array(1.0, -1.0, 1.0)(i) else 0.0)
    SpectralLaws.isLeftOf(
      Eigen.eigSymmetricGeneralized(a, b, EigenSelection.All),
      { case _: LinAlgError.NotPositiveDefinite => () }
    )
  }

  // ===========================================================================
  // Generalized SVD
  // ===========================================================================

  property("gsvd: reconstruction, C²+S²=I, descending ratio, well-determined orthonormal") {
    forAll(gsvdGen) { (pencil: (DMat, DMat)) =>
      val (a, b) = pencil
      Svds.gsvd(a, b) match
        case Left(_: LinAlgError.RankDeficient) => () // rare rank-deficient draw: not in scope
        case Left(other)                        => fail(s"unexpected Left: $other")
        case Right(g) =>
          SpectralLaws.gsvdReconstruction(a, b, g, 1e-6)
          SpectralLaws.csIdentity(g, 1e-10)
          SpectralLaws.gsvdDescendingRatio(g)
          SpectralLaws.gsvdWellDeterminedOrthonormal(g, 1e-7)
    }
  }

  test("gsvd: rank-deficient stacked pencil → Left(RankDeficient)") {
    // A and B share an identical column pair, so [A;B] has rank 2 < 3.
    val a = Matrix.dense(3, 3)(
      1.0, 2.0, 1.0,
      3.0, 1.0, 3.0,
      0.0, 4.0, 0.0
    )
    val b = Matrix.dense(3, 3)(
      2.0, 1.0, 2.0,
      0.0, 5.0, 0.0,
      1.0, 1.0, 1.0
    )
    SpectralLaws.isLeftOf(Svds.gsvd(a, b), { case _: LinAlgError.RankDeficient => () })
  }
