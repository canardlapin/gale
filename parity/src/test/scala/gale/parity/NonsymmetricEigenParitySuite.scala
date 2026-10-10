package gale.parity

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.eig as breezeEig
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.spectral.*

/** Nonsymmetric dense eigen parity versus `breeze.linalg.eig`.
  *
  * Breeze does not document eigenvalue order; gale returns a canonical
  * magnitude-ordered layout with conjugate pairs adjacent (positive-imaginary
  * first). Tests therefore compare '''eigenvalue multisets''' and '''eigenpair
  * residuals''' `‖A v − λ v‖`, not raw factor columns. For matrices with
  * distinct real eigenvalues, sign-normalized right eigenvectors are compared
  * after sorting both spectra by λ.
  */
class NonsymmetricEigenParitySuite extends munit.FunSuite:

  private val eigTol = 1e-8
  private val resTol = 1e-8
  private val vecTol = 1e-7

  /** Same frozen seed-2200 LCG fixture as DenseDecompositionBreezeJmh, kept
    * local so the parity module does not depend on the JMH module.
    */
  private def benchmarkFixture(n: Int): Array[Array[Double]] =
    Array.tabulate(n) { i =>
      var state = (2200L + i * 0x9e3779b9L) * 6364136223846793005L + 1442695040888963407L
      Array.tabulate(n) { _ =>
        state = state * 6364136223846793005L + 1442695040888963407L
        ((state >>> 11).toDouble / (1L << 53).toDouble) * 2.0 - 1.0
      }
    }

  /** Compute on A/max|A| so the residual stays meaningful at 1e+/-200. */
  private def normalizedResidual(a: Array[Array[Double]], re: Double, im: Double, vr: Array[Double], vi: Array[Double]): Double =
    val scale = a.iterator.flatMap(_.iterator).map(math.abs).max
    val divisor = if scale == 0.0 then 1.0 else scale
    val lr = re / divisor
    val li = im / divisor
    var residual = 0.0
    var aNorm = 0.0
    var vNorm = 0.0
    for i <- a.indices do
      var ar = 0.0
      var ai = 0.0
      for j <- a.indices do
        val entry = a(i)(j) / divisor
        ar += entry * vr(j)
        ai += entry * vi(j)
        aNorm = math.hypot(aNorm, entry)
      residual = math.hypot(residual, ar - lr * vr(i) + li * vi(i))
      residual = math.hypot(residual, ai - lr * vi(i) - li * vr(i))
      vNorm = math.hypot(vNorm, math.hypot(vr(i), vi(i)))
    assert(vNorm > 0.0 && vNorm.isFinite, "finite nonzero eigenvector")
    residual / ((aNorm + math.hypot(lr, li)) * vNorm)

  private def checkNormalized(data: Array[Array[Double]], clue: String): Unit =
    val n = data.length
    val (wr, wi, v) = breezeParts(data)
    assertEquals(wr.length, n, clue)
    for k <- 0 until n do
      val rc = if wi(k) < 0.0 then k - 1 else k
      val ic = if wi(k) > 0.0 then k + 1 else k
      val br = normalizedResidual(data, wr(k), wi(k), Array.tabulate(n)(v(_, rc)), Array.tabulate(n)(i => if wi(k) == 0.0 then 0.0 else math.signum(wi(k)) * v(i, ic)))
      assert(br <= 64.0 * n * Eps, s"$clue Breeze pair[$k] residual=$br")
    val g = galeAll(data)
    assert(g.diagnostics.allConverged, s"$clue Gale convergence")
    assertEquals(g.size, n, clue)
    for k <- 0 until n do
      val lambda = g.eigenvalue(k)
      val (vr, vi) = g.eigenvector(k)
      val gr = normalizedResidual(data, lambda.re, lambda.im, vr.toArray, vi.toArray)
      assert(gr <= 64.0 * n * Eps, s"$clue Gale pair[$k] residual=$gr")

  test("timed nonsymmetric fixtures n=16/64/256 have normalized residual evidence on both sides") {
    for n <- Seq(16, 64, 256) do checkNormalized(benchmarkFixture(n), s"seed2200 n=$n")
  }

  test("strongly nonnormal and nearly defective systems have small backward error") {
    for coupling <- Seq(1.0, 1e8); gap <- Seq(0.0, 1e-8) do
      val a = Array.tabulate(6, 6)((i, j) => if i == j then 1.0 + i * gap else if j == i + 1 then coupling else 0.0)
      // Forward roots/vectors are ill-conditioned here; compare backward error.
      checkNormalized(a, s"Jordan-like coupling=$coupling gap=$gap")
  }

  test("normal real/complex spectra retain eigenpairs at tiny and huge scales") {
    val n = 6
    val v = Array.tabulate(n)(i => (i + 1).toDouble / math.sqrt(91.0))
    val q = Array.tabulate(n, n)((i, j) => (if i == j then 1.0 else 0.0) - 2.0 * v(i) * v(j))
    val t = Array(Array(1.0, -2.0, 0.0, 0.0, 0.0, 0.0), Array(2.0, 1.0, 0.0, 0.0, 0.0, 0.0), Array(0.0, 0.0, -3.0, 0.0, 0.0, 0.0), Array(0.0, 0.0, 0.0, 5.0, -0.5, 0.0), Array(0.0, 0.0, 0.0, 0.5, 5.0, 0.0), Array(0.0, 0.0, 0.0, 0.0, 0.0, 7.0))
    val normal = Array.tabulate(n, n)((i, j) => (0 until n).flatMap(k => (0 until n).map(l => q(i)(k) * t(k)(l) * q(j)(l))).sum)
    for scale <- Seq(1e-200, 1.0, 1e200) do
      val a = normal.map(_.map(_ * scale))
      checkNormalized(a, s"normal scale=$scale")
      val g = galeValuesOnly(a)
      val (br, bi, _) = breezeParts(a)
      val expectedRe = IndexedSeq(1.0, 1.0, -3.0, 5.0, 5.0, 7.0)
      val expectedIm = IndexedSeq(2.0, -2.0, 0.0, 0.5, -0.5, 0.0)
      matchEigenvalueMultisets((0 until n).map(k => g.eigenvalue(k).re / scale), (0 until n).map(k => g.eigenvalue(k).im / scale), expectedRe, expectedIm, 64.0 * n * Eps, s"analytic Gale scale=$scale")
      matchEigenvalueMultisets(br.map(_ / scale), bi.map(_ / scale), expectedRe, expectedIm, 64.0 * n * Eps, s"analytic Breeze scale=$scale")
  }

  private def galeAll(data: Array[Array[Double]]): NonsymmetricEigenDecomposition =
    Eigen.eigNonsymmetric(galeMatrix(data), EigenSelection.All, EigenVectors.Right).orThrow

  private def galeValuesOnly(data: Array[Array[Double]]): NonsymmetricEigenDecomposition =
    Eigen.eigNonsymmetric(galeMatrix(data), EigenSelection.All, EigenVectors.ValuesOnly).orThrow

  /** Breeze returns `(wr, wi, vr)` with LAPACK real-Schur packing for `vr`. */
  private def breezeParts(data: Array[Array[Double]]): (IndexedSeq[Double], IndexedSeq[Double], BDM[Double]) =
    val e = breezeEig(breezeMatrix(data))
    val wr = (0 until e.eigenvalues.length).map(e.eigenvalues(_))
    val wi = (0 until e.eigenvaluesComplex.length).map(e.eigenvaluesComplex(_))
    (wr, wi, e.eigenvectors)

  private def complexClose(aRe: Double, aIm: Double, bRe: Double, bIm: Double, tol: Double): Boolean =
    val scale = math.max(1.0, math.max(math.hypot(aRe, aIm), math.hypot(bRe, bIm)))
    math.hypot(aRe - bRe, aIm - bIm) <= tol * scale

  /** Greedy multiset match of complex eigenvalues within `tol`. */
  private def matchEigenvalueMultisets(
      gRe: IndexedSeq[Double],
      gIm: IndexedSeq[Double],
      bRe: IndexedSeq[Double],
      bIm: IndexedSeq[Double],
      tol: Double,
      clue: String
  ): Unit =
    assertEquals(gRe.length, bRe.length, clue)
    val used = Array.fill(bRe.length)(false)
    var i = 0
    while i < gRe.length do
      var found = -1
      var j = 0
      while j < bRe.length && found < 0 do
        if !used(j) && complexClose(gRe(i), gIm(i), bRe(j), bIm(j), tol) then found = j
        j += 1
      assert(found >= 0, s"$clue: no Breeze match for gale λ=$i (${gRe(i)}+${gIm(i)}i)")
      used(found) = true
      i += 1

  /** ‖A (vr + i·vi) − λ (vr + i·vi)‖ for a Gale eigenpair. */
  private def galeResidual(a: DMat, d: NonsymmetricEigenDecomposition, i: Int): Double =
    val lambda = d.eigenvalue(i)
    val (vr, vi) = d.eigenvector(i)
    val realPart = (a * vr) - (vr * lambda.re - vi * lambda.im)
    val imagPart = (a * vi) - (vi * lambda.re + vr * lambda.im)
    math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))

  /** ‖A (vr + i·vi) − λ (vr + i·vi)‖ for a Breeze eigenpair in real-Schur packing. */
  private def breezeResidual(a: BDM[Double], wr: IndexedSeq[Double], wi: IndexedSeq[Double], vr: BDM[Double], i: Int): Double =
    val n = a.rows
    val imag = wi(i)
    val (reCol, imCol) =
      if imag == 0.0 then (i, -1)
      else if imag > 0.0 then (i, i + 1)
      else (i - 1, i) // conjugate: packed under the positive-imaginary column
    var realSq = 0.0
    var imagSq = 0.0
    var row = 0
    while row < n do
      var axRe = 0.0
      var axIm = 0.0
      var col = 0
      while col < n do
        val vre = vr(col, reCol)
        val vim = if imCol >= 0 then (if imag > 0.0 then vr(col, imCol) else -vr(col, imCol)) else 0.0
        axRe += a(row, col) * vre
        axIm += a(row, col) * vim
        col += 1
      val vre = vr(row, reCol)
      val vim = if imCol >= 0 then (if imag > 0.0 then vr(row, imCol) else -vr(row, imCol)) else 0.0
      val lamRe = wr(i)
      val lamIm = wi(i)
      val rRe = axRe - (vre * lamRe - vim * lamIm)
      val rIm = axIm - (vim * lamRe + vre * lamIm)
      realSq += rRe * rRe
      imagSq += rIm * rIm
      row += 1
    math.sqrt(realSq + imagSq)

  private def signNormalized(col: IndexedSeq[Double]): IndexedSeq[Double] =
    var idx = 0
    var i = 1
    while i < col.length do
      if math.abs(col(i)) > math.abs(col(idx)) then idx = i
      i += 1
    if col(idx) < 0.0 then col.map(-_) else col

  /** A real diagonalizable matrix `P diag(λ) P⁻¹` with distinct real eigenvalues. */
  private def realDiagonalizable(n: Int, seed: Long): (Array[Array[Double]], Array[Double]) =
    val spectrum = Array.tabulate(n)(i => (i + 1).toDouble * (if i % 2 == 0 then 1.0 else -1.0))
    val p = galeMatrix(diagonallyDominant(n, seed))
    val pInv = p.solve(Matrix.eye(n)).orThrow
    val d = Matrix.tabulate(n, n)((i, j) => if i == j then spectrum(i) else 0.0)
    val a = (p * d) * pInv
    (Array.tabulate(n, n)((i, j) => a(i, j)), spectrum)

  /** A real matrix with at least one conjugate eigenvalue pair (2-D rotation block). */
  private def withComplexPair(seed: Long): Array[Array[Double]] =
    // Block diagonal: 2×2 rotation-scale for 1±2i, plus a real 3.
    val base = Array(
      Array(1.0, -2.0, 0.0),
      Array(2.0, 1.0, 0.0),
      Array(0.0, 0.0, 3.0)
    )
    // Similarity-scramble so eigenvectors are non-trivial.
    val p = galeMatrix(diagonallyDominant(3, seed))
    val pInv = p.solve(Matrix.eye(3)).orThrow
    val a = (p * galeMatrix(base)) * pInv
    Array.tabulate(3, 3)((i, j) => a(i, j))

  // ---------------------------------------------------------------------------
  // Eigenvalue multisets + residuals
  // ---------------------------------------------------------------------------

  test("dense eig All: eigenvalue multisets match; both sides have small residuals") {
    for n <- List(4, 7, 12); seed <- List(1L, 2L, 3L) do
      val data = matrixData(n, n, seed)
      val gd = galeAll(data)
      val (bRe, bIm, bV) = breezeParts(data)
      val gRe = (0 until gd.size).map(gd.eigenvalue(_).re)
      val gIm = (0 until gd.size).map(gd.eigenvalue(_).im)
      matchEigenvalueMultisets(gRe, gIm, bRe, bIm, eigTol, s"eig multiset n=$n seed=$seed")

      val ga = galeMatrix(data)
      val ba = breezeMatrix(data)
      var i = 0
      while i < n do
        val gr = galeResidual(ga, gd, i)
        assert(gr < resTol, s"gale residual n=$n seed=$seed i=$i res=$gr")
        val br = breezeResidual(ba, bRe, bIm, bV, i)
        assert(br < resTol * 10, s"breeze residual n=$n seed=$seed i=$i res=$br")
        i += 1
  }

  test("dense eig ValuesOnly: eigenvalue multiset matches Breeze") {
    for n <- List(5, 9); seed <- List(4L, 5L) do
      val data = matrixData(n, n, seed)
      val gd = galeValuesOnly(data)
      val (bRe, bIm, _) = breezeParts(data)
      val gRe = (0 until gd.size).map(gd.eigenvalue(_).re)
      val gIm = (0 until gd.size).map(gd.eigenvalue(_).im)
      matchEigenvalueMultisets(gRe, gIm, bRe, bIm, eigTol, s"values-only n=$n seed=$seed")
  }

  test("dense eig: conjugate-pair matrix — multiset and residuals") {
    for seed <- List(1L, 2L, 3L) do
      val data = withComplexPair(seed)
      val gd = galeAll(data)
      val (bRe, bIm, bV) = breezeParts(data)
      val gRe = (0 until gd.size).map(gd.eigenvalue(_).re)
      val gIm = (0 until gd.size).map(gd.eigenvalue(_).im)
      matchEigenvalueMultisets(gRe, gIm, bRe, bIm, eigTol, s"complex-pair seed=$seed")
      assert(gIm.exists(_ != 0.0), s"fixture lost its complex pair seed=$seed")

      val ga = galeMatrix(data)
      val ba = breezeMatrix(data)
      var i = 0
      while i < gd.size do
        assert(galeResidual(ga, gd, i) < resTol, s"gale residual complex-pair seed=$seed i=$i")
        assert(breezeResidual(ba, bRe, bIm, bV, i) < resTol * 10, s"breeze residual complex-pair seed=$seed i=$i")
        i += 1
  }

  // ---------------------------------------------------------------------------
  // Distinct real spectrum: sign-aware vectors
  // ---------------------------------------------------------------------------

  test("dense eig: distinct real eigenvalues — sign-normalized eigenvectors match") {
    for n <- List(5, 8); seed <- List(1L, 2L) do
      val (data, spectrum) = realDiagonalizable(n, seed)
      val gd = galeAll(data)
      val (bRe, bIm, bV) = breezeParts(data)

      // All eigenvalues must be real (fixture guarantee).
      assert(bIm.forall(math.abs(_) < 1e-10), s"breeze produced complex λ n=$n seed=$seed")
      assert((0 until gd.size).forall(i => math.abs(gd.eigenvalue(i).im) < 1e-10), s"gale produced complex λ n=$n seed=$seed")

      val gPairs = (0 until n).map(i => (gd.eigenvalue(i).re, i)).sortBy(_._1)
      val bPairs = (0 until n).map(i => (bRe(i), i)).sortBy(_._1)
      gPairs.zip(bPairs).zipWithIndex.foreach { case (((gLam, gIdx), (bLam, bIdx)), k) =>
        assertScalarClose(gLam, bLam, eigTol, s"sorted λ[$k] n=$n seed=$seed")
        assertScalarClose(gLam, spectrum.sorted.apply(k), eigTol * 10, s"vs prescribed λ[$k] n=$n seed=$seed")
        val (gVr, _) = gd.eigenvector(gIdx)
        val gCol = signNormalized((0 until n).map(gVr(_)))
        val bCol = signNormalized((0 until n).map(bV(_, bIdx)))
        var row = 0
        while row < n do
          assertScalarClose(gCol(row), bCol(row), vecTol, s"eigvec n=$n seed=$seed λ[$k] row=$row")
          row += 1
      }
  }
