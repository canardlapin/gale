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
