package gale.laws

import gale.linalg.*
import gale.spectral.*
import munit.Assertions

/** Reusable laws for the frozen v0.3.5 spectral surface (`gale.spectral`),
  * expressed against the public API and the parity-doc guarantees
  * (`docs/spectral-parity.md`). Each law is a throwing/`Unit` assertion so it can
  * be driven from both munit tests and ScalaCheck property bodies, matching the
  * [[MatrixLaws]] / [[SolverLaws]] convention.
  *
  * The six families: residual (`A v = λ v` and friends), orthogonality, the
  * canonical ordering guarantees, top/bottom membership, rank-deficiency, and the
  * generalized-problem identities.
  */
object SpectralLaws extends Assertions:

  // ===========================================================================
  // Norm helpers
  // ===========================================================================

  /** Frobenius norm — a convenient upper bound on the spectral norm, used to scale
    * the residual tolerances.
    */
  def frobenius(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  private def scaleOf(a: DMat): Double =
    math.max(1.0, frobenius(a))

  // ===========================================================================
  // Residual laws
  // ===========================================================================

  /** Symmetric (dense or the converged iterative pairs): `‖A vᵢ − λᵢ vᵢ‖ ≤
    * relTol·max(1, ‖A‖)` for every returned pair. `a` must be symmetric.
    */
  def symmetricResidual(a: DMat, d: EigenDecomposition, relTol: Double): Unit =
    val bound = relTol * scaleOf(a)
    var i = 0
    while i < d.eigenvectors.cols do
      val v = d.eigenvectors.col(i)
      val res = ((a * v) - v * d.eigenvalues(i)).norm2
      assert(res <= bound, s"symmetric residual $i = $res exceeds $bound")
      i += 1

  /** Nonsymmetric (dense or the converged Arnoldi pairs): the complex residual
    * `‖A(v_re + i·v_im) − λ(v_re + i·v_im)‖`, computed in real arithmetic, is
    * `≤ relTol·max(1, ‖A‖)` for every returned eigenvalue.
    */
  def nonsymmetricResidual(a: DMat, d: NonsymmetricEigenDecomposition, relTol: Double): Unit =
    val bound = relTol * scaleOf(a)
    var i = 0
    while i < d.size do
      val lambda = d.eigenvalue(i)
      val (vr, vi) = d.eigenvector(i)
      val realPart = (a * vr) - (vr * lambda.re - vi * lambda.im)
      val imagPart = (a * vi) - (vi * lambda.re + vr * lambda.im)
      val res = math.sqrt(realPart.dot(realPart) + imagPart.dot(imagPart))
      assert(res <= bound, s"nonsymmetric residual $i = $res exceeds $bound")
      i += 1

  /** SVD (dense-`A` or operator): the two-sided singular-triplet residuals
    * `‖A vᵢ − σᵢ uᵢ‖` and `‖Aᵀ uᵢ − σᵢ vᵢ‖` are both `≤ relTol·max(1, ‖A‖)`.
    */
  def svdResidual(a: DMat, svd: SVD, relTol: Double): Unit =
    val bound = relTol * scaleOf(a)
    val at = a.t
    var i = 0
    while i < svd.singularValues.length do
      val u = svd.u.col(i)
      val v = svd.vt.row(i)
      val sigma = svd.singularValues(i)
      val r1 = ((a * v) - u * sigma).norm2
      val r2 = ((at * u) - v * sigma).norm2
      assert(r1 <= bound, s"‖A v − σ u‖ $i = $r1 exceeds $bound")
      assert(r2 <= bound, s"‖Aᵀ u − σ v‖ $i = $r2 exceeds $bound")
      i += 1

  /** Generalized symmetric-definite: `‖A xᵢ − λᵢ B xᵢ‖ ≤ relTol·max(1, ‖A‖, ‖B‖)`
    * for every returned pair. `a`, `b` must be symmetric.
    */
  def generalizedResidual(a: DMat, b: DMat, d: EigenDecomposition, relTol: Double): Unit =
    val bound = relTol * math.max(scaleOf(a), frobenius(b))
    var i = 0
    while i < d.eigenvectors.cols do
      val x = d.eigenvectors.col(i)
      val res = ((a * x) - (b * x) * d.eigenvalues(i)).norm2
      assert(res <= bound, s"generalized residual $i = $res exceeds $bound")
      i += 1

  /** GSVD reconstruction: `‖A − U C Xᵀ‖ ≤ relTol·max(1, ‖A‖)` and the analogous
    * `‖B − V S Xᵀ‖`.
    */
  def gsvdReconstruction(a: DMat, b: DMat, g: GeneralizedSVD, relTol: Double): Unit =
    val reconA = g.u * diagOf(g.c) * g.x.t
    val reconB = g.v * diagOf(g.s) * g.x.t
    assert(frobenius(a - reconA) <= relTol * scaleOf(a), s"A reconstruction ${frobenius(a - reconA)}")
    assert(frobenius(b - reconB) <= relTol * scaleOf(b), s"B reconstruction ${frobenius(b - reconB)}")

  private def diagOf(d: DVec): DMat =
    DMat.tabulate(d.length, d.length)((i, j) => if i == j then d(i) else 0.0)

  // ===========================================================================
  // Orthogonality laws
  // ===========================================================================

  /** The columns of `m` are orthonormal: `‖MᵀM − I‖_F ≤ tol`. */
  def orthonormalColumns(m: DMat, tol: Double): Unit =
    val err = gramError((0 until m.cols).map(m.col))
    assert(err <= tol, s"columns not orthonormal: ‖MᵀM − I‖ = $err")

  /** The columns of `x` are `B`-orthonormal: `‖Xᵀ B X − I‖_F ≤ tol` (`b`
    * symmetric).
    */
  def bOrthonormal(x: DMat, b: DMat, tol: Double): Unit =
    val g = x.t * (b * x)
    val err = frobenius(g - Matrix.eye(g.rows))
    assert(err <= tol, s"columns not B-orthonormal: ‖XᵀBX − I‖ = $err")

  /** GSVD `U`/`V` orthonormality '''on the well-determined columns only''' — a
    * [[GeneralizedSingularValue.Zero]] leaves its `U` column undetermined (zeroed),
    * an [[GeneralizedSingularValue.Infinite]] its `V` column (the documented
    * contract).
    */
  def gsvdWellDeterminedOrthonormal(g: GeneralizedSVD, tol: Double): Unit =
    val uCols = (0 until g.size).filter(k => g.values(k) != GeneralizedSingularValue.Zero).map(g.u.col)
    val vCols = (0 until g.size).filter(k => g.values(k) != GeneralizedSingularValue.Infinite).map(g.v.col)
    assert(gramError(uCols) <= tol, s"U not orthonormal on well-determined columns: ${gramError(uCols)}")
    assert(gramError(vCols) <= tol, s"V not orthonormal on well-determined columns: ${gramError(vCols)}")

  private def gramError(cols: IndexedSeq[DVec]): Double =
    val k = cols.length
    if k == 0 then 0.0
    else
      var sum = 0.0
      var i = 0
      while i < k do
        var j = 0
        while j < k do
          val d = cols(i).dot(cols(j)) - (if i == j then 1.0 else 0.0)
          sum += d * d
          j += 1
        i += 1
      math.sqrt(sum)

  // ===========================================================================
  // Ordering laws
  // ===========================================================================

  /** Values are ascending within `tol` (symmetric eigenvalue guarantee). */
  def ascending(values: DVec, tol: Double): Unit =
    var i = 1
    while i < values.length do
      assert(values(i) >= values(i - 1) - tol, s"not ascending at $i: ${values(i - 1)} > ${values(i)}")
      i += 1

  /** Values are descending within `tol` (singular value guarantee). */
  def descending(values: DVec, tol: Double): Unit =
    var i = 1
    while i < values.length do
      assert(values(i) <= values(i - 1) + tol, s"not descending at $i: ${values(i - 1)} < ${values(i)}")
      i += 1

  /** Nonsymmetric canonical order: the selection criterion is monotonic in the
    * order's direction (ties allowed), and every conjugate pair is '''adjacent
    * with the positive-imaginary member first''' and exact conjugate symmetry.
    */
  def nonsymmetricOrdering(d: NonsymmetricEigenDecomposition, order: EigenOrder, tol: Double): Unit =
    val largest = isLargest(order)
    var i = 0
    while i < d.size do
      val lam = d.eigenvalue(i)
      // Pair structure.
      if lam.im > 0.0 then
        assert(i + 1 < d.size, s"positive-imag eigenvalue $i has no successor")
        val nxt = d.eigenvalue(i + 1)
        assertEquals(nxt.im, -lam.im, s"pair $i imaginary parts not exact negatives")
        assertEquals(nxt.re, lam.re, s"pair $i real parts differ")
      // Criterion monotonicity.
      if i > 0 then
        val prev = criterion(d.eigenvalue(i - 1), order)
        val cur = criterion(lam, order)
        if largest then assert(prev >= cur - tol, s"criterion not non-increasing at $i: $prev < $cur")
        else assert(prev <= cur + tol, s"criterion not non-decreasing at $i: $prev > $cur")
      i += 1

  /** GSVD ratios are non-increasing — which, since [[GeneralizedSingularValue]]'s
    * value is `+∞` / finite / `0`, is exactly "Infinite first, descending, Zero
    * last".
    */
  def gsvdDescendingRatio(g: GeneralizedSVD): Unit =
    var i = 1
    while i < g.size do
      assert(g.ratio(i - 1) >= g.ratio(i), s"ratio not descending at $i: ${g.ratio(i - 1)} < ${g.ratio(i)}")
      i += 1

  private def isLargest(order: EigenOrder): Boolean =
    order match
      case EigenOrder.LargestMagnitude | EigenOrder.LargestRealPart | EigenOrder.LargestAlgebraic => true
      case _                                                                                       => false

  private def criterion(lambda: Complex, order: EigenOrder): Double =
    order match
      case EigenOrder.LargestMagnitude | EigenOrder.SmallestMagnitude => lambda.magnitude
      case EigenOrder.LargestRealPart | EigenOrder.SmallestRealPart   => lambda.re
      case _                                                          => lambda.re // algebraic (symmetric)

  // ===========================================================================
  // Top/bottom membership laws
  // ===========================================================================

  /** Two real spectra agree as multisets (sorted, within `tol`). */
  def sameSortedValues(actual: Seq[Double], expected: Seq[Double], tol: Double): Unit =
    assertEquals(actual.length, expected.length, s"count mismatch: $actual vs $expected")
    actual.sorted.zip(expected.sorted).foreach: (x, y) =>
      assert(math.abs(x - y) <= tol, s"membership mismatch: $x != $y")

  /** Symmetric `Count(k, order)` selects exactly the `k` extremes of the full
    * ascending spectrum named by `order` (compared against the dense `All` solve).
    */
  def symmetricMembership(full: DVec, subsetValues: Seq[Double], k: Int, order: EigenOrder, tol: Double): Unit =
    val all = (0 until full.length).map(full(_))
    val expected =
      order match
        case EigenOrder.SmallestAlgebraic => all.sorted.take(k)
        case EigenOrder.LargestAlgebraic  => all.sorted.takeRight(k)
        case EigenOrder.SmallestMagnitude => all.sortBy(math.abs).take(k)
        case EigenOrder.LargestMagnitude  => all.sortBy(v => -math.abs(v)).take(k)
        case _                            => all.take(k)
    sameSortedValues(subsetValues, expected, tol)

  /** Singular `Count(k, order)` selects exactly the `k` largest/smallest singular
    * values of the full descending spectrum.
    */
  def singularMembership(full: DVec, subsetValues: Seq[Double], k: Int, order: SingularOrder, tol: Double): Unit =
    val all = (0 until full.length).map(full(_))
    val expected =
      order match
        case SingularOrder.Largest  => all.sorted.reverse.take(k)
        case SingularOrder.Smallest => all.sorted.take(k)
    sameSortedValues(subsetValues, expected, tol)

  /** Nonsymmetric `Count(k, order)`: the '''never-split''' contract — the result
    * holds `k` or `k+1` eigenvalues (a boundary conjugate pair is kept whole), and
    * they are exactly the top-`size` of the full spectrum by the criterion.
    */
  def nonsymmetricMembership(
      full: NonsymmetricEigenDecomposition,
      subset: NonsymmetricEigenDecomposition,
      k: Int,
      order: EigenOrder,
      tol: Double
  ): Unit =
    assert(subset.size == k || subset.size == k + 1, s"expected k or k+1 (=$k/${k + 1}), got ${subset.size}")
    val largest = isLargest(order)
    val fullCrit = (0 until full.size).map(i => criterion(full.eigenvalue(i), order))
    val ranked = if largest then fullCrit.sortBy(c => -c) else fullCrit.sorted
    val expected = ranked.take(subset.size)
    val actual = (0 until subset.size).map(i => criterion(subset.eigenvalue(i), order))
    sameSortedValues(actual, expected, tol)

  /** Each returned value lies within `tol` of some value of the full spectrum — a
    * weaker membership check for the iterative paths, robust to partial
    * convergence.
    */
  def subsetOfSpectrum(subsetValues: Seq[Double], full: DVec, tol: Double): Unit =
    val all = (0 until full.length).map(full(_))
    subsetValues.foreach: v =>
      assert(all.exists(f => math.abs(f - v) <= tol), s"value $v not in the full spectrum")

  // ===========================================================================
  // Rank-deficiency + generalized identities
  // ===========================================================================

  /** `SVD.rank` counts exactly the returned singular values above `tol·σ_max`. */
  def singularRankConsistent(svd: SVD, tol: Double): Unit =
    if svd.size > 0 then
      val sigmaMax = svd.singularValues(0)
      var expected = 0
      var i = 0
      while i < svd.size do
        if svd.singularValues(i) > tol * sigmaMax then expected += 1
        i += 1
      assertEquals(svd.rank, expected, s"rank ${svd.rank} != counted $expected")

  /** GSVD `CᵀC + SᵀS = I`: `c_i² + s_i² = 1` within `tol` for every value. */
  def csIdentity(g: GeneralizedSVD, tol: Double): Unit =
    var i = 0
    while i < g.size do
      val d = g.c(i) * g.c(i) + g.s(i) * g.s(i) - 1.0
      assert(math.abs(d) <= tol, s"c²+s²≠1 at $i: ${g.c(i)}² + ${g.s(i)}²")
      i += 1

  /** A structural precondition failure is the expected `Left`. */
  def isLeftOf(result: Either[LinAlgError, ?], expected: PartialFunction[LinAlgError, Unit]): Unit =
    result match
      case Left(error) if expected.isDefinedAt(error) => ()
      case other                                      => fail(s"expected a specific Left, got $other")
