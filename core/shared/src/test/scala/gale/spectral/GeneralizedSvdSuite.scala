package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix

/** Tests for the pure generalized SVD `Svds.gsvd(A, B)` — `A = U C Xᵀ`,
  * `B = V S Xᵀ`, `CᵀC + SᵀS = I` — over full-column-rank pencils. Covers the
  * analytic Infinite/Zero construction, reconstruction and the CS identity on
  * random pencils, descending-ratio ordering, the B = I cross-check against the
  * ordinary SVD, and the rank-deficiency / shape `Left`s.
  */
class GeneralizedSvdSuite extends munit.FunSuite:

  // --- helpers ---------------------------------------------------------------

  private def randomMatrix(rows: Int, cols: Int, seed: Long): DMat =
    val rng = new scala.util.Random(seed)
    Matrix.tabulate(rows, cols)((_, _) => rng.nextDouble() * 2.0 - 1.0)

  private def frobenius(a: DMat): Double =
    var sum = 0.0
    var i = 0
    while i < a.rows do
      var j = 0
      while j < a.cols do
        sum += a(i, j) * a(i, j)
        j += 1
      i += 1
    math.sqrt(sum)

  /** `U diag(d) Xᵀ` — the GSVD reconstruction of A (with d = c) or B (with d = s). */
  private def reconstruct(u: DMat, d: DVec, x: DMat): DMat =
    val n = d.length
    val cMat = DMat.tabulate(n, n)((i, j) => if i == j then d(i) else 0.0)
    u * cMat * x.t

  private def assertDescendingRatio(g: GeneralizedSVD): Unit =
    var k = 1
    while k < g.size do
      assert(g.ratio(k - 1) >= g.ratio(k), s"ratios not descending at $k: ${g.ratio(k - 1)} < ${g.ratio(k)}")
      k += 1

  private def columnOrthoError(mat: DMat): Double =
    frobenius(mat.t * mat - Matrix.eye(mat.cols))

  private def diag(d: Array[Double]): DMat =
    DMat.tabulate(d.length, d.length)((i, j) => if i == j then d(i) else 0.0)

  // --- analytic Infinite + Zero ----------------------------------------------

  test("gsvd: constructed pencil recovers Infinite, Finite, and Zero values") {
    // A, B with X = I: value 0 has (c,s)=(1,0) → Infinite; value 1 (0.8,0.6) →
    // Finite(4/3); value 2 (0,1) → Zero. [A;B] stays full column rank (n=3).
    val a = Matrix.dense(4, 3)(
      1.0, 0.0, 0.0,
      0.0, 0.8, 0.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0
    )
    val b = Matrix.dense(4, 3)(
      0.0, 0.0, 0.0,
      0.0, 0.6, 0.0,
      0.0, 0.0, 1.0,
      0.0, 0.0, 0.0
    )
    val g = Svds.gsvd(a, b).toOption.get
    assertEquals(g.size, 3)
    assertEquals(g.values(0), GeneralizedSingularValue.Infinite)
    assertEquals(g.values(2), GeneralizedSingularValue.Zero)
    g.values(1) match
      case GeneralizedSingularValue.Finite(r) => assert(math.abs(r - 4.0 / 3.0) < 1e-8, s"ratio $r")
      case other                              => fail(s"expected Finite, got $other")
    // Cosine/sine diagonals in the descending-ratio order.
    assert(math.abs(g.c(0) - 1.0) < 1e-8 && math.abs(g.s(0)) < 1e-8)
    assert(math.abs(g.c(1) - 0.8) < 1e-8 && math.abs(g.s(1) - 0.6) < 1e-8)
    assert(math.abs(g.c(2)) < 1e-8 && math.abs(g.s(2) - 1.0) < 1e-8)
    assertDescendingRatio(g)
    // Reconstruction holds even with the undetermined (zeroed) U/V columns, since
    // their C/S diagonal entries are 0.
    assert(frobenius(reconstruct(g.u, g.c, g.x) - a) < 1e-8, "A reconstruction")
    assert(frobenius(reconstruct(g.v, g.s, g.x) - b) < 1e-8, "B reconstruction")
  }

  test("gsvd: nontrivial X with an Infinite and a Zero — reordering + reconstruction") {
    // Build from factors with a NON-orthogonal X (so the stacked columns are
    // non-orthonormal and R1 ≠ I, unlike the X = I case above): U = V = [I₃; 0],
    // c = (1, 0.6, 0), s = (0, 0.8, 1) ⇒ value 0 Infinite, value 1 Finite(0.75),
    // value 2 Zero. Degenerate columns must reorder to head/last, and A/B must
    // reconstruct through the nontrivial X with the undetermined columns zeroed.
    val u = Matrix.dense(4, 3)(
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
      0.0, 0.0, 0.0
    )
    val x = Matrix.dense(3, 3)(
      1.0, 0.0, 0.0,
      0.5, 1.0, 0.0,
      0.3, 0.2, 1.0
    )
    val a = u * diag(Array(1.0, 0.6, 0.0)) * x.t
    val b = u * diag(Array(0.0, 0.8, 1.0)) * x.t
    val g = Svds.gsvd(a, b).toOption.get
    assertEquals(g.size, 3)
    assertEquals(g.values.head, GeneralizedSingularValue.Infinite)
    assertEquals(g.values.last, GeneralizedSingularValue.Zero)
    g.values(1) match
      case GeneralizedSingularValue.Finite(r) => assert(math.abs(r - 0.75) < 1e-6, s"ratio $r")
      case other                              => fail(s"expected Finite, got $other")
    // Reconstruction through the nontrivial X (zeroed undetermined columns included).
    val scale = math.max(frobenius(a), frobenius(b))
    assert(frobenius(reconstruct(g.u, g.c, g.x) - a) < 1e-8 * scale, "A reconstruction")
    assert(frobenius(reconstruct(g.v, g.s, g.x) - b) < 1e-8 * scale, "B reconstruction")
  }

  test("gsvd: repeated Infinite values order deterministically (index tie-break)") {
    // Two directions with s = 0 (c = 1) give a repeated ∞ ratio; the finite value
    // sits last. The (−ratio, index) key is a total order, so the layout is fixed.
    val u = Matrix.dense(4, 3)(
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
      0.0, 0.0, 0.0
    )
    val a = u * diag(Array(1.0, 1.0, 0.6)) * Matrix.eye(3).t
    val b = u * diag(Array(0.0, 0.0, 0.8)) * Matrix.eye(3).t
    val g1 = Svds.gsvd(a, b).toOption.get
    assertEquals(g1.values(0), GeneralizedSingularValue.Infinite)
    assertEquals(g1.values(1), GeneralizedSingularValue.Infinite)
    g1.values(2) match
      case GeneralizedSingularValue.Finite(r) => assert(math.abs(r - 0.75) < 1e-6, s"ratio $r")
      case other                              => fail(s"expected Finite, got $other")
    // Deterministic: a second solve yields identical c/s in the same order.
    val g2 = Svds.gsvd(a, b).toOption.get
    var i = 0
    while i < 3 do
      assertEquals(g2.c(i), g1.c(i), s"c$i")
      assertEquals(g2.s(i), g1.s(i), s"s$i")
      i += 1
  }

  // --- random full-rank pencils ----------------------------------------------

  test("gsvd: random pencils reconstruct, satisfy C²+S²=I, and are orthonormal") {
    val configs = Seq((6, 4, 5, 11L), (8, 5, 8, 22L), (5, 3, 4, 33L))
    for (m, n, p, seed) <- configs do
      val a = randomMatrix(m, n, seed)
      val b = randomMatrix(p, n, seed + 1)
      val g = Svds.gsvd(a, b).toOption.get
      assertEquals(g.size, n)
      val scale = math.max(frobenius(a), frobenius(b))
      assert(frobenius(reconstruct(g.u, g.c, g.x) - a) < 1e-8 * scale, s"A recon seed $seed")
      assert(frobenius(reconstruct(g.v, g.s, g.x) - b) < 1e-8 * scale, s"B recon seed $seed")
      // CᵀC + SᵀS = I.
      var i = 0
      while i < n do
        assert(math.abs(g.c(i) * g.c(i) + g.s(i) * g.s(i) - 1.0) < 1e-10, s"c²+s²≠1 at $i")
        i += 1
      assertDescendingRatio(g)
      // Random pencils have all c_i, s_i > 0, so U and V are fully orthonormal.
      assert(columnOrthoError(g.u) < 1e-9, s"U not orthonormal seed $seed")
      assert(columnOrthoError(g.v) < 1e-9, s"V not orthonormal seed $seed")
      assert(g.diagnostics.allConverged)
      assertEquals(g.diagnostics.rank, Some(n))
  }

  // --- B = I cross-check against the ordinary SVD ----------------------------

  test("gsvd: with B = I the generalized values are A's singular values") {
    val m = 6
    val n = 4
    val a = randomMatrix(m, n, 4242L)
    val identity = Matrix.eye(n)
    val g = Svds.gsvd(a, identity, EigenVectors.ValuesOnly).toOption.get
    // All finite (B = I is full rank), descending.
    assertDescendingRatio(g)

    val k = 3
    val svd = Svds.svd(a, SingularSelection.Count(k, SingularOrder.Largest)).toOption.get
    assert(svd.diagnostics.allConverged, s"reference SVD not converged: ${svd.diagnostics}")
    var i = 0
    while i < k do
      assert(math.abs(g.ratio(i) - svd.singularValues(i)) < 1e-6, s"σ$i: ${g.ratio(i)} vs ${svd.singularValues(i)}")
      i += 1
  }

  // --- values-only -----------------------------------------------------------

  test("gsvd: values-only omits U, V, X but keeps the values") {
    val a = randomMatrix(6, 4, 7L)
    val b = randomMatrix(5, 4, 8L)
    val full = Svds.gsvd(a, b).toOption.get
    val vo = Svds.gsvd(a, b, EigenVectors.ValuesOnly).toOption.get
    assertEquals(vo.u.cols, 0)
    assertEquals(vo.v.cols, 0)
    assertEquals(vo.x.cols, 0)
    assertEquals(vo.diagnostics.orthogonalityError, 0.0)
    var i = 0
    while i < full.size do
      assert(math.abs(vo.c(i) - full.c(i)) < 1e-12, s"c$i")
      assert(math.abs(vo.s(i) - full.s(i)) < 1e-12, s"s$i")
      i += 1
  }

  // --- structural Lefts ------------------------------------------------------

  test("gsvd: rank-deficient stacked pencil returns Left(RankDeficient)") {
    // A and B share an identical column pair, so the stacked matrix has rank 2 < 3.
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
    Svds.gsvd(a, b) match
      case Left(_: LinAlgError.RankDeficient) => ()
      case other                              => fail(s"expected RankDeficient, got $other")
  }

  test("gsvd: m + p < n cannot be full column rank") {
    val a = randomMatrix(1, 3, 1L)
    val b = randomMatrix(1, 3, 2L)
    Svds.gsvd(a, b) match
      case Left(_: LinAlgError.RankDeficient) => ()
      case other                              => fail(s"expected RankDeficient, got $other")
  }

  test("gsvd: column-count mismatch returns Left(DimensionMismatch)") {
    val a = randomMatrix(3, 4, 1L)
    val b = randomMatrix(3, 5, 2L)
    Svds.gsvd(a, b) match
      case Left(_: LinAlgError.DimensionMismatch) => ()
      case other                                  => fail(s"expected DimensionMismatch, got $other")
  }

  test("gsvd: one-sided vector flags are rejected") {
    val a = randomMatrix(5, 3, 1L)
    val b = randomMatrix(4, 3, 2L)
    for flag <- Seq(EigenVectors.Left, EigenVectors.LeftAndRight) do
      assert(Svds.gsvd(a, b, flag).isLeft, s"$flag should be rejected")
  }
