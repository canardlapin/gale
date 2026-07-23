package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix

/** The S8 dense-symmetric-eigen dispatch seam on `Eigen.eigSymmetric` (seam S8 of
  * `docs/spectral-backend-boundary.md`), in the `BackendSeamSuite` style: with no
  * acceleration import the `given` resolves to `SpectralBackend.none` and the pure
  * kernel runs — the entire existing symmetric-eigen suite is the byte-identical
  * witness for that; these tests add the ROUTING behaviour. A doubling provider
  * makes a routed call unmistakable; recording counters make "the provider was
  * never consulted" observable; a declining provider isolates the fallback; and a
  * malformed provider proves the loud conformance failure.
  */
class EigSymmetricBackendSeamSuite extends munit.FunSuite:

  private val a: DMat = Matrix.dense(3, 3)(
    4.0, 1.0, 0.5,
    1.0, 3.0, -0.25,
    0.5, -0.25, 2.0
  )

  private def pureEig(
      m: DMat,
      selection: EigenSelection,
      vectors: EigenVectors = EigenVectors.Right
  ): EigenDecomposition =
    Eigen.eigSymmetric(m, selection, vectors)(using SpectralBackend.none).toOption.get

  private def assertSameDecomposition(actual: EigenDecomposition, expected: EigenDecomposition): Unit =
    assertEquals(actual.size, expected.size)
    var i = 0
    while i < expected.size do
      assertEquals(actual.eigenvalues(i), expected.eigenvalues(i), s"eigenvalue $i")
      i += 1
    assertEquals(actual.eigenvectors.shape, expected.eigenvectors.shape)
    var r = 0
    while r < expected.eigenvectors.rows do
      var c = 0
      while c < expected.eigenvectors.cols do
        assertEquals(actual.eigenvectors(r, c), expected.eigenvectors(r, c), s"vector ($r,$c)")
        c += 1
      r += 1

  /** Computes the true spectrum with the pure kernel and DOUBLES the eigenvalues —
    * a routed call is unmistakably distinguishable from the pure result — while
    * recording invocations.
    */
  private final class DoublingEigenProvider(minSize: Int) extends SpectralBackend:
    var calls: Int = 0
    val name: String = "doubling-eigen"
    val capabilities: Set[SpectralCapability] = Set(SpectralCapability.DenseSymmetricEigen)
    override val denseSymmetricEigenMinSize: Int = minSize
    override def denseSymmetricEigen(m: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
      calls += 1
      val kernel = DenseSpectralKernels.symmetricEigen(m, wantVectors).toOption.get
      Right(
        RawSymmetricEigen(
          DVec.tabulate(kernel.values.length)(i => 2.0 * kernel.values(i)),
          kernel.vectors.getOrElse(DMat.zeros(m.rows, 0))
        )
      )

  /** Computes the true spectrum but hands it back in DESCENDING order (values and
    * vector columns in lockstep) — a legal raw carrier whose order the facade must
    * silently correct.
    */
  private final class ReversingEigenProvider extends SpectralBackend:
    val name: String = "reversing-eigen"
    val capabilities: Set[SpectralCapability] = Set(SpectralCapability.DenseSymmetricEigen)
    override def denseSymmetricEigen(m: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
      val kernel = DenseSpectralKernels.symmetricEigen(m, wantVectors).toOption.get
      val n = kernel.values.length
      Right(
        RawSymmetricEigen(
          DVec.tabulate(n)(i => kernel.values(n - 1 - i)),
          kernel.vectors match
            case Some(v) => DMat.tabulate(n, n)((r, c) => v(r, n - 1 - c))
            case None    => DMat.zeros(m.rows, 0)
        )
      )

  /** Advertises the capability but declines every input, so the facade's
    * fallback-to-pure is observable in isolation.
    */
  private final class DecliningEigenProvider extends SpectralBackend:
    var calls: Int = 0
    val name: String = "declining-eigen"
    val capabilities: Set[SpectralCapability] = Set(SpectralCapability.DenseSymmetricEigen)
    override def denseSymmetricEigen(m: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
      calls += 1
      Left(LinAlgError.UnsupportedOperation("declined dense symmetric eigen"))

  /** Returns structurally malformed factors (a conformance violation). */
  private final class MalformedEigenProvider(raw: (DMat, Boolean) => RawSymmetricEigen) extends SpectralBackend:
    val name: String = "malformed-eigen"
    val capabilities: Set[SpectralCapability] = Set(SpectralCapability.DenseSymmetricEigen)
    override def denseSymmetricEigen(m: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
      Right(raw(m, wantVectors))

  // --- resolution ------------------------------------------------------------

  test("no backend import → the companion `none` resolves and the pure kernel computes") {
    val resolved = Eigen.eigSymmetric(a, EigenSelection.All).toOption.get // companion `given none`
    assertSameDecomposition(resolved, pureEig(a, EigenSelection.All))
    assert(!SpectralBackend.none.routesDenseSymmetricEigen(a.rows), "none must never route")
  }

  test("a `given` capable backend in scope routes a plain call (the ergonomic pattern)") {
    val provider = DoublingEigenProvider(minSize = 0)
    given SpectralBackend = provider // shadows the companion `given none`
    val routed = Eigen.eigSymmetric(a, EigenSelection.All).toOption.get
    val pure = pureEig(a, EigenSelection.All)
    assertEquals(provider.calls, 1)
    var i = 0
    while i < pure.size do
      assertEqualsDouble(routed.eigenvalues(i), 2.0 * pure.eigenvalues(i), 1e-12, s"eigenvalue $i")
      i += 1
  }

  // --- threshold gate ---------------------------------------------------------

  test("above threshold the facade routes to the provider (doubled ⇒ it routed)") {
    val provider = DoublingEigenProvider(minSize = a.rows)
    val routed = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly)(using provider).toOption.get
    val pure = pureEig(a, EigenSelection.All, EigenVectors.ValuesOnly)
    assertEquals(provider.calls, 1)
    var i = 0
    while i < pure.size do
      assertEqualsDouble(routed.eigenvalues(i), 2.0 * pure.eigenvalues(i), 1e-12, s"eigenvalue $i")
      i += 1
  }

  test("below threshold the facade stays on the pure path and never consults the provider") {
    val provider = DoublingEigenProvider(minSize = a.rows + 1)
    val notRouted = Eigen.eigSymmetric(a, EigenSelection.All)(using provider).toOption.get
    assertEquals(provider.calls, 0)
    assertSameDecomposition(notRouted, pureEig(a, EigenSelection.All)) // NOT doubled ⇒ pure
  }

  test("routesDenseSymmetricEigen is the single gate: capability AND size") {
    val provider = DoublingEigenProvider(minSize = 5)
    assert(!provider.routesDenseSymmetricEigen(4))
    assert(provider.routesDenseSymmetricEigen(5))
    assert(!SpectralBackend.none.routesDenseSymmetricEigen(Int.MaxValue))
  }

  test("compose: the routing threshold follows the part that serves the capability") {
    val provider = DoublingEigenProvider(minSize = 5)
    val composite = SpectralBackend.compose(SpectralBackend.none, provider)
    assertEquals(composite.denseSymmetricEigenMinSize, 5)
    assert(composite.routesDenseSymmetricEigen(5))
    assert(!composite.routesDenseSymmetricEigen(4))
    // No capable part: the gate is closed regardless of n.
    val empty = SpectralBackend.compose(SpectralBackend.none)
    assert(!empty.routesDenseSymmetricEigen(Int.MaxValue))
  }

  // --- validation precedes the provider ---------------------------------------

  test("validation runs BEFORE the gate: structural Lefts never reach the provider") {
    val provider = DoublingEigenProvider(minSize = 0)
    val rect = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    Eigen.eigSymmetric(rect, EigenSelection.All)(using provider) match
      case Left(_: LinAlgError.NonSquareMatrix) => ()
      case other                                => fail(s"expected Left(NonSquareMatrix), got $other")
    Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.LeftAndRight)(using provider) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected Left(InvalidArgument), got $other")
    Eigen.eigSymmetric(a, EigenSelection.Count(0, EigenOrder.LargestAlgebraic))(using provider) match
      case Left(_: LinAlgError.InvalidArgument) => ()
      case other                                => fail(s"expected Left(InvalidArgument), got $other")
    assertEquals(provider.calls, 0)
  }

  // --- fallback and loud failure ----------------------------------------------

  test("a provider Left is a decline: the pure kernel computes the identical answer") {
    val provider = DecliningEigenProvider()
    val viaDeclining = Eigen.eigSymmetric(a, EigenSelection.All)(using provider)
    assertEquals(provider.calls, 1) // consulted, declined
    assert(viaDeclining.isRight, "a decline must not surface as a failure")
    assertSameDecomposition(viaDeclining.toOption.get, pureEig(a, EigenSelection.All))
  }

  test("malformed provider factors fail loudly (conformance violation, never silent)") {
    // Wrong eigenvalue count.
    val wrongCount = MalformedEigenProvider((m, _) => RawSymmetricEigen(DVec.zeros(m.rows - 1), DMat.zeros(m.rows, m.rows)))
    intercept[LinAlgError.InvalidArgument] {
      Eigen.eigSymmetric(a, EigenSelection.All)(using wrongCount)
    }
    // Vectors requested but too few columns returned.
    val tooFewColumns = MalformedEigenProvider((m, _) => RawSymmetricEigen(DVec.zeros(m.rows), DMat.zeros(m.rows, m.rows - 1)))
    intercept[LinAlgError.InvalidArgument] {
      Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.Right)(using tooFewColumns)
    }
    // Wrong vector row count.
    val wrongRows = MalformedEigenProvider((m, _) => RawSymmetricEigen(DVec.zeros(m.rows), DMat.zeros(m.rows - 1, m.rows)))
    intercept[LinAlgError.InvalidArgument] {
      Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.Right)(using wrongRows)
    }
  }

  // --- facade-owned canonicalization -------------------------------------------

  test("the facade re-imposes the ascending order on a backend's raw output") {
    val routed = Eigen.eigSymmetric(a, EigenSelection.All)(using ReversingEigenProvider()).toOption.get
    assertSameDecomposition(routed, pureEig(a, EigenSelection.All)) // silently corrected, not passed through
    assert(routed.diagnostics.worstResidual < 1e-10, "facade-re-derived residuals must be honest")
    assert(routed.diagnostics.orthogonalityError < 1e-10)
  }

  test("selection is realized on the canonicalized spectrum, not the provider's order") {
    val provider = ReversingEigenProvider()
    assertSameDecomposition(
      Eigen.eigSymmetric(a, EigenSelection.Count(1, EigenOrder.SmallestAlgebraic))(using provider).toOption.get,
      pureEig(a, EigenSelection.Count(1, EigenOrder.SmallestAlgebraic))
    )
    assertSameDecomposition(
      Eigen.eigSymmetric(a, EigenSelection.IndexRange(1, 2))(using provider).toOption.get,
      pureEig(a, EigenSelection.IndexRange(1, 2))
    )
  }

  test("excess provider vector columns are tolerated: the leading n columns are the factors") {
    // The RawSymmetricEigen contract mirrors RawSvd's slack: AT LEAST the leading
    // n columns; a wider matrix (extra workspace columns) is sliced, not rejected.
    val padded = new SpectralBackend:
      val name: String = "padded-eigen"
      val capabilities: Set[SpectralCapability] = Set(SpectralCapability.DenseSymmetricEigen)
      override def denseSymmetricEigen(m: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
        val kernel = DenseSpectralKernels.symmetricEigen(m, wantVectors).toOption.get
        val n = m.rows
        val v = kernel.vectors.get
        Right(
          RawSymmetricEigen(
            kernel.values,
            DMat.tabulate(n, n + 2)((r, c) => if c < n then v(r, c) else 99.0)
          )
        )
    val routed = Eigen.eigSymmetric(a, EigenSelection.All)(using padded).toOption.get
    assertSameDecomposition(routed, pureEig(a, EigenSelection.All))
  }

  test("values-only routing carries an empty vector matrix through cleanly") {
    val provider = DoublingEigenProvider(minSize = 0)
    val routed = Eigen.eigSymmetric(a, EigenSelection.All, EigenVectors.ValuesOnly)(using provider).toOption.get
    assertEquals(routed.eigenvectors.cols, 0)
    assertEquals(provider.calls, 1)
  }
