package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.solvers.Preconditioner

/** Tests for the `SpectralBackend` contract skeleton: the `none` fallback (no capabilities, every op
  * `Left(UnsupportedOperation)`) and `compose` (capability union, first-capable-wins dispatch). Test doubles are
  * stateless, matching the thread-safety requirement (G1) — nothing here relies on mutable state.
  */
class SpectralBackendSuite extends munit.FunSuite:

  private val a: DMat = Matrix.tabulate(2, 2)((i, j) => if i == j then 1.0 else 0.0)
  private val b: DMat = Matrix.tabulate(2, 2)((i, j) => if i == j then 2.0 else 0.0)

  // --- none ------------------------------------------------------------------

  test("none: no capabilities and every operation is Left(UnsupportedOperation)") {
    val none = SpectralBackend.none
    assertEquals(none.name, "none")
    assert(none.capabilities.isEmpty)
    assert(none.denseSymmetricEigen(a, true).isLeft)
    assert(none.denseNonsymmetricEigen(a, EigenVectors.Right).isLeft)
    assert(none.denseSvd(a, true).isLeft)
    assert(none.generalizedNonsymmetricEigen(a, b, EigenVectors.Right).isLeft)
    assert(none.rankDeficientGsvd(a, b, true).isLeft)
    assert(
      none
        .iterativeGeneralizedEigen(
          a,
          b,
          2,
          1,
          EigenOrder.SmallestAlgebraic,
          GeneralizedSpectralOptions(),
          Preconditioner.Identity
        )
        .isLeft
    )
    assert(none.shiftInvertSolve(a, None, 0.0).isLeft)
    none.generalizedNonsymmetricEigen(a, b, EigenVectors.Right) match
      case Left(_: LinAlgError.UnsupportedOperation) => ()
      case other                                     => fail(s"expected UnsupportedOperation, got $other")
  }

  test("none is the resolved given when no backend is imported") {
    // summon finds the companion `given none`.
    val resolved = summon[SpectralBackend]
    assertEquals(resolved.name, "none")
    assert(resolved.capabilities.isEmpty)
  }

  // --- stateless test doubles ------------------------------------------------

  /** A backend that "supports" one capability, returning a marker raw value so dispatch can be observed. Stateless.
    */
  private def fakeQz(tag: String): SpectralBackend =
    new SpectralBackend:
      def name: String = tag
      def capabilities: Set[SpectralCapability] = Set(SpectralCapability.GeneralizedNonsymmetricEigen)
      override def generalizedNonsymmetricEigen(
          x: DMat,
          y: DMat,
          v: EigenVectors
      ): Either[LinAlgError, RawGeneralizedEigen] =
        Right(
          RawGeneralizedEigen(
            DVec.tabulate(1)(_ => tag.length.toDouble), // marker: encodes which backend answered
            DVec.zeros(1),
            DVec.tabulate(1)(_ => 1.0),
            DMat.zeros(x.rows, 0),
            None,
            None
          )
        )

  private def fakeSvd(tag: String): SpectralBackend =
    new SpectralBackend:
      def name: String = tag
      def capabilities: Set[SpectralCapability] = Set(SpectralCapability.RankDeficientGsvd)
      override def rankDeficientGsvd(x: DMat, y: DMat, wantVectors: Boolean): Either[LinAlgError, RawGsvd] =
        Right(
          RawGsvd(
            DMat.zeros(x.rows, 0),
            DMat.zeros(y.rows, 0),
            DMat.zeros(x.cols, 0),
            DVec.tabulate(1)(_ => tag.length.toDouble),
            DVec.tabulate(1)(_ => 1.0)
          )
        )

  private def fakeIterative(tag: String): SpectralBackend =
    new SpectralBackend:
      def name: String = tag
      def capabilities: Set[SpectralCapability] =
        Set(SpectralCapability.IterativeGeneralized)
      override def iterativeGeneralizedEigen(
          x: gale.linalg.DoubleLinearOperator,
          y: gale.linalg.DoubleLinearOperator,
          n: Int,
          k: Int,
          order: EigenOrder,
          options: GeneralizedSpectralOptions,
          preconditioner: Preconditioner
      ): Either[LinAlgError, RawIterativeGeneralizedEigen] =
        Right(
          RawIterativeGeneralizedEigen(
            DVec.tabulate(k)(_ => tag.length.toDouble),
            DMat.tabulate(n, k)((row, col) => if row == col then 1.0 else 0.0),
            BackendConvergence(k, k, iterations = tag.length)
          )
        )

  private def fakeShiftSolve(tag: String): SpectralBackend =
    new SpectralBackend:
      def name: String = tag
      def capabilities: Set[SpectralCapability] =
        Set(SpectralCapability.ShiftInvertSolve)
      override def shiftInvertSolve(
          x: DMat,
          y: Option[DMat],
          sigma: Double
      ): Either[LinAlgError, LinearSolveOperator] =
        LinearSolveOperator.backendProvided(x.rows): rhs =>
          Right(
            LinearSolveResult(
              DVec.tabulate(rhs.length)(i => rhs(i) + tag.length.toDouble),
              LinearSolveDiagnostics(true, 0, None, 0L)
            )
          )

  // --- compose ---------------------------------------------------------------

  test("compose: capabilities are the union of the parts") {
    val composite = SpectralBackend.compose(fakeQz("qz"), fakeSvd("svd"))
    assertEquals(
      composite.capabilities,
      Set(SpectralCapability.GeneralizedNonsymmetricEigen, SpectralCapability.RankDeficientGsvd)
    )
    assert(composite.name.contains("qz") && composite.name.contains("svd"))
  }

  test("compose: dispatch routes each op to the first capable part") {
    val composite = SpectralBackend.compose(fakeQz("qz"), fakeSvd("svd"))
    // QZ goes to the QZ part; rank-deficient GSVD to the SVD part.
    assert(composite.generalizedNonsymmetricEigen(a, b, EigenVectors.Right).isRight)
    assert(composite.rankDeficientGsvd(a, b, true).isRight)
    // An unsupported op falls through to Left(UnsupportedOperation).
    assert(composite.denseSvd(a, true).isLeft)
  }

  test("compose: earlier part wins on capability overlap (documented precedence)") {
    // Both parts support QZ; the FIRST (marker length 5) must answer.
    val first = fakeQz("first") // length 5
    val second = fakeQz("secondxxx") // length 9
    val composite = SpectralBackend.compose(first, second)
    val raw = composite.generalizedNonsymmetricEigen(a, b, EigenVectors.Right).toOption.get
    assertEquals(raw.alphaRe(0), 5.0, "first part should have answered (marker 5)")
  }

  test("compose: IterativeGeneralized dispatches to the first capable part") {
    val composite =
      SpectralBackend.compose(fakeIterative("first"), fakeIterative("secondxxx"))
    val raw = composite
      .iterativeGeneralizedEigen(
        a,
        b,
        2,
        1,
        EigenOrder.SmallestAlgebraic,
        GeneralizedSpectralOptions(),
        Preconditioner.Identity
      )
      .toOption
      .get

    assert(composite.capabilities.contains(SpectralCapability.IterativeGeneralized))
    assertEquals(raw.values(0), 5.0)
    assertEquals(raw.convergence.iterations, 5)
  }

  test("compose: ShiftInvertSolve dispatches a typed solve to the first capable part") {
    val composite =
      SpectralBackend.compose(fakeShiftSolve("first"), fakeShiftSolve("secondxxx"))
    val solve = composite.shiftInvertSolve(a, None, 0.5).toOption.get
    val result = solve.solve(DVec.tabulate(2)(i => i.toDouble)).toOption.get

    assertEquals(solve.kind, LinearSolveKind.BackendProvided)
    assertEquals(result.solution(0), 5.0)
    assertEquals(result.solution(1), 6.0)
  }
