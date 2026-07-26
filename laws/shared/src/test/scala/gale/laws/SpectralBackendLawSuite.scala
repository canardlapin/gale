package gale.laws

import gale.linalg.*
import gale.solvers.Preconditioner
import gale.spectral.*

class SpectralBackendLawSuite extends munit.FunSuite:

  private def provider(
      raw: RawIterativeGeneralizedEigen
  ): SpectralBackend =
    new SpectralBackend:
      def name: String = "law-fixture"
      def capabilities: Set[SpectralCapability] =
        Set(SpectralCapability.IterativeGeneralized)
      override def iterativeGeneralizedEigen(
          a: DoubleLinearOperator,
          b: DoubleLinearOperator,
          n: Int,
          k: Int,
          order: EigenOrder,
          options: GeneralizedSpectralOptions,
          preconditioner: Preconditioner
      ): Either[LinAlgError, RawIterativeGeneralizedEigen] =
        Right(raw)

  test("conforming IterativeGeneralized provider satisfies the reusable law") {
    val raw = RawIterativeGeneralizedEigen(
      DVec.fromSeq(Seq(2.0, 1.0)),
      DMat.tabulate(4, 2): (row, col) =>
        if (row == 1 && col == 0) || (row == 0 && col == 1) then
          (col + 2).toDouble
        else 0.0,
      BackendConvergence(requested = 2, converged = 2, iterations = 3)
    )
    SpectralBackendLaws.iterativeGeneralizedCapability(provider(raw))
  }

  test("advertised but unimplemented IterativeGeneralized capability is rejected") {
    val unimplemented = new SpectralBackend:
      def name: String = "unimplemented"
      def capabilities: Set[SpectralCapability] =
        Set(SpectralCapability.IterativeGeneralized)

    intercept[AssertionError]:
      SpectralBackendLaws.iterativeGeneralizedCapability(unimplemented)
  }

  test("malformed iterative generalized carrier is rejected") {
    val malformed = RawIterativeGeneralizedEigen(
      DVec.fromSeq(Seq(1.0)),
      DMat.zeros(3, 1),
      BackendConvergence(requested = 2, converged = 2, iterations = 1)
    )
    intercept[AssertionError]:
      SpectralBackendLaws.iterativeGeneralizedCapability(provider(malformed))
  }

  test("false provider convergence is rejected by facade-derived residuals") {
    val falseConvergence = RawIterativeGeneralizedEigen(
      DVec.fromSeq(Seq(999.0)),
      DMat.tabulate(4, 1)((row, _) => if row == 0 then 1.0 else 0.0),
      BackendConvergence(requested = 2, converged = 1, iterations = 1)
    )
    intercept[AssertionError]:
      SpectralBackendLaws.iterativeGeneralizedCapability(
        provider(falseConvergence)
      )
  }
