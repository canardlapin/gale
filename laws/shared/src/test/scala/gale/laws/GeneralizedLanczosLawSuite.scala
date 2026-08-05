package gale.laws

import gale.linalg.*
import gale.spectral.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

class GeneralizedLanczosLawSuite extends ScalaCheckSuite:

  override def scalaCheckInitialSeed =
    "3IpKYLWqvse9f3GvOj9DgO4pqDgfpJ3mSfRFzU9i4yL="

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(8).withWorkers(1)

  private final case class DiagonalPencil(
      values: IndexedSeq[Double],
      metric: IndexedSeq[Double]
  ):
    val operatorDiagonal: IndexedSeq[Double] =
      values.zip(metric).map((value, weight) => value * weight)

  private val pencilGen: Gen[DiagonalPencil] =
    for
      n <- Gen.choose(4, 9)
      seed <- Gen.choose(1L, 1000000L)
    yield
      val random = new scala.util.Random(seed)
      val values =
        IndexedSeq.tabulate(n)(i => i.toDouble + 1.0 + 0.2 * random.nextDouble())
      val metric =
        IndexedSeq.tabulate(n)(_ => math.exp(random.nextDouble() * 4.0 - 2.0))
      DiagonalPencil(values, metric)

  private def diagonalMatrix(diagonal: IndexedSeq[Double]): DMat =
    DMat.tabulate(diagonal.length, diagonal.length): (row, col) =>
      if row == col then diagonal(row) else 0.0

  private def solve(
      operatorDiagonal: IndexedSeq[Double],
      metricDiagonal: IndexedSeq[Double],
      order: EigenOrder = EigenOrder.SmallestAlgebraic
  ): EigenDecomposition =
    val a = diagonalMatrix(operatorDiagonal)
    val b = diagonalMatrix(metricDiagonal)
    val positiveMetric = b.assumePositiveDefinite
    val metricSolve =
      MetricSolveOperator
        .bind(
          positiveMetric,
          LinearSolveOperator.direct(b.cholesky.orThrow)
        )
        .orThrow
    Eigen
      .eigSymmetricGeneralizedLanczos(
        a.assumeSymmetric,
        metricSolve,
        a.rows,
        EigenSelection.Count(2, order),
        GeneralizedLanczosOptions(tolerance = 1e-8, maxIterations = 4)
      )
      .fold(error => fail(s"generalized Lanczos solve failed: $error"), identity)

  property("analytic diagonal oracle and true diagnostic invariants") {
    forAll(pencilGen): pencil =>
      val result = solve(pencil.operatorDiagonal, pencil.metric)
      val a = diagonalMatrix(pencil.operatorDiagonal)
      val b = diagonalMatrix(pencil.metric)

      assert(result.diagnostics.allConverged, result.diagnostics.toString)
      assertEqualsDouble(result.eigenvalues(0), pencil.values(0), 1e-7)
      assertEqualsDouble(result.eigenvalues(1), pencil.values(1), 1e-7)
      assert(result.diagnostics.innerSolve.exists(_.solves > 0))
      GeneralizedOperatorLaws.trueResiduals(a, b, result, 1e-8)
      GeneralizedOperatorLaws.bOrthonormal(b, result, 1e-8)
  }

  property("permutation and diagonal congruence preserve the spectrum") {
    forAll(pencilGen): pencil =>
      val base = solve(pencil.operatorDiagonal, pencil.metric)
      val permutation = pencil.values.indices.reverse
      val permuted = solve(
        permutation.map(pencil.operatorDiagonal),
        permutation.map(pencil.metric)
      )
      val congruence =
        IndexedSeq.tabulate(pencil.values.length)(i => 0.5 + i.toDouble)
      val congruent = solve(
        pencil.operatorDiagonal.indices.map(i => pencil.operatorDiagonal(i) * congruence(i) * congruence(i)),
        pencil.metric.indices.map(i => pencil.metric(i) * congruence(i) * congruence(i))
      )

      GeneralizedOperatorLaws.scaledSpectrum(
        permuted,
        base,
        1.0,
        1e-7,
        1e-7
      )
      GeneralizedOperatorLaws.scaledSpectrum(
        congruent,
        base,
        1.0,
        1e-7,
        1e-7
      )
  }

  property("A, B, and common scaling obey generalized spectral laws") {
    forAll(pencilGen): pencil =>
      val scale = 3.25
      val base = solve(pencil.operatorDiagonal, pencil.metric)
      val scaledA =
        solve(pencil.operatorDiagonal.map(_ * scale), pencil.metric)
      val scaledB =
        solve(pencil.operatorDiagonal, pencil.metric.map(_ * scale))
      val scaledTogether =
        solve(
          pencil.operatorDiagonal.map(_ * scale),
          pencil.metric.map(_ * scale)
        )

      GeneralizedOperatorLaws.scaledSpectrum(
        scaledA,
        base,
        scale,
        1e-7,
        1e-7
      )
      GeneralizedOperatorLaws.scaledSpectrum(
        scaledB,
        base,
        1.0 / scale,
        1e-7,
        1e-7
      )
      GeneralizedOperatorLaws.scaledSpectrum(
        scaledTogether,
        base,
        1.0,
        1e-7,
        1e-7
      )
  }
