package gale.laws

import gale.linalg.*
import gale.solvers.Preconditioner
import gale.spectral.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

class GeneralizedOperatorLawSuite extends ScalaCheckSuite:

  override def scalaCheckInitialSeed =
    "_Xa8Y-7Jwvj_jHt0GpX-fggSVAD0CKUbGaCwwWQsCzN="

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(8).withWorkers(1)

  private final case class DiagonalPencil(
      values: IndexedSeq[Double],
      metric: IndexedSeq[Double]
  ):
    val operatorDiagonal: IndexedSeq[Double] =
      values.zip(metric).map:
        case (value, weight) => value * weight

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

  private def diagonalOperator(diagonal: IndexedSeq[Double]): DoubleLinearOperator =
    LinearOperator.fromFunction(diagonal.length, diagonal.length): (x, into) =>
      var i = 0
      while i < diagonal.length do
        into(i) = diagonal(i) * x(i)
        i += 1

  private def solve(
      operatorDiagonal: IndexedSeq[Double],
      metricDiagonal: IndexedSeq[Double],
      order: EigenOrder = EigenOrder.SmallestAlgebraic
  ): EigenDecomposition =
    val n = operatorDiagonal.length
    val ratios =
      operatorDiagonal.indices.map(i => operatorDiagonal(i) / metricDiagonal(i))
    val selected =
      order match
        case EigenOrder.SmallestAlgebraic => ratios.indices.sortBy(ratios)
        case EigenOrder.LargestAlgebraic  => ratios.indices.sortBy(ratios).reverse
        case other                        => fail(s"unsupported diagonal fixture order: $other")
    val initial = DMat.tabulate(n, 2): (row, col) =>
      if row == selected(col) then 1.0 / math.sqrt(metricDiagonal(row)) else 0.0
    val preconditioner = new Preconditioner:
      def solve(r: DVec, into: MutableVec[Double]): Unit =
        var i = 0
        while i < r.length do
          into(i) = r(i) / operatorDiagonal(i)
          i += 1
    Eigen
      .eigSymmetricGeneralized(
        diagonalOperator(operatorDiagonal).assumeSymmetricOperator,
        diagonalOperator(metricDiagonal).assumePositiveDefiniteOperator,
        n,
        EigenSelection.Count(2, order),
        GeneralizedSpectralOptions(
          tolerance = 1e-8,
          maxIterations = 300,
          initialSubspace = Some(initial)
        ),
        preconditioner
      )
      .fold(error => fail(s"generalized solve failed: $error"), identity)

  property("analytic diagonal oracle and true diagnostic invariants") {
    forAll(pencilGen): pencil =>
      val result = solve(pencil.operatorDiagonal, pencil.metric)
      assert(result.diagnostics.allConverged, result.diagnostics.toString)
      assertEqualsDouble(result.eigenvalues(0), pencil.values(0), 1e-7)
      assertEqualsDouble(result.eigenvalues(1), pencil.values(1), 1e-7)
      GeneralizedOperatorLaws.trueResiduals(
        diagonalOperator(pencil.operatorDiagonal),
        diagonalOperator(pencil.metric),
        result,
        1e-8
      )
      GeneralizedOperatorLaws.bOrthonormal(
        diagonalOperator(pencil.metric),
        result,
        1e-8
      )
  }

  property("permutation and diagonal congruence preserve the spectrum") {
    forAll(pencilGen): pencil =>
      val base = solve(pencil.operatorDiagonal, pencil.metric)
      val permutation = pencil.values.indices.reverse
      val permuted = solve(
        permutation.map(pencil.operatorDiagonal),
        permutation.map(pencil.metric)
      )
      val congruence = IndexedSeq.tabulate(pencil.values.length)(i => 0.5 + i.toDouble)
      val congruent = solve(
        pencil.operatorDiagonal.indices.map(i =>
          pencil.operatorDiagonal(i) * congruence(i) * congruence(i)
        ),
        pencil.metric.indices.map(i =>
          pencil.metric(i) * congruence(i) * congruence(i)
        )
      )

      GeneralizedOperatorLaws.scaledSpectrum(permuted, base, 1.0, 1e-7, 1e-7)
      GeneralizedOperatorLaws.scaledSpectrum(congruent, base, 1.0, 1e-7, 1e-7)
  }

  property("A scaling, B scaling, and common-pencil scaling obey their laws") {
    forAll(pencilGen): pencil =>
      val scale = 3.25
      val base = solve(pencil.operatorDiagonal, pencil.metric)
      val scaledA = solve(pencil.operatorDiagonal.map(_ * scale), pencil.metric)
      val scaledB = solve(pencil.operatorDiagonal, pencil.metric.map(_ * scale))
      val scaledTogether =
        solve(
          pencil.operatorDiagonal.map(_ * scale),
          pencil.metric.map(_ * scale)
        )

      GeneralizedOperatorLaws.scaledSpectrum(scaledA, base, scale, 1e-7, 1e-7)
      GeneralizedOperatorLaws.scaledSpectrum(scaledB, base, 1.0 / scale, 1e-7, 1e-7)
      GeneralizedOperatorLaws.scaledSpectrum(scaledTogether, base, 1.0, 1e-7, 1e-7)
  }

  property("identity metric agrees with ordinary partial symmetric eigen") {
    forAll(pencilGen): pencil =>
      val identity = IndexedSeq.fill(pencil.values.length)(1.0)
      val operator = diagonalOperator(pencil.values)
      val generalized = solve(pencil.values, identity)
      val ordinary = Eigen
        .eigSymmetric(
          operator,
          pencil.values.length,
          EigenSelection.Count(2, EigenOrder.SmallestAlgebraic),
          SpectralOptions(tolerance = 1e-8, maxIterations = 100)
        )
        .toOption
        .get

      assert(generalized.diagnostics.allConverged)
      assert(ordinary.diagnostics.allConverged)
      GeneralizedOperatorLaws.scaledSpectrum(
        generalized,
        ordinary,
        1.0,
        1e-7,
        1e-7
      )
  }
