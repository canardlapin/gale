package gale.optim

import gale.linalg.{DMat, LinAlgError, Matrix}
import gale.spectral.{Eigen, EigenSelection, EigenVectors, SpectralBackend}

/** Numerical rotation of an already centered, population-whitened matrix. Preparation, statistical rank selection and
  * artifact identification belong to the caller. Only the portable symmetric eigensolver is admitted here.
  */
final case class TanhIcaPolicy(
    maxIterations: Int = 500,
    tolerance: Double = 1e-7,
    tileSamples: Int = 1024,
    whiteningTolerance: Double = 1e-8,
    orthogonalityTolerance: Double = 1e-8
)

enum TanhIcaFailure:
  case Invalid(detail: String)
  case Budget(requiredBytes: BigInt, availableBytes: BigInt)
  case Cancelled(completedIterations: Int)
  case Numerical(detail: String)
  case Spectral(error: LinAlgError)
  case SingularRotation(iteration: Int, minimumEigenvalue: Double, maximumEigenvalue: Double)
  case NonConverged(iterations: Int, change: Double)

/** Only a converged, independently orthogonality-checked rotation is publishable. */
final class TanhIcaRotation private[optim] (
    val unmixing: DMat,
    val iterations: Int,
    val change: Double,
    val orthogonalityError: Double,
    val policy: TanhIcaPolicy,
    val numericalBytes: BigInt
):
  val spectralBackend: String = "portable-pure"
  val decorrelationRelativeCutoff: Double = 1e-12
  // One input validation pass, followed by one tiled pass per complete update.
  val observationPasses: Long = iterations.toLong + 1

object SymmetricTanhIca:
  /** Conservative primitive payload for the portable path, excluding immutable caller inputs, object headers and GC
    * retention. The 96 square-array slots cover iteration, Gram/product staging, owned output, full symmetric
    * eigendecomposition/canonicalization and its residual certificates. No full samples-by-components
    * projection/nonlinearity arrays are retained.
    */
  def numericalBytes(rank: Int, samples: Int, tileSamples: Int): BigInt =
    val r = BigInt(rank)
    8 * (96 * r * r + 2 * r * math.min(samples, tileSamples) + 12 * r)

  /** Pure preflight for a caller which must reserve before acquiring inputs. */
  def inspect(rank: Int, samples: Int, policy: TanhIcaPolicy): Either[TanhIcaFailure, BigInt] =
    if rank < 2 || samples <= rank || BigInt(rank) * samples > Int.MaxValue then
      Left(TanhIcaFailure.Invalid("bounded rank >= 2 and samples > rank required"))
    else if policy.maxIterations <= 0 || policy.tileSamples <= 0 ||
      !policy.tolerance.isFinite || policy.tolerance <= 0 || policy.tolerance > 0.01 ||
      !policy.whiteningTolerance.isFinite || policy.whiteningTolerance <= 0 || policy.whiteningTolerance > 0.01 ||
      !policy.orthogonalityTolerance.isFinite || policy.orthogonalityTolerance <= 0 || policy.orthogonalityTolerance > 0.01
    then Left(TanhIcaFailure.Invalid("positive iteration/tile limits and finite tolerances in (0,0.01] required"))
    else Right(numericalBytes(rank, samples, policy.tileSamples))

  def fitWhitened(
      observations: DMat,
      initial: DMat,
      numericBudgetBytes: BigInt,
      policy: TanhIcaPolicy = TanhIcaPolicy(),
      cancelled: () => Boolean = () => false
  ): Either[TanhIcaFailure, TanhIcaRotation] =
    val rank = observations.rows
    val samples = observations.cols
    val bytes = inspect(rank, samples, policy) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    if initial.rows != rank || initial.cols != rank then
      return Left(TanhIcaFailure.Invalid("square matched initialization required"))
    if bytes > numericBudgetBytes then return Left(TanhIcaFailure.Budget(bytes, numericBudgetBytes))
    if cancelled() then return Left(TanhIcaFailure.Cancelled(0))
    val sums = new Array[Double](rank)
    val gram = new Array[Double](rank * rank)
    val tile = math.min(samples, policy.tileSamples)
    var start = 0
    while start < samples do
      if cancelled() then return Left(TanhIcaFailure.Cancelled(0))
      var t = start
      val stop = math.min(samples.toLong, start.toLong + tile).toInt
      while t < stop do
        var i = 0
        while i < rank do
          val value = observations(i, t)
          if !value.isFinite then return Left(TanhIcaFailure.Numerical("nonfinite observation"))
          sums(i) += value
          var j = 0
          while j <= i do
            gram(i * rank + j) += value * observations(j, t)
            j += 1
          i += 1
        t += 1
      start = stop
    var i = 0
    while i < rank do
      if !sums(i).isFinite || math.abs(sums(i) / samples) > policy.whiteningTolerance then
        return Left(TanhIcaFailure.Invalid("observations are not centered within the declared tolerance"))
      var j = 0
      while j <= i do
        val covariance = gram(i * rank + j) / samples
        if !covariance.isFinite || math.abs(covariance - (if i == j then 1.0 else 0.0)) > policy.whiteningTolerance then
          return Left(TanhIcaFailure.Invalid("population covariance is not identity within the declared tolerance"))
        j += 1
      i += 1
    var current = decorrelate(initial, 0, policy.orthogonalityTolerance) match
      case Left(error)  => return Left(error)
      case Right(value) => value._1
    val projected = new Array[Double](rank * tile)
    val nonlinear = new Array[Double](rank * tile)
    val products = new Array[Double](rank * rank)
    val derivatives = new Array[Double](rank)
    var iteration = 0
    var change = Double.PositiveInfinity
    while iteration < policy.maxIterations do
      java.util.Arrays.fill(products, 0.0)
      java.util.Arrays.fill(derivatives, 0.0)
      start = 0
      while start < samples do
        if cancelled() then return Left(TanhIcaFailure.Cancelled(iteration))
        val count = math.min(tile, samples - start)
        i = 0
        while i < rank do
          var t = 0
          while t < count do
            var value = 0.0
            var j = 0
            while j < rank do
              value += current(i, j) * observations(j, start + t)
              j += 1
            if !value.isFinite then return Left(TanhIcaFailure.Numerical("nonfinite projection"))
            projected(i * tile + t) = value
            nonlinear(i * tile + t) = math.tanh(projected(i * tile + t))
            t += 1
          i += 1
        // Each scalar accumulator sees observations in the same order for all
        // tile sizes; partitions never change the reduction tree.
        var t = 0
        while t < count do
          i = 0
          while i < rank do
            val g = nonlinear(i * tile + t)
            derivatives(i) += 1.0 - g * g
            var j = 0
            while j < rank do
              products(i * rank + j) += g * observations(j, start + t)
              j += 1
            i += 1
          t += 1
        start += count
      val update = Matrix.tabulate(rank, rank)((row, col) =>
        products(row * rank + col) / samples - derivatives(row) / samples * current(row, col)
      )
      val (next, orthogonality) = decorrelate(update, iteration + 1, policy.orthogonalityTolerance) match
        case Left(error)  => return Left(error)
        case Right(value) => value
      change = 0.0
      i = 0
      while i < rank do
        var dot = 0.0
        var j = 0
        while j < rank do
          dot += next(i, j) * current(i, j)
          j += 1
        change = math.max(change, math.abs(math.abs(dot) - 1.0))
        i += 1
      if !change.isFinite then return Left(TanhIcaFailure.Numerical("nonfinite convergence statistic"))
      iteration += 1
      if cancelled() then return Left(TanhIcaFailure.Cancelled(iteration))
      if change <= policy.tolerance then
        return Right(new TanhIcaRotation(next, iteration, change, orthogonality, policy, bytes))
      current = next
    Left(TanhIcaFailure.NonConverged(iteration, change))

  private def decorrelate(matrix: DMat, iteration: Int, tolerance: Double): Either[TanhIcaFailure, (DMat, Double)] =
    val rank = matrix.rows
    if !(0 until rank).forall(i => (0 until rank).forall(j => matrix(i, j).isFinite)) then
      return Left(TanhIcaFailure.Numerical("nonfinite initialization or fixed-point update"))
    val gram = Matrix.tabulate(rank, rank): (i, j) =>
      var sum = 0.0
      var k = 0
      while k < rank do
        sum += matrix(i, k) * matrix(j, k)
        k += 1
      sum
    if !(0 until rank).forall(i => (0 until rank).forall(j => gram(i, j).isFinite)) then
      return Left(TanhIcaFailure.Numerical("nonfinite rotation Gram matrix"))
    given SpectralBackend = SpectralBackend.none
    Eigen
      .eigSymmetric(gram, EigenSelection.All, EigenVectors.Right)
      .flatMap(_.requireConverged)
      .left
      .map(TanhIcaFailure.Spectral.apply)
      .flatMap: eigen =>
        val minimum = eigen.eigenvalues(0)
        val maximum = eigen.eigenvalues(rank - 1)
        if !minimum.isFinite || !maximum.isFinite || maximum <= 0 || minimum <= 1e-12 * maximum then
          Left(TanhIcaFailure.SingularRotation(iteration, minimum, maximum))
        else
          val inverseRoot = Matrix.tabulate(rank, rank): (i, j) =>
            var sum = 0.0
            var k = 0
            while k < rank do
              sum += eigen.eigenvectors(i, k) / math.sqrt(eigen.eigenvalues(k)) * eigen.eigenvectors(j, k)
              k += 1
            sum
          val result = Matrix.tabulate(rank, rank): (i, j) =>
            var sum = 0.0
            var k = 0
            while k < rank do
              sum += inverseRoot(i, k) * matrix(k, j)
              k += 1
            sum
          var error = 0.0
          var i = 0
          while i < rank do
            var j = 0
            while j <= i do
              var dot = 0.0
              var k = 0
              while k < rank do
                dot += result(i, k) * result(j, k)
                k += 1
              error = math.max(error, math.abs(dot - (if i == j then 1.0 else 0.0)))
              j += 1
            i += 1
          if !error.isFinite || error > tolerance then
            Left(TanhIcaFailure.Numerical("rotation failed orthogonality certificate"))
          else Right((result, error))
