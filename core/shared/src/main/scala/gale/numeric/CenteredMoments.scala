package gale.numeric

import gale.linalg.{DMat, DVec, Matrix, Vec}

/** Worst-case numerical limb/array capacities, not JVM/JS object, allocator, GC or RSS bounds. */
final case class CenteredMomentsResources(
    stateBytes: BigInt,
    copyAdditionalBytes: BigInt,
    addScratchBytes: BigInt,
    mergeScratchBytes: BigInt,
    readoutScratchBytes: BigInt,
    ownedResultBytes: BigInt,
    maximumSumBits: Int,
    maximumProductSumBits: Int,
    maximumReadoutOperandBits: Int
)
enum CenteredMomentsError:
  case InvalidDimension(dimension: Int)
  case ShapeOverflow(dimension: Int)
  case MissingInput
  case LengthMismatch(expected: Int, actual: Int)
  case NonFiniteInput(index: Int)
  case InvalidFrequency(frequency: Long)
  case CapacityExceeded(current: Long, incoming: Long)
  case DimensionMismatch(expected: Int, actual: Int)
  case SelfMerge
  case UnsupportedDivisor(degreesRemoved: Int)
  case InsufficientFrequency(frequency: Long, degreesRemoved: Int)
  case NonFiniteReadout(row: Int, column: Int)
  def message: String = this match
    case InvalidDimension(n)         => s"positive dimension required, got $n"
    case ShapeOverflow(n)            => s"dense output/triangular moment shape exceeds Int storage: $n"
    case MissingInput                => "observation is absent"
    case LengthMismatch(e, a)        => s"observation length $a differs from dimension $e"
    case NonFiniteInput(i)           => s"observation is nonfinite at $i"
    case InvalidFrequency(f)         => s"nonnegative integral frequency required, got $f"
    case CapacityExceeded(a, b)      => s"total frequency exceeds ${CenteredMoments.MaxFrequency}: $a + $b"
    case DimensionMismatch(e, a)     => s"moment dimension $a differs from $e"
    case SelfMerge                   => "self-merge duplicates this exact state"
    case UnsupportedDivisor(d)       => s"only population/sample degrees removed 0 or 1 are admitted, got $d"
    case InsufficientFrequency(n, d) => s"frequency $n does not exceed degrees removed $d"
    case NonFiniteReadout(r, c)      => s"normalized mean/covariance is outside finite Double range at ($r,$c)"

/** Owned immutable mean/covariance for integer frequency replication, not an independence claim. */
final class CenteredMomentsSummary private[numeric] (
    val frequency: Long,
    val positiveOccurrences: Long,
    val degreesRemoved: Int,
    val divisor: Long,
    val mean: DVec,
    val covariance: DMat
)

/** Worker-local exact finite-input moments. Inputs/merge sources are not retained or mutated. Each x is stored in
  * integer units 2^-1074 and products in units 2^-2148. Exact centered covariance is
  * (W*sum(w*x*y)-sum(w*x)*sum(w*y))/(W*(W-ddof)); the apparently subtractive formula uses exact integers, so
  * cancellation never discards bits. Mean/covariance round once, ties to even. Frequency <= MaxFrequency bounds every
  * limb, independently of traversal/partitions. This is a correctness primitive, not a cheap-loop/RSS claim.
  */
final class CenteredMoments private (
    val dimension: Int,
    private var sums: Array[BigInt],
    private var products: Array[BigInt],
    private var mass: Long,
    private var occurrences: Long
):
  import CenteredMomentsError.*
  def frequency: Long = mass
  def positiveOccurrences: Long = occurrences
  def resources: CenteredMomentsResources = CenteredMoments.resources(dimension).toOption.get
  def add(values: Array[Double], frequency: Long = 1L): Either[CenteredMomentsError, Unit] =
    if values == null then Left(MissingInput)
    else if values.length != dimension then Left(LengthMismatch(dimension, values.length))
    else if frequency < 0 then Left(InvalidFrequency(frequency))
    else if frequency > CenteredMoments.MaxFrequency - mass then Left(CapacityExceeded(mass, frequency))
    else
      var row = 0
      while row < dimension do
        if !values(row).isFinite then return Left(NonFiniteInput(row))
        row += 1
      if frequency == 0 then Right(())
      else
        val exact = values.map(ExactDyadic.units)
        val weight = BigInt(frequency)
        val nextSums = new Array[BigInt](dimension)
        val nextProducts = new Array[BigInt](products.length)
        row = 0; var cell = 0
        while row < dimension do
          nextSums(row) = sums(row) + weight * exact(row)
          var column = row
          while column < dimension do
            nextProducts(cell) = products(cell) + weight * exact(row) * exact(column)
            cell += 1; column += 1
          row += 1
        sums = nextSums; products = nextProducts; mass += frequency; occurrences += 1
        Right(())
  def copy(): CenteredMoments = new CenteredMoments(dimension, sums.clone(), products.clone(), mass, occurrences)
  def merge(other: CenteredMoments): Either[CenteredMomentsError, Unit] =
    if this eq other then Left(SelfMerge)
    else if dimension != other.dimension then Left(DimensionMismatch(dimension, other.dimension))
    else if other.mass > CenteredMoments.MaxFrequency - mass then Left(CapacityExceeded(mass, other.mass))
    else
      val nextSums = sums.indices.map(i => sums(i) + other.sums(i)).toArray
      val nextProducts = products.indices.map(i => products(i) + other.products(i)).toArray
      sums = nextSums; products = nextProducts; mass += other.mass; occurrences += other.occurrences
      Right(())
  def finish(degreesRemoved: Int = 1): Either[CenteredMomentsError, CenteredMomentsSummary] =
    if degreesRemoved != 0 && degreesRemoved != 1 then Left(UnsupportedDivisor(degreesRemoved))
    else if mass <= degreesRemoved then Left(InsufficientFrequency(mass, degreesRemoved))
    else
      val weight = BigInt(mass); val meanDenominator = weight << 1074
      val covarianceDenominator = (weight * (mass - degreesRemoved)) << 2148
      val mean = new Array[Double](dimension); val covariance = new Array[Double](dimension * dimension)
      var row = 0; var cell = 0
      while row < dimension do
        mean(row) = ExactDyadic.roundRatio(sums(row), meanDenominator)
        if !mean(row).isFinite then return Left(NonFiniteReadout(row, -1))
        var column = row
        while column < dimension do
          val centered = weight * products(cell) - sums(row) * sums(column)
          val value = ExactDyadic.roundRatio(centered, covarianceDenominator)
          if !value.isFinite then return Left(NonFiniteReadout(row, column))
          covariance(row * dimension + column) = value; covariance(column * dimension + row) = value
          cell += 1; column += 1
        row += 1
      Right(
        new CenteredMomentsSummary(
          mass,
          occurrences,
          degreesRemoved,
          mass - degreesRemoved,
          Vec(mean*),
          Matrix.tabulate(dimension, dimension)((r, c) => covariance(r * dimension + c))
        )
      )

object CenteredMoments:
  val MaxFrequency: Long = 1L << 60
  private val SumBits = 2098 + 60
  private val ProductBits = 4196 + 60
  // Centered numerator < 2^(4196+120+1); ratio rounding may left-shift it by1074.
  private val ReadoutBits = 4196 + 120 + 1 + 1074
  def resources(dimension: Int): Either[CenteredMomentsError, CenteredMomentsResources] =
    if dimension <= 0 then Left(CenteredMomentsError.InvalidDimension(dimension))
    else if dimension.toLong * dimension > Int.MaxValue then Left(CenteredMomentsError.ShapeOverflow(dimension))
    else
      val n = BigInt(dimension); val triangle = n * (n + 1) / 2
      def limbs(bits: Int): BigInt = BigInt((bits + 31) / 32) * 4
      val state = n * limbs(SumBits) + triangle * limbs(ProductBits) + 16
      val scalar = limbs(ReadoutBits)
      val output = 8 * (n + n * n)
      Right(
        CenteredMomentsResources(
          state,
          state,
          state + n * limbs(2098) + 4096 * scalar,
          state + 4096 * scalar,
          4096 * scalar + 2 * output,
          output,
          SumBits,
          ProductBits,
          ReadoutBits
        )
      )
  def make(dimension: Int): Either[CenteredMomentsError, CenteredMoments] =
    resources(dimension).map(_ =>
      new CenteredMoments(
        dimension,
        Array.fill(dimension)(BigInt(0)),
        Array.fill((dimension.toLong * (dimension + 1) / 2).toInt)(BigInt(0)),
        0L,
        0L
      )
    )
