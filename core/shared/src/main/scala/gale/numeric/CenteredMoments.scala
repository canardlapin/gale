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

/** Worker-local exact finite-input moments. Fixed signed radix30 limbs represent absolute units2^-1074/2^-2148. Total
  * frequency<=2^60 bounds sums below2^2158 and products below2^4256, with spare sign bits for exact modular
  * carry/borrow. Exact covariance remains (W*sum(w*x*y)-sum(w*x)*sum(w*y))/(W*(W-ddof)). Only final readout constructs
  * BigInts and rounds once, ties to even.
  */
final class CenteredMoments private (
    val dimension: Int,
    private val sums: FixedMomentLimbs,
    private val products: FixedMomentLimbs,
    private var mass: Long,
    private var occurrences: Long
):
  import CenteredMomentsError.*
  private val mantissas = new Array[Long](dimension)
  private val shifts = new Array[Int](dimension)
  private val negative = new Array[Boolean](dimension)
  private val significand = new Array[Int](2)
  private val weight = new Array[Int](3)
  private val weighted = new Array[Int](5)
  private val product = new Array[Int](7)
  def frequency: Long = mass
  def positiveOccurrences: Long = occurrences
  def resources: CenteredMomentsResources = CenteredMoments.resources(dimension).toOption.get
  private def split(value: Long): Unit =
    significand(0) = (value & FixedMomentLimbs.Mask).toInt
    significand(1) = (value >>> 30).toInt
  def add(values: Array[Double], frequency: Long = 1L): Either[CenteredMomentsError, Unit] =
    if values == null then Left(MissingInput)
    else if values.length != dimension then Left(LengthMismatch(dimension, values.length))
    else if frequency < 0 then Left(InvalidFrequency(frequency))
    else if frequency > CenteredMoments.MaxFrequency - mass then Left(CapacityExceeded(mass, frequency))
    else
      // Validate/decompose one caller value at a time. Scratch may change on
      // refusal, but no observable moment, frequency or occurrence changes.
      var row = 0
      while row < dimension do
        val value = values(row)
        if !value.isFinite then return Left(NonFiniteInput(row))
        val bits = java.lang.Double.doubleToRawLongBits(value)
        val exponent = ((bits >>> 52) & 0x7ffL).toInt
        val fraction = bits & 0xfffffffffffffL
        mantissas(row) = if exponent == 0 then fraction else fraction | (1L << 52)
        shifts(row) = if exponent == 0 then 0 else exponent - 1
        negative(row) = bits < 0
        row += 1
      if frequency == 0 then CenteredMoments.success
      else
        weight(0) = (frequency & FixedMomentLimbs.Mask).toInt
        weight(1) = ((frequency >>> 30) & FixedMomentLimbs.Mask).toInt
        weight(2) = (frequency >>> 60).toInt
        val weightLength = FixedMomentLimbs.length(weight)
        row = 0; var cell = 0
        while row < dimension do
          split(mantissas(row))
          val weightedLength = FixedMomentLimbs.multiply(significand, 2, weight, weightLength, weighted)
          sums.addShifted(row, weighted, weightedLength, shifts(row), negative(row))
          var column = row
          while column < dimension do
            split(mantissas(column))
            val productLength = FixedMomentLimbs.multiply(weighted, weightedLength, significand, 2, product)
            products.addShifted(
              cell,
              product,
              productLength,
              shifts(row) + shifts(column),
              negative(row) != negative(column)
            )
            cell += 1; column += 1
          row += 1
        mass += frequency; occurrences += 1
        CenteredMoments.success
  def copy(): CenteredMoments = new CenteredMoments(dimension, sums.copy(), products.copy(), mass, occurrences)
  def merge(other: CenteredMoments): Either[CenteredMomentsError, Unit] =
    if this eq other then Left(SelfMerge)
    else if dimension != other.dimension then Left(DimensionMismatch(dimension, other.dimension))
    else if other.mass > CenteredMoments.MaxFrequency - mass then Left(CapacityExceeded(mass, other.mass))
    else
      sums.merge(other.sums); products.merge(other.products)
      mass += other.mass; occurrences += other.occurrences
      CenteredMoments.success
  def finish(degreesRemoved: Int = 1): Either[CenteredMomentsError, CenteredMomentsSummary] =
    if degreesRemoved != 0 && degreesRemoved != 1 then Left(UnsupportedDivisor(degreesRemoved))
    else if mass <= degreesRemoved then Left(InsufficientFrequency(mass, degreesRemoved))
    else
      val weight = BigInt(mass); val meanDenominator = weight << 1074
      val covarianceDenominator = (weight * (mass - degreesRemoved)) << 2148
      val absoluteSums = Array.tabulate(dimension)(sums.integer)
      val mean = new Array[Double](dimension); val covariance = new Array[Double](dimension * dimension)
      var row = 0; var cell = 0
      while row < dimension do
        mean(row) = ExactDyadic.roundRatio(absoluteSums(row), meanDenominator)
        if !mean(row).isFinite then return Left(NonFiniteReadout(row, -1))
        var column = row
        while column < dimension do
          val centered = weight * products.integer(cell) - absoluteSums(row) * absoluteSums(column)
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
  private val ReadoutBits = 4196 + 120 + 1 + 1074
  private val SumDigits = (SumBits + 29) / 30
  private val ProductDigits = (ProductBits + 29) / 30
  private val success: Either[CenteredMomentsError, Unit] = Right(())
  def resources(dimension: Int): Either[CenteredMomentsError, CenteredMomentsResources] =
    if dimension <= 0 then Left(CenteredMomentsError.InvalidDimension(dimension))
    else if dimension.toLong * dimension > Int.MaxValue then Left(CenteredMomentsError.ShapeOverflow(dimension))
    else
      val n = BigInt(dimension); val triangle = n * (n + 1) / 2
      def limbs(bits: Int): BigInt = BigInt((bits + 31) / 32) * 4
      // Padded signed physical radix30 arrays, primitive decode/multiply
      // scratch retained by this worker-local state, and two Long counters.
      val state = (n * SumDigits + triangle * ProductDigits) * 4 + n * 13 + 68 + 16
      val scalar = limbs(ReadoutBits); val output = 8 * (n + n * n)
      Right(
        CenteredMomentsResources(
          state,
          state,
          0,
          0,
          4096 * scalar + n * limbs(SumBits) + 2 * output,
          output,
          SumBits,
          ProductBits,
          ReadoutBits
        )
      )
  def make(dimension: Int): Either[CenteredMomentsError, CenteredMoments] =
    resources(dimension).map: _ =>
      val cells = (dimension.toLong * (dimension + 1) / 2).toInt
      new CenteredMoments(
        dimension,
        FixedMomentLimbs.zero(dimension, SumDigits),
        FixedMomentLimbs.zero(cells, ProductDigits),
        0L,
        0L
      )
