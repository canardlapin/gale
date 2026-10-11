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

/** Worker-local exact finite-input moments. Inputs/merge sources are not retained or mutated. Exact dyadics retain a
  * separate power-of-two shift instead of repeatedly allocating low zero-padding. Their absolute units remain 2^-1074
  * for sums and 2^-2148 for products. Exact centered covariance is (W*sum(w*x*y)-sum(w*x)*sum(w*y))/(W*(W-ddof)); the
  * apparently subtractive formula uses exact integers, so cancellation never discards bits. Mean/covariance round once,
  * ties to even. Frequency <= MaxFrequency bounds every limb, independently of traversal/partitions. This is a
  * correctness primitive, not a cheap-loop/RSS claim.
  */
final class CenteredMoments private (
    val dimension: Int,
    private var sums: Array[BigInt],
    private var products: Array[BigInt],
    private var sumShifts: Array[Int],
    private var productShifts: Array[Int],
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
        val exact = new Array[BigInt](dimension)
        val shifts = new Array[Int](dimension)
        row = 0
        while row < dimension do
          val bits = java.lang.Double.doubleToRawLongBits(values(row))
          val exponent = ((bits >>> 52) & 0x7ffL).toInt
          val fraction = bits & 0xfffffffffffffL
          val mantissa = if exponent == 0 then fraction else fraction | (1L << 52)
          val trailing = if mantissa == 0 then 0 else java.lang.Long.numberOfTrailingZeros(mantissa)
          val compact = mantissa >>> trailing
          exact(row) = BigInt(if bits < 0L then -compact else compact)
          shifts(row) = (if exponent == 0 then 0 else exponent - 1) + trailing
          row += 1
        val weight = BigInt(frequency)
        val nextSums = new Array[BigInt](dimension)
        val nextProducts = new Array[BigInt](products.length)
        val nextSumShifts = sumShifts.clone()
        val nextProductShifts = productShifts.clone()
        row = 0; var cell = 0
        while row < dimension do
          val weighted = weight * exact(row)
          nextSums(row) = combine(sums(row), sumShifts(row), weighted, shifts(row), nextSumShifts, row)
          var column = row
          while column < dimension do
            nextProducts(cell) = combine(
              products(cell),
              productShifts(cell),
              weighted * exact(column),
              shifts(row) + shifts(column),
              nextProductShifts,
              cell
            )
            cell += 1; column += 1
          row += 1
        sums = nextSums; products = nextProducts; sumShifts = nextSumShifts; productShifts = nextProductShifts
        mass += frequency; occurrences += 1
        Right(())

  /** Align only the meaningful bits. A zero operand does not impose an unnecessarily low scale. */
  private def combine(a: BigInt, aShift: Int, b: BigInt, bShift: Int, outputShifts: Array[Int], index: Int): BigInt =
    if b == 0 then a
    else if a == 0 then
      outputShifts(index) = bShift
      b
    else
      val shift = math.min(aShift, bShift)
      outputShifts(index) = shift
      (a << (aShift - shift)) + (b << (bShift - shift))
  def copy(): CenteredMoments =
    new CenteredMoments(
      dimension,
      sums.clone(),
      products.clone(),
      sumShifts.clone(),
      productShifts.clone(),
      mass,
      occurrences
    )
  def merge(other: CenteredMoments): Either[CenteredMomentsError, Unit] =
    if this eq other then Left(SelfMerge)
    else if dimension != other.dimension then Left(DimensionMismatch(dimension, other.dimension))
    else if other.mass > CenteredMoments.MaxFrequency - mass then Left(CapacityExceeded(mass, other.mass))
    else
      val nextSumShifts = sumShifts.clone()
      val nextProductShifts = productShifts.clone()
      val nextSums = sums.indices
        .map(i => combine(sums(i), sumShifts(i), other.sums(i), other.sumShifts(i), nextSumShifts, i))
        .toArray
      val nextProducts = products.indices
        .map(i =>
          combine(products(i), productShifts(i), other.products(i), other.productShifts(i), nextProductShifts, i)
        )
        .toArray
      sums = nextSums; products = nextProducts; sumShifts = nextSumShifts; productShifts = nextProductShifts
      mass += other.mass; occurrences += other.occurrences
      Right(())
  def finish(degreesRemoved: Int = 1): Either[CenteredMomentsError, CenteredMomentsSummary] =
    if degreesRemoved != 0 && degreesRemoved != 1 then Left(UnsupportedDivisor(degreesRemoved))
    else if mass <= degreesRemoved then Left(InsufficientFrequency(mass, degreesRemoved))
    else
      val weight = BigInt(mass); val meanDenominator = weight << 1074
      val covarianceDenominator = (weight * (mass - degreesRemoved)) << 2148
      val absoluteSums = sums.indices.map(i => sums(i) << sumShifts(i)).toArray
      val mean = new Array[Double](dimension); val covariance = new Array[Double](dimension * dimension)
      var row = 0; var cell = 0
      while row < dimension do
        mean(row) = ExactDyadic.roundRatio(absoluteSums(row), meanDenominator)
        if !mean(row).isFinite then return Left(NonFiniteReadout(row, -1))
        var column = row
        while column < dimension do
          val centered = weight * (products(cell) << productShifts(cell)) - absoluteSums(row) * absoluteSums(column)
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
      val shifts = (n + triangle) * 4
      val state = n * limbs(SumBits) + triangle * limbs(ProductBits) + shifts + 16
      val scalar = limbs(ReadoutBits)
      val output = 8 * (n + n * n)
      Right(
        CenteredMomentsResources(
          state,
          state,
          state + n * (limbs(2098) + 4) + 4096 * scalar,
          state + 4096 * scalar,
          4096 * scalar + n * limbs(SumBits) + 2 * output,
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
        new Array[Int](dimension),
        new Array[Int]((dimension.toLong * (dimension + 1) / 2).toInt),
        0L,
        0L
      )
    )
