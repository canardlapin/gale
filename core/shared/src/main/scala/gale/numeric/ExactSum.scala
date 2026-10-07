package gale.numeric

/** Source-level numerical payload estimates, excluding object/host/GC overhead. State/copy/value capacities are exact
  * Long-array payloads. Ratio is a deliberately conservative estimate for normalized arrays and bounded BigInt limb
  * scratch; it is not a JVM heap/RSS or JavaScript engine-memory guarantee.
  */
final case class ExactSumResources(
    stateBytes: Long,
    copyAdditionalBytes: Long,
    valueScratchBytes: Long,
    ratioScratchBytes: Long,
    maximumExactIntegerBits: Int,
    maximumRatioOperandBits: Int
)

/** Exact, order-independent accumulation of IEEE doubles.
  *
  * Every finite double is an integer multiple of `2^-1074`, so the running total is held exactly as a fixed-point
  * integer in carry-save base-`2^32` digits. Addition and [[addAll]] are therefore associative and commutative, and
  * [[value]] rounds the exact total to nearest (ties to even) once. Finite-input totals agree bit-for-bit across input
  * permutations and arbitrary merge trees within [[ExactSum.MaxTerms]] nonzero finite inputs.
  *
  * Non-finite inputs follow IEEE addition in a separate channel, so a total that saw `NaN` or an infinity reports it
  * exactly as a naive sum would. NaN payload bits and their ordering are unspecified. Zero and non-finite inputs do not
  * consume the finite-input budget. Capacity refusals leave the accumulator unchanged.
  *
  * Instances are mutable builders and are not thread-safe. [[copy]] owns independent storage, [[addAll]] leaves its
  * source unchanged unless merging an accumulator with itself, and [[value]] never changes the builder.
  */
final class ExactSum private (
    private val digits: Array[Long],
    private var pending: Int,
    private var special: Double,
    private var hasSpecial: Boolean,
    private var termCount: Long
):
  import ExactSum.*

  /** Immutable facts shared by every state, including empty/special/capacity states. Inspection neither normalizes nor
    * rounds the builder and allocates no workspace.
    */
  def resources: ExactSumResources = ExactSum.resources

  /** Number of nonzero finite inputs represented, including inputs inherited through merges. */
  def finiteTerms: Long = termCount

  /** Add one input. Successful additions reuse a shared result and allocate no result objects. */
  def add(value: Double): Either[ExactSumError.CapacityExceeded, Unit] =
    if value != 0.0 then
      if !value.isFinite then
        special = if hasSpecial then special + value else value
        hasSpecial = true
      else
        if termCount == MaxTerms then return Left(ExactSumError.CapacityExceeded(termCount, 1L))
        val bits = java.lang.Double.doubleToRawLongBits(value)
        val exponent = ((bits >>> 52) & 0x7ffL).toInt
        val fraction = bits & 0xfffffffffffffL
        val mantissa = if exponent == 0 then fraction else fraction | (1L << 52)
        val position = if exponent == 0 then 0 else exponent - 1
        val index = position >>> 5
        val shift = position & 31
        val low = (mantissa << shift) & DigitMask
        val middle = (mantissa >>> (32 - shift)) & DigitMask
        val high = if shift == 0 then 0L else mantissa >>> (64 - shift)
        if bits < 0L then
          digits(index) -= low
          digits(index + 1) -= middle
          digits(index + 2) -= high
        else
          digits(index) += low
          digits(index + 1) += middle
          digits(index + 2) += high
        termCount += 1L
        pending += 1
        if pending >= NormalizeEvery then normalize()
    Added

  /** Merge exact state without rounding; self-merge doubles both the total and finite-input count. */
  def addAll(other: ExactSum): Either[ExactSumError.CapacityExceeded, Unit] =
    if MaxTerms - termCount < other.termCount then
      return Left(ExactSumError.CapacityExceeded(termCount, other.termCount))
    normalize()
    // This side is now below 2^32 per digit and the other side, whatever its pending count, below 2^62 + 2^32,
    // so adding its carry-save digits directly cannot overflow and leaves `other` untouched.
    var index = 0
    while index < DigitCount do
      digits(index) += other.digits(index)
      index += 1
    pending = 1
    normalize()
    if other.hasSpecial then
      special = if hasSpecial then special + other.special else other.special
      hasSpecial = true
    termCount += other.termCount
    Added

  def copy(): ExactSum =
    new ExactSum(digits.clone(), pending, special, hasSpecial, termCount)

  /** The exact total rounded once to the nearest double. */
  def value: Double =
    if hasSpecial then special
    else
      val magnitude = digits.clone()
      normalizeDigits(magnitude)
      val negative = magnitude(DigitCount - 1) < 0L
      if negative then
        var index = 0
        while index < DigitCount do
          magnitude(index) = -magnitude(index)
          index += 1
        normalizeDigits(magnitude)
      val rounded = roundMagnitude(magnitude)
      if negative then -rounded else rounded

  /** Divide this exact finite total by another exact finite total before rounding once to nearest, ties to even.
    * Neither builder is mutated; each retains its independent [[ExactSum.MaxTerms]] capacity. Finite totals may exceed
    * Double range. A genuinely overflowing ratio returns signed infinity; underflow returns a subnormal or signed zero.
    * An exactly zero numerator returns positive zero. Non-finite input channels and an exactly zero denominator are
    * refused, even when their rounded readouts would suggest a usable ratio. This readout allocates arbitrary-precision
    * integer scratch storage; addition and merging retain their existing allocation behavior.
    */
  def ratio(normalizer: ExactSum): Either[ExactSumError, Double] =
    if hasSpecial || normalizer.hasSpecial then Left(ExactSumError.NonFiniteRatio)
    else
      val denominator = normalizer.exactInteger
      if denominator == 0 then Left(ExactSumError.ZeroNormalizer)
      else Right(roundRatio(exactInteger, denominator))

  private def exactInteger: BigInt =
    val normalized = digits.clone()
    normalizeDigits(normalized)
    var result = BigInt(normalized(DigitCount - 1))
    var index = DigitCount - 2
    while index >= 0 do
      result = (result << 32) + normalized(index)
      index -= 1
    result

  private def normalize(): Unit =
    if pending > 0 then
      normalizeDigits(digits)
      pending = 0

enum ExactSumError:
  case CapacityExceeded(currentTerms: Long, incomingTerms: Long)
  case NonFiniteRatio
  case ZeroNormalizer

  def message: String = this match
    case CapacityExceeded(currentTerms, incomingTerms) =>
      s"Exact sum capacity ${ExactSum.MaxTerms} nonzero finite inputs exceeded: $currentTerms + $incomingTerms"
    case NonFiniteRatio => "Exact sum ratio requires finite input channels in both accumulators"
    case ZeroNormalizer => "Exact sum ratio requires a nonzero exact normalizer"

object ExactSum:
  /** Supported number of nonzero finite inputs, counted across merges without resetting after cancellation. */
  val MaxTerms: Long = 1L << 60

  private val Added: Either[ExactSumError.CapacityExceeded, Unit] = Right(())
  // Every finite double has magnitude below 2^1024 = 2^2098 units, so bits 0 through 2097 (digits 0 through 65)
  // cover one term; the remaining digits absorb carries from up to 2^60 terms.
  private val DigitCount = 70
  private val DigitMask = 0xffffffffL
  // Each addition moves a digit by less than 2^32, so 2^30 additions keep every carry-save digit below 2^63.
  private val NormalizeEvery = 1 << 30
  private val SignificandWindow = 62

  /** Qualified numerical capacity model for this implementation, not host allocation. A finite state is below 2^(2098 +
    * log2(MaxTerms)) fixed-point units. A ratio operand can additionally shift by at most 1074 bits. Each readout
    * converts two DigitCount arrays with repeated bounded shift/add scratch, then divides bounded operands. The limb
    * allowance deliberately exceeds the ordinary live temporaries and cumulative primitive allocations of the qualified
    * JVM path. Hosts needing total allocator/native bounds must supply separate evidence.
    */
  val resources: ExactSumResources =
    val state = DigitCount.toLong * java.lang.Long.BYTES
    val exactBits = 2098 + (63 - java.lang.Long.numberOfLeadingZeros(MaxTerms))
    val operandBits = exactBits + 1074
    val limbBytes = ((operandBits.toLong + 31) / 32) * java.lang.Integer.BYTES
    ExactSumResources(
      state,
      state,
      state,
      2 * state + (16L * DigitCount + 4096L) * limbBytes,
      exactBits,
      operandBits
    )

  def zero(): ExactSum =
    new ExactSum(new Array[Long](DigitCount), 0, 0.0, false, 0L)

  private def roundRatio(numerator: BigInt, denominator: BigInt): Double =
    if numerator == 0 then 0.0
    else
      val negative = numerator.signum != denominator.signum
      val n = numerator.abs
      val d = denominator.abs
      // Locate the leading binary exponent using exact comparisons, including ratios below one.
      var exponent = n.bitLength - d.bitLength
      val below = if exponent >= 0 then n < (d << exponent) else (n << -exponent) < d
      if below then exponent -= 1
      val magnitude =
        if exponent >= 1024 then Double.PositiveInfinity
        else
          var quantum = math.max(exponent - 52, -1074)
          val (whole, remainder) =
            if quantum >= 0 then n /% (d << quantum)
            else (n << -quantum) /% d
          val divisor = if quantum >= 0 then d << quantum else d
          val comparison = (remainder << 1).compare(divisor)
          var rounded = whole
          if comparison > 0 || (comparison == 0 && whole.testBit(0)) then rounded += 1
          if rounded.bitLength > 53 then
            rounded >>= 1
            quantum += 1
          val significand = rounded.toLong
          val bits =
            if significand < (1L << 52) then significand // Subnormal or zero, in units of 2^-1074.
            else
              val field = quantum + 52 + 1023
              if field >= 2047 then 0x7ff0000000000000L
              else (field.toLong << 52) | (significand & 0xfffffffffffffL)
          java.lang.Double.longBitsToDouble(bits)
      if negative then -magnitude else magnitude

  private def normalizeDigits(digits: Array[Long]): Unit =
    var index = 0
    while index < DigitCount - 1 do
      val carry = digits(index) >> 32
      digits(index) -= carry << 32
      digits(index + 1) += carry
      index += 1

  /** Round a normalized non-negative fixed-point integer, in units of `2^-1074`, to the nearest double. */
  private def roundMagnitude(digits: Array[Long]): Double =
    var top = DigitCount - 1
    while top >= 0 && digits(top) == 0L do top -= 1
    if top < 0 then 0.0
    else
      val bitLength = 32 * top + (64 - java.lang.Long.numberOfLeadingZeros(digits(top)))
      if bitLength <= SignificandWindow then
        var integer = 0L
        var index = top
        while index >= 0 do
          integer = (integer << 32) | digits(index)
          index -= 1
        // Either the integer is below 2^53 and exact, or the scaled result is normal; one rounding either way.
        scale(integer.toDouble, -1074)
      else
        val dropped = bitLength - SignificandWindow
        var window = 0L
        var bit = bitLength - 1
        while bit >= dropped do
          window = (window << 1) | bitAt(digits, bit)
          bit -= 1
        // Bit zero of the 62-bit window sits below the rounding position, so it can carry the sticky bit.
        val sticky = if anyBitBelow(digits, dropped) then 1L else 0L
        scale((window | sticky).toDouble, dropped - 1074)

  private def bitAt(digits: Array[Long], bit: Int): Long =
    val index = math.min(bit >>> 5, DigitCount - 1)
    (digits(index) >>> (bit - 32 * index)) & 1L

  private def anyBitBelow(digits: Array[Long], bitCount: Int): Boolean =
    val whole = bitCount >>> 5
    var index = 0
    while index < whole do
      if digits(index) != 0L then return true
      index += 1
    val remainder = bitCount & 31
    remainder > 0 && (digits(whole) & ((1L << remainder) - 1L)) != 0L

  /** `value * 2^exponent` for `value` an integer below 2^62. An exponent above 1023 only arises from a window of at
    * least 2^61 units scaled past 2^1084, which overflows; it is answered directly because `powerOfTwo` would otherwise
    * encode an out-of-range exponent field.
    */
  private def scale(value: Double, exponent: Int): Double =
    if exponent > 1023 then Double.PositiveInfinity
    else value * powerOfTwo(exponent)

  private def powerOfTwo(exponent: Int): Double =
    if exponent >= -1022 then java.lang.Double.longBitsToDouble((exponent + 1023).toLong << 52)
    else java.lang.Double.longBitsToDouble(1L << (exponent + 1074))
