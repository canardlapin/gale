package gale.numeric

/** Shared exact binary64 conversion/readout; no finite arithmetic is rounded until the final ratio. */
private[numeric] object ExactDyadic:
  def units(value: Double): BigInt =
    val bits = java.lang.Double.doubleToRawLongBits(value)
    val exponent = ((bits >>> 52) & 0x7ffL).toInt
    val fraction = bits & 0xfffffffffffffL
    val mantissa = if exponent == 0 then fraction else fraction | (1L << 52)
    val result = BigInt(mantissa) << (if exponent == 0 then 0 else exponent - 1)
    if bits < 0L then -result else result

  def roundRatio(numerator: BigInt, denominator: BigInt): Double =
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
