package gale.spectral

import gale.platform.DoubleArray

/** Exact binary normalization for dense kernels that form squared intermediates. */
private[spectral] object DenseInputScaling:
  /** SVD and Francis QR recurrences form products at input scale
    * squared. Normalize extreme inputs before either stage; ordinary inputs
    * receive no scaling multiplications, preserving their arithmetic order.
    * The +/-256 interval leaves headroom for dimension growth and epsilon-scale
    * active bidiagonal entries. This is normwise protection, not a guarantee of
    * relative accuracy for tiny singular values in a mixed-scale matrix.
    * Returns e such that the owned copy was multiplied by 2^e.
    */
  def normalize(a: DoubleArray, length: Int): Int =
    var norm = 0.0
    var i = 0
    while i < length do
      norm = math.max(norm, math.abs(a(i)))
      i += 1
    if !norm.isFinite || norm == 0.0 || (norm >= powerOfTwo(-256) && norm <= powerOfTwo(256)) then 0
    else
      var exponent = 0
      if norm < java.lang.Double.longBitsToDouble(1L << 52) then
        // Lift a subnormal maximum before extracting its binary exponent.
        scaleBy(a, length, 600)
        norm *= powerOfTwo(600)
        exponent = 600
      val normExponent = ((java.lang.Double.doubleToRawLongBits(norm) >>> 52) & 0x7ffL).toInt - 1023
      scaleBy(a, length, -normExponent)
      exponent - normExponent

  def powerOfTwo(k: Int): Double = java.lang.Double.longBitsToDouble((1023L + k) << 52)

  /** Two normal factors also cover exponents needed for subnormal inputs;
    * never form an infinite reciprocal scaling factor.
    */
  def scaleBy(a: DoubleArray, length: Int, exponent: Int): Unit =
    val first = powerOfTwo(exponent / 2)
    val second = powerOfTwo(exponent - exponent / 2)
    var i = 0
    while i < length do
      a(i) = a(i) * first * second
      i += 1

