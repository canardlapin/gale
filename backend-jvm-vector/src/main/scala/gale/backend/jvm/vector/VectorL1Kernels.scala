package gale.backend.jvm.vector

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray

import jdk.incubator.vector.DoubleVector
import jdk.incubator.vector.VectorOperators
import jdk.incubator.vector.VectorSpecies

/** W2.1 spike: explicit Vector API versions of the level-1 and reduction kernels.
  *
  * Each kernel has the strided signature of its [[gale.kernel.DoubleKernels]]
  * twin. The unit-stride path runs whole `SPECIES_PREFERRED` vectors plus a
  * scalar tail; any other stride, or a runtime narrower than two doubles,
  * forwards to the pure kernel. Nothing routes here yet: these exist to be
  * measured against the pure kernels and Breeze VectorBLAS (ADR A-2b gate,
  * `docs/verification/w21-simd-spike/README.md`).
  *
  * Semantics relative to the pure kernels:
  *   - `daxpy` is bit-identical (one fma per element, same operand order).
  *   - `ddot` and `dsum` reassociate (lane order, four vector accumulators), so they
  *     agree within the usual `n·ε·Σ|terms|` summation bound, not bitwise. The
  *     horizontal combine has a fixed order, so a result is reproducible for a given
  *     species width (but differs between, say, 2-lane NEON and 4-lane AVX2).
  *   - `dnrm2` keeps the pure contract (overflow/underflow-safe; NaN iff an entry is
  *     NaN, else `+Inf` iff an entry is infinite) and is accurate to a few ulp.
  *   - `dmaxIndex` is exact: first maximum, first NaN wins, `±0` compare equal.
  *   - `dexpInto` uses `lanewise(EXP)`. Measured (JDK 25, aarch64, 2 lanes) within
  *     1 ulp of `StrictMath.exp` (so within 2 ulp of the true value) over
  *     `[-745, 710]`, but neither bit-identical to the
  *     pure `math.exp` nor stable across JIT tiers: about 1% of results change by
  *     1 ulp once the kernel is C2-compiled. The scalar tail uses `Math.exp`, so a
  *     value's result also depends on its position. Recommendation: keep it out of
  *     default routing (see `docs/verification/w21-simd-spike/README.md`).
  */
private[gale] object VectorL1Kernels:
  private final val Species: VectorSpecies[java.lang.Double] = DoubleVector.SPECIES_PREFERRED
  private final val Lanes: Int = Species.length()
  private final val Simd: Boolean = Lanes >= 2

  // Mirrors DoubleKernels.dnrmFrobenius: an optimistic unscaled sum of squares is
  // trusted when finite and at least this large.
  private final val TrustedSumsqMin = 1e-280

  def lanes: Int = Lanes

  /** Horizontal sum in a fixed pairwise order. `reduceLanes(ADD)` leaves the lane
    * order unspecified, so its rounding may differ between the interpreter, C1 and
    * the C2 intrinsic on wide species. Lane extraction is exact, so this combine is
    * reproducible for a given species width. The width itself is part of the result:
    * a 2-lane and an 8-lane machine reassociate differently.
    */
  private def sumLanes(v: DoubleVector): Double =
    Lanes match
      case 2 => v.lane(0) + v.lane(1)
      case 4 => (v.lane(0) + v.lane(1)) + (v.lane(2) + v.lane(3))
      case 8 =>
        ((v.lane(0) + v.lane(1)) + (v.lane(2) + v.lane(3))) +
          ((v.lane(4) + v.lane(5)) + (v.lane(6) + v.lane(7)))
      case _ =>
        var r = v.lane(0)
        var i = 1
        while i < Lanes do
          r += v.lane(i)
          i += 1
        r

  def ddot(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Double =
    if Simd && xStride == 1 && yStride == 1 then
      dotContiguous(n, DoubleArray.asArray(x), xOffset, DoubleArray.asArray(y), yOffset)
    else DoubleKernels.ddot(n, x, xOffset, xStride, y, yOffset, yStride)

  private def dotContiguous(n: Int, x: Array[Double], xOffset: Int, y: Array[Double], yOffset: Int): Double =
    val species = Species
    val step = Lanes
    val bound4 = n - n % (4 * step)
    var s0 = DoubleVector.zero(species)
    var s1 = DoubleVector.zero(species)
    var s2 = DoubleVector.zero(species)
    var s3 = DoubleVector.zero(species)
    var i = 0
    while i < bound4 do
      s0 = DoubleVector.fromArray(species, x, xOffset + i).fma(DoubleVector.fromArray(species, y, yOffset + i), s0)
      s1 = DoubleVector
        .fromArray(species, x, xOffset + i + step)
        .fma(DoubleVector.fromArray(species, y, yOffset + i + step), s1)
      s2 = DoubleVector
        .fromArray(species, x, xOffset + i + 2 * step)
        .fma(DoubleVector.fromArray(species, y, yOffset + i + 2 * step), s2)
      s3 = DoubleVector
        .fromArray(species, x, xOffset + i + 3 * step)
        .fma(DoubleVector.fromArray(species, y, yOffset + i + 3 * step), s3)
      i += 4 * step
    val bound = species.loopBound(n)
    while i < bound do
      s0 = DoubleVector.fromArray(species, x, xOffset + i).fma(DoubleVector.fromArray(species, y, yOffset + i), s0)
      i += step
    var acc = sumLanes(s0.add(s1).add(s2.add(s3)))
    while i < n do
      acc = Math.fma(x(xOffset + i), y(yOffset + i), acc)
      i += 1
    acc

  def daxpy(
      n: Int,
      alpha: Double,
      x: DoubleArray,
      xOffset: Int,
      xStride: Int,
      y: DoubleArray,
      yOffset: Int,
      yStride: Int
  ): Unit =
    if Simd && xStride == 1 && yStride == 1 then
      val species = Species
      val xs = DoubleArray.asArray(x)
      val ys = DoubleArray.asArray(y)
      val a = DoubleVector.broadcast(species, alpha)
      val bound = species.loopBound(n)
      var i = 0
      while i < bound do
        // x·alpha + y with one rounding: the same value as the pure fma(alpha, x, y).
        DoubleVector
          .fromArray(species, xs, xOffset + i)
          .fma(a, DoubleVector.fromArray(species, ys, yOffset + i))
          .intoArray(ys, yOffset + i)
        i += Lanes
      while i < n do
        ys(yOffset + i) = Math.fma(alpha, xs(xOffset + i), ys(yOffset + i))
        i += 1
    else DoubleKernels.daxpy(n, alpha, x, xOffset, xStride, y, yOffset, yStride)

  def dsum(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    if Simd && xStride == 1 then
      val species = Species
      val xs = DoubleArray.asArray(x)
      val step = Lanes
      val bound4 = n - n % (4 * step)
      var s0 = DoubleVector.zero(species)
      var s1 = DoubleVector.zero(species)
      var s2 = DoubleVector.zero(species)
      var s3 = DoubleVector.zero(species)
      var i = 0
      while i < bound4 do
        s0 = s0.add(DoubleVector.fromArray(species, xs, xOffset + i))
        s1 = s1.add(DoubleVector.fromArray(species, xs, xOffset + i + step))
        s2 = s2.add(DoubleVector.fromArray(species, xs, xOffset + i + 2 * step))
        s3 = s3.add(DoubleVector.fromArray(species, xs, xOffset + i + 3 * step))
        i += 4 * step
      val bound = species.loopBound(n)
      while i < bound do
        s0 = s0.add(DoubleVector.fromArray(species, xs, xOffset + i))
        i += step
      var acc = sumLanes(s0.add(s1).add(s2.add(s3)))
      while i < n do
        acc += xs(xOffset + i)
        i += 1
      acc
    else DoubleKernels.dsum(n, x, xOffset, xStride)

  /** Euclidean norm with the [[gale.kernel.DoubleKernels.dnrm2]] contract.
    *
    * One optimistic SIMD pass forms the fma sum of squares. Only when that total
    * is non-finite, zero, or below a safe floor do two more SIMD passes find the
    * largest magnitude `m` and sum `(x/m)^2`, giving `m * sqrt(...)` (the
    * [[gale.kernel.DoubleKernels.dnrmFrobenius]] strategy).
    */
  def dnrm2(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    if !(Simd && xStride == 1) then DoubleKernels.dnrm2(n, x, xOffset, xStride)
    else if n < 1 then 0.0
    else
      val xs = DoubleArray.asArray(x)
      val ssq = dotContiguous(n, xs, xOffset, xs, xOffset)
      if ssq.isFinite && ssq >= TrustedSumsqMin then math.sqrt(ssq)
      else
        val m = maxAbsContiguous(n, xs, xOffset)
        // NaN anywhere makes `m` NaN (MAX propagates NaN); otherwise an infinite
        // entry makes it +Inf. Zero means every entry is ±0.
        if m == 0.0 || m.isNaN || m.isInfinite then m
        else m * math.sqrt(scaledSumsqContiguous(n, xs, xOffset, m))

  private def maxAbsContiguous(n: Int, x: Array[Double], xOffset: Int): Double =
    val species = Species
    val bound = species.loopBound(n)
    var m = DoubleVector.zero(species)
    var i = 0
    while i < bound do
      m = m.lanewise(VectorOperators.MAX, DoubleVector.fromArray(species, x, xOffset + i).abs())
      i += Lanes
    var r = m.reduceLanes(VectorOperators.MAX)
    while i < n do
      r = math.max(r, math.abs(x(xOffset + i)))
      i += 1
    r

  private def scaledSumsqContiguous(n: Int, x: Array[Double], xOffset: Int, scale: Double): Double =
    val species = Species
    val s = DoubleVector.broadcast(species, scale)
    val bound = species.loopBound(n)
    var acc = DoubleVector.zero(species)
    var i = 0
    while i < bound do
      val v = DoubleVector.fromArray(species, x, xOffset + i).div(s)
      acc = v.fma(v, acc)
      i += Lanes
    var r = sumLanes(acc)
    while i < n do
      val v = x(xOffset + i) / scale
      r = Math.fma(v, v, r)
      i += 1
    r

  /** Index of the first maximum, `-1` when `n == 0`; the first NaN wins. Identical
    * to [[gale.kernel.DoubleKernels.dmaxIndex]].
    *
    * One blocked pass. Each block's maximum is formed with lane-wise `MAX`, which
    * propagates NaN: a NaN block maximum means the block holds the first NaN, which
    * a scalar scan of that block returns. Otherwise a block replaces the running
    * winner only when its maximum is strictly larger, so ties keep the earliest
    * block. Finally the winning block is scanned for the first element equal to the
    * maximum. Equality treats `-0.0 == 0.0`, matching the pure kernel's strict `>`
    * (a later `±0` never displaces an earlier one); an all-`-Inf` input reports 0.
    */
  def dmaxIndex(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Int =
    if !(Simd && xStride == 1) || n < 2 * Lanes then DoubleKernels.dmaxIndex(n, x, xOffset, xStride)
    else
      val xs = DoubleArray.asArray(x)
      var best = Double.NegativeInfinity
      var bestBlock = 0
      var start = 0
      while start < n do
        val length = math.min(ArgBlock, n - start)
        val m = maxContiguous(length, xs, xOffset + start)
        if m.isNaN then return start + firstNaN(length, xs, xOffset + start)
        if m > best then
          best = m
          bestBlock = start
        start += length
      var i = bestBlock
      while xs(xOffset + i) != best do i += 1
      i

  // Elements per argmax block: large enough to amortize the block reduction, small
  // enough that the final rescan of the winning block is cheap.
  private final val ArgBlock = 512

  /** Lane-wise `MAX` over four accumulators; NaN if any element is NaN. */
  private def maxContiguous(n: Int, x: Array[Double], xOffset: Int): Double =
    val species = Species
    val step = Lanes
    val bound4 = n - n % (4 * step)
    var m0 = DoubleVector.broadcast(species, Double.NegativeInfinity)
    var m1 = m0
    var m2 = m0
    var m3 = m0
    var i = 0
    while i < bound4 do
      m0 = m0.lanewise(VectorOperators.MAX, DoubleVector.fromArray(species, x, xOffset + i))
      m1 = m1.lanewise(VectorOperators.MAX, DoubleVector.fromArray(species, x, xOffset + i + step))
      m2 = m2.lanewise(VectorOperators.MAX, DoubleVector.fromArray(species, x, xOffset + i + 2 * step))
      m3 = m3.lanewise(VectorOperators.MAX, DoubleVector.fromArray(species, x, xOffset + i + 3 * step))
      i += 4 * step
    val bound = species.loopBound(n)
    while i < bound do
      m0 = m0.lanewise(VectorOperators.MAX, DoubleVector.fromArray(species, x, xOffset + i))
      i += step
    var r = m0.lanewise(VectorOperators.MAX, m1).lanewise(VectorOperators.MAX, m2.lanewise(VectorOperators.MAX, m3))
      .reduceLanes(VectorOperators.MAX)
    while i < n do
      r = math.max(r, x(xOffset + i))
      i += 1
    r

  private def firstNaN(n: Int, x: Array[Double], xOffset: Int): Int =
    var i = 0
    while i < n && !x(xOffset + i).isNaN do i += 1
    i

  /** `y_i := exp(x_i)` via `lanewise(EXP)`; `y` may be `x` itself. */
  def dexpInto(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    if Simd && xStride == 1 && yStride == 1 then
      val species = Species
      val xs = DoubleArray.asArray(x)
      val ys = DoubleArray.asArray(y)
      val bound = species.loopBound(n)
      var i = 0
      while i < bound do
        DoubleVector.fromArray(species, xs, xOffset + i).lanewise(VectorOperators.EXP).intoArray(ys, yOffset + i)
        i += Lanes
      while i < n do
        ys(yOffset + i) = Math.exp(xs(xOffset + i))
        i += 1
    else DoubleKernels.dexpInto(n, x, xOffset, xStride, y, yOffset, yStride)

end VectorL1Kernels
