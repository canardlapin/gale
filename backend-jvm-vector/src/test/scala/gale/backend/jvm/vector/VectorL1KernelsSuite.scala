package gale.backend.jvm.vector

import gale.kernel.DoubleKernels
import gale.numeric.ExactSum
import gale.platform.DoubleArray

/** W2.1 spike correctness: the SIMD level-1/reduction kernels against the pure
  * [[gale.kernel.DoubleKernels]] reference. Exact ops (axpy, argmax) must match
  * bitwise; reassociating sums must stay within the summation error bound; `exp`
  * reports its measured worst ulp error over a dense sweep of its finite range.
  */
class VectorL1KernelsSuite extends munit.FunSuite:
  private val Eps = Math.ulp(1.0) / 2 // unit roundoff
  // Every length class: below one vector, ragged tails, the 4×-unrolled body, large.
  private val Lengths = Seq(0, 1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 33, 64, 127, 1000, 4096, 65537)

  private def lcg(n: Int, seed: Long): Array[Double] =
    var state = seed * 6364136223846793005L + 1442695040888963407L
    Array.fill(n) {
      state = state * 6364136223846793005L + 1442695040888963407L
      ((state >>> 11).toDouble / (1L << 53).toDouble) * 2.0 - 1.0
    }

  // Shares storage (no copy): the in-place kernels must write through to `values`.
  private def da(values: Array[Double]): DoubleArray = DoubleArray.adopt(values)

  /** Reassociated sum bound: γ_n·Σ|terms| with γ_n ≈ n·u (plus a 2u slack for fma). */
  private def sumBound(n: Int, absSum: Double): Double = (n + 2) * Eps * absSum

  test("simd lane count is reported") {
    assert(VectorL1Kernels.lanes >= 1)
  }

  test("ddot agrees with the pure kernel within the reassociation bound, at offsets") {
    for n <- Lengths; offset <- Seq(0, 3) do
      val x = lcg(n + offset, 11L + n)
      val y = lcg(n + offset, 29L + n)
      val simd = VectorL1Kernels.ddot(n, da(x), offset, 1, da(y), offset, 1)
      val pure = DoubleKernels.ddot(n, da(x), offset, 1, da(y), offset, 1)
      var absSum = 0.0
      for i <- 0 until n do absSum += math.abs(x(offset + i) * y(offset + i))
      assert(math.abs(simd - pure) <= 2 * sumBound(n, absSum), s"n=$n offset=$offset simd=$simd pure=$pure")
  }

  /** Correctly rounded `Σ x_i y_i` via exact `BigDecimal` products. */
  private def exactDot(n: Int, x: Array[Double], y: Array[Double], offset: Int): Double =
    var acc = java.math.BigDecimal.ZERO
    for i <- 0 until n do
      acc = acc.add(new java.math.BigDecimal(x(offset + i)).multiply(new java.math.BigDecimal(y(offset + i))))
    acc.doubleValue()

  private def exactSum(n: Int, x: Array[Double], offset: Int): Double =
    val acc = ExactSum.zero()
    for i <- 0 until n do assert(acc.add(x(offset + i)).isRight)
    acc.value

  test("ddot and dsum stay within the summation bound of an exact reference") {
    for n <- Lengths; offset <- Seq(0, 3) do
      val x = lcg(n + offset, 101L + n)
      val y = lcg(n + offset, 103L + n)
      val dotExact = exactDot(n, x, y, offset)
      var absDot = 0.0
      var absSum = 0.0
      for i <- 0 until n do
        absDot += math.abs(x(offset + i) * y(offset + i))
        absSum += math.abs(x(offset + i))
      val dot = VectorL1Kernels.ddot(n, da(x), offset, 1, da(y), offset, 1)
      assert(
        math.abs(dot - dotExact) <= sumBound(n, absDot) + Eps * math.abs(dotExact),
        s"dot n=$n simd=$dot exact=$dotExact"
      )
      val sumExact = exactSum(n, x, offset)
      val sum = VectorL1Kernels.dsum(n, da(x), offset, 1)
      assert(
        math.abs(sum - sumExact) <= sumBound(n, absSum) + Eps * math.abs(sumExact),
        s"sum n=$n simd=$sum exact=$sumExact"
      )
  }

  test("ddot, dsum and dnrm2 are bit-stable across JIT tiers") {
    // The first call runs interpreted; 20,000 more push the kernels through C1 and
    // C2. Every rerun must reproduce the first result's bits (fixed-order combine).
    val n = 1027 // a full unrolled body, a single-vector remainder and a scalar tail
    val x = lcg(n, 211L)
    val y = lcg(n, 223L)
    val tiny = x.map(_ * 1e-300) // forces the nrm2 rescan path
    def bits(v: Double): Long = java.lang.Double.doubleToRawLongBits(v)
    val dot0 = bits(VectorL1Kernels.ddot(n, da(x), 0, 1, da(y), 0, 1))
    val sum0 = bits(VectorL1Kernels.dsum(n, da(x), 0, 1))
    val nrm0 = bits(VectorL1Kernels.dnrm2(n, da(y), 0, 1))
    val tiny0 = bits(VectorL1Kernels.dnrm2(n, da(tiny), 0, 1))
    var diffs = 0
    var rep = 0
    while rep < 20000 do
      if bits(VectorL1Kernels.ddot(n, da(x), 0, 1, da(y), 0, 1)) != dot0 then diffs += 1
      if bits(VectorL1Kernels.dsum(n, da(x), 0, 1)) != sum0 then diffs += 1
      if bits(VectorL1Kernels.dnrm2(n, da(y), 0, 1)) != nrm0 then diffs += 1
      if bits(VectorL1Kernels.dnrm2(n, da(tiny), 0, 1)) != tiny0 then diffs += 1
      rep += 1
    assertEquals(diffs, 0, "results changed bits across reruns")
  }

  test("ddot strided input falls back to the pure kernel bitwise") {
    val x = lcg(300, 3L)
    val y = lcg(300, 4L)
    val simd = VectorL1Kernels.ddot(100, da(x), 1, 3, da(y), 0, 2)
    val pure = DoubleKernels.ddot(100, da(x), 1, 3, da(y), 0, 2)
    assertEquals(java.lang.Double.doubleToRawLongBits(simd), java.lang.Double.doubleToRawLongBits(pure))
  }

  test("daxpy is bit-identical to the pure kernel (contiguous and strided)") {
    for n <- Lengths; offset <- Seq(0, 5) do
      val x = lcg(n + offset, 7L + n)
      val y0 = lcg(n + offset, 13L + n)
      val ySimd = y0.clone()
      val yPure = y0.clone()
      VectorL1Kernels.daxpy(n, -1.37, da(x), offset, 1, da(ySimd), offset, 1)
      DoubleKernels.daxpy(n, -1.37, da(x), offset, 1, da(yPure), offset, 1)
      assert(ySimd.sameElements(yPure), s"n=$n offset=$offset")
    val x = lcg(90, 1L)
    val ySimd = lcg(60, 2L)
    val yPure = ySimd.clone()
    VectorL1Kernels.daxpy(30, 0.5, da(x), 0, 3, da(ySimd), 1, 2)
    DoubleKernels.daxpy(30, 0.5, da(x), 0, 3, da(yPure), 1, 2)
    assert(ySimd.sameElements(yPure))
  }

  test("dsum agrees with the pure kernel within the reassociation bound") {
    for n <- Lengths; offset <- Seq(0, 1) do
      val x = lcg(n + offset, 17L + n).map(_ * 1e3)
      val simd = VectorL1Kernels.dsum(n, da(x), offset, 1)
      val pure = DoubleKernels.dsum(n, da(x), offset, 1)
      val absSum = x.iterator.drop(offset).map(math.abs).sum
      assert(math.abs(simd - pure) <= 2 * sumBound(n, absSum), s"n=$n simd=$simd pure=$pure")
  }

  test("dsum propagates NaN and infinities like the pure kernel") {
    val x = lcg(37, 5L)
    x(20) = Double.PositiveInfinity
    assertEquals(VectorL1Kernels.dsum(37, da(x), 0, 1), Double.PositiveInfinity)
    x(3) = Double.NegativeInfinity
    assert(VectorL1Kernels.dsum(37, da(x), 0, 1).isNaN)
    val y = lcg(37, 6L)
    y(36) = Double.NaN
    assert(VectorL1Kernels.dsum(37, da(y), 0, 1).isNaN)
  }

  private def assertNrm2Close(x: Array[Double], label: String)(using munit.Location): Unit =
    val simd = VectorL1Kernels.dnrm2(x.length, da(x), 0, 1)
    val pure = DoubleKernels.dnrm2(x.length, da(x), 0, 1)
    if pure.isNaN then assert(simd.isNaN, s"$label: simd=$simd pure=NaN")
    else if pure.isInfinite || pure == 0.0 then assertEquals(simd, pure, label)
    else
      // Both are O(n·u) accurate (sum of squares, then sqrt halves the relative error).
      val tolerance = ((x.length + 2) * Eps + 8 * Eps) * pure
      assert(math.abs(simd - pure) <= tolerance, s"$label: simd=$simd pure=$pure")

  test("dnrm2 agrees with the pure kernel on ordinary, huge and tiny inputs") {
    for n <- Lengths do
      val x = lcg(n, 23L + n)
      assertNrm2Close(x, s"n=$n")
      assertNrm2Close(x.map(_ * 1e300), s"n=$n x1e300")
      assertNrm2Close(x.map(_ * 1e-300), s"n=$n x1e-300")
      assertNrm2Close(x.map(_ * 1e-170), s"n=$n x1e-170")
    // A finite result whose sum of squares overflows: sqrt(4096)·1e300 must be exact-ish.
    val big = Array.fill(4096)(1e300)
    assert(math.abs(VectorL1Kernels.dnrm2(4096, da(big), 0, 1) - 64e300) <= 4 * Eps * 64e300)
    val tiny = Array.fill(4096)(1e-300)
    assert(math.abs(VectorL1Kernels.dnrm2(4096, da(tiny), 0, 1) - 64e-300) <= 4 * Eps * 64e-300)
  }

  test("dnrm2 keeps the NaN/Inf/zero contract") {
    assertEquals(VectorL1Kernels.dnrm2(0, da(Array.empty), 0, 1), 0.0)
    assertEquals(VectorL1Kernels.dnrm2(40, da(Array.fill(40)(-0.0)), 0, 1), 0.0)
    for at <- Seq(0, 13, 39) do
      val inf = lcg(40, 1L)
      inf(at) = Double.NegativeInfinity
      assertEquals(VectorL1Kernels.dnrm2(40, da(inf), 0, 1), Double.PositiveInfinity, s"inf at $at")
      inf((at + 7) % 40) = Double.PositiveInfinity
      assertEquals(VectorL1Kernels.dnrm2(40, da(inf), 0, 1), Double.PositiveInfinity, s"two infs at $at")
      val nan = inf.clone()
      nan((at + 20) % 40) = Double.NaN
      assert(VectorL1Kernels.dnrm2(40, da(nan), 0, 1).isNaN, s"nan with infs at $at")
    // Strided input forwards to the pure recurrence.
    val x = lcg(99, 9L)
    assertEquals(VectorL1Kernels.dnrm2(33, da(x), 0, 3), DoubleKernels.dnrm2(33, da(x), 0, 3))
  }

  private def assertArgmax(x: Array[Double], label: String)(using munit.Location): Unit =
    for offset <- Seq(0, 1) if offset <= x.length do
      val n = x.length - offset
      assertEquals(
        VectorL1Kernels.dmaxIndex(n, da(x), offset, 1),
        DoubleKernels.dmaxIndex(n, da(x), offset, 1),
        s"$label offset=$offset"
      )

  test("dmaxIndex matches the pure kernel exactly: random, ties, NaN, ±0, ±Inf") {
    for n <- Lengths do
      val x = lcg(n, 31L + n)
      assertArgmax(x, s"random n=$n")
      // Coarse values: many ties, first occurrence must win.
      assertArgmax(x.map(v => math.rint(v * 2)), s"ties n=$n")
      assertArgmax(Array.fill(n)(Double.NegativeInfinity), s"all -Inf n=$n")
      assertArgmax(Array.tabulate(n)(i => if i % 2 == 0 then -0.0 else 0.0), s"±0 n=$n")
      assertArgmax(Array.tabulate(n)(i => if i % 3 == 0 then 0.0 else -0.0), s"0/-0 n=$n")
      if n > 0 then
        for at <- Seq(0, n / 3, n / 2, n - 1).distinct do
          val withNaN = x.clone()
          withNaN(at) = Double.NaN
          assertArgmax(withNaN, s"NaN at $at n=$n")
          val second = withNaN.clone()
          second((at + n / 2 + 1) % n) = Double.NaN
          assertArgmax(second, s"two NaNs n=$n")
          val withInf = x.clone()
          withInf(at) = Double.PositiveInfinity
          withInf((at + 1) % n) = Double.PositiveInfinity
          assertArgmax(withInf, s"tied +Inf n=$n")
    // Ties and NaNs at and across the kernel's 512-element block boundaries.
    for (a, b) <- Seq((100, 600), (511, 512), (512, 1023), (1023, 1024), (5, 1500)) do
      val x = lcg(1601, 41L)
      x(b) = 9.0
      x(a) = 9.0
      assertArgmax(x, s"tie $a/$b")
      x(b) = Double.NaN
      assertArgmax(x, s"max $a, NaN $b")
      x(a) = Double.NaN
      assertArgmax(x, s"NaNs $a/$b")
    // Ties straddling lanes and the scalar tail.
    for n <- Seq(17, 33, 65) do
      val x = Array.fill(n)(1.0)
      assertArgmax(x, s"all equal n=$n")
      x(n - 1) = 2.0
      x(n - 2) = 2.0
      assertArgmax(x, s"tail max n=$n")
    val strided = lcg(120, 77L)
    assertEquals(
      VectorL1Kernels.dmaxIndex(40, da(strided), 2, 3),
      DoubleKernels.dmaxIndex(40, da(strided), 2, 3)
    )
    assertEquals(VectorL1Kernels.dmaxIndex(0, da(Array.empty), 0, 1), -1)
  }

  private def ulpError(actual: Double, expected: Double): Double =
    if actual == expected then 0.0
    else if expected.isNaN then (if actual.isNaN then 0.0 else Double.PositiveInfinity)
    else if expected.isInfinite || actual.isInfinite || actual.isNaN then Double.PositiveInfinity
    else math.abs(actual - expected) / Math.ulp(expected)

  test("dexpInto: specials exact; within 2 ulp of StrictMath over [-745, 710] (measured value printed)") {
    val specials = Array(
      0.0, -0.0, 1.0, -1.0, Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity,
      709.782712893384, 709.79, 710.0, -745.1332191019411, -745.14, -746.0, -708.4, -720.0,
      Double.MinPositiveValue, -Double.MinPositiveValue, 1e-300, Double.MaxValue, -Double.MaxValue
    )
    val sOut = new Array[Double](specials.length)
    VectorL1Kernels.dexpInto(specials.length, da(specials), 0, 1, da(sOut), 0, 1)
    for i <- specials.indices do
      val e = StrictMath.exp(specials(i))
      if e.isNaN || e.isInfinite || e == 0.0 || e == 1.0 then
        // Boxed equality: NaN equals NaN, and +0.0 is distinguished from -0.0.
        assert(java.lang.Double.valueOf(sOut(i)).equals(e), s"exp(${specials(i)}) = ${sOut(i)} vs $e")
      else assert(ulpError(sOut(i), e) <= 1.0, s"exp(${specials(i)}) = ${sOut(i)} vs $e")

    val n = 2_000_003
    val lo = -745.0
    val hi = 710.0
    val x = Array.tabulate(n)(i => lo + (hi - lo) * i / (n - 1))
    val y = new Array[Double](n)
    VectorL1Kernels.dexpInto(n, da(x), 0, 1, da(y), 0, 1)
    var maxUlp = 0.0
    var maxUlpMath = 0.0
    var worstAt = 0.0
    var i = 0
    while i < n do
      val reference = StrictMath.exp(x(i))
      val err = ulpError(y(i), reference)
      if err > maxUlp then
        maxUlp = err
        worstAt = x(i)
      maxUlpMath = math.max(maxUlpMath, ulpError(Math.exp(x(i)), reference))
      i += 1
    println(
      f"[w21] dexpInto lanes=${VectorL1Kernels.lanes} max ulp vs StrictMath.exp over $n points " +
        f"in [$lo, $hi]: $maxUlp%.3f at x=$worstAt (Math.exp: $maxUlpMath%.3f)"
    )
    // The documented bound is <= 2 ulp from the true value. Against StrictMath (fdlibm,
    // itself within 1 ulp) allow 2 ulp: x86 SVML may legitimately sit 2 ulp from fdlibm.
    // The measured value is printed above (1.000 on aarch64/JDK 25).
    assert(maxUlp <= 2.0, s"SIMD exp max ulp $maxUlp at $worstAt")

    // `lanewise(EXP)` is not tier-stable: the interpreter/C1 path and the C2
    // intrinsic stub may round differently. Repeat the sweep (warming the kernel)
    // and require every rerun to stay within the same ulp bound of the first.
    var rerunDiffs = 0
    var rerunMaxUlp = 0.0
    var rep = 0
    while rep < 5 do
      val again = x.clone() // in place: y aliases x
      VectorL1Kernels.dexpInto(n, da(again), 0, 1, da(again), 0, 1)
      var j = 0
      while j < n do
        if !java.lang.Double.valueOf(again(j)).equals(y(j)) then
          rerunDiffs += 1
          val err = ulpError(again(j), StrictMath.exp(x(j)))
          rerunMaxUlp = math.max(rerunMaxUlp, err)
          assert(err <= 2.0, s"rerun $rep at x=${x(j)}")
        j += 1
      rep += 1
    println(f"[w21] dexpInto reruns differing bitwise from the first run: $rerunDiffs of ${5L * n}; max ulp among them $rerunMaxUlp%.3f")
  }
