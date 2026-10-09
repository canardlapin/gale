package gale.bench

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Private copies of the pure kernels as they stood before the L1 hill climb
  * (breeze/integration d102321), shared by the JVM JMH and Scala.js runners so
  * each run pairs the old kernel with the live one on the same data.
  */
object L1HillClimbRef:
  import gale.platform.PlatformMath.fma

  def refDaxpy(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Unit =
    if xStride == 1 && yStride == 1 then
      val limit = n - (n & 3)
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < limit do
        y(yi) = fma(alpha, x(xi), y(yi))
        y(yi + 1) = fma(alpha, x(xi + 1), y(yi + 1))
        y(yi + 2) = fma(alpha, x(xi + 2), y(yi + 2))
        y(yi + 3) = fma(alpha, x(xi + 3), y(yi + 3))
        xi += 4
        yi += 4
        i += 4
      while i < n do
        y(yi) = fma(alpha, x(xi), y(yi))
        xi += 1
        yi += 1
        i += 1
    else
      var i = 0
      var xi = xOffset
      var yi = yOffset
      while i < n do
        y(yi) = fma(alpha, x(xi), y(yi))
        xi += xStride
        yi += yStride
        i += 1

  def refDdot(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int): Double =
    var acc0 = 0.0
    var acc1 = 0.0
    var acc2 = 0.0
    var acc3 = 0.0
    val limit = n - (n & 3)
    var i = 0
    var xi = xOffset
    var yi = yOffset
    while i < limit do
      acc0 = fma(x(xi), y(yi), acc0)
      acc1 = fma(x(xi + 1), y(yi + 1), acc1)
      acc2 = fma(x(xi + 2), y(yi + 2), acc2)
      acc3 = fma(x(xi + 3), y(yi + 3), acc3)
      xi += 4
      yi += 4
      i += 4
    var acc = (acc0 + acc1) + (acc2 + acc3)
    while i < n do
      acc = fma(x(xi), y(yi), acc)
      xi += 1
      yi += 1
      i += 1
    acc

  def refDnrm2(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Double =
    if n < 1 then 0.0
    else if n == 1 then math.abs(x(xOffset))
    else
      var scale = 0.0
      var ssq = 1.0
      var i = 0
      var xi = xOffset
      while i < n do
        val value = x(xi)
        if value != 0.0 then
          val abs = math.abs(value)
          if scale < abs then
            val ratio = scale / abs
            ssq = 1.0 + ssq * ratio * ratio
            scale = abs
          else
            val ratio = abs / scale
            ssq += ratio * ratio
        xi += xStride
        i += 1
      scale * math.sqrt(ssq)

  def refMaxIndex(n: Int, x: DoubleArray, xOffset: Int, xStride: Int): Int =
    var b0 = Double.NegativeInfinity
    var b1 = Double.NegativeInfinity
    var b2 = Double.NegativeInfinity
    var b3 = Double.NegativeInfinity
    var i0 = Int.MaxValue
    var i1 = Int.MaxValue
    var i2 = Int.MaxValue
    var i3 = Int.MaxValue
    val limit = n - (n & 3)
    var i = 0
    var xi = xOffset
    var nan = -1
    while i < limit && nan < 0 do
      val v0 = x(xi)
      val v1 = x(xi + 1)
      val v2 = x(xi + 2)
      val v3 = x(xi + 3)
      if v0 > b0 then
        b0 = v0; i0 = i
      else if v0 != v0 then nan = i
      if nan < 0 then
        if v1 > b1 then
          b1 = v1; i1 = i + 1
        else if v1 != v1 then nan = i + 1
      if nan < 0 then
        if v2 > b2 then
          b2 = v2; i2 = i + 2
        else if v2 != v2 then nan = i + 2
      if nan < 0 then
        if v3 > b3 then
          b3 = v3; i3 = i + 3
        else if v3 != v3 then nan = i + 3
      xi += 4
      i += 4
    if nan >= 0 then nan
    else
      var best = b0
      var bestIndex = i0
      if b1 > best || (b1 == best && i1 < bestIndex) then
        best = b1; bestIndex = i1
      if b2 > best || (b2 == best && i2 < bestIndex) then
        best = b2; bestIndex = i2
      if b3 > best || (b3 == best && i3 < bestIndex) then
        best = b3; bestIndex = i3
      while i < n && nan < 0 do
        val v = x(xi)
        if v > best then
          best = v; bestIndex = i
        else if v != v then nan = i
        xi += 1
        i += 1
      if nan >= 0 then nan
      else if bestIndex == Int.MaxValue then 0
      else bestIndex

  def refStreamMax(rows: Int, cols: Int, data: DoubleArray, out: DoubleArray): Unit =
    var line = 0
    while line < cols do
      out(line) = data(line)
      line += 1
    var k = 1
    while k < rows do
      val start = k * cols
      line = 0
      while line < cols do
        val current = out(line)
        val value = data(start + line)
        if value > current || (value != value && current == current) then out(line) = value
        line += 1
      k += 1

  def refDadd(n: Int, x: DoubleArray, xOffset: Int, xStride: Int, y: DoubleArray, yOffset: Int, yStride: Int, out: DoubleArray, outOffset: Int, outStride: Int): Unit =
    var i = 0
    var xi = xOffset
    var yi = yOffset
    var oi = outOffset
    while i < n do
      out(oi) = x(xi) + y(yi)
      xi += xStride
      yi += yStride
      oi += outStride
      i += 1

  private inline def refMap(n: Int, x: DoubleArray, y: DoubleArray)(inline f: Double => Double): Unit =
    val limit = n - (n & 3)
    var i = 0
    var xi = 0
    var yi = 0
    while i < limit do
      y(yi) = f(x(xi))
      y(yi + 1) = f(x(xi + 1))
      y(yi + 2) = f(x(xi + 2))
      y(yi + 3) = f(x(xi + 3))
      xi += 4
      yi += 4
      i += 4
    while i < n do
      y(yi) = f(x(xi))
      xi += 1
      yi += 1
      i += 1

  def refSigmoidInto(n: Int, x: DoubleArray, y: DoubleArray): Unit =
    refMap(n, x, y) { v =>
      val t = math.exp(-math.abs(v))
      if v >= 0.0 then 1.0 / (1.0 + t) else t / (1.0 + t)
    }

  def refExpInto(n: Int, x: DoubleArray, y: DoubleArray): Unit =
    refMap(n, x, y)(math.exp)

  def refLogSumExp(n: Int, x: DoubleArray): Double =
    val m = x(refMaxIndex(n, x, 0, 1))
    if m.isNaN || m.isInfinite then m
    else
      var acc0 = 0.0
      var acc1 = 0.0
      var acc2 = 0.0
      var acc3 = 0.0
      val limit = n - (n & 3)
      var i = 0
      var xi = 0
      while i < limit do
        acc0 += math.exp(x(xi) - m)
        acc1 += math.exp(x(xi + 1) - m)
        acc2 += math.exp(x(xi + 2) - m)
        acc3 += math.exp(x(xi + 3) - m)
        xi += 4
        i += 4
      var acc = (acc0 + acc1) + (acc2 + acc3)
      while i < n do
        acc += math.exp(x(xi) - m)
        xi += 1
        i += 1
      m + math.log(acc)
