package gale.backend.jvm.vector

import jdk.incubator.vector.DoubleVector

/** Benchmark-only frozen copy of the row-major GEMV kernel as it stood before the
  * HC-4 hill climb (main 90c74a3): a four-row tile with two accumulators per row and
  * a scalar column tail. `gale.bench.GemvHillClimbJmh` times it against the live
  * kernel in the same fork. Nothing in the library calls it; it lives here, not in
  * the benchmarks module, because the JMH annotation generator verifies every bench
  * class in a JVM without `jdk.incubator.vector`. Do not edit; delete it once the
  * HC-4 verdict is recorded.
  */
private[vector] object VectorGemvBaseline:
  def gemvRowMajor(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      x: Array[Double],
      xOffset: Int,
      beta: Double,
      y: Array[Double],
      yOffset: Int
  ): Unit =
    val species = DoubleVector.SPECIES_PREFERRED
    val lanes = species.length()
    val colBound = species.loopBound(cols)
    val pairBound = cols - cols % (2 * lanes)
    val rowBound = rows - (rows & 3)
    var i = 0
    while i < rowBound do
      val aRow0 = aOffset + i * aRowStride
      val aRow1 = aRow0 + aRowStride
      val aRow2 = aRow1 + aRowStride
      val aRow3 = aRow2 + aRowStride
      var s0 = DoubleVector.zero(species)
      var s1 = DoubleVector.zero(species)
      var s2 = DoubleVector.zero(species)
      var s3 = DoubleVector.zero(species)
      // Two accumulators per row: eight independent vector FMA chains.
      var t0 = DoubleVector.zero(species)
      var t1 = DoubleVector.zero(species)
      var t2 = DoubleVector.zero(species)
      var t3 = DoubleVector.zero(species)
      var j = 0
      while j < pairBound do
        val xv = DoubleVector.fromArray(species, x, xOffset + j)
        val xw = DoubleVector.fromArray(species, x, xOffset + j + lanes)
        s0 = DoubleVector.fromArray(species, a, aRow0 + j).fma(xv, s0)
        s1 = DoubleVector.fromArray(species, a, aRow1 + j).fma(xv, s1)
        s2 = DoubleVector.fromArray(species, a, aRow2 + j).fma(xv, s2)
        s3 = DoubleVector.fromArray(species, a, aRow3 + j).fma(xv, s3)
        t0 = DoubleVector.fromArray(species, a, aRow0 + j + lanes).fma(xw, t0)
        t1 = DoubleVector.fromArray(species, a, aRow1 + j + lanes).fma(xw, t1)
        t2 = DoubleVector.fromArray(species, a, aRow2 + j + lanes).fma(xw, t2)
        t3 = DoubleVector.fromArray(species, a, aRow3 + j + lanes).fma(xw, t3)
        j += 2 * lanes
      s0 = s0.add(t0)
      s1 = s1.add(t1)
      s2 = s2.add(t2)
      s3 = s3.add(t3)
      while j < colBound do
        val xv = DoubleVector.fromArray(species, x, xOffset + j)
        s0 = DoubleVector.fromArray(species, a, aRow0 + j).fma(xv, s0)
        s1 = DoubleVector.fromArray(species, a, aRow1 + j).fma(xv, s1)
        s2 = DoubleVector.fromArray(species, a, aRow2 + j).fma(xv, s2)
        s3 = DoubleVector.fromArray(species, a, aRow3 + j).fma(xv, s3)
        j += species.length()
      var sum0 = s0.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      var sum1 = s1.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      var sum2 = s2.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      var sum3 = s3.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      while j < cols do
        val xj = x(xOffset + j)
        sum0 = Math.fma(a(aRow0 + j), xj, sum0)
        sum1 = Math.fma(a(aRow1 + j), xj, sum1)
        sum2 = Math.fma(a(aRow2 + j), xj, sum2)
        sum3 = Math.fma(a(aRow3 + j), xj, sum3)
        j += 1
      y(yOffset + i) = combine(alpha, sum0, beta, y(yOffset + i))
      y(yOffset + i + 1) = combine(alpha, sum1, beta, y(yOffset + i + 1))
      y(yOffset + i + 2) = combine(alpha, sum2, beta, y(yOffset + i + 2))
      y(yOffset + i + 3) = combine(alpha, sum3, beta, y(yOffset + i + 3))
      i += 4

    while i < rows do
      val aRow = aOffset + i * aRowStride
      var sum = DoubleVector.zero(species)
      var j = 0
      while j < colBound do
        val av = DoubleVector.fromArray(species, a, aRow + j)
        val xv = DoubleVector.fromArray(species, x, xOffset + j)
        sum = av.fma(xv, sum)
        j += species.length()
      var scalar = sum.reduceLanes(jdk.incubator.vector.VectorOperators.ADD)
      while j < cols do
        scalar = Math.fma(a(aRow + j), x(xOffset + j), scalar)
        j += 1
      y(yOffset + i) = combine(alpha, scalar, beta, y(yOffset + i))
      i += 1


  private inline def combine(alpha: Double, sum: Double, beta: Double, prior: Double): Double =
    if beta == 0.0 then alpha * sum
    else Math.fma(beta, prior, alpha * sum)

  /** The live row-major tile without the eight-lane row-sum tree (per-row
    * `reduceLanes` and scalar stores instead), so a benchmark can attribute the
    * tree's share of the HC-4 change. Identical to the live kernel off eight lanes.
    */
  def gemvRowMajorNoTree(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      x: Array[Double],
      xOffset: Int,
      beta: Double,
      y: Array[Double],
      yOffset: Int
  ): Unit =
    VectorDenseDoubleKernel.gemvRowMajorTiles(
      DoubleVector.SPECIES_PREFERRED, false, rows, cols, alpha, a, aOffset, aRowStride, x, xOffset, beta, y, yOffset
    )
