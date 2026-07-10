package gale.bench

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.DenseVector as BDV
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Matrix
import gale.linalg.Vec

/** Seeded input generators shared by the paired gale-vs-Breeze JMH benchmarks.
  *
  * Every input is generated once as a plain array and then handed to '''both'''
  * libraries, so a gale benchmark and its Breeze twin run on identical data at the
  * same `@Param` size. Generation happens in `@Setup` (never in a timed method), so
  * these allocations do not count against per-invocation cost.
  *
  * Self-contained on purpose (it does not reach into the parity module): the two
  * harnesses evolve independently.
  */
object BreezeBenchData:

  /** Deterministic pseudo-random `[-1, 1)` sequence — a tiny LCG so a given seed
    * yields the same data on every fork.
    */
  private def fill(out: Array[Double], seed: Long): Unit =
    var state = seed * 6364136223846793005L + 1442695040888963407L
    var i = 0
    while i < out.length do
      state = state * 6364136223846793005L + 1442695040888963407L
      out(i) = ((state >>> 11).toDouble / (1L << 53).toDouble) * 2.0 - 1.0
      i += 1

  def vectorData(n: Int, seed: Long): Array[Double] =
    val out = new Array[Double](n)
    fill(out, seed)
    out

  def matrixData(rows: Int, cols: Int, seed: Long): Array[Array[Double]] =
    Array.tabulate(rows)(i => vectorData(cols, seed + i * 0x9e3779b9L))

  /** Strictly diagonally dominant (well-conditioned, nonsingular): each diagonal is
    * set to its row's absolute off-diagonal sum plus one.
    */
  def diagonallyDominant(n: Int, seed: Long): Array[Array[Double]] =
    val a = matrixData(n, n, seed)
    var i = 0
    while i < n do
      var sum = 0.0
      var j = 0
      while j < n do
        if j != i then sum += math.abs(a(i)(j))
        j += 1
      a(i)(i) = sum + 1.0
      i += 1
    a

  /** Symmetric positive-definite `B Bᵀ + n·I`, built exactly symmetric. */
  def spd(n: Int, seed: Long): Array[Array[Double]] =
    val b = matrixData(n, n, seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        var s = 0.0
        var k = 0
        while k < n do
          s += b(i)(k) * b(j)(k)
          k += 1
        if i == j then s += n.toDouble
        a(i)(j) = s
        a(j)(i) = s
        j += 1
      i += 1
    a

  /** Symmetric matrix with entries in `[-1, 1)` (lower triangle mirrored). */
  def symmetric(n: Int, seed: Long): Array[Array[Double]] =
    val src = matrixData(n, n, seed)
    val a = Array.ofDim[Double](n, n)
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        a(i)(j) = src(i)(j)
        a(j)(i) = src(i)(j)
        j += 1
      i += 1
    a

  def galeMatrix(data: Array[Array[Double]]): DMat =
    Matrix.tabulate(data.length, if data.isEmpty then 0 else data(0).length)((i, j) => data(i)(j))

  def breezeMatrix(data: Array[Array[Double]]): BDM[Double] =
    BDM.tabulate(data.length, if data.isEmpty then 0 else data(0).length)((i, j) => data(i)(j))

  def galeVector(data: Array[Double]): DVec =
    Vec.tabulate(data.length)(i => data(i))

  def breezeVector(data: Array[Double]): BDV[Double] =
    BDV(data.clone())
