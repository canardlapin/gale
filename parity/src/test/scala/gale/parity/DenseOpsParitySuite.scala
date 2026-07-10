package gale.parity

import breeze.linalg.*
import gale.linalg.*
import gale.parity.ParitySupport.*

/** Dense BLAS-level parity: gale's `DMat`/`DVec` products and combinations versus
  * `breeze.linalg.DenseMatrix`/`DenseVector` on bit-identical random data.
  *
  * These are same-arithmetic operations (both ultimately a sum of products), so
  * they must agree to a tight elementwise relative tolerance; the only expected
  * difference is summation order (gale's straight loops vs Breeze's netlib BLAS),
  * which stays far below `1e-12` at these sizes.
  */
class DenseOpsParitySuite extends munit.FunSuite:

  private val tol  = 1e-12
  private val sizes = List(5, 13, 27, 50)
  private val seeds = List(1L, 2L, 3L)

  test("matrix–vector product A·x") {
    for size <- sizes; seed <- seeds do
      val aData = matrixData(size, size, seed)
      val xData = vectorData(size, seed * 7 + 1)
      val ga    = galeMatrix(aData)
      val ba    = breezeMatrix(aData)
      val gx    = galeVector(xData)
      val bx    = breezeVector(xData)
      assertVecClose(ga * gx, ba * bx, tol, s"A·x size=$size seed=$seed")
  }

  test("matrix–matrix product A·B (rectangular)") {
    for size <- sizes; seed <- seeds do
      val m = size
      val k = size + 3
      val n = size - 1
      val aData = matrixData(m, k, seed)
      val bData = matrixData(k, n, seed * 13 + 5)
      val ga    = galeMatrix(aData)
      val gb    = galeMatrix(bData)
      val ba    = breezeMatrix(aData)
      val bb    = breezeMatrix(bData)
      assertMatClose(ga * gb, ba * bb, tol, s"A·B ${m}x${k}·${k}x$n seed=$seed")
  }

  test("transpose–matrix product Aᵀ·B") {
    for size <- sizes; seed <- seeds do
      // A is m×k so Aᵀ is k×m; B is m×n; Aᵀ·B is k×n.
      val m = size
      val k = size + 2
      val n = size + 1
      val aData = matrixData(m, k, seed)
      val bData = matrixData(m, n, seed * 17 + 2)
      val ga    = galeMatrix(aData)
      val gb    = galeMatrix(bData)
      val ba    = breezeMatrix(aData)
      val bb    = breezeMatrix(bData)
      assertMatClose(ga.t * gb, ba.t * bb, tol, s"Aᵀ·B size=$size seed=$seed")
  }

  test("transpose–vector product Aᵀ·x") {
    for size <- sizes; seed <- seeds do
      val m = size + 4
      val n = size
      val aData = matrixData(m, n, seed)
      val xData = vectorData(m, seed * 19 + 3)
      val ga    = galeMatrix(aData)
      val ba    = breezeMatrix(aData)
      val gx    = galeVector(xData)
      val bx    = breezeVector(xData)
      assertVecClose(ga.t * gx, ba.t * bx, tol, s"Aᵀ·x size=$size seed=$seed")
  }

  test("vector dot product x·y") {
    for size <- sizes; seed <- seeds do
      val xData = vectorData(size, seed)
      val yData = vectorData(size, seed * 23 + 9)
      val gx    = galeVector(xData)
      val gy    = galeVector(yData)
      val bx    = breezeVector(xData)
      val by    = breezeVector(yData)
      assertScalarClose(gx.dot(gy), bx.dot(by), tol, s"x·y size=$size seed=$seed")
  }

  test("axpy y := α·x + y") {
    for size <- sizes; seed <- seeds do
      val alpha = 2.5
      val xData = vectorData(size, seed)
      val yData = vectorData(size, seed * 29 + 4)
      val gy    = galeVector(yData).mutableCopy
      gy.axpyInPlace(alpha, galeVector(xData))
      val expected = breezeVector(yData) + (breezeVector(xData) * alpha)
      assertVecClose(gy.asVec, expected, tol, s"axpy size=$size seed=$seed")
  }

  test("matrix add/subtract A ± B") {
    for size <- sizes; seed <- seeds do
      val aData = matrixData(size, size, seed)
      val bData = matrixData(size, size, seed * 31 + 6)
      val ga    = galeMatrix(aData)
      val gb    = galeMatrix(bData)
      val ba    = breezeMatrix(aData)
      val bb    = breezeMatrix(bData)
      assertMatClose(ga + gb, ba + bb, tol, s"A+B size=$size seed=$seed")
      assertMatClose(ga - gb, ba - bb, tol, s"A-B size=$size seed=$seed")
  }

  test("vector scale, add, subtract") {
    for size <- sizes; seed <- seeds do
      val xData = vectorData(size, seed)
      val yData = vectorData(size, seed * 37 + 8)
      val gx    = galeVector(xData)
      val gy    = galeVector(yData)
      val bx    = breezeVector(xData)
      val by    = breezeVector(yData)
      assertVecClose(gx * 3.0, bx * 3.0, tol, s"3·x size=$size seed=$seed")
      assertVecClose(gx + gy, bx + by, tol, s"x+y size=$size seed=$seed")
      assertVecClose(gx - gy, bx - by, tol, s"x-y size=$size seed=$seed")
  }
