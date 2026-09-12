package gale.parity

import breeze.linalg.cholesky
import gale.linalg.*
import gale.parity.ParitySupport.*

class BandedCholeskyParitySuite extends munit.FunSuite:
  test("packed SPD factors, multi-RHS solves and logdet agree with Breeze LAPACK"):
    for n <- List(5, 23, 80); b <- List(0, 1, 4); seed <- List(11L, 91L) do
      val random = new scala.util.Random(seed)
      val data = Array.fill(n, n)(0.0)
      for i <- 0 until n do
        data(i)(i) = 2.0 * b + 1.0
        for j <- math.max(0, i - b) until i do
          val value = random.nextDouble() * 2.0 - 1.0
          data(i)(j) = value
          data(j)(i) = value
      val bands = Matrix.tabulate(n, b + 1)((i, d) => if d <= i then data(i)(i - d) else 0.0)
      val f = BandedCholesky.factorLower(bands).orThrow
      val ba = breezeMatrix(data)
      val lower = cholesky(ba)
      for i <- 0 until n; d <- 0 to math.min(b, i) do
        assertEqualsDouble(f.lowerBands(i, d), lower(i, i - d), 3e-13)
      val rhs = Array.tabulate(n, 4)((i, j) => math.sin(i * 2.0 + j))
      assertMatClose(f.solve(galeMatrix(rhs)).orThrow, ba \ breezeMatrix(rhs), 2e-12, s"n=$n b=$b seed=$seed")
      assertEqualsDouble(f.logDet, 2.0 * (0 until n).map(i => math.log(lower(i, i))).sum, 3e-12)
