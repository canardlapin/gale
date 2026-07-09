package gale.laws

import gale.linalg.*
import gale.sparse.*
import org.scalacheck.Prop.forAll

/** Sparse round-trip and matvec laws, driven through [[SparseLaws]]. Lives in the
  * laws module: public API only.
  */
class MatrixMarketLawSuite extends munit.ScalaCheckSuite:
  property("CSR matvec equals dense matvec for small generated matrices") {
    forAll { (
        a00: Int,
        a01: Int,
        a02: Int,
        a10: Int,
        a11: Int,
        a12: Int,
        a20: Int,
        a21: Int,
        a22: Int,
        x0: Int,
        x1: Int,
        x2: Int
    ) =>
      val raw = Array(a00, a01, a02, a10, a11, a12, a20, a21, a22).map(v => (v % 5).toDouble)
      val x = Vec((x0 % 5).toDouble, (x1 % 5).toDouble, (x2 % 5).toDouble)
      val dense = Matrix.dense(3, 3, raw.toIndexedSeq)
      val builder = Sparse.coo(3, 3)
      var row = 0
      while row < 3 do
        var col = 0
        while col < 3 do
          val value = raw(row * 3 + col)
          if value != 0.0 then builder.add(row, col, value)
          col += 1
        row += 1
      val sparse = builder.toCSR()

      SparseLaws.matvecMatchesDense(sparse, dense, x)
    }
  }
