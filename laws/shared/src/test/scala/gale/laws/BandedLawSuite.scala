package gale.laws

import gale.linalg.*
import gale.sparse.*

/** Banded matvec / transpose-matvec laws against the dense reference, via the
  * [[SparseLaws]] bundle.
  */
class BandedLawSuite extends munit.FunSuite:
  private val dense = Matrix.dense(4, 4)(
    2.0, -1.0, 0.0, 0.0,
    -1.0, 2.0, -1.0, 0.0,
    0.0, -1.0, 2.0, -1.0,
    0.0, 0.0, -1.0, 2.0
  )
  private val banded = Sparse.banded(dense)
  private val x = Vec(1.0, -2.0, 3.0, -4.0)

  test("banded matvec matches the dense reference") {
    SparseLaws.matvecMatchesDense(banded, dense, x)
  }

  test("banded transpose-matvec matches the dense transpose") {
    SparseLaws.transposeMatvecMatchesDense(banded, dense, x)
  }
