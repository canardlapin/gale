package userland

import gale.linalg.*
import gale.sparse.*

class PropertyWrapperInteropSuite extends munit.FunSuite:
  private def matrixRows(matrix: Matrix[Double]): Int = matrix.rows
  private def sparseRows(matrix: SparseMatrix[Double]): Int = matrix.rows

  test("property wrappers and extensions are available through public wildcard imports") {
    val matrix = Matrix.dense(2, 2)(
      2.0, 0.5,
      0.5, 1.0
    )
    val symmetric: Symmetric[DMat] = matrix.verifySymmetric().orThrow
    val spd: PositiveDefinite[DMat] = matrix.verifyPositiveDefinite.orThrow
    val canonical: CanonicalSparse[CSR] =
      Sparse.coo(2, 2).add(0, 0, 1.0).add(1, 1, 2.0).toCSR().verifyCanonicalSparse.orThrow

    assertEquals(matrixRows(symmetric), 2)
    assertEquals(matrixRows(spd), 2)
    assertEquals(sparseRows(canonical), 2)
    assertEquals((canonical * Vec(3.0, 4.0)).toSeq, Seq(3.0, 8.0))
  }
