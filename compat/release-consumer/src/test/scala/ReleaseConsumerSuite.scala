import gale.linalg.*
import gale.laws.*
import gale.solvers.*

class ReleaseConsumerSuite extends munit.FunSuite:
  test("published core and laws solve and validate a dense workflow") {
    val design = Matrix(4, 2)(1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0)
    val fit = design.leastSquares(Vec(1.0, 3.0, 5.0, 7.0)).orThrow
    VecLaws.assertClose(fit, Vec(1.0, 2.0), 1e-12)
    MatrixLaws.transposeIsInvolution(design)
    val a = Matrix(1, 1)(1e12)
    val result = gmres(a, Vec(1.0), preconditioner = Preconditioner.Jacobi(a))
    SolverLaws.solvesSystem(a, result.orThrow, Vec(1.0), 1e-10)
  }
