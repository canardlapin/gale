package externalconsumer
import gale.linalg.*
import gale.solvers.*
class GmresContractProbe extends munit.FunSuite:
  for mode <- List(ToleranceMode.Absolute, ToleranceMode.RelativeToRhs) do
    test(s"GMRES must satisfy documented true residual under Jacobi scaling: $mode"):
      val a = Matrix.dense(1, 1)(1e12)
      val b = Vec(1.0)
      val result = gmres(a, b, SolverConfig(tolerance = 1e-10), Preconditioner.Jacobi(a), toleranceMode = mode)
      val trueResidual = (b - a * result.x).norm2
      println(s"AUDIT GMRES mode=$mode converged=${result.converged} iterations=${result.iterations} reportedResidual=${result.residual} x=${result.x(0)} trueResidual=$trueResidual")
      assert(!result.converged || trueResidual <= 1e-10, s"False convergence: true residual $trueResidual exceeds 1e-10")
