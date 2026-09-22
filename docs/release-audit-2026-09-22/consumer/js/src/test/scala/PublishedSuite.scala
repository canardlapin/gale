package externalconsumer
import gale.backend.{Backend, PureBackend}
import gale.linalg.*
import gale.laws.{BackendConformanceSuite, VecLaws}
import gale.sparse.*
import gale.sparse.direct.*
class PublishedPureBackendLaws extends BackendConformanceSuite:
  def backend: Backend = PureBackend
class PublishedSmokeSuite extends munit.FunSuite:
  test("banded analytic solve from published core"):
    val n = 7
    val bands = Matrix.tabulate(n, 2)((_, d) => if d == 0 then 2.0 else -1.0)
    val f = BandedCholesky.factorLower(bands).orThrow
    val x = f.solve(Vec.fill(n)(1.0)).orThrow
    for i <- 0 until n do assertEqualsDouble(x(i), (i + 1.0) * (n - i) / 2.0, 1e-12)
  test("sparse Cholesky from public artifact and explicit provider"):
    import gale.sparse.direct.pure.given
    val b = Sparse.coo(2, 2)
    b.add(0, 0, 4.0); b.add(0, 1, 1.0); b.add(1, 0, 1.0); b.add(1, 1, 3.0)
    val a = b.toCSR()
    val ws = SparseDirect.newWorkspace().orThrow
    val analysis = SparseDirect.analyze(a.pattern, SparseDirectFactorization.Cholesky, ws).orThrow
    val factor = SparseDirect.factor(analysis, a, ws).orThrow
    try
      val x = SparseDirect.solve(factor, Vec(6.0, 7.0), ws).orThrow.solution
      VecLaws.assertCloseRel(x, Vec(1.0, 2.0), 1e-12)
    finally
      factor.close(); analysis.close(); ws.close()
