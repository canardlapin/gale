package gale.bench

import breeze.linalg.DenseMatrix as BDM
import breeze.linalg.eigSym
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import gale.spectral.*
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Dense symmetric eigendecomposition (values + vectors) paired benchmark:
  * gale `Eigen.eigSymmetric(A, All)` (tridiagonal QL/QR) vs breeze `eigSym(A)`
  * (LAPACK `dsyev`). Both compute the full spectrum with eigenvectors.
  *
  * Sizes stay modest because dense eigen has a larger constant than the pure BLAS-3
  * kernels above.
  */
@BenchmarkMode(Array(Mode.Throughput))
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Fork(2)
@State(Scope.Thread)
class SymEigenBreezeJmh:
  @Param(Array("16", "64", "128"))
  var n: Int = 0

  private var gA: DMat        = _
  private var bA: BDM[Double] = _

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    val data = symmetric(n, 600L)
    gA = galeMatrix(data)
    bA = breezeMatrix(data)

  @Benchmark def galeEigSym(bh: Blackhole): Unit =
    bh.consume(Eigen.eigSymmetric(gA, EigenSelection.All, EigenVectors.Right))

  @Benchmark def breezeEigSym(bh: Blackhole): Unit =
    bh.consume(eigSym(bA))
