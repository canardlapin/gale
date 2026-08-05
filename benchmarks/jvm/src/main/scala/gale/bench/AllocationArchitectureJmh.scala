package gale.bench

import scala.compiletime.uninitialized

import gale.backend.PureBackend
import gale.linalg.DMat
import gale.linalg.DMatBuilder
import gale.linalg.DVec
import gale.linalg.DenseWorkspace
import gale.linalg.MutableDVec
import gale.linalg.Vec
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*
import gale.sparse.CSR
import gale.sparse.CSRProductPlan
import gale.sparse.CSRUnionPlan
import gale.sparse.CSRValuesDestination
import gale.sparse.Sparse
import gale.spectral.Eigen
import gale.spectral.EigenSelection
import gale.spectral.EigenVectors
import java.util.concurrent.TimeUnit
import org.openjdk.jmh.annotations.*
import org.openjdk.jmh.infra.Blackhole

/** Allocation baseline for the allocation-control architecture epic.
  *
  * Run with JMH's `-prof gc`. The paired allocating/reuse methods deliberately consume the same logical work so
  * normalized bytes/op identify storage that a destination or workspace contract can actually remove. Input
  * construction is confined to `@Setup` and is never part of a timed invocation.
  */
@BenchmarkMode(Array(Mode.AverageTime))
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 3, time = 300, timeUnit = TimeUnit.MILLISECONDS)
@Fork(1)
@State(Scope.Thread)
class AllocationArchitectureJmh:
  @Param(Array("64", "128"))
  var n: Int = 0

  private var denseA: DMat = uninitialized
  private var denseB: DMat = uninitialized
  private var symmetric: DMat = uninitialized
  private var denseGemmDestination: DMatBuilder = uninitialized
  private var denseCombinationDestination: DMatBuilder = uninitialized
  private var qrWorkspace: DenseWorkspace = uninitialized
  private var symmetricWorkspace: DenseWorkspace = uninitialized
  private var sparseA: CSR = uninitialized
  private var sparseB: CSR = uninitialized
  private var sparseCanonicalInput: CSR = uninitialized
  private var sparseCanonicalWorkspace: DenseWorkspace = uninitialized
  private var sparseUnionPlan: CSRUnionPlan = uninitialized
  private var sparseUnionDestination: CSRValuesDestination = uninitialized
  private var sparseProductPlan: CSRProductPlan = uninitialized
  private var sparseProductDestination: CSRValuesDestination = uninitialized
  private var sparseX: DVec = uninitialized
  private var sparseDestination: MutableDVec = uninitialized

  @Setup(Level.Trial)
  def setupTrial(): Unit =
    denseA = DMat.tabulate(n, n): (row, col) =>
      if row == col then n.toDouble + 1.0
      else ((row * 17 + col * 13) % 23).toDouble / 23.0
    denseB = DMat.tabulate(n, n): (row, col) =>
      if row == col then 2.0
      else ((row * 7 - col * 11) % 19).toDouble / 19.0
    symmetric = DMat.tabulate(n, n): (row, col) =>
      val lo = math.min(row, col)
      val hi = math.max(row, col)
      if lo == hi then n.toDouble + 2.0
      else ((lo * 5 + hi * 3) % 17).toDouble / 34.0
    denseGemmDestination = DMatBuilder.zeros(n, n)
    denseCombinationDestination = DMatBuilder.zeros(n, n)
    qrWorkspace = DenseWorkspace.forQR(n, n)
    symmetricWorkspace = DenseWorkspace.forRequirement(
      Eigen.symmetricScratchRequirement(n, EigenVectors.ValuesOnly).toOption.get
    )

    val aBuilder = Sparse.coo(n, n)
    val bBuilder = Sparse.coo(n, n)
    var row = 0
    while row < n do
      val _ = aBuilder.add(row, row, 3.0)
      val _ = aBuilder.add(row, (row + 1) % n, -0.5)
      val _ = bBuilder.add(row, row, -1.0)
      val _ = bBuilder.add(row, (row + n - 1) % n, 0.25)
      row += 1
    sparseA = aBuilder.toCSR()
    sparseB = bBuilder.toCSR()
    sparseUnionPlan = CSRUnionPlan.analyze(sparseA.pattern, sparseB.pattern).toOption.get
    sparseUnionDestination = sparseUnionPlan.newDestination()
    sparseProductPlan = CSRProductPlan.analyze(sparseA.pattern, sparseB.pattern).toOption.get
    sparseProductDestination = sparseProductPlan.newDestination()
    val rowPtr = IndexArray.alloc(n + 1)
    val colIdx = IndexArray.alloc(4 * n)
    val sparseValues = DoubleArray.alloc(4 * n)
    row = 0
    while row < n do
      val start = 4 * row
      rowPtr(row) = start
      colIdx(start) = (row + 1) % n
      sparseValues(start) = 0.75
      colIdx(start + 1) = row
      sparseValues(start + 1) = 3.0
      colIdx(start + 2) = (row + 1) % n
      sparseValues(start + 2) = -0.25
      colIdx(start + 3) = (row + n - 1) % n
      sparseValues(start + 3) = -0.5
      row += 1
    rowPtr(n) = 4 * n
    sparseCanonicalInput = new CSR(n, n, rowPtr, colIdx, sparseValues)
    sparseCanonicalWorkspace = DenseWorkspace.forRequirement(
      sparseCanonicalInput.canonicalizeScratchRequirement.toOption.get
    )
    sparseX = Vec.tabulate(n)(i => (i % 11).toDouble - 5.0)
    sparseDestination = MutableDVec.zeros(n)

  /** One owned output matrix is unavoidable through the current high-level API. */
  @Benchmark
  def denseGemmAllocating(blackhole: Blackhole): Unit =
    blackhole.consume(denseA.*(denseB)(using PureBackend))

  @Benchmark
  def denseGemmReusedDestination(): Double =
    denseA.gemmInto(denseB, denseGemmDestination)(using PureBackend)
    denseGemmDestination(0, 0)

  /** Current public construction route for `0.5 * (A + B)`: the sum and scaled result are separately owned. This is the
    * baseline for a future fused AXPBY destination path, not a proposed public implementation.
    */
  @Benchmark
  def denseAddScalePipeline(blackhole: Blackhole): Unit =
    val sum = denseA + denseB
    val scaled = DMatBuilder.zeros(n, n)
    var row = 0
    while row < n do
      var col = 0
      while col < n do
        scaled(row, col) = 0.5 * sum(row, col)
        col += 1
      row += 1
    blackhole.consume(scaled.result())

  @Benchmark
  def denseLinearCombinationReusedDestination(): Double =
    denseA.linearCombinationInto(denseB, denseCombinationDestination, alpha = 0.5, beta = 0.5)
    denseCombinationDestination(0, 0)

  @Benchmark
  def qrFreshWorkspace(blackhole: Blackhole): Unit =
    blackhole.consume(denseA.qr(using PureBackend))

  @Benchmark
  def qrReusedWorkspace(blackhole: Blackhole): Unit =
    blackhole.consume(denseA.qrWith(qrWorkspace)(using PureBackend))

  @Benchmark
  def denseSymmetricEigenvalues(blackhole: Blackhole): Unit =
    blackhole.consume(Eigen.eigSymmetric(symmetric, EigenSelection.All, EigenVectors.ValuesOnly))

  @Benchmark
  def denseSymmetricEigenvaluesReusedWorkspace(blackhole: Blackhole): Unit =
    blackhole.consume(
      Eigen.eigSymmetricWith(symmetric, EigenSelection.All, EigenVectors.ValuesOnly, symmetricWorkspace)
    )

  @Benchmark
  def sparseAddAllocating(blackhole: Blackhole): Unit =
    blackhole.consume(sparseA + sparseB)

  @Benchmark
  def sparseUnionPlanAnalysis(blackhole: Blackhole): Unit =
    blackhole.consume(CSRUnionPlan.analyze(sparseA.pattern, sparseB.pattern))

  @Benchmark
  def sparseUnionPlannedReplay(): Double =
    sparseUnionPlan.evaluateInto(sparseA, sparseB, sparseUnionDestination) match
      case Left(error) => throw error
      case Right(())   => sparseUnionDestination(0, 0)

  @Benchmark
  def sparseProductPlanAnalysis(blackhole: Blackhole): Unit =
    blackhole.consume(CSRProductPlan.analyze(sparseA.pattern, sparseB.pattern))

  @Benchmark
  def sparseProductPlannedReplay(): Double =
    sparseProductPlan.evaluateInto(sparseA, sparseB, sparseProductDestination) match
      case Left(error) => throw error
      case Right(())   => sparseProductDestination(0, 0)

  @Benchmark
  def sparseMapAllocating(blackhole: Blackhole): Unit =
    blackhole.consume(sparseA.mapValues(_ * 1.0001))

  @Benchmark
  def sparseCanonicalizeFreshWorkspace(blackhole: Blackhole): Unit =
    blackhole.consume(sparseCanonicalInput.canonicalize)

  @Benchmark
  def sparseCanonicalizeReusedWorkspace(blackhole: Blackhole): Unit =
    blackhole.consume(sparseCanonicalInput.canonicalizeWith(sparseCanonicalWorkspace))

  @Benchmark
  def sparseMatvecAllocating(blackhole: Blackhole): Unit =
    blackhole.consume(sparseA.*(sparseX)(using PureBackend))

  @Benchmark
  def sparseMatvecReused(): Double =
    sparseA.mulInto(sparseX, sparseDestination)
    sparseDestination(0)
