package gale.bench

import gale.backend.PureBackend
import gale.linalg.*
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*
import gale.sparse.CSR
import gale.sparse.CSRProductPlan
import gale.sparse.CSRUnionPlan
import gale.sparse.Sparse
import gale.spectral.Eigen
import gale.spectral.EigenSelection
import gale.spectral.EigenVectors
import scala.scalajs.js

object DenseKernelJsBench:
  def main(args: Array[String]): Unit =
    val x = Vec.tabulate(65536)(i => ((i % 17) - 8).toDouble * 0.125)
    val y = Vec.tabulate(65536)(i => ((i % 11) - 5).toDouble * 0.25)
    val a = Matrix.tabulate(512, 512)((i, j) => ((i * 3 + j * 5) % 19).toDouble * 0.03125)
    val gemvX = Vec.tabulate(512)(i => ((i % 7) - 3).toDouble * 0.5)
    val gemmA = Matrix.tabulate(128, 128)((i, j) => ((i * 3 + j * 5) % 19).toDouble * 0.03125)
    val gemmB = Matrix.tabulate(128, 128)((i, j) => ((i * 7 - j * 2) % 23).toDouble * 0.015625)
    val gemmDestination = DMatBuilder.zeros(128, 128)
    val combinationDestination = DMatBuilder.zeros(128, 128)
    val symmetric = Matrix.tabulate(64, 64): (i, j) =>
      val lo = math.min(i, j)
      val hi = math.max(i, j)
      if lo == hi then 66.0 else ((lo * 5 + hi * 3) % 17).toDouble / 34.0
    val qrInput = gemmA.slice(0, 64, 0, 64)
    val qrWorkspace = DenseWorkspace.forQR(64, 64)
    val symmetricWorkspace = DenseWorkspace.forRequirement(
      Eigen.symmetricScratchRequirement(64, EigenVectors.ValuesOnly).toOption.get
    )
    val sparseBuilderA = Sparse.coo(256, 256)
    val sparseBuilderB = Sparse.coo(256, 256)
    var sparseRow = 0
    while sparseRow < 256 do
      sparseBuilderA.add(sparseRow, sparseRow, 3.0)
      sparseBuilderA.add(sparseRow, (sparseRow + 1) % 256, -0.5)
      sparseBuilderB.add(sparseRow, sparseRow, -1.0)
      sparseBuilderB.add(sparseRow, (sparseRow + 255) % 256, 0.25)
      sparseRow += 1
    val sparseA: CSR = sparseBuilderA.toCSR()
    val sparseB: CSR = sparseBuilderB.toCSR()
    val sparseUnionPlan = CSRUnionPlan.analyze(sparseA.pattern, sparseB.pattern).toOption.get
    val sparseUnionDestination = sparseUnionPlan.newDestination()
    val sparseProductPlan = CSRProductPlan.analyze(sparseA.pattern, sparseB.pattern).toOption.get
    val sparseProductDestination = sparseProductPlan.newDestination()
    val canonicalRowPtr = IndexArray.alloc(257)
    val canonicalColIdx = IndexArray.alloc(1024)
    val canonicalValues = DoubleArray.alloc(1024)
    sparseRow = 0
    while sparseRow < 256 do
      val start = 4 * sparseRow
      canonicalRowPtr(sparseRow) = start
      canonicalColIdx(start) = (sparseRow + 1) % 256
      canonicalValues(start) = 0.75
      canonicalColIdx(start + 1) = sparseRow
      canonicalValues(start + 1) = 3.0
      canonicalColIdx(start + 2) = (sparseRow + 1) % 256
      canonicalValues(start + 2) = -0.25
      canonicalColIdx(start + 3) = (sparseRow + 255) % 256
      canonicalValues(start + 3) = -0.5
      sparseRow += 1
    canonicalRowPtr(256) = 1024
    val sparseCanonicalInput = new CSR(256, 256, canonicalRowPtr, canonicalColIdx, canonicalValues)
    val sparseCanonicalWorkspace = DenseWorkspace.forRequirement(
      sparseCanonicalInput.canonicalizeScratchRequirement.toOption.get
    )
    val sparseX = Vec.tabulate(256)(i => (i % 11).toDouble - 5.0)
    val sparseDestination = MutableDVec.zeros(256)

    profile("dot-65536", iterations = 500)((_, _) => x.dot(y))

    val gemvOut = MutableVec.zeros(512)
    profile("gemv-512-reuse", iterations = 200) { (counters, _) =>
      a.mulInto(gemvX, gemvOut)
      counters.destinationReuses += 1
      gemvOut(0)
    }

    profile("gemm-128-alloc", iterations = 40) { (counters, _) =>
      val result = gemmA.*(gemmB)(using PureBackend)
      counters.ownedResults += 1
      result(7, 11)
    }

    profile("gemm-128-reuse", iterations = 40) { (counters, _) =>
      gemmA.gemmInto(gemmB, gemmDestination)(using PureBackend)
      counters.destinationReuses += 1
      gemmDestination(7, 11)
    }

    profile("add-scale-128", iterations = 15) { (counters, _) =>
      val sum = gemmA + gemmB
      counters.ownedResults += 1
      val scaled = DMatBuilder.zeros(128, 128)
      var row = 0
      while row < 128 do
        var col = 0
        while col < 128 do
          scaled(row, col) = 0.5 * sum(row, col)
          col += 1
        row += 1
      val result = scaled.result()
      counters.ownedResults += 1
      result(7, 11)
    }

    profile("add-scale-128-reuse", iterations = 500) { (counters, _) =>
      gemmA.linearCombinationInto(gemmB, combinationDestination, alpha = 0.5, beta = 0.5)
      counters.destinationReuses += 1
      combinationDestination(7, 11)
    }

    profile("qr-64-fresh", iterations = 8) { (counters, _) =>
      val result = qrInput.qr(using PureBackend)
      counters.ownedResults += 1
      result.r(0, 0)
    }

    profile("qr-64-reuse", iterations = 8) { (counters, _) =>
      val result = qrInput.qrWith(qrWorkspace)(using PureBackend)
      counters.ownedResults += 1
      counters.workspaceReuses += 1
      result.r(0, 0)
    }

    profile("eigen-64-values", iterations = 5) { (counters, _) =>
      val result = Eigen.eigSymmetric(symmetric, EigenSelection.All, EigenVectors.ValuesOnly).toOption.get
      counters.ownedResults += 1
      result.eigenvalues(0)
    }

    profile("eigen-64-values-reuse", iterations = 5) { (counters, _) =>
      val result =
        Eigen
          .eigSymmetricWith(symmetric, EigenSelection.All, EigenVectors.ValuesOnly, symmetricWorkspace)
          .toOption
          .get
      counters.ownedResults += 1
      counters.workspaceReuses += 1
      result.eigenvalues(0)
    }

    profile("sparse-add-256", iterations = 80) { (counters, _) =>
      val result = sparseA + sparseB
      counters.ownedResults += 1
      result(0, 0)
    }

    profile("sparse-union-plan", iterations = 500) { (counters, _) =>
      val result = CSRUnionPlan.analyze(sparseA.pattern, sparseB.pattern).toOption.get
      counters.planAnalyses += 1
      result.resultPattern.nnz.toDouble
    }

    profile("sparse-union-replay", iterations = 5000) { (counters, _) =>
      sparseUnionPlan.evaluateInto(sparseA, sparseB, sparseUnionDestination).fold(throw _, identity)
      counters.destinationReuses += 1
      counters.planReuses += 1
      sparseUnionDestination(0, 0)
    }

    profile("sparse-product-plan", iterations = 200) { (counters, _) =>
      val result = CSRProductPlan.analyze(sparseA.pattern, sparseB.pattern).toOption.get
      counters.planAnalyses += 1
      result.contributionCount.toDouble
    }

    profile("sparse-product-replay", iterations = 2000) { (counters, _) =>
      sparseProductPlan.evaluateInto(sparseA, sparseB, sparseProductDestination).fold(throw _, identity)
      counters.destinationReuses += 1
      counters.planReuses += 1
      sparseProductDestination(0, 0)
    }

    profile("sparse-map-256", iterations = 5000) { (counters, _) =>
      val result = sparseA.mapValues(_ * 1.0001)
      counters.ownedResults += 1
      counters.patternReuses += 1
      result(0, 0)
    }

    profile("sparse-canon-256-fresh", iterations = 1000) { (counters, _) =>
      val result = sparseCanonicalInput.canonicalize
      counters.ownedResults += 1
      result(0, 0)
    }

    profile("sparse-canon-256-reuse", iterations = 1000) { (counters, _) =>
      val result = sparseCanonicalInput.canonicalizeWith(sparseCanonicalWorkspace)
      counters.ownedResults += 1
      counters.workspaceReuses += 1
      result(0, 0)
    }

    profile("spmv-256-alloc", iterations = 5000) { (counters, _) =>
      val result = sparseA.*(sparseX)(using PureBackend)
      counters.ownedResults += 1
      result(0)
    }

    profile("spmv-256-reuse", iterations = 5000) { (counters, _) =>
      sparseA.mulInto(sparseX, sparseDestination)
      counters.destinationReuses += 1
      sparseDestination(0)
    }

  /** A small wall-clock profile rather than a statistical benchmark framework:
    * five complete warmup batches, then the median of nine batches. Both JS and
    * Wasm run this exact linked program with the checksum consumed and printed.
    */
  private def profile(name: String, iterations: Int)(body: (ScenarioCounters, Int) => Double): Unit =
    val counters = new ScenarioCounters
    var warmup = 0
    var checksum = 0.0
    while warmup < 5 do
      var i = 0
      while i < iterations do
        checksum += body(counters, i)
        i += 1
      warmup += 1

    val samples = new Array[Double](9)
    var sample = 0
    while sample < samples.length do
      val start = js.Date.now()
      var i = 0
      while i < iterations do
        checksum += body(counters, i)
        i += 1
      samples(sample) = (js.Date.now() - start) * 1000000.0 / iterations.toDouble
      sample += 1
    scala.util.Sorting.quickSort(samples)
    val medianNs = samples(samples.length / 2)
    val operations = 14L * iterations.toLong
    println(
      f"$name%-20s median=$medianNs%.1f ns/op samples=9 " +
        f"owned=${counters.ownedResults.toDouble / operations.toDouble}%.1f/op " +
        f"dest-reuse=${counters.destinationReuses.toDouble / operations.toDouble}%.1f/op " +
        f"work-reuse=${counters.workspaceReuses.toDouble / operations.toDouble}%.1f/op " +
        f"pattern-reuse=${counters.patternReuses.toDouble / operations.toDouble}%.1f/op " +
        f"plan-analysis=${counters.planAnalyses.toDouble / operations.toDouble}%.1f/op " +
        f"plan-reuse=${counters.planReuses.toDouble / operations.toDouble}%.1f/op " +
        f"checksum=$checksum%.6f"
    )

  /** Contract-level construction counters for Scala.js, where JMH's JVM GC
    * profiler is unavailable. They count public owned results and explicit
    * destination/workspace reuse in the executed scenario; they are not a claim
    * about all internal JavaScript engine allocations.
    */
  private final class ScenarioCounters:
    var ownedResults = 0L
    var destinationReuses = 0L
    var workspaceReuses = 0L
    var patternReuses = 0L
    var planAnalyses = 0L
    var planReuses = 0L
