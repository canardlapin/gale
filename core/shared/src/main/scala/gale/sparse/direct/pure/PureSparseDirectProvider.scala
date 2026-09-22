package gale.sparse.direct.pure

import gale.backend.BackendConfig
import gale.linalg.*
import gale.sparse.CSR
import gale.sparse.CSRPattern
import gale.sparse.Permutation
import gale.sparse.Sparse
import gale.sparse.direct.*

/** Portable sparse Cholesky provider. Advertise only Cholesky and the optional
  * solve/ordering features this implementation actually has. Import
  * `gale.sparse.direct.pure.given` or pass [[PureSparseDirectProvider]]
  * explicitly; the default `given` remains [[SparseDirectProvider.none]].
  */
object PureSparseDirectProvider extends SparseDirectProvider:
  val name: String = "pure"
  val capabilities: Set[SparseDirectCapability] = Set(
    SparseDirectCapability.Cholesky,
    SparseDirectCapability.UserOrdering,
    SparseDirectCapability.TransposeSolve,
    SparseDirectCapability.MultipleRhs
  )
  val config: BackendConfig = BackendConfig.singleThreaded

  def createWorkspace(): Either[LinAlgError, SparseDirectWorkspace] =
    Right(new PureSparseDirectWorkspace(this))

  def analyze(
      pattern: CSRPattern,
      factorization: SparseDirectFactorization,
      ordering: SparseDirectOrdering,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseDirectSymbolicAnalysis] =
    if factorization != SparseDirectFactorization.Cholesky then
      Left(LinAlgError.UnsupportedOperation(s"sparse direct $factorization is not provided by '$name'"))
    else
      requireWorkspace(workspace).flatMap: _ =>
        val n = pattern.rows
        val adj = SparseCholeskyKernels.undirectedLowerAdj(pattern)
        val perm = ordering match
          case SparseDirectOrdering.Natural         => SparseCholeskyKernels.identityPerm(n)
          case SparseDirectOrdering.ProviderDefault => SparseCholeskyKernels.minimumDegree(adj)
          case SparseDirectOrdering.User(user)      => user.toArray
        val inv = SparseCholeskyKernels.invertPerm(perm)
        val adjPerm = SparseCholeskyKernels.permuteAdj(adj, perm, inv)
        val parent = SparseCholeskyKernels.eliminationTree(n, adjPerm)
        SparseCholeskyKernels.symbolicColumns(n, adjPerm, parent, pattern.nnz).map: columns =>
          val (colPtr, rowIdx) = SparseCholeskyKernels.packCsc(columns)
          val (rowPtr, hitCol, hitOff) = SparseCholeskyKernels.rowHits(n, colPtr, rowIdx)
          val predicted = colPtr(n).toLong
          new PureSparseCholeskyAnalysis(
            this,
            pattern,
            ordering,
            copyPermutation(perm),
            colPtr,
            rowIdx,
            rowPtr,
            hitCol,
            hitOff,
            predicted
          )

  private[pure] def requireWorkspace(workspace: SparseDirectWorkspace): Either[LinAlgError, PureSparseDirectWorkspace] =
    workspace match
      case w: PureSparseDirectWorkspace if w.provider eq this =>
        if w.isClosed then Left(LinAlgError.UnsupportedOperation("closed sparse-direct workspace"))
        else Right(w)
      case _ =>
        Left(LinAlgError.InvalidArgument("sparse-direct workspace belongs to a different provider"))

  private[pure] def copyPermutation(perm: Array[Int]): Permutation =
    Sparse.permutation(perm.toIndexedSeq*)

given pureSparseDirect: SparseDirectProvider = PureSparseDirectProvider

private final class PureSparseDirectWorkspace(val provider: SparseDirectProvider) extends SparseDirectWorkspace:
  private var closed = false
  private var work: Array[Double] = Array.empty

  def isClosed: Boolean = closed

  def close(): Unit =
    if !closed then
      closed = true
      work = Array.empty

  def scratch(n: Int): Array[Double] =
    if work.length < n then work = new Array[Double](n)
    work

private final class PureSparseCholeskyAnalysis(
    val provider: PureSparseDirectProvider.type,
    val pattern: CSRPattern,
    val ordering: SparseDirectOrdering,
    val columnPermutation: Permutation,
    colPtr: Array[Int],
    rowIdx: Array[Int],
    rowPtr: Array[Int],
    hitCol: Array[Int],
    hitOff: Array[Int],
    predictedFactorNnz: Long
) extends SparseDirectSymbolicAnalysis:
  private var closed = false

  val factorization: SparseDirectFactorization = SparseDirectFactorization.Cholesky
  val diagnostics: SparseSymbolicDiagnostics =
    SparseSymbolicDiagnostics(
      provider.name,
      factorization,
      pattern.nnz,
      predictedFactorNnz = Some(predictedFactorNnz),
      ordering,
      deterministic = true
    )

  def isClosed: Boolean = closed

  def close(): Unit =
    closed = true

  def factorNumeric(
      matrix: CSR,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseDirectNumericFactor] =
    if closed then Left(LinAlgError.UnsupportedOperation("closed sparse-direct symbolic analysis"))
    else
      provider.requireWorkspace(workspace).flatMap: ws =>
        val n = pattern.rows
        val perm = columnPermutation.toArray
        val work = ws.scratch(n)
        SparseCholeskyKernels
          .factorValues(matrix, perm, colPtr, rowIdx, rowPtr, hitCol, hitOff, work)
          .map: values =>
            new PureSparseCholeskyFactor(
              provider,
              n,
              PureSparseDirectProvider.copyPermutation(perm),
              colPtr.clone(),
              rowIdx.clone(),
              values,
              predictedFactorNnz
            )

private final class PureSparseCholeskyFactor(
    val provider: PureSparseDirectProvider.type,
    n: Int,
    val columnPermutation: Permutation,
    colPtr: Array[Int],
    rowIdx: Array[Int],
    values: Array[Double],
    factorNnz: Long
) extends SparseDirectNumericFactor:
  private var closed = false

  val factorization: SparseDirectFactorization = SparseDirectFactorization.Cholesky
  val inputRows: Int = n
  val inputCols: Int = n
  val rowPermutation: Permutation = PureSparseDirectProvider.copyPermutation(SparseCholeskyKernels.identityPerm(n))
  val diagnostics: SparseNumericDiagnostics =
    SparseNumericDiagnostics(
      provider.name,
      factorization,
      factorNnz = Some(factorNnz),
      rank = Some(n),
      pivotCount = Some(0),
      reciprocalConditionEstimate = None
    )

  def rhsRows(operation: SparseSolveOperation): Int = n
  def solutionRows(operation: SparseSolveOperation): Int = n

  def isClosed: Boolean = closed

  def close(): Unit =
    closed = true

  def solveVectorInto(
      rhs: DVec,
      destination: MutableDVec,
      operation: SparseSolveOperation,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseSolveDiagnostics] =
    if closed then Left(LinAlgError.UnsupportedOperation("closed sparse-direct numeric factor"))
    else
      provider.requireWorkspace(workspace).map: ws =>
        val perm = columnPermutation.toArray
        val work = ws.scratch(n)
        var i = 0
        while i < n do
          work(i) = rhs(perm(i))
          i += 1
        SparseCholeskyKernels.solveInPlace(colPtr, rowIdx, values, work)
        i = 0
        while i < n do
          destination(perm(i)) = work(i)
          i += 1
        SparseSolveDiagnostics(provider.name, operation, 1, residualNorm = None, refinementSteps = 0)

  def solveMatrixInto(
      rhs: DMat,
      destination: DMatBuilder,
      operation: SparseSolveOperation,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseSolveDiagnostics] =
    if closed then Left(LinAlgError.UnsupportedOperation("closed sparse-direct numeric factor"))
    else
      provider.requireWorkspace(workspace).map: ws =>
        val perm = columnPermutation.toArray
        val work = ws.scratch(n)
        var col = 0
        while col < rhs.cols do
          var i = 0
          while i < n do
            work(i) = rhs(perm(i), col)
            i += 1
          SparseCholeskyKernels.solveInPlace(colPtr, rowIdx, values, work)
          i = 0
          while i < n do
            destination(perm(i), col) = work(i)
            i += 1
          col += 1
        SparseSolveDiagnostics(provider.name, operation, rhs.cols, residualNorm = None, refinementSteps = 0)
