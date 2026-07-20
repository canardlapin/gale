package gale.sparse.direct

import gale.backend.BackendConfig
import gale.linalg.*
import gale.sparse.CSR
import gale.sparse.CSRPattern
import gale.sparse.Permutation
import scala.util.control.NonFatal

/** Sparse factorization families a JVM provider may explicitly implement. */
enum SparseDirectFactorization:
  case LU, Cholesky, QR

  def requiredCapability: SparseDirectCapability =
    this match
      case LU       => SparseDirectCapability.LU
      case Cholesky => SparseDirectCapability.Cholesky
      case QR       => SparseDirectCapability.QR

/** Discoverable operations and optional solve/ordering features. An empty set
  * means that no sparse-direct implementation is available.
  */
enum SparseDirectCapability:
  case LU, Cholesky, QR
  case UserOrdering, TransposeSolve, MultipleRhs

/** Column-ordering request for symbolic analysis.
  *
  * A user permutation maps new column position to original column through
  * `columnPermutation.toIndexSeq`. `Natural` is identity; `ProviderDefault`
  * lets the provider choose and report an owned permutation.
  */
enum SparseDirectOrdering:
  case Natural, ProviderDefault
  case User(columnPermutation: Permutation)

enum SparseSolveOperation:
  case Normal, Transpose

/** Immutable symbolic-stage diagnostics. No wall-clock time is stored, keeping
  * deterministic providers deterministic.
  */
final case class SparseSymbolicDiagnostics(
    providerName: String,
    factorization: SparseDirectFactorization,
    inputNnz: Int,
    predictedFactorNnz: Option[Long],
    ordering: SparseDirectOrdering,
    deterministic: Boolean
)

final case class SparseNumericDiagnostics(
    providerName: String,
    factorization: SparseDirectFactorization,
    factorNnz: Option[Long],
    rank: Option[Int],
    pivotCount: Option[Int],
    reciprocalConditionEstimate: Option[Double]
)

final case class SparseSolveDiagnostics(
    providerName: String,
    operation: SparseSolveOperation,
    rightHandSides: Int,
    residualNorm: Option[Double],
    refinementSteps: Int
)

final case class SparseVectorSolve(solution: DVec, diagnostics: SparseSolveDiagnostics)
final case class SparseMatrixSolve(solution: DMat, diagnostics: SparseSolveDiagnostics)

/** Explicit lifecycle shared by provider workspaces and native-capable handles.
  * `close()` must be idempotent. A closed resource is rejected by the Gale
  * facade before provider code is entered.
  */
trait SparseDirectResource extends AutoCloseable:
  def isClosed: Boolean
  override def close(): Unit

/** Mutable provider-owned scratch. A workspace belongs to exactly one provider,
  * is single-user (not safe for concurrent calls), and may be reused
  * sequentially for analysis, factorization, and solves.
  */
trait SparseDirectWorkspace extends SparseDirectResource:
  def provider: SparseDirectProvider

/** Reusable symbolic analysis tied to the exact immutable CSR storage analyzed.
  * Implementations may wrap native resources. The analysis itself may be used
  * concurrently only with distinct workspaces.
  */
trait SparseDirectSymbolicAnalysis extends SparseDirectResource:
  def provider: SparseDirectProvider
  def pattern: CSRPattern
  def factorization: SparseDirectFactorization
  def ordering: SparseDirectOrdering
  def columnPermutation: Permutation
  def diagnostics: SparseSymbolicDiagnostics

  /** Low-level provider hook. Call [[SparseDirect.factor]] for Gale's pattern,
    * workspace, lifecycle, and returned-handle validation.
    */
  def factorNumeric(
      matrix: CSR,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseDirectNumericFactor]

/** Numeric factor handle. Permutations are Gale-owned immutable values: a
  * provider must copy native index buffers before returning them.
  *
  * A factor may be shared across concurrent solves when each invocation uses a
  * distinct workspace. The same workspace must never be used concurrently.
  */
trait SparseDirectNumericFactor extends SparseDirectResource:
  def provider: SparseDirectProvider
  def factorization: SparseDirectFactorization
  def inputRows: Int
  def inputCols: Int
  def rowPermutation: Permutation
  def columnPermutation: Permutation
  def diagnostics: SparseNumericDiagnostics

  /** Provider-declared solve dimensions. For normal square LU/Cholesky both are
    * the matrix order; a rectangular QR provider may report rows -> columns.
    */
  def rhsRows(operation: SparseSolveOperation): Int
  def solutionRows(operation: SparseSolveOperation): Int

  /** Low-level provider hooks. Use the [[SparseDirect.solve]] / `solveInto`
    * facades for capability, dimension, workspace, and lifecycle validation.
    */
  def solveVectorInto(
      rhs: DVec,
      destination: MutableDVec,
      operation: SparseSolveOperation,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseSolveDiagnostics]

  def solveMatrixInto(
      rhs: DMat,
      destination: DMatBuilder,
      operation: SparseSolveOperation,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseSolveDiagnostics]

/** JVM-only provider contract. Provider objects are resolved explicitly through
  * a `given`, are safe for concurrent invocation, and fix their thread policy in
  * [[config]] at construction time. Stateful native scratch belongs in distinct
  * [[SparseDirectWorkspace]] instances, never in the shared provider singleton.
  */
trait SparseDirectProvider:
  def name: String
  def capabilities: Set[SparseDirectCapability]
  def config: BackendConfig

  def createWorkspace(): Either[LinAlgError, SparseDirectWorkspace]

  /** Low-level symbolic hook. Use [[SparseDirect.analyze]] for Gale validation. */
  def analyze(
      pattern: CSRPattern,
      factorization: SparseDirectFactorization,
      ordering: SparseDirectOrdering,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseDirectSymbolicAnalysis]

final case class SparseDirectProviderReport(
    name: String,
    capabilities: Set[SparseDirectCapability],
    config: BackendConfig,
    available: Boolean
)

object SparseDirectProvider:
  private object NoSparseDirectProvider extends SparseDirectProvider:
    val name: String = "none"
    val capabilities: Set[SparseDirectCapability] = Set.empty
    val config: BackendConfig = BackendConfig.singleThreaded

    def createWorkspace(): Either[LinAlgError, SparseDirectWorkspace] =
      Left(LinAlgError.UnsupportedOperation("sparse direct factorization: no JVM provider is installed"))

    def analyze(
        pattern: CSRPattern,
        factorization: SparseDirectFactorization,
        ordering: SparseDirectOrdering,
        workspace: SparseDirectWorkspace
    ): Either[LinAlgError, SparseDirectSymbolicAnalysis] =
      Left(LinAlgError.UnsupportedOperation("sparse direct factorization: no JVM provider is installed"))

  /** Capability-less default. The existence of this seam does not imply that LU,
    * Cholesky, or QR is implemented.
    */
  val none: SparseDirectProvider = NoSparseDirectProvider
  given default: SparseDirectProvider = none

  def current(using provider: SparseDirectProvider): SparseDirectProviderReport =
    SparseDirectProviderReport(provider.name, provider.capabilities, provider.config, provider.capabilities.nonEmpty)

  def validationErrors(provider: SparseDirectProvider): List[String] =
    val errors = List.newBuilder[String]
    if provider.name.trim.isEmpty then errors += "provider name must be non-empty"
    val factorCapabilities = Set(
      SparseDirectCapability.LU,
      SparseDirectCapability.Cholesky,
      SparseDirectCapability.QR
    )
    if provider.capabilities.intersect(factorCapabilities).isEmpty && provider.capabilities.nonEmpty then
      errors += "feature capabilities require at least one LU, Cholesky, or QR capability"
    errors.result()

  def requireValid(provider: SparseDirectProvider): provider.type =
    val errors = validationErrors(provider)
    require(errors.isEmpty, errors.mkString(s"invalid sparse-direct provider '${provider.name}': ", "; ", ""))
    provider

/** Validated sparse-direct facade. Every method is JVM-only and explicitly
  * provider/workspace based; nothing routes through Gale's portable sparse APIs.
  */
object SparseDirect:
  def capabilities(using provider: SparseDirectProvider): Set[SparseDirectCapability] =
    provider.capabilities

  def newWorkspace()(using provider: SparseDirectProvider): Either[LinAlgError, SparseDirectWorkspace] =
    val validation = SparseDirectProvider.validationErrors(provider)
    if validation.nonEmpty then
      Left(LinAlgError.InvalidArgument(validation.mkString("invalid sparse-direct provider: ", "; ", "")))
    else
      provider.createWorkspace().flatMap: workspace =>
        acceptReturnedResource(workspace):
          if !(workspace.provider eq provider) then
            Left(LinAlgError.InvalidArgument("sparse-direct workspace belongs to a different provider"))
          else if workspace.isClosed then
            Left(LinAlgError.InvalidArgument("sparse-direct provider returned a closed workspace"))
          else Right(())

  def analyze(
      pattern: CSRPattern,
      factorization: SparseDirectFactorization,
      workspace: SparseDirectWorkspace,
      ordering: SparseDirectOrdering = SparseDirectOrdering.ProviderDefault
  )(using provider: SparseDirectProvider): Either[LinAlgError, SparseDirectSymbolicAnalysis] =
    validateAnalysisRequest(pattern, factorization, ordering, workspace, provider).flatMap: _ =>
      provider.analyze(pattern, factorization, ordering, workspace).flatMap: analysis =>
        acceptReturnedResource(analysis):
          validateAnalysisResult(analysis, pattern, factorization, provider)

  def factor(
      analysis: SparseDirectSymbolicAnalysis,
      matrix: CSR,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseDirectNumericFactor] =
    if analysis.isClosed then closed("symbolic analysis")
    else if !matrix.sharesPatternStorage(analysis.pattern) then
      Left(LinAlgError.InvalidArgument("numeric CSR pattern differs from the exact pattern used for symbolic analysis"))
    else
      validateWorkspace(workspace, analysis.provider).flatMap: _ =>
        analysis.factorNumeric(matrix, workspace).flatMap: factor =>
          acceptReturnedResource(factor):
            validateNumericResult(factor, analysis)

  def solve(
      factor: SparseDirectNumericFactor,
      rhs: DVec,
      workspace: SparseDirectWorkspace,
      operation: SparseSolveOperation = SparseSolveOperation.Normal
  ): Either[LinAlgError, SparseVectorSolve] =
    validateVectorSolve(factor, rhs, workspace, operation).flatMap: _ =>
      val destination = MutableDVec.zeros(factor.solutionRows(operation))
      factor.solveVectorInto(rhs, destination, operation, workspace).map: diagnostics =>
        SparseVectorSolve(destination.asVec, diagnostics)

  def solveInto(
      factor: SparseDirectNumericFactor,
      rhs: DVec,
      destination: MutableDVec,
      workspace: SparseDirectWorkspace,
      operation: SparseSolveOperation
  ): Either[LinAlgError, SparseSolveDiagnostics] =
    validateVectorSolve(factor, rhs, workspace, operation).flatMap: _ =>
      val expected = factor.solutionRows(operation)
      if destination.length != expected then Left(LinAlgError.VectorLengthMismatch(expected, destination.length))
      else factor.solveVectorInto(rhs, destination, operation, workspace)

  def solveInto(
      factor: SparseDirectNumericFactor,
      rhs: DVec,
      destination: MutableDVec,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseSolveDiagnostics] =
    solveInto(factor, rhs, destination, workspace, SparseSolveOperation.Normal)

  def solve(
      factor: SparseDirectNumericFactor,
      rhs: DMat,
      workspace: SparseDirectWorkspace,
      operation: SparseSolveOperation
  ): Either[LinAlgError, SparseMatrixSolve] =
    validateMatrixSolve(factor, rhs, workspace, operation).flatMap: _ =>
      val destination = DMatBuilder.zeros(factor.solutionRows(operation), rhs.cols)
      factor.solveMatrixInto(rhs, destination, operation, workspace).map: diagnostics =>
        SparseMatrixSolve(destination.result(), diagnostics)

  def solve(
      factor: SparseDirectNumericFactor,
      rhs: DMat,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseMatrixSolve] =
    solve(factor, rhs, workspace, SparseSolveOperation.Normal)

  def solveInto(
      factor: SparseDirectNumericFactor,
      rhs: DMat,
      destination: DMatBuilder,
      workspace: SparseDirectWorkspace,
      operation: SparseSolveOperation
  ): Either[LinAlgError, SparseSolveDiagnostics] =
    validateMatrixSolve(factor, rhs, workspace, operation).flatMap: _ =>
      val expectedRows = factor.solutionRows(operation)
      if destination.rows != expectedRows || destination.cols != rhs.cols then
        Left(
          LinAlgError.InvalidArgument(
            s"sparse solve destination must be ${expectedRows}x${rhs.cols}, got ${destination.rows}x${destination.cols}"
          )
        )
      else factor.solveMatrixInto(rhs, destination, operation, workspace)

  def solveInto(
      factor: SparseDirectNumericFactor,
      rhs: DMat,
      destination: DMatBuilder,
      workspace: SparseDirectWorkspace
  ): Either[LinAlgError, SparseSolveDiagnostics] =
    solveInto(factor, rhs, destination, workspace, SparseSolveOperation.Normal)

  private def validateAnalysisRequest(
      pattern: CSRPattern,
      factorization: SparseDirectFactorization,
      ordering: SparseDirectOrdering,
      workspace: SparseDirectWorkspace,
      provider: SparseDirectProvider
  ): Either[LinAlgError, Unit] =
    val providerErrors = SparseDirectProvider.validationErrors(provider)
    if providerErrors.nonEmpty then
      Left(LinAlgError.InvalidArgument(providerErrors.mkString("invalid sparse-direct provider: ", "; ", "")))
    else if !provider.capabilities.contains(factorization.requiredCapability) then
      Left(
        LinAlgError.UnsupportedOperation(
          s"sparse direct $factorization is not provided by '${provider.name}'"
        )
      )
    else if !pattern.hasCanonicalFormat then
      Left(LinAlgError.InvalidArgument("sparse-direct symbolic analysis requires a canonical CSR pattern"))
    else if factorization != SparseDirectFactorization.QR && pattern.rows != pattern.cols then
      Left(LinAlgError.NonSquareMatrix(Shape(Rows(pattern.rows), Cols(pattern.cols))))
    else
      ordering match
        case SparseDirectOrdering.User(permutation) =>
          if !provider.capabilities.contains(SparseDirectCapability.UserOrdering) then
            Left(LinAlgError.UnsupportedOperation(s"user sparse ordering is not provided by '${provider.name}'"))
          else if permutation.rows != pattern.cols || permutation.cols != pattern.cols then
            Left(
              LinAlgError.InvalidArgument(
                s"column ordering must have size ${pattern.cols}, got ${permutation.rows}x${permutation.cols}"
              )
            )
          else validateWorkspace(workspace, provider)
        case _ => validateWorkspace(workspace, provider)

  private def validateAnalysisResult(
      analysis: SparseDirectSymbolicAnalysis,
      pattern: CSRPattern,
      factorization: SparseDirectFactorization,
      provider: SparseDirectProvider
  ): Either[LinAlgError, Unit] =
    if analysis.isClosed then Left(LinAlgError.InvalidArgument("provider returned a closed symbolic analysis"))
    else if !(analysis.provider eq provider) then
      Left(LinAlgError.InvalidArgument("symbolic analysis belongs to a different provider"))
    else if analysis.factorization != factorization then
      Left(LinAlgError.InvalidArgument("symbolic analysis reports a different factorization family"))
    else if !analysis.pattern.sharesStorageWith(pattern) then
      Left(LinAlgError.InvalidArgument("symbolic analysis reports a different input pattern"))
    else if analysis.columnPermutation.rows != pattern.cols || analysis.columnPermutation.cols != pattern.cols then
      Left(LinAlgError.InvalidArgument("symbolic analysis returned a malformed column permutation"))
    else if analysis.diagnostics.providerName != provider.name || analysis.diagnostics.factorization != factorization then
      Left(LinAlgError.InvalidArgument("symbolic analysis returned inconsistent diagnostics"))
    else Right(())

  private def validateNumericResult(
      factor: SparseDirectNumericFactor,
      analysis: SparseDirectSymbolicAnalysis
  ): Either[LinAlgError, Unit] =
    if factor.isClosed then Left(LinAlgError.InvalidArgument("provider returned a closed numeric factor"))
    else if !(factor.provider eq analysis.provider) then
      Left(LinAlgError.InvalidArgument("numeric factor belongs to a different provider"))
    else if factor.factorization != analysis.factorization then
      Left(LinAlgError.InvalidArgument("numeric factor reports a different factorization family"))
    else if factor.inputRows != analysis.pattern.rows || factor.inputCols != analysis.pattern.cols then
      Left(LinAlgError.InvalidArgument("numeric factor reports dimensions different from symbolic analysis"))
    else if factor.rowPermutation.rows != factor.inputRows || factor.rowPermutation.cols != factor.inputRows then
      Left(LinAlgError.InvalidArgument("numeric factor returned a malformed row permutation"))
    else if factor.columnPermutation.rows != factor.inputCols || factor.columnPermutation.cols != factor.inputCols then
      Left(LinAlgError.InvalidArgument("numeric factor returned a malformed column permutation"))
    else if
      factor.diagnostics.providerName != factor.provider.name ||
        factor.diagnostics.factorization != factor.factorization
    then Left(LinAlgError.InvalidArgument("numeric factor returned inconsistent diagnostics"))
    else Right(())

  private def validateVectorSolve(
      factor: SparseDirectNumericFactor,
      rhs: DVec,
      workspace: SparseDirectWorkspace,
      operation: SparseSolveOperation
  ): Either[LinAlgError, Unit] =
    validateSolve(factor, workspace, operation).flatMap: _ =>
      val expected = factor.rhsRows(operation)
      if expected < 0 || factor.solutionRows(operation) < 0 then
        Left(LinAlgError.InvalidArgument("provider reported negative sparse solve dimensions"))
      else if rhs.length != expected then Left(LinAlgError.VectorLengthMismatch(expected, rhs.length))
      else Right(())

  private def validateMatrixSolve(
      factor: SparseDirectNumericFactor,
      rhs: DMat,
      workspace: SparseDirectWorkspace,
      operation: SparseSolveOperation
  ): Either[LinAlgError, Unit] =
    validateSolve(factor, workspace, operation).flatMap: _ =>
      val expected = factor.rhsRows(operation)
      if expected < 0 || factor.solutionRows(operation) < 0 then
        Left(LinAlgError.InvalidArgument("provider reported negative sparse solve dimensions"))
      else if rhs.cols > 1 && !factor.provider.capabilities.contains(SparseDirectCapability.MultipleRhs) then
        Left(
          LinAlgError.UnsupportedOperation(
            s"multiple-right-hand-side sparse solve is not provided by '${factor.provider.name}'"
          )
        )
      else if rhs.rows != expected then
        Left(
          LinAlgError.InvalidArgument(
            s"sparse solve right-hand side must have $expected rows, got ${rhs.rows}"
          )
        )
      else Right(())

  private def validateSolve(
      factor: SparseDirectNumericFactor,
      workspace: SparseDirectWorkspace,
      operation: SparseSolveOperation
  ): Either[LinAlgError, Unit] =
    if factor.isClosed then closed("numeric factor")
    else if operation == SparseSolveOperation.Transpose &&
      !factor.provider.capabilities.contains(SparseDirectCapability.TransposeSolve)
    then
      Left(
        LinAlgError.UnsupportedOperation(
          s"transpose sparse solve is not provided by '${factor.provider.name}'"
        )
      )
    else validateWorkspace(workspace, factor.provider)

  private def validateWorkspace(
      workspace: SparseDirectWorkspace,
      provider: SparseDirectProvider
  ): Either[LinAlgError, Unit] =
    if !(workspace.provider eq provider) then
      Left(LinAlgError.InvalidArgument("sparse-direct workspace belongs to a different provider"))
    else if workspace.isClosed then closed("workspace")
    else Right(())

  /** A provider-created resource is provisionally owned by the facade until its
    * returned-handle contract passes. Rejected resources must not escape or leak.
    */
  private def acceptReturnedResource[R <: SparseDirectResource](
      resource: R
  )(validation: => Either[LinAlgError, Unit]): Either[LinAlgError, R] =
    validation match
      case Right(()) => Right(resource)
      case Left(error) =>
        try resource.close()
        catch case NonFatal(closeFailure) => error.addSuppressed(closeFailure)
        Left(error)

  private def closed(resource: String): Either[LinAlgError, Nothing] =
    Left(LinAlgError.UnsupportedOperation(s"closed sparse-direct $resource"))
