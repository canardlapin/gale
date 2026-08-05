package gale.spectral

import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.PositiveDefinite
import gale.linalg.Rows
import gale.linalg.Shape
import scala.collection.mutable.ArrayBuffer

/** Portable generalized block-Lanczos engine for `A x = lambda B x`.
  *
  * The Krylov action is `T = B^-1 A`, but the inverse is never formed: `A` is applied and each resulting column is
  * passed to the caller's explicit [[MetricSolveOperator]]. Since `T` is self-adjoint in the B inner product, the basis
  * is maintained with twice-reorthogonalized B geometry and the dense Rayleigh-Ritz projection remains symmetric.
  */
private[spectral] object GeneralizedLanczos:

  private final case class RitzState(
      block: GeneralizedBlockKernels.MetricBlock,
      operatorImages: DMat,
      values: DVec,
      residualNorms: DVec,
      converged: Array[Boolean]
  )

  private final case class ExpandedBasis(
      block: GeneralizedBlockKernels.MetricBlock,
      operatorImages: DMat,
      innerWork: LinearSolveSummary,
      fullSpace: Boolean
  )

  def solve[B <: DoubleLinearOperator](
      operator: DoubleLinearOperator,
      metric: PositiveDefinite[B],
      metricSolve: MetricSolveOperator[B],
      n: Int,
      k: Int,
      order: EigenOrder,
      options: GeneralizedLanczosOptions
  ): Either[LinAlgError, EigenDecomposition] =
    validate(operator, metric, metricSolve, n, k, order, options).flatMap: _ =>
      val initial =
        options.initialSubspace
          .map(matrix => DMat.tabulate(matrix.rows, matrix.cols)((row, col) => matrix(row, col)))
          .getOrElse(deterministicInitial(n, k))
      val orthogonalityTolerance =
        math.max(1e-12, math.min(1e-8, options.tolerance * 1e-2))
      val targetDimension =
        options.subspaceDimension.getOrElse(math.min(n, math.max(4 * k, 20)))

      for
        initialBlock <- GeneralizedBlockKernels.bOrthonormalizeAndReplenish(
          initial,
          metric,
          k,
          orthogonalityTolerance
        )
        initialImages <- GeneralizedBlockKernels.applyBlock(operator, initialBlock.vectors)
        initialState <- rayleighRitz(
          initialBlock,
          initialImages,
          k,
          order,
          options.tolerance
        )
        result <- iterate(
          operator,
          metric,
          metricSolve,
          n,
          k,
          order,
          options,
          targetDimension,
          orthogonalityTolerance,
          initialState
        )
      yield result

  private def iterate[B <: DoubleLinearOperator](
      operator: DoubleLinearOperator,
      metric: PositiveDefinite[B],
      metricSolve: MetricSolveOperator[B],
      n: Int,
      k: Int,
      order: EigenOrder,
      options: GeneralizedLanczosOptions,
      targetDimension: Int,
      orthogonalityTolerance: Double,
      initialState: RitzState
  ): Either[LinAlgError, EigenDecomposition] =
    var state = initialState
    var innerWork = LinearSolveSummary.Empty
    var iterations = 0
    var extremalityCertified = false
    var invariantStartExplorationPending =
      convergedCount(initialState) == k && k < n
    var failure: Option[LinAlgError] = None

    while iterations < options.maxIterations &&
      (convergedCount(state) < k || invariantStartExplorationPending) &&
      failure.isEmpty
    do
      expandBasis(
        operator,
        metric,
        metricSolve,
        state,
        targetDimension,
        k,
        orthogonalityTolerance,
        iterations + 1,
        innerWork
      ) match
        case Left(error) =>
          failure = Some(error)
        case Right(expanded) =>
          innerWork = expanded.innerWork
          rayleighRitz(
            expanded.block,
            expanded.operatorImages,
            k,
            order,
            options.tolerance
          ) match
            case Left(error) =>
              failure = Some(error)
            case Right(next) =>
              state = next
              extremalityCertified = expanded.fullSpace
              invariantStartExplorationPending = false
      iterations += 1

    failure match
      case Some(error) => Left(error)
      case None        =>
        Right(
          assemble(
            state,
            n,
            k,
            options.returnVectors == EigenVectors.Right,
            iterations,
            extremalityCertified,
            innerWork
          )
        )

  /** Thick restart begins with all wanted Ritz vectors. Converged vectors remain in the basis but are soft-locked: they
    * are not used as Krylov frontiers.
    */
  private def expandBasis[B <: DoubleLinearOperator](
      operator: DoubleLinearOperator,
      metric: PositiveDefinite[B],
      metricSolve: MetricSolveOperator[B],
      state: RitzState,
      targetDimension: Int,
      blockWidth: Int,
      tolerance: Double,
      outerIteration: Int,
      initialWork: LinearSolveSummary
  ): Either[LinAlgError, ExpandedBasis] =
    var basis = state.block
    var operatorImages = state.operatorImages
    var innerWork = initialWork
    val frontier = ArrayBuffer.empty[Int]
    var index = 0
    while index < state.block.cols do
      if !state.converged(index) then frontier += index
      index += 1
    var frontierOffset = 0
    var streamOffset = outerIteration * math.max(blockWidth, 1)

    def appendBlock(
        directions: GeneralizedBlockKernels.MetricBlock
    ): Either[LinAlgError, Unit] =
      if directions.cols == 0 then Right(())
      else
        GeneralizedBlockKernels
          .applyBlock(operator, directions.vectors)
          .flatMap: images =>
            val oldColumns = basis.cols
            GeneralizedBlockKernels
              .concatenate(basis, directions)
              .map: combined =>
                basis = combined
                operatorImages = concatenateMatrices(operatorImages, images)
                var column = 0
                while column < directions.cols do
                  frontier += oldColumns + column
                  column += 1

    while basis.cols < targetDimension do
      if frontierOffset >= frontier.length then
        val width = math.min(blockWidth, targetDimension - basis.cols)
        val empty =
          GeneralizedBlockKernels.MetricBlock(
            DMat.zeros(basis.rows, 0),
            DMat.zeros(basis.rows, 0)
          )
        GeneralizedBlockKernels
          .bOrthonormalizeAgainstAndReplenish(
            empty,
            metric,
            basis,
            width,
            tolerance,
            streamOffset
          ) match
          case Left(error)        => return Left(error)
          case Right(replenished) =>
            appendBlock(replenished) match
              case Left(error) => return Left(error)
              case Right(_)    => ()
            streamOffset += width
      else
        val take =
          math.min(
            targetDimension - basis.cols,
            math.min(blockWidth, frontier.length - frontierOffset)
          )
        val indices = Array.tabulate(take)(i => frontier(frontierOffset + i))
        frontierOffset += take
        val rightHandSides = selectColumns(operatorImages, indices)
        applyMetricSolveBlock(
          metricSolve,
          rightHandSides,
          outerIteration,
          innerWork
        ) match
          case Left(error)                      => return Left(error)
          case Right((directions, updatedWork)) =>
            innerWork = updatedWork
            GeneralizedBlockKernels.applyBlock(metric, directions) match
              case Left(error)         => return Left(error)
              case Right(metricImages) =>
                val raw =
                  GeneralizedBlockKernels.MetricBlock(directions, metricImages)
                GeneralizedBlockKernels
                  .bOrthonormalizeAgainst(
                    raw,
                    metric,
                    basis,
                    tolerance
                  ) match
                  case Left(error)        => return Left(error)
                  case Right(independent) =>
                    appendBlock(independent) match
                      case Left(error) => return Left(error)
                      case Right(_)    => ()

    // Incremental MGS keeps the Krylov expansion cheap, but a long clustered
    // basis can still accumulate enough drift for its projected B Gram to lose
    // a Cholesky pivot. Stabilize the complete basis once before Rayleigh-Ritz,
    // deterministically replenishing any cancellation-dominated directions.
    GeneralizedBlockKernels
      .bOrthonormalizeAndReplenish(
        basis.vectors,
        metric,
        targetDimension,
        tolerance,
        streamOffset
      )
      .flatMap: stabilized =>
        GeneralizedBlockKernels
          .applyBlock(operator, stabilized.vectors)
          .map: stabilizedImages =>
            ExpandedBasis(
              stabilized,
              stabilizedImages,
              innerWork,
              fullSpace = stabilized.cols == stabilized.rows
            )

  private def applyMetricSolveBlock[B <: DoubleLinearOperator](
      metricSolve: MetricSolveOperator[B],
      rightHandSides: DMat,
      outerIteration: Int,
      initialWork: LinearSolveSummary
  ): Either[LinAlgError, (DMat, LinearSolveSummary)] =
    val output = DMat.newBuilder(metricSolve.size, rightHandSides.cols)
    var work = initialWork
    var column = 0
    while column < rightHandSides.cols do
      metricSolve.solve(rightHandSides.col(column)) match
        case Left(error) =>
          return Left(
            LinAlgError.InnerSolveFailed(
              outerIteration = outerIteration,
              completedSolves = work.solves,
              innerIterations = work.iterations,
              operatorApplications = work.operatorApplications,
              failure = error
            )
          )
        case Right(result) =>
          work = work.append(result.diagnostics)
          if !result.diagnostics.converged then
            return Left(
              LinAlgError.InnerSolveDidNotConverge(
                outerIteration = outerIteration,
                completedSolves = work.solves,
                innerIterations = work.iterations,
                residual = result.diagnostics.residualNorm.getOrElse(Double.PositiveInfinity),
                operatorApplications = work.operatorApplications
              )
            )
          val destination = output.mutableColumn(column)
          var row = 0
          while row < result.solution.length do
            destination(row) = result.solution(row)
            row += 1
      column += 1
    Right((output.result(), work))

  private def rayleighRitz(
      block: GeneralizedBlockKernels.MetricBlock,
      operatorImages: DMat,
      k: Int,
      order: EigenOrder,
      tolerance: Double
  ): Either[LinAlgError, RitzState] =
    for
      projectedA <- GeneralizedBlockKernels.symmetricProjection(
        block.vectors,
        operatorImages
      )
      projectedB <- GeneralizedBlockKernels.symmetricProjection(
        block.vectors,
        block.metricImages
      )
      projected <- GeneralizedBlockKernels.projectedGeneralizedEigen(
        projectedA,
        projectedB,
        k,
        order
      )
      transformed <- GeneralizedBlockKernels.transform(
        block,
        projected.eigenvectors
      )
    yield
      val transformedImages = operatorImages * projected.eigenvectors
      val residualNorms = DVec.tabulate(k): column =>
        val residual = transformedImages.col(column).mutableCopy
        residual.axpyInPlace(
          -projected.eigenvalues(column),
          transformed.metricImages.col(column)
        )
        residual.asVec.norm2
      RitzState(
        transformed,
        transformedImages,
        projected.eigenvalues,
        residualNorms,
        Array.tabulate(k)(column => residualNorms(column) <= tolerance)
      )

  private[spectral] def validate[B <: DoubleLinearOperator](
      operator: DoubleLinearOperator,
      metric: PositiveDefinite[B],
      metricSolve: MetricSolveOperator[B],
      n: Int,
      k: Int,
      order: EigenOrder,
      options: GeneralizedLanczosOptions
  ): Either[LinAlgError, Unit] =
    val targetDimension =
      options.subspaceDimension.getOrElse(math.min(n, math.max(4 * k, 20)))
    if n <= 0 then Left(LinAlgError.InvalidArgument(s"dimension must be positive, got $n"))
    else if operator.rows != operator.cols then
      Left(LinAlgError.NonSquareMatrix(Shape(Rows(operator.rows), Cols(operator.cols))))
    else if metric.rows != metric.cols then
      Left(LinAlgError.NonSquareMatrix(Shape(Rows(metric.rows), Cols(metric.cols))))
    else if operator.rows != n then
      Left(
        LinAlgError.DimensionMismatch(
          Shape(Rows(n), Cols(n)),
          Shape(Rows(operator.rows), Cols(operator.cols))
        )
      )
    else if metric.rows != n then
      Left(
        LinAlgError.DimensionMismatch(
          Shape(Rows(n), Cols(n)),
          Shape(Rows(metric.rows), Cols(metric.cols))
        )
      )
    else if metricSolve.size != n then Left(LinAlgError.VectorLengthMismatch(n, metricSolve.size))
    else if k <= 0 || k >= n then
      Left(
        LinAlgError.InvalidArgument(
          s"k=$k must be in [1, ${n - 1}] for generalized block Lanczos"
        )
      )
    else if order != EigenOrder.SmallestAlgebraic &&
      order != EigenOrder.LargestAlgebraic
    then
      Left(
        LinAlgError.InvalidArgument(
          s"generalized block Lanczos supports smallest or largest algebraic selection, got $order"
        )
      )
    else if !options.tolerance.isFinite || options.tolerance < 0.0 then
      Left(
        LinAlgError.InvalidArgument(
          s"tolerance must be finite and non-negative, got ${options.tolerance}"
        )
      )
    else if options.maxIterations < 0 then
      Left(
        LinAlgError.InvalidArgument(
          s"maxIterations must be non-negative, got ${options.maxIterations}"
        )
      )
    else if targetDimension <= k || targetDimension > n then
      Left(
        LinAlgError.InvalidArgument(
          s"subspaceDimension=$targetDimension must be in [${k + 1}, $n]"
        )
      )
    else if options.returnVectors != EigenVectors.ValuesOnly &&
      options.returnVectors != EigenVectors.Right
    then
      Left(
        LinAlgError.InvalidArgument(
          s"generalized block Lanczos supports ValuesOnly or Right vectors, got ${options.returnVectors}"
        )
      )
    else
      options.initialSubspace match
        case Some(initial) if initial.rows != n || initial.cols != k =>
          Left(
            LinAlgError.InvalidArgument(
              s"initialSubspace must have shape ${n}x$k, got ${initial.rows}x${initial.cols}"
            )
          )
        case _ => Right(())

  private def assemble(
      state: RitzState,
      n: Int,
      requested: Int,
      wantVectors: Boolean,
      iterations: Int,
      extremalityCertified: Boolean,
      innerWork: LinearSolveSummary
  ): EigenDecomposition =
    val selected = state.converged.indices.filter(state.converged).toArray
    val values = DVec.tabulate(selected.length)(i => state.values(selected(i)))
    val residuals =
      DVec.tabulate(selected.length)(i => state.residualNorms(selected(i)))
    val selectedBlock =
      GeneralizedBlockKernels.MetricBlock(
        DMat.tabulate(n, selected.length)((row, col) => state.block.vectors(row, selected(col))),
        DMat.tabulate(n, selected.length)((row, col) => state.block.metricImages(row, selected(col)))
      )
    val vectors =
      if wantVectors then selectedBlock.vectors
      else DMat.zeros(n, 0)
    EigenDecomposition(
      values,
      vectors,
      SpectralDiagnostics(
        requested = requested,
        converged = selected.length,
        residuals = residuals,
        orthogonalityError =
          if wantVectors then GeneralizedBlockKernels.metricOrthogonalityError(selectedBlock)
          else 0.0,
        iterations = iterations,
        rank = None,
        extremalityCertified = extremalityCertified,
        innerSolve = Some(innerWork)
      )
    )

  private def convergedCount(state: RitzState): Int =
    state.converged.count(identity)

  private def selectColumns(matrix: DMat, columns: Array[Int]): DMat =
    DMat.tabulate(matrix.rows, columns.length)((row, col) => matrix(row, columns(col)))

  private def concatenateMatrices(left: DMat, right: DMat): DMat =
    DMat.tabulate(left.rows, left.cols + right.cols): (row, col) =>
      if col < left.cols then left(row, col)
      else right(row, col - left.cols)

  private def deterministicInitial(n: Int, k: Int): DMat =
    val output = DMat.newBuilder(n, k)
    var column = 0
    while column < k do
      var state = 0x2c9277b5 ^ ((column + 1) * 0x9e3779b9)
      val destination = output.mutableColumn(column)
      var row = 0
      while row < n do
        state = state * 1103515245 + 12345
        destination(row) = ((state >>> 9) & 0x7fffff).toDouble /
          0x800000.toDouble * 2.0 - 1.0
        row += 1
      column += 1
    output.result()
