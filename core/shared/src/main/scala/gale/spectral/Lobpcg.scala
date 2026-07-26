package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Cols
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import gale.linalg.Rows
import gale.linalg.Shape
import gale.solvers.Preconditioner

/** Portable locally optimal block preconditioned conjugate-gradient kernel for
  * the symmetric-definite generalized eigenproblem `A X = B X Λ`.
  *
  * Ambient operators are used only through block applications. Dense
  * generalized eigendecompositions are confined to the current trial subspace.
  * The current Ritz block and every retained history block carry cached metric
  * images; `A X` is likewise transformed with Ritz coefficients instead of
  * being recomputed.
  */
private[spectral] object Lobpcg:

  private final case class RitzState(
      block: GeneralizedBlockKernels.MetricBlock,
      operatorImages: DMat,
      values: DVec,
      residuals: DMat,
      residualNorms: DVec,
      converged: Array[Boolean]
  )

  def solve(
      operator: DoubleLinearOperator,
      metric: DoubleLinearOperator,
      n: Int,
      k: Int,
      order: EigenOrder,
      options: GeneralizedSpectralOptions,
      preconditioner: Preconditioner
  ): Either[LinAlgError, EigenDecomposition] =
    validate(operator, metric, n, k, order, options).flatMap: _ =>
      val initial =
        options.initialSubspace
          .map(matrix => DMat.tabulate(matrix.rows, matrix.cols)((row, col) => matrix(row, col)))
          .getOrElse(deterministicInitial(n, k))
      val orthogonalityTolerance =
        math.max(1e-12, math.min(1e-8, options.tolerance * 1e-2))

      GeneralizedBlockKernels
        .bOrthonormalizeAndReplenish(initial, metric, k, orthogonalityTolerance)
        .flatMap: initialBlock =>
          GeneralizedBlockKernels.applyBlock(operator, initialBlock.vectors).flatMap: initialImages =>
            rayleighRitz(initialBlock, initialImages, k, order, options.tolerance).flatMap: initialState =>
              var state = initialState
              var history = emptyBlock(n)
              var iterations = 0
              var extremalityCertified = false
              var invariantStartExplorationPending =
                convergedCount(initialState) == k && k < n
              var failure: Option[LinAlgError] = None

              while iterations < options.maxIterations &&
                  (convergedCount(state) < k || invariantStartExplorationPending) &&
                  failure.isEmpty
              do
                val active = unconvergedIndices(state)
                val residualBlock = selectColumns(state.residuals, active)

                applyPreconditioner(preconditioner, residualBlock) match
                  case Left(error) =>
                    failure = Some(error)
                  case Right(preconditioned) =>
                    GeneralizedBlockKernels.applyBlock(metric, preconditioned) match
                      case Left(error) =>
                        failure = Some(error)
                      case Right(preconditionedMetricImages) =>
                        val rawDirections =
                          GeneralizedBlockKernels.MetricBlock(preconditioned, preconditionedMetricImages)
                        // An exactly invariant initial block can be the wrong
                        // spectral end. Force one deterministic complement
                        // probe instead of accepting residual zero as target
                        // identity.
                        val requestedDirections =
                          if active.nonEmpty then active.length else k
                        val directionCount = math.min(requestedDirections, n - k)
                        GeneralizedBlockKernels.bOrthonormalizeAgainstAndReplenish(
                          rawDirections,
                          metric,
                          state.block,
                          directionCount,
                          orthogonalityTolerance,
                          streamOffset = (iterations + 1) * math.max(k, 1)
                        ) match
                          case Left(error) =>
                            failure = Some(error)
                          case Right(directions) =>
                            GeneralizedBlockKernels.concatenate(state.block, directions) match
                              case Left(error) =>
                                failure = Some(error)
                              case Right(stateAndDirections) =>
                                GeneralizedBlockKernels
                                  .bOrthonormalizeAgainst(history, stateAndDirections, orthogonalityTolerance)
                                  .flatMap(candidate =>
                                    retainStableHistory(
                                      candidate,
                                      stateAndDirections,
                                      orthogonalityTolerance
                                    )
                                  )
                                match
                                  case Left(error) =>
                                    failure = Some(error)
                                  case Right(usableHistory) =>
                                    GeneralizedBlockKernels.concatenate(
                                      state.block,
                                      directions,
                                      usableHistory
                                    ) match
                                      case Left(error) =>
                                        failure = Some(error)
                                      case Right(trial) =>
                                        val directionImagesEither =
                                          GeneralizedBlockKernels.applyBlock(operator, directions.vectors)
                                        val historyImagesEither =
                                          GeneralizedBlockKernels.applyBlock(operator, usableHistory.vectors)

                                        (directionImagesEither, historyImagesEither) match
                                          case (Left(error), _) =>
                                            failure = Some(error)
                                          case (_, Left(error)) =>
                                            failure = Some(error)
                                          case (Right(directionImages), Right(historyImages)) =>
                                            val trialImages =
                                              concatenateMatrices(
                                                state.operatorImages,
                                                directionImages,
                                                historyImages
                                              )
                                            rayleighRitz(
                                              trial,
                                              trialImages,
                                              k,
                                              order,
                                              options.tolerance
                                            ) match
                                              case Left(error) =>
                                                failure = Some(error)
                                              case Right(nextState) =>
                                                locallyOptimalHistory(
                                                  previous = state.block,
                                                  current = nextState.block,
                                                  tolerance = orthogonalityTolerance
                                                ) match
                                                  case Left(error) =>
                                                    failure = Some(error)
                                                  case Right(nextHistory) =>
                                                    extremalityCertified =
                                                      extremalityCertified || trial.cols == n
                                                    invariantStartExplorationPending = false
                                                    history = nextHistory
                                                    state = nextState
                                                    iterations += 1

              failure match
                case Some(error) => Left(error)
                case None =>
                  Right(
                    assemble(
                      state,
                      n,
                      k,
                      options.returnVectors == EigenVectors.Right,
                      iterations,
                      extremalityCertified
                    )
                  )

  private[spectral] def validate(
      operator: DoubleLinearOperator,
      metric: DoubleLinearOperator,
      n: Int,
      k: Int,
      order: EigenOrder,
      options: GeneralizedSpectralOptions
  ): Either[LinAlgError, Unit] =
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
    else if k <= 0 || k >= n then
      Left(
        LinAlgError.InvalidArgument(
          s"k=$k must be in [1, ${n - 1}] for matrix-free LOBPCG"
        )
      )
    else if order != EigenOrder.SmallestAlgebraic && order != EigenOrder.LargestAlgebraic then
      Left(LinAlgError.InvalidArgument(s"LOBPCG supports smallest or largest algebraic selection, got $order"))
    else if options.tolerance < 0.0 || !options.tolerance.isFinite then
      Left(LinAlgError.InvalidArgument(s"tolerance must be finite and non-negative, got ${options.tolerance}"))
    else if options.maxIterations < 0 then
      Left(LinAlgError.InvalidArgument(s"maxIterations must be non-negative, got ${options.maxIterations}"))
    else if options.returnVectors != EigenVectors.ValuesOnly && options.returnVectors != EigenVectors.Right then
      Left(
        LinAlgError.InvalidArgument(
          s"generalized symmetric LOBPCG supports ValuesOnly or Right vectors, got ${options.returnVectors}"
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

  private def rayleighRitz(
      block: GeneralizedBlockKernels.MetricBlock,
      operatorImages: DMat,
      k: Int,
      order: EigenOrder,
      tolerance: Double
  ): Either[LinAlgError, RitzState] =
    for
      projectedA <- GeneralizedBlockKernels.symmetricProjection(block.vectors, operatorImages)
      projectedB <- GeneralizedBlockKernels.symmetricProjection(block.vectors, block.metricImages)
      projected <- GeneralizedBlockKernels.projectedGeneralizedEigen(projectedA, projectedB, k, order)
      transformed <- GeneralizedBlockKernels.transform(block, projected.eigenvectors)
    yield
      val transformedImages = operatorImages * projected.eigenvectors
      val residuals = DMat.tabulate(transformed.rows, k): (row, col) =>
        transformedImages(row, col) - transformed.metricImages(row, col) * projected.eigenvalues(col)
      val residualNorms = DVec.tabulate(k)(col => residuals.col(col).norm2)
      val converged = Array.tabulate(k)(col => residualNorms(col) <= tolerance)
      RitzState(
        transformed,
        transformedImages,
        projected.eigenvalues,
        residuals,
        residualNorms,
        converged
      )

  private def locallyOptimalHistory(
      previous: GeneralizedBlockKernels.MetricBlock,
      current: GeneralizedBlockKernels.MetricBlock,
      tolerance: Double
  ): Either[LinAlgError, GeneralizedBlockKernels.MetricBlock] =
    val overlap = previous.vectors.t * current.metricImages
    val raw = GeneralizedBlockKernels.MetricBlock(
      current.vectors - previous.vectors * overlap,
      current.metricImages - previous.metricImages * overlap
    )
    GeneralizedBlockKernels.bOrthonormalizeAgainst(raw, current, tolerance)

  private def retainStableHistory(
      candidate: GeneralizedBlockKernels.MetricBlock,
      stateAndDirections: GeneralizedBlockKernels.MetricBlock,
      tolerance: Double
  ): Either[LinAlgError, GeneralizedBlockKernels.MetricBlock] =
    GeneralizedBlockKernels.concatenate(stateAndDirections, candidate).map: trial =>
      val guard =
        math.max(1e-4, 10.0 * tolerance * math.sqrt(math.max(1, trial.cols).toDouble))
      if GeneralizedBlockKernels.metricOrthogonalityError(trial) <= guard then candidate
      else emptyBlock(candidate.rows)

  private def applyPreconditioner(
      preconditioner: Preconditioner,
      residuals: DMat
  ): Either[LinAlgError, DMat] =
    try
      val output = DMat.newBuilder(residuals.rows, residuals.cols)
      var column = 0
      while column < residuals.cols do
        // A user preconditioner may retain its destination. Copy out of an
        // isolated temporary so no mutable alias reaches the returned block.
        val temporary = MutableDVec.zeros(residuals.rows)
        preconditioner.solve(residuals.col(column), temporary)
        val destination = output.mutableColumn(column)
        var row = 0
        while row < residuals.rows do
          destination(row) = temporary(row)
          row += 1
        column += 1
      Right(output.result())
    catch case error: LinAlgError => Left(error)

  private def assemble(
      state: RitzState,
      n: Int,
      requested: Int,
      wantVectors: Boolean,
      iterations: Int,
      extremalityCertified: Boolean
  ): EigenDecomposition =
    val selected = state.converged.indices.filter(state.converged).toArray
    val values = DVec.tabulate(selected.length)(i => state.values(selected(i)))
    val residuals = DVec.tabulate(selected.length)(i => state.residualNorms(selected(i)))
    val vectors =
      if wantVectors then DMat.tabulate(n, selected.length)((row, col) => state.block.vectors(row, selected(col)))
      else DMat.zeros(n, 0)
    val selectedMetricImages =
      GeneralizedBlockKernels.MetricBlock(
        DMat.tabulate(n, selected.length)((row, col) => state.block.vectors(row, selected(col))),
        DMat.tabulate(n, selected.length)((row, col) => state.block.metricImages(row, selected(col)))
      )
    EigenDecomposition(
      values,
      vectors,
      SpectralDiagnostics(
        requested = requested,
        converged = selected.length,
        residuals = residuals,
        orthogonalityError =
          if wantVectors then GeneralizedBlockKernels.metricOrthogonalityError(selectedMetricImages)
          else 0.0,
        iterations = iterations,
        rank = None,
        extremalityCertified = extremalityCertified
      )
    )

  private def convergedCount(state: RitzState): Int =
    state.converged.count(identity)

  private def unconvergedIndices(state: RitzState): Array[Int] =
    state.converged.indices.filterNot(state.converged).toArray

  private def selectColumns(matrix: DMat, columns: Array[Int]): DMat =
    DMat.tabulate(matrix.rows, columns.length)((row, col) => matrix(row, columns(col)))

  private def concatenateMatrices(blocks: DMat*): DMat =
    val rows = blocks.head.rows
    val offsets = new Array[Int](blocks.length + 1)
    var block = 0
    while block < blocks.length do
      offsets(block + 1) = offsets(block) + blocks(block).cols
      block += 1
    DMat.tabulate(rows, offsets.last): (row, col) =>
      var index = 0
      while index + 1 < offsets.length && col >= offsets(index + 1) do
        index += 1
      blocks(index)(row, col - offsets(index))

  private def emptyBlock(n: Int): GeneralizedBlockKernels.MetricBlock =
    GeneralizedBlockKernels.MetricBlock(DMat.zeros(n, 0), DMat.zeros(n, 0))

  private def deterministicInitial(n: Int, k: Int): DMat =
    val output = DMat.newBuilder(n, k)
    var column = 0
    while column < k do
      var state = 0x4f1bbcdc ^ ((column + 1) * 0x9e3779b9)
      val destination = output.mutableColumn(column)
      var row = 0
      while row < n do
        state = state * 1103515245 + 12345
        destination(row) =
          ((state >>> 9) & 0x7fffff).toDouble / 0x800000.toDouble * 2.0 - 1.0
        row += 1
      column += 1
    output.result()
