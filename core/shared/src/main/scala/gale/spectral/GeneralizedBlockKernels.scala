package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DenseDecompositions
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import gale.linalg.TriangularSolve
import scala.collection.mutable.ArrayBuffer

/** Portable block primitives shared by matrix-free generalized
  * symmetric-definite eigensolvers.
  *
  * All returned matrices own their storage. Mutable destinations and work
  * vectors are confined to this object because user operators may retain a
  * destination passed to `applyTo`.
  */
private[spectral] object GeneralizedBlockKernels:

  /** A block `X` and its aligned metric image `B X`. */
  final case class MetricBlock(vectors: DMat, metricImages: DMat):
    require(
      vectors.rows == metricImages.rows && vectors.cols == metricImages.cols,
      s"metric block shape mismatch: X=${vectors.rows}x${vectors.cols}, BX=${metricImages.rows}x${metricImages.cols}"
    )

    def rows: Int = vectors.rows
    def cols: Int = vectors.cols

  /** Apply an operator to every block column. The operator overload owns the
    * result and totalizes `LinAlgError` thrown by a primitive implementation.
    */
  def applyBlock(operator: DoubleLinearOperator, input: DMat): Either[LinAlgError, DMat] =
    if input.rows != operator.cols then
      Left(
        LinAlgError.InvalidArgument(
          s"operator block input has ${input.rows} rows but operator expects ${operator.cols}"
        )
      )
    else
      try
        val output = DMat.newBuilder(operator.rows, input.cols)
        var column = 0
        while column < input.cols do
          // Never expose the result builder to a user operator: it may retain
          // the mutable destination and mutate it after returning.
          val temporary = MutableDVec.zeros(operator.rows)
          operator.applyTo(input.col(column), temporary)
          val destination = output.mutableColumn(column)
          var row = 0
          while row < operator.rows do
            destination(row) = temporary(row)
            row += 1
          column += 1
        requireFinite(output.result(), "operator block image")
      catch case error: LinAlgError => Left(error)

  /** Symmetric Rayleigh/Gram projection `0.5 * (XᵀY + YᵀX)`. */
  def symmetricProjection(left: DMat, appliedRight: DMat): Either[LinAlgError, DMat] =
    if left.rows != appliedRight.rows || left.cols != appliedRight.cols then
      Left(
        LinAlgError.InvalidArgument(
          s"symmetric projection requires aligned blocks, got ${left.rows}x${left.cols} and " +
            s"${appliedRight.rows}x${appliedRight.cols}"
        )
      )
    else if left.cols == 0 then Right(DMat.zeros(0, 0))
    else
      requireFinite(left, "projection basis").flatMap: finiteLeft =>
        requireFinite(appliedRight, "projected operator image").map: finiteImage =>
          val raw = finiteLeft.t * finiteImage
          raw.symmetrizedAverage

  /** `‖Xᵀ B X - I‖_F` for an aligned metric block. */
  def metricOrthogonalityError(block: MetricBlock): Double =
    val gram = block.vectors.t * block.metricImages
    var sum = 0.0
    var i = 0
    while i < gram.rows do
      var j = 0
      while j < gram.cols do
        val delta = if i == j then gram(i, j) - 1.0 else gram(i, j)
        sum += delta * delta
        j += 1
      i += 1
    math.sqrt(sum)

  /** B-orthonormalize the independent part of `candidates`.
    *
    * A full-rank block uses a Cholesky whitening transform. If its B-Gram is
    * rank-deficient, twice-reorthogonalized modified Gram-Schmidt retains the
    * independent columns instead. An encountered negative B norm returns
    * [[LinAlgError.NotPositiveDefinite]] rather than silently treating an
    * indefinite direction as dependence.
    */
  def bOrthonormalize(
      candidates: DMat,
      metric: DoubleLinearOperator,
      tolerance: Double = 1e-12
  ): Either[LinAlgError, MetricBlock] =
    validate(candidates, metric, candidates.cols, tolerance).flatMap: _ =>
      applyBlock(metric, candidates).flatMap: metricImages =>
        choleskyWhiten(candidates, metricImages, tolerance) match
          case Some(result) => Right(result)
          case None         => rankReveal(candidates, metricImages, metric, candidates.cols, tolerance, replenish = false)

  /** B-orthonormalize and deterministically replenish a block to exactly
    * `targetColumns`.
    *
    * Dependent caller columns are not failures. They are removed, then portable
    * deterministic probes and finally canonical directions fill the missing
    * span. Failure to reach the target means the encountered metric geometry is
    * not positive-definite on enough independent directions.
    */
  def bOrthonormalizeAndReplenish(
      candidates: DMat,
      metric: DoubleLinearOperator,
      targetColumns: Int,
      tolerance: Double = 1e-12,
      streamOffset: Int = 0
  ): Either[LinAlgError, MetricBlock] =
    validate(candidates, metric, targetColumns, tolerance).flatMap: _ =>
      applyBlock(metric, candidates).flatMap: metricImages =>
        if candidates.cols == targetColumns then
          choleskyWhiten(candidates, metricImages, tolerance) match
            case Some(result) => Right(result)
            case None =>
              rankReveal(
                candidates,
                metricImages,
                metric,
                targetColumns,
                tolerance,
                replenish = true,
                streamOffset
              )
        else
          rankReveal(
            candidates,
            metricImages,
            metric,
            targetColumns,
            tolerance,
            replenish = true,
            streamOffset
          )

  /** Apply a coefficient matrix to both aligned halves of a metric block. */
  def transform(block: MetricBlock, coefficients: DMat): Either[LinAlgError, MetricBlock] =
    if coefficients.rows != block.cols then
      Left(
        LinAlgError.InvalidArgument(
          s"block transform requires ${block.cols} coefficient rows, got ${coefficients.rows}"
        )
      )
    else
      requireFinite(coefficients, "block coefficients").map: finite =>
        MetricBlock(block.vectors * finite, block.metricImages * finite)

  /** Concatenate aligned metric blocks without reapplying the metric. */
  def concatenate(blocks: MetricBlock*): Either[LinAlgError, MetricBlock] =
    if blocks.isEmpty then
      Left(LinAlgError.InvalidArgument("at least one metric block is required"))
    else
      val rows = blocks.head.rows
      if blocks.exists(_.rows != rows) then
        Left(LinAlgError.InvalidArgument("all concatenated metric blocks must have the same row count"))
      else
        val offsets = new Array[Int](blocks.length + 1)
        var blockIndex = 0
        while blockIndex < blocks.length do
          offsets(blockIndex + 1) = offsets(blockIndex) + blocks(blockIndex).cols
          blockIndex += 1

        def elementAt(metricImage: Boolean, row: Int, column: Int): Double =
          var index = 0
          while index + 1 < offsets.length && column >= offsets(index + 1) do
            index += 1
          val block = blocks(index)
          val localColumn = column - offsets(index)
          if metricImage then block.metricImages(row, localColumn)
          else block.vectors(row, localColumn)

        val totalColumns = offsets.last
        Right(
          MetricBlock(
            DMat.tabulate(rows, totalColumns)((row, col) => elementAt(false, row, col)),
            DMat.tabulate(rows, totalColumns)((row, col) => elementAt(true, row, col))
          )
        )

  /** B-orthonormalize the independent part of an already imaged block against
    * an existing B-orthonormal block. No metric application is repeated.
    */
  def bOrthonormalizeAgainst(
      candidates: MetricBlock,
      against: MetricBlock,
      tolerance: Double = 1e-12
  ): Either[LinAlgError, MetricBlock] =
    rankRevealAgainst(
      candidates,
      against,
      targetColumns = candidates.cols,
      metric = None,
      tolerance,
      replenish = false,
      streamOffset = 0
    )

  /** B-orthonormalize against an existing block with access to the metric.
    *
    * The metric is normally not reapplied. It is used only if separately
    * updated `x` and `B x` work vectors produce a negative squared norm, so the
    * implementation can distinguish accumulated cancellation from genuinely
    * indefinite geometry.
    */
  def bOrthonormalizeAgainst(
      candidates: MetricBlock,
      metric: DoubleLinearOperator,
      against: MetricBlock,
      tolerance: Double
  ): Either[LinAlgError, MetricBlock] =
    rankRevealAgainst(
      candidates,
      against,
      targetColumns = candidates.cols,
      metric = Some(metric),
      tolerance,
      replenish = false,
      streamOffset = 0
    )

  /** B-orthonormalize against an existing block and deterministically replenish
    * the independent complement to `targetColumns`.
    *
    * Supplied candidate metric images are consumed directly. The metric is
    * applied only to replacement probes.
    */
  def bOrthonormalizeAgainstAndReplenish(
      candidates: MetricBlock,
      metric: DoubleLinearOperator,
      against: MetricBlock,
      targetColumns: Int,
      tolerance: Double = 1e-12,
      streamOffset: Int = 0
  ): Either[LinAlgError, MetricBlock] =
    rankRevealAgainst(
      candidates,
      against,
      targetColumns,
      metric = Some(metric),
      tolerance,
      replenish = true,
      streamOffset
    )

  /** Solve a small dense projected symmetric-definite problem. */
  def projectedGeneralizedEigen(
      projectedA: DMat,
      projectedB: DMat,
      count: Int,
      order: EigenOrder
  ): Either[LinAlgError, EigenDecomposition] =
    Eigen.eigSymmetricGeneralized(
      projectedA,
      projectedB,
      EigenSelection.Count(count, order),
      EigenVectors.Right
    )

  private def validate(
      candidates: DMat,
      metric: DoubleLinearOperator,
      targetColumns: Int,
      tolerance: Double
  ): Either[LinAlgError, Unit] =
    if metric.rows != metric.cols then
      Left(LinAlgError.InvalidArgument(s"metric operator must be square, got ${metric.rows}x${metric.cols}"))
    else if candidates.rows != metric.cols then
      Left(
        LinAlgError.InvalidArgument(
          s"candidate block row count ${candidates.rows} does not match metric order ${metric.cols}"
        )
      )
    else if targetColumns <= 0 || targetColumns > metric.cols then
      Left(
        LinAlgError.InvalidArgument(
          s"target block width $targetColumns must be in [1, ${metric.cols}]"
        )
      )
    else if candidates.cols > targetColumns then
      Left(
        LinAlgError.InvalidArgument(
          s"candidate block has ${candidates.cols} columns but target width is $targetColumns"
        )
      )
    else if tolerance < 0.0 || !tolerance.isFinite then
      Left(LinAlgError.InvalidArgument(s"B-orthonormalization tolerance must be finite and non-negative, got $tolerance"))
    else requireFinite(candidates, "candidate block").map(_ => ())

  /** Fast full-rank path: if `G = Xᵀ B X = L Lᵀ`, then
    * `Q = X L⁻ᵀ` satisfies `Qᵀ B Q = I`.
    */
  private def choleskyWhiten(
      candidates: DMat,
      metricImages: DMat,
      tolerance: Double
  ): Option[MetricBlock] =
    if candidates.cols == 0 then Some(MetricBlock(candidates, metricImages))
    else
      val gram = (candidates.t * metricImages).symmetrizedAverage
      DenseDecompositions.cholesky(gram).toOption.flatMap: factor =>
        inverseUpper(factor.lower.t).toOption.flatMap: transform =>
          val result = MetricBlock(candidates * transform, metricImages * transform)
          val error = metricOrthogonalityError(result)
          val guard = math.max(1e-10, 100.0 * tolerance * math.sqrt(candidates.cols.toDouble))
          Option.when(error.isFinite && error <= guard)(result)

  private def inverseUpper(upper: DMat): Either[LinAlgError, DMat] =
    val n = upper.rows
    val columns = new Array[DVec](n)
    var j = 0
    while j < n do
      val unit = DVec.tabulate(n)(i => if i == j then 1.0 else 0.0)
      TriangularSolve.upper(upper, unit) match
        case Left(error) => return Left(error)
        case Right(col)  => columns(j) = col
      j += 1
    Right(DMat.tabulate(n, n)((i, col) => columns(col)(i)))

  private def rankReveal(
      candidates: DMat,
      initialMetricImages: DMat,
      metric: DoubleLinearOperator,
      targetColumns: Int,
      tolerance: Double,
      replenish: Boolean,
      streamOffset: Int = 0
  ): Either[LinAlgError, MetricBlock] =
    val vectors = ArrayBuffer.empty[DVec]
    val metricImages = ArrayBuffer.empty[DVec]

    def add(candidate: DVec, suppliedMetricImage: Option[DVec]): Either[LinAlgError, Boolean] =
      requireFinite(candidate, "B-orthogonalization candidate").flatMap: finiteCandidate =>
        val imageEither =
          suppliedMetricImage match
            case Some(image) => requireFinite(image, "candidate metric image")
            case None        => applyVector(metric, finiteCandidate)
        imageEither.flatMap: finiteImage =>
          val work = finiteCandidate.mutableCopy
          val metricWork = finiteImage.mutableCopy
          var pass = 0
          while pass < 2 do
            var i = 0
            while i < vectors.length do
              val coefficient = vectors(i).dot(metricWork.asVec)
              work.axpyInPlace(-coefficient, vectors(i))
              metricWork.axpyInPlace(-coefficient, metricImages(i))
              i += 1
            pass += 1

          val initialNormSquared = work.asVec.dot(metricWork.asVec)
          val checkedMetricImage =
            if initialNormSquared < 0.0 then
              applyVector(metric, work.asVec)
            else Right(metricWork.asVec)
          checkedMetricImage.flatMap: refreshedMetricWork =>
            if initialNormSquared < 0.0 then
              metricWork := refreshedMetricWork
            val normSquared = work.asVec.dot(metricWork.asVec)
            val scale = math.max(1.0, finiteCandidate.norm2 * finiteImage.norm2)
            val dependenceThreshold = tolerance * tolerance * scale
            val negativeThreshold =
              negativeBNormThreshold(
                dependenceThreshold,
                scale,
                vectors.length
              )
            val effectiveDependenceThreshold =
              if initialNormSquared < 0.0 then
                math.max(dependenceThreshold, tolerance * scale)
              else dependenceThreshold
            if !normSquared.isFinite then
              Left(LinAlgError.InvalidArgument("B-orthogonalization produced a non-finite squared norm"))
            else if normSquared < -negativeThreshold then
              Left(LinAlgError.NotPositiveDefinite(vectors.length))
            else if normSquared <= effectiveDependenceThreshold || normSquared < 0.0 then
              Right(false)
            else
              val inverseNorm = 1.0 / math.sqrt(normSquared)
              val vector = (work.asVec * inverseNorm).copy
              val metricImage = (metricWork.asVec * inverseNorm).copy
              vectors += vector
              metricImages += metricImage
              Right(true)

    var col = 0
    while col < candidates.cols do
      add(candidates.col(col), Some(initialMetricImages.col(col))) match
        case Left(error) => return Left(error)
        case Right(_)    => ()
      col += 1

    if replenish then
      var stream = streamOffset
      val probeLimit = math.max(8, 4 * metric.cols + targetColumns)
      while vectors.length < targetColumns && stream < streamOffset + probeLimit do
        add(deterministicVector(metric.cols, stream), None) match
          case Left(error) => return Left(error)
          case Right(_)    => ()
        stream += 1

      var coordinate = 0
      while vectors.length < targetColumns && coordinate < metric.cols do
        val unit = DVec.tabulate(metric.cols)(i => if i == coordinate then 1.0 else 0.0)
        add(unit, None) match
          case Left(error) => return Left(error)
          case Right(_)    => ()
        coordinate += 1

    if replenish && vectors.length < targetColumns then
      Left(LinAlgError.NotPositiveDefinite(vectors.length))
    else
      Right(
        MetricBlock(
          DMat.tabulate(metric.cols, vectors.length)((i, j) => vectors(j)(i)),
          DMat.tabulate(metric.rows, metricImages.length)((i, j) => metricImages(j)(i))
        )
      )

  private def rankRevealAgainst(
      candidates: MetricBlock,
      against: MetricBlock,
      targetColumns: Int,
      metric: Option[DoubleLinearOperator],
      tolerance: Double,
      replenish: Boolean,
      streamOffset: Int
  ): Either[LinAlgError, MetricBlock] =
    if candidates.rows != against.rows then
      Left(
        LinAlgError.InvalidArgument(
          s"candidate and existing blocks must have the same row count, got ${candidates.rows} and ${against.rows}"
        )
      )
    else if tolerance < 0.0 || !tolerance.isFinite then
      Left(LinAlgError.InvalidArgument(s"B-orthonormalization tolerance must be finite and non-negative, got $tolerance"))
    else if targetColumns < 0 || (replenish && targetColumns > candidates.rows - against.cols) then
      Left(
        LinAlgError.InvalidArgument(
          s"target complement width $targetColumns must be in [0, ${candidates.rows - against.cols}]"
        )
      )
    else if replenish && metric.isEmpty then
      Left(LinAlgError.InvalidArgument("deterministic complement replenishment requires a metric operator"))
    else if replenish && metric.exists(op => op.rows != candidates.rows || op.cols != candidates.rows) then
      Left(
        LinAlgError.InvalidArgument(
          s"metric operator must match block order ${candidates.rows}"
        )
      )
    else
      for
        _ <- requireFinite(candidates.vectors, "candidate block")
        _ <- requireFinite(candidates.metricImages, "candidate metric image")
        _ <- requireFinite(against.vectors, "existing basis")
        _ <- requireFinite(against.metricImages, "existing metric image")
        result <- rankRevealAgainstFinite(
          candidates,
          against,
          targetColumns,
          metric,
          tolerance,
          replenish,
          streamOffset
        )
      yield result

  private def rankRevealAgainstFinite(
      candidates: MetricBlock,
      against: MetricBlock,
      targetColumns: Int,
      metric: Option[DoubleLinearOperator],
      tolerance: Double,
      replenish: Boolean,
      streamOffset: Int
  ): Either[LinAlgError, MetricBlock] =
    val newVectors = ArrayBuffer.empty[DVec]
    val newMetricImages = ArrayBuffer.empty[DVec]

    def add(candidate: DVec, candidateMetricImage: DVec): Either[LinAlgError, Boolean] =
      val work = candidate.mutableCopy
      val metricWork = candidateMetricImage.mutableCopy
      var pass = 0
      while pass < 2 do
        var i = 0
        while i < against.cols do
          val coefficient = against.vectors.col(i).dot(metricWork.asVec)
          work.axpyInPlace(-coefficient, against.vectors.col(i))
          metricWork.axpyInPlace(-coefficient, against.metricImages.col(i))
          i += 1
        i = 0
        while i < newVectors.length do
          val coefficient = newVectors(i).dot(metricWork.asVec)
          work.axpyInPlace(-coefficient, newVectors(i))
          metricWork.axpyInPlace(-coefficient, newMetricImages(i))
          i += 1
        pass += 1

      var normSquared = work.asVec.dot(metricWork.asVec)
      var refreshed = false
      if normSquared < 0.0 then
        metric match
          case Some(operator) =>
            applyVector(operator, work.asVec) match
              case Left(error) => return Left(error)
              case Right(refreshedMetricWork) =>
                metricWork := refreshedMetricWork
                normSquared = work.asVec.dot(metricWork.asVec)
                refreshed = true
          case None => ()
      val scale = math.max(1.0, candidate.norm2 * candidateMetricImage.norm2)
      val dependenceThreshold = tolerance * tolerance * scale
      val negativeThreshold =
        negativeBNormThreshold(
          dependenceThreshold,
          scale,
          against.cols + newVectors.length
        )
      val effectiveDependenceThreshold =
        if refreshed then
          math.max(dependenceThreshold, tolerance * scale)
        else dependenceThreshold
      if !normSquared.isFinite then
        Left(LinAlgError.InvalidArgument("B-orthogonalization produced a non-finite squared norm"))
      else if normSquared < -negativeThreshold then
        Left(LinAlgError.NotPositiveDefinite(against.cols + newVectors.length))
      else if normSquared <= effectiveDependenceThreshold || normSquared < 0.0 then
        Right(false)
      else
        val inverseNorm = 1.0 / math.sqrt(normSquared)
        newVectors += (work.asVec * inverseNorm).copy
        newMetricImages += (metricWork.asVec * inverseNorm).copy
        Right(true)

    var column = 0
    while column < candidates.cols && newVectors.length < targetColumns do
      add(candidates.vectors.col(column), candidates.metricImages.col(column)) match
        case Left(error) => return Left(error)
        case Right(_)    => ()
      column += 1

    if replenish then
      val operator =
        metric match
          case Some(value) => value
          case None =>
            return Left(
              LinAlgError.InvalidArgument(
                "deterministic complement replenishment requires a metric operator"
              )
            )
      var stream = streamOffset
      val probeLimit = math.max(8, 4 * candidates.rows + targetColumns)
      while newVectors.length < targetColumns && stream < streamOffset + probeLimit do
        val probe = deterministicVector(candidates.rows, stream)
        applyVector(operator, probe) match
          case Left(error) => return Left(error)
          case Right(image) =>
            add(probe, image) match
              case Left(error) => return Left(error)
              case Right(_)    => ()
        stream += 1

      var coordinate = 0
      while newVectors.length < targetColumns && coordinate < candidates.rows do
        val unit = DVec.tabulate(candidates.rows)(i => if i == coordinate then 1.0 else 0.0)
        applyVector(operator, unit) match
          case Left(error) => return Left(error)
          case Right(image) =>
            add(unit, image) match
              case Left(error) => return Left(error)
              case Right(_)    => ()
        coordinate += 1

    if replenish && newVectors.length < targetColumns then
      Left(LinAlgError.NotPositiveDefinite(against.cols + newVectors.length))
    else
      Right(
        MetricBlock(
          DMat.tabulate(candidates.rows, newVectors.length)((i, j) => newVectors(j)(i)),
          DMat.tabulate(candidates.rows, newMetricImages.length)((i, j) => newMetricImages(j)(i))
        )
      )

  /** Dependence/indefiniteness threshold for a twice-reorthogonalized B norm.
    *
    * Subtracting an almost represented direction updates `x` and `B x`
    * separately. Their final dot product can therefore be a few accumulated
    * ulps negative even for an exactly SPD metric. Negative values within this
    * basis-size-scaled roundoff envelope are dependence, while the ordinary
    * positive dependence threshold remains caller-controlled and a materially
    * negative norm still reports `NotPositiveDefinite`. A direction whose
    * metric image had to be refreshed must also clear the unsquared requested
    * tolerance before normalization; this avoids amplifying a
    * cancellation-dominated remainder into the basis.
    */
  private def negativeBNormThreshold(
      dependenceThreshold: Double,
      scale: Double,
      basisSize: Int
  ): Double =
    val roundoff =
      8.0 * 2.220446049250313e-16 * math.max(1, basisSize).toDouble * scale
    math.max(dependenceThreshold, roundoff)

  private def applyVector(operator: DoubleLinearOperator, vector: DVec): Either[LinAlgError, DVec] =
    if vector.length != operator.cols then
      Left(LinAlgError.VectorLengthMismatch(operator.cols, vector.length))
    else
      try
        val destination = MutableDVec.zeros(operator.rows)
        operator.applyTo(vector, destination)
        requireFinite(destination.toVec, "operator vector image")
      catch case error: LinAlgError => Left(error)

  private def requireFinite(matrix: DMat, label: String): Either[LinAlgError, DMat] =
    var i = 0
    while i < matrix.rows do
      var j = 0
      while j < matrix.cols do
        if !matrix(i, j).isFinite then
          return Left(LinAlgError.InvalidArgument(s"$label contains a non-finite value at ($i, $j)"))
        j += 1
      i += 1
    Right(matrix)

  private def requireFinite(vector: DVec, label: String): Either[LinAlgError, DVec] =
    var i = 0
    while i < vector.length do
      if !vector(i).isFinite then
        return Left(LinAlgError.InvalidArgument(s"$label contains a non-finite value at index $i"))
      i += 1
    Right(vector)

  /** Portable deterministic probe. Integer overflow is defined identically on
    * JVM and Scala.js.
    */
  private def deterministicVector(n: Int, stream: Int): DVec =
    var state = 0x6d2b79f5 ^ (stream * 0x9e3779b9)
    DVec.tabulate(n): _ =>
      state = state * 1103515245 + 12345
      ((state >>> 9) & 0x7fffff).toDouble / 0x800000.toDouble * 2.0 - 1.0
