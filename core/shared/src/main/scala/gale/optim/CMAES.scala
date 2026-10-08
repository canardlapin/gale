package gale.optim

import gale.linalg.{DMat, DVec}
import gale.spectral.{Eigen, EigenSelection}
import scala.util.control.NonFatal

/** Full-covariance CMA-ES for continuous black-box objectives. Sessions are single-owner and not thread-safe;
  * evaluating an immutable batch in parallel is safe. No objective calls are made by ask/tell itself.
  */
object CMAES:
  def start(
      initialMean: DVec,
      config: CMAESConfig = CMAESConfig(),
      bounds: Option[MatrixBoxBounds] = None
  ): Either[CMAESError, CMAESSession] =
    import CMAESError.InvalidConfiguration
    val size = initialMean.length
    if size == 0 || !config.initialStepSize.isFinite || config.initialStepSize <= 0.0 ||
      config.populationSize < 0 || config.populationSize == 1 || config.maxEvaluations <= 0 ||
      config.maxGenerations < 0 || !config.stepTolerance.isFinite || config.stepTolerance < 0.0 ||
      config.stagnationGenerations < 0 || !config.maxCondition.isFinite || config.maxCondition <= 1.0 ||
      config.maxSamplingAttempts <= 0 || config.eigenUpdatePeriod < 0 || config.targetValue.exists(!_.isFinite)
    then Left(InvalidConfiguration("invalid CMA-ES dimensions, scale, population, limits, or stopping settings"))
    else if bounds.exists(b => b.lower.rows != size || b.lower.cols != 1) then
      Left(InvalidConfiguration("CMA-ES bounds must be a column matching the initial vector"))
    else
      val lower = Array.tabulate(size)(i => bounds.fold(Double.NegativeInfinity)(_.lower(i, 0)))
      val upper = Array.tabulate(size)(i => bounds.fold(Double.PositiveInfinity)(_.upper(i, 0)))
      if (0 until size).exists(i => !initialMean(i).isFinite || initialMean(i) < lower(i) || initialMean(i) > upper(i))
      then Left(InvalidConfiguration("the initial mean must be finite and feasible"))
      else
        val free = (0 until size).filter(i => lower(i) != upper(i)).toArray
        val n = math.max(1, free.length)
        val population =
          if config.populationSize == 0 then 4 + (3.0 * math.log(n.toDouble)).toInt else config.populationSize
        if n.toLong * n > Int.MaxValue || population.toLong * size > Int.MaxValue then
          Left(InvalidConfiguration("CMA-ES dense storage exceeds indexed array capacity"))
        else Right(new CMAESSession(initialMean.copy, lower, upper, free, config, population))

  /** Evaluate owned candidates serially. Callback failure preserves all preceding finite observations and charges the
    * failed call. Cancellation is checked before every evaluation. Partial generations never adapt the distribution.
    */
  def minimize(
      objective: DVec => Double,
      initialMean: DVec,
      config: CMAESConfig = CMAESConfig(),
      bounds: Option[MatrixBoxBounds] = None,
      control: CMAESControl = CMAESControl(),
      restarts: CMAESRestarts = CMAESRestarts()
  ): Either[CMAESError, CMAESResult] =
    if control.traceCapacity < 0 || restarts.maxRestarts < 0 || restarts.populationMultiplier < 2 then
      return Left(CMAESError.InvalidConfiguration("invalid trace capacity or restart policy"))
    var session = start(initialMean, config, bounds) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    val seedStream = new CMAESRandom(config.seed)
    var runs = Vector.empty[CMAESRun]
    var trace = Vector.empty[CMAESProgress]
    var best = Option.empty[CMAESPoint]
    var total = CMAESWork(0, 0, 0L, 0L)
    var restart = 0
    var finished = false
    def errorMessage(error: Throwable): String = Option(error.getMessage).getOrElse(error.getClass.getName)
    while !finished do
      while session.stopReason.isEmpty do
        try if control.cancelled() then session.halt(CMAESStop.Cancelled)
        catch case NonFatal(error) => session.halt(CMAESStop.ControlFailure(errorMessage(error)))
        if session.stopReason.isEmpty then
          session.ask() match
            case Left(error)        => return Left(error)
            case Right(None)        => ()
            case Right(Some(batch)) =>
              val values = Vector.newBuilder[Double]
              var index = 0
              var aborted = false
              while index < batch.points.size && !aborted do
                var stop = Option.empty[CMAESStop]
                try if control.cancelled() then stop = Some(CMAESStop.Cancelled)
                catch case NonFatal(error) => stop = Some(CMAESStop.ControlFailure(errorMessage(error)))
                stop match
                  case Some(reason) =>
                    session.abort(batch, values.result(), failedAttempt = false, reason)
                    aborted = true
                  case None =>
                    try
                      val value = objective(batch.points(index).copy)
                      if !value.isFinite then
                        session.abort(
                          batch,
                          values.result(),
                          failedAttempt = true,
                          CMAESStop.ObjectiveFailure("objective returned a non-finite value")
                        )
                        aborted = true
                      else values += value
                    catch
                      case NonFatal(error) =>
                        session.abort(
                          batch,
                          values.result(),
                          failedAttempt = true,
                          CMAESStop.ObjectiveFailure(errorMessage(error))
                        )
                        aborted = true
                    index += 1
              if !aborted then
                session.tell(batch, values.result()) match
                  case Left(error) => return Left(error)
                  case Right(_)    => ()
        val currentBest = session.bestPoint
        if currentBest.exists(p => best.forall(_.value > p.value)) then best = currentBest
        if control.progress.nonEmpty || control.traceCapacity > 0 then
          val progress = CMAESProgress(
            total.generations + session.work.generations,
            total.evaluations + session.work.evaluations,
            best.map(_.value),
            session.stepSize,
            restart
          )
          if control.traceCapacity > 0 then trace = (trace :+ progress).takeRight(control.traceCapacity)
          try if control.progress.exists(callback => !callback(progress)) then session.halt(CMAESStop.Cancelled)
          catch case NonFatal(error) => session.halt(CMAESStop.ControlFailure(errorMessage(error)))
      val reason = session.stopReason.get
      val work = session.work
      total = CMAESWork(
        total.evaluations + work.evaluations,
        total.generations + work.generations,
        total.sampled + work.sampled,
        total.rejected + work.rejected
      )
      runs :+= CMAESRun(session.config.seed, session.populationSize, work, reason)
      val mayRestart = reason match
        case CMAESStop.StepTolerance | CMAESStop.FitnessStagnation | CMAESStop.ConditionLimit |
            CMAESStop.NoRepresentableStep | CMAESStop.SamplingLimit =>
          true
        case _ => false
      val nextPopulation = session.populationSize.toLong * restarts.populationMultiplier
      if mayRestart && restart < restarts.maxRestarts && total.evaluations < config.maxEvaluations &&
        total.generations < config.maxGenerations && nextPopulation <= Int.MaxValue.toLong / initialMean.length
      then
        restart += 1
        val next = config.copy(
          seed = seedStream.nextLong(),
          populationSize = nextPopulation.toInt,
          maxEvaluations = config.maxEvaluations - total.evaluations,
          maxGenerations = config.maxGenerations - total.generations
        )
        session = start(initialMean, next, bounds) match
          case Left(error)  => return Left(error)
          case Right(value) => value
      else finished = true
    Right(CMAESResult(best, session.snapshot.distribution, total, session.stopReason.get, runs, trace))

/** One outstanding generation at a time. Invalid tell/cancel calls leave it available for correction. A final
  * budget-limited batch can contain fewer than populationSize points; its values are recorded without adaptation.
  */
final class CMAESSession private[optim] (
    initial: DVec,
    lower: Array[Double],
    upper: Array[Double],
    free: Array[Int],
    val config: CMAESConfig,
    val populationSize: Int
):
  private val n = free.length
  private val random = new CMAESRandom(config.seed)
  private var mean = Array.tabulate(initial.length)(initial(_))
  private var sigma = config.initialStepSize
  private var covariance = Array.tabulate(n * n)(i => if i / n == i % n then 1.0 else 0.0)
  private var basis = covariance.clone()
  private var axes = Array.fill(n)(1.0)
  private var pc = Array.fill(n)(0.0)
  private var ps = Array.fill(n)(0.0)
  private var evaluations = 0
  private var generations = 0
  private var sampled = 0L
  private var rejected = 0L
  private var stagnant = 0
  private var best = Option.empty[CMAESPoint]
  private var stopped: Option[CMAESStop] = if config.maxGenerations == 0 then Some(CMAESStop.GenerationLimit) else None
  private var pending = Option.empty[CMAESBatch]
  private var steps = Array.empty[Array[Double]]
  private val mu = populationSize / 2
  private val weights = Array.tabulate(mu)(i => math.log(mu + 0.5) - math.log(i + 1.0))
  private val weightSum = weights.sum
  private var wi = 0
  while wi < mu do
    weights(wi) /= weightSum
    wi += 1
  private val mueff = 1.0 / weights.iterator.map(w => w * w).sum
  private val dimension = math.max(n, 1).toDouble
  private val cc = (4.0 + mueff / dimension) / (dimension + 4.0 + 2.0 * mueff / dimension)
  private val cs = (mueff + 2.0) / (dimension + mueff + 5.0)
  private val c1 = 2.0 / (math.pow(dimension + 1.3, 2.0) + mueff)
  private val cmu = math.min(1.0 - c1, 2.0 * (mueff - 2.0 + 1.0 / mueff) / (math.pow(dimension + 2.0, 2.0) + mueff))
  private val damping = 1.0 + 2.0 * math.max(0.0, math.sqrt((mueff - 1.0) / (dimension + 1.0)) - 1.0) + cs
  private val expectedNorm =
    math.sqrt(dimension) * (1.0 - 1.0 / (4.0 * dimension) + 1.0 / (21.0 * dimension * dimension))
  private val refreshPeriod =
    if config.eigenUpdatePeriod > 0 then config.eigenUpdatePeriod
    else math.max(1, (1.0 / ((c1 + cmu) * dimension * 10.0)).toInt)

  def stopReason: Option[CMAESStop] = stopped
  def work: CMAESWork = CMAESWork(evaluations, generations, sampled, rejected)
  def stepSize: Double = sigma
  def bestPoint: Option[CMAESPoint] = best.map(p => p.copy(coordinates = p.coordinates.copy))
  def snapshot: CMAESSnapshot = CMAESSnapshot(
    bestPoint,
    CMAESDistribution(
      DVec.fromArray(mean.clone()),
      sigma,
      DMat.tabulate(n, n)((i, j) => covariance(i * n + j)),
      free.toVector
    ),
    work,
    stopped
  )

  /** Returns None when stopped. Rejected bound samples never consume objective evaluations. */
  def ask(): Either[CMAESError, Option[CMAESBatch]] =
    if pending.nonEmpty then
      Left(CMAESError.ProtocolViolation("tell or cancel the outstanding batch before asking again"))
    else if stopped.nonEmpty then Right(None)
    else
      val count = if n == 0 then 1 else math.min(populationSize, config.maxEvaluations - evaluations)
      val candidates = Vector.newBuilder[DVec]
      val proposedSteps = new Array[Array[Double]](count)
      var index = 0
      var attempts = 0
      var anyMoved = false
      while index < count && stopped.isEmpty do
        if n > 0 && attempts >= config.maxSamplingAttempts then stopped = Some(CMAESStop.SamplingLimit)
        else
          attempts += 1
          sampled += 1
          val z = Array.fill(n)(random.gaussian())
          val y = new Array[Double](n)
          val x = mean.clone()
          var i = 0
          var feasible = true
          var finite = true
          while i < n do
            var sum = 0.0
            var j = 0
            while j < n do
              sum += basis(i * n + j) * axes(j) * z(j)
              j += 1
            val coordinate = free(i)
            x(coordinate) = mean(coordinate) + sigma * sum
            y(i) = (x(coordinate) - mean(coordinate)) / sigma
            finite = finite && x(coordinate).isFinite && y(i).isFinite
            feasible = feasible && x(coordinate) >= lower(coordinate) && x(coordinate) <= upper(coordinate)
            i += 1
          if !finite then stopped = Some(CMAESStop.NumericalFailure("sample is not representable"))
          else if !feasible then rejected += 1
          else
            anyMoved = anyMoved || y.exists(_ != 0.0)
            candidates += DVec.fromArray(x)
            proposedSteps(index) = y
            index += 1
      if stopped.isEmpty && n > 0 && !anyMoved then stopped = Some(CMAESStop.NoRepresentableStep)
      if stopped.nonEmpty then Right(None)
      else
        val batch = new CMAESBatch(candidates.result())
        pending = Some(batch)
        steps = proposedSteps
        Right(Some(batch))

  def tell(batch: CMAESBatch, fitness: IndexedSeq[Double]): Either[CMAESError, Unit] =
    validateBatch(batch, fitness, partial = false).map: _ =>
      val previous = best.map(_.value)
      receive(batch, fitness)
      pending = None
      val complete = n > 0 && fitness.size == populationSize
      if complete then generations += 1
      if config.targetValue.exists(t => best.exists(_.value <= t)) then stopped = Some(CMAESStop.TargetReached)
      else if n == 0 then stopped = Some(CMAESStop.AllCoordinatesFixed)
      else if complete then
        stagnant = if best.exists(p => previous.forall(_ > p.value)) then 0 else stagnant + 1
        update(fitness)
      if stopped.isEmpty && evaluations >= config.maxEvaluations then stopped = Some(CMAESStop.EvaluationLimit)
      if stopped.isEmpty && generations >= config.maxGenerations then stopped = Some(CMAESStop.GenerationLimit)
      steps = Array.empty
      ()

  /** Stop with finite scores for an evaluated prefix; an incomplete population is never used to update covariance. */
  def cancel(
      batch: CMAESBatch,
      completedFitness: IndexedSeq[Double] = Vector.empty
  ): Either[CMAESError, Unit] =
    validateBatch(batch, completedFitness, partial = true).map: _ =>
      abort(batch, completedFitness, failedAttempt = false, CMAESStop.Cancelled)
      ()

  private def validateBatch(
      batch: CMAESBatch,
      fitness: IndexedSeq[Double],
      partial: Boolean
  ): Either[CMAESError, Unit] =
    if !pending.exists(_ eq batch) then
      Left(CMAESError.ProtocolViolation("batch is foreign, stale, or already consumed"))
    else if (if partial then fitness.size > batch.points.size else fitness.size != batch.points.size) then
      Left(
        CMAESError.InvalidFitness("provide one finite score per candidate in batch order (or a prefix when cancelling)")
      )
    else if fitness.exists(!_.isFinite) then Left(CMAESError.InvalidFitness("scores must be finite"))
    else Right(())

  private def receive(batch: CMAESBatch, fitness: IndexedSeq[Double]): Unit =
    var i = 0
    while i < fitness.size do
      if best.forall(_.value > fitness(i)) then best = Some(CMAESPoint(batch.points(i).copy, fitness(i)))
      i += 1
    evaluations += fitness.size

  private[optim] def halt(reason: CMAESStop): Unit =
    stopped = Some(reason)

  private[optim] def abort(
      batch: CMAESBatch,
      fitness: IndexedSeq[Double],
      failedAttempt: Boolean,
      reason: CMAESStop
  ): Unit =
    receive(batch, fitness)
    if failedAttempt then evaluations += 1
    pending = None
    steps = Array.empty
    stopped = Some(reason)

  private def update(fitness: IndexedSeq[Double]): Unit =
    val order = fitness.indices.sortBy(i => (fitness(i), i))
    val yw = new Array[Double](n)
    var k = 0
    while k < mu do
      var i = 0
      while i < n do
        yw(i) += weights(k) * steps(order(k))(i)
        i += 1
      k += 1
    val whitened = new Array[Double](n)
    var j = 0
    while j < n do
      var projection = 0.0
      var i = 0
      while i < n do
        projection += basis(i * n + j) * yw(i)
        i += 1
      projection /= axes(j)
      i = 0
      while i < n do
        whitened(i) += basis(i * n + j) * projection
        i += 1
      j += 1
    val nextPs = Array.tabulate(n)(i => (1.0 - cs) * ps(i) + math.sqrt(cs * (2.0 - cs) * mueff) * whitened(i))
    val psNorm = OptimNumerics.norm(nextPs)
    val correction = math.sqrt(-math.expm1(2.0 * generations * math.log1p(-cs)))
    val hsig = psNorm / correction / expectedNorm < 1.4 + 2.0 / (dimension + 1.0)
    val nextPc = Array.tabulate(n)(i =>
      (1.0 - cc) * pc(i) +
        (if hsig then math.sqrt(cc * (2.0 - cc) * mueff) * yw(i) else 0.0)
    )
    val nextMean = mean.clone()
    var i = 0
    while i < n do
      nextMean(free(i)) += sigma * yw(i)
      i += 1
    val nextCovariance = new Array[Double](n * n)
    val decay = 1.0 - c1 - cmu + (if hsig then 0.0 else c1 * cc * (2.0 - cc))
    i = 0
    while i < n do
      j = 0
      while j <= i do
        var rankMu = 0.0
        k = 0
        while k < mu do
          rankMu += weights(k) * steps(order(k))(i) * steps(order(k))(j)
          k += 1
        val value = decay * covariance(i * n + j) + c1 * nextPc(i) * nextPc(j) + cmu * rankMu
        nextCovariance(i * n + j) = value
        nextCovariance(j * n + i) = value
        j += 1
      i += 1
    val nextSigma = sigma * math.exp(cs / damping * (psNorm / expectedNorm - 1.0))
    if !nextSigma.isFinite || nextSigma <= 0.0 ||
      !OptimNumerics.allFinite(nextMean) || !OptimNumerics.allFinite(nextCovariance) ||
      !OptimNumerics.allFinite(nextPs) || !OptimNumerics.allFinite(nextPc)
    then stopped = Some(CMAESStop.NumericalFailure("non-finite adaptation state or vanished step size"))
    else if free.exists(i => nextMean(i) < lower(i) || nextMean(i) > upper(i)) then
      stopped = Some(CMAESStop.NumericalFailure("rounded recombination mean left the box"))
    else
      mean = nextMean
      covariance = nextCovariance
      sigma = nextSigma
      ps = nextPs
      pc = nextPc
      if generations % refreshPeriod == 0 then refresh()
      var largestVariance = 0.0
      i = 0
      while i < n do
        largestVariance = math.max(largestVariance, covariance(i * n + i))
        i += 1
      if stopped.isEmpty && config.stepTolerance > 0.0 && sigma * math.sqrt(largestVariance) <= config.stepTolerance
      then stopped = Some(CMAESStop.StepTolerance)
      if stopped.isEmpty && config.stagnationGenerations > 0 && stagnant >= config.stagnationGenerations then
        stopped = Some(CMAESStop.FitnessStagnation)

  private def refresh(): Unit =
    val matrix = DMat.tabulate(n, n)((i, j) => covariance(i * n + j))
    Eigen.eigSymmetric(matrix, EigenSelection.All) match
      case Left(error)          => stopped = Some(CMAESStop.NumericalFailure(error.toString))
      case Right(decomposition) =>
        val values = decomposition.eigenvalues
        if (0 until n).exists(i => !values(i).isFinite || values(i) <= 0.0) then
          stopped = Some(CMAESStop.NumericalFailure("covariance lost positive definiteness"))
        else if values(n - 1) / values(0) > config.maxCondition then stopped = Some(CMAESStop.ConditionLimit)
        else
          axes = Array.tabulate(n)(i => math.sqrt(values(i)))
          basis = Array.tabulate(n * n)(i => decomposition.eigenvectors(i / n, i % n))

/** SplitMix64 and Box-Muller. Reproducible within a runtime; transcendental/eigensolver rounding may differ across
  * platforms.
  */
private[optim] final class CMAESRandom(seed: Long):
  private var state = seed
  private var spare = Option.empty[Double]
  def nextLong(): Long =
    state += 0x9e3779b97f4a7c15L
    var z = state
    z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L
    z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL
    z ^ (z >>> 31)
  private def uniform(): Double = ((nextLong() >>> 12).toDouble + 0.5) / 4503599627370496.0
  def gaussian(): Double = spare match
    case Some(value) => spare = None; value
    case None        =>
      val radius = math.sqrt(-2.0 * math.log(uniform()))
      val angle = 2.0 * math.Pi * uniform()
      spare = Some(radius * math.sin(angle))
      radius * math.cos(angle)
