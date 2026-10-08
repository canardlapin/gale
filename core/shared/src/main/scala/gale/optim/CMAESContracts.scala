package gale.optim

import gale.linalg.{DMat, DVec}

/** Invalid configuration or ask/tell usage. Invalid tell calls do not consume a batch. */
enum CMAESError:
  case InvalidConfiguration(detail: String)
  case ProtocolViolation(detail: String)
  case InvalidFitness(detail: String)

/** A termination reason, never a stationarity or global-optimality certificate. */
enum CMAESStop:
  case TargetReached, EvaluationLimit, GenerationLimit, StepTolerance, FitnessStagnation
  case ConditionLimit, NoRepresentableStep, SamplingLimit, AllCoordinatesFixed, Cancelled
  case ObjectiveFailure(detail: String)
  case ControlFailure(detail: String)
  case NumericalFailure(detail: String)

/** Positive-weight rank-one/rank-mu CMA-ES with cumulative step-size adaptation. `populationSize = 0` selects 4 +
  * floor(3 log(number of free coordinates)). `maxSamplingAttempts` caps whole-vector draws per generation, including
  * rejected draws. `eigenUpdatePeriod = 0` selects the standard amortized refresh interval. `stagnationGenerations = 0`
  * disables the no-improvement stop; zero step tolerance disables that stop.
  */
final case class CMAESConfig(
    initialStepSize: Double = 0.5,
    populationSize: Int = 0,
    seed: Long = 0L,
    maxEvaluations: Int = 100000,
    maxGenerations: Int = 10000,
    targetValue: Option[Double] = None,
    stepTolerance: Double = 1e-12,
    stagnationGenerations: Int = 100,
    maxCondition: Double = 1e14,
    maxSamplingAttempts: Int = 100000,
    eigenUpdatePeriod: Int = 0
)

/** Restarts begin at the original mean and step size, with identity covariance and fresh derived seeds. Population size
  * grows geometrically. Evaluation and generation limits remain global.
  */
final case class CMAESRestarts(maxRestarts: Int = 0, populationMultiplier: Int = 2)

final case class CMAESProgress(
    generation: Int,
    evaluations: Int,
    bestValue: Option[Double],
    stepSize: Double,
    restart: Int
)

final case class CMAESControl(
    cancelled: () => Boolean = () => false,
    progress: Option[CMAESProgress => Boolean] = None,
    traceCapacity: Int = 0
)

final case class CMAESPoint(coordinates: DVec, value: Double)

/** The dimensionless adapted, unconstrained covariance model on `freeCoordinates`. The sampling factor can lag this
  * model until its next eigendecomposition; bounds condition accepted samples.
  */
final case class CMAESDistribution(mean: DVec, stepSize: Double, covariance: DMat, freeCoordinates: Vector[Int])

/** Minimize counts attempted calls; ask/tell counts finite scores in accepted receipts. Sampling rejections do not call
  * the objective. Generations count complete received populations (a final partial population is not a generation).
  */
final case class CMAESWork(evaluations: Int, generations: Int, sampled: Long, rejected: Long)

final case class CMAESSnapshot(
    best: Option[CMAESPoint],
    distribution: CMAESDistribution,
    work: CMAESWork,
    stopReason: Option[CMAESStop]
)

final case class CMAESRun(seed: Long, populationSize: Int, work: CMAESWork, stopReason: CMAESStop)

final case class CMAESResult(
    best: Option[CMAESPoint],
    distribution: CMAESDistribution,
    work: CMAESWork,
    stopReason: CMAESStop,
    runs: Vector[CMAESRun],
    trace: Vector[CMAESProgress]
)

/** Owned immutable candidates, in the order expected by tell. A batch belongs to one session and can be consumed once.
  */
final class CMAESBatch private[optim] (val points: Vector[DVec])
