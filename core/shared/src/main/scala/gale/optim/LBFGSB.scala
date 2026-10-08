package gale.optim

import gale.linalg.DMat

final case class LBFGSBConfig(
    maxIterations: Int = 1000,
    tolerance: FirstOrderTolerance = FirstOrderTolerance.strict,
    historySize: Int = 10,
    c1: Double = 1e-4,
    c2: Double = 0.9,
    maxLineSearch: Int = 40,
    control: SolverControl = SolverControl()
)

/** L-BFGS-B with a generalized Cauchy point and an iterative free-variable subspace solve. Storage is O(history * n).
  * Rebuilding the limited-memory Hessian factors costs O(history² * n); no dense n-by-n Hessian is constructed. Initial
  * points must be feasible. Convergence checks the unit projected-gradient mapping at the returned point.
  */
object LBFGSB:
  import OptimNumerics.*
  private[optim] final case class Pair(s: Array[Double], y: Array[Double])
  private[optim] final case class Model(theta: Double, factors: Vector[(Array[Double], Double)]):
    def apply(v: Array[Double]): Array[Double] =
      val out = new Array[Double](v.length)
      var i = 0
      while i < out.length do
        out(i) = v(i) * theta
        i += 1
      var k = 0
      while k < factors.size do
        val (a, sign) = factors(k)
        val weight = sign * dot(a, v)
        i = 0
        while i < out.length do
          out(i) += weight * a(i)
          i += 1
        k += 1
      out

    /** The CG direction is zero on active coordinates. Skip their zero products and reuse scratch; free indices retain
      * ascending order, so nonzero summands have the same accumulation order.
      */
    def restricted(v: Array[Double], free: Array[Int], out: Array[Double]): Unit =
      var j = 0
      while j < free.length do
        val i = free(j)
        out(i) = theta * v(i)
        j += 1
      var k = 0
      while k < factors.size do
        val (a, sign) = factors(k)
        var product = 0.0
        j = 0
        while j < free.length do
          val i = free(j)
          product += a(i) * v(i)
          j += 1
        val weight = sign * product
        j = 0
        while j < free.length do
          val i = free(j)
          out(i) += weight * a(i)
          j += 1
        k += 1

  private[optim] def model(history: Vector[Pair], initialTheta: Double = 1.0): Model =
    val theta = history.lastOption.map(p => dot(p.y, p.y) / dot(p.s, p.y)).getOrElse(initialTheta)
    if !theta.isFinite || theta <= 0.0 then return Model(1.0, Vector.empty)
    var result = Model(theta, Vector.empty)
    history.foreach: pair =>
      val bs = result(pair.s)
      val sbs = dot(pair.s, bs)
      val sy = dot(pair.s, pair.y)
      if sbs.isFinite && sbs > 0.0 && sy.isFinite && sy > 0.0 then
        val rootB = Math.sqrt(sbs)
        val rootY = Math.sqrt(sy)
        val v = new Array[Double](bs.length)
        val u = new Array[Double](bs.length)
        var i = 0
        while i < bs.length do
          v(i) = bs(i) / rootB
          u(i) = pair.y(i) / rootY
          i += 1
        if allFinite(v) && allFinite(u) then result = Model(theta, result.factors ++ Vector(v -> -1.0, u -> 1.0))
    result

  private[optim] def projected(
      x: Array[Double],
      g: Array[Double],
      lower: Array[Double],
      upper: Array[Double]
  ): Array[Double] =
    tabulate(x.length): i =>
      // Algebraically x - clamp(x-g); avoids rounded-away subtraction at translated points.
      if g(i) > 0.0 then Math.min(g(i), x(i) - lower(i))
      else Math.max(g(i), x(i) - upper(i))

  private def infinity(v: Array[Double]): Double =
    var largest = 0.0
    var i = 0
    while i < v.length do
      largest = Math.max(largest, Math.abs(v(i)))
      i += 1
    largest

  /** First minimum of the quadratic along the projected negative-gradient path. Each breakpoint updates O(history)
    * scalars; coordinates are materialized only at the final path time.
    */
  private[optim] def cauchy(
      x: Array[Double],
      g: Array[Double],
      l: Array[Double],
      u: Array[Double],
      b: Model
  ): Array[Double] =
    val n = x.length
    val direction = tabulate(n)(i =>
      if l(i) == u(i) || (x(i) == l(i) && g(i) >= 0.0) ||
        (x(i) == u(i) && g(i) <= 0.0)
      then 0.0
      else -g(i)
    )
    val breaks = (0 until n)
      .flatMap: i =>
        val t =
          if direction(i) < 0.0 then (l(i) - x(i)) / direction(i)
          else if direction(i) > 0.0 then (u(i) - x(i)) / direction(i)
          else Double.PositiveInfinity
        if t.isFinite then Some((Math.max(0.0, t), i)) else None
      .sortBy(_._1)
    val fixed = Array.fill(n)(false)
    var a = dot(direction, direction)
    var h = dot(g, direction)
    var dz = 0.0
    val projections = b.factors.map((v, _) => dot(v, direction)).toArray
    val displacementProjections = new Array[Double](b.factors.size)
    var time = 0.0
    var index = 0
    var finished = false
    while !finished do
      var slope = h + b.theta * dz
      var curvature = b.theta * a
      var k = 0
      while k < b.factors.size do
        val sign = b.factors(k)._2
        slope += sign * projections(k) * displacementProjections(k)
        curvature += sign * projections(k) * projections(k)
        k += 1
      val end = if index < breaks.size then breaks(index)._1 else Double.PositiveInfinity
      if !slope.isFinite || !curvature.isFinite || curvature <= 0.0 || slope >= 0.0 then finished = true
      else
        val minimizer = -slope / curvature
        val interval = end - time
        if minimizer < interval || !end.isFinite then
          time += minimizer
          finished = true
        else
          k = 0
          while k < b.factors.size do
            displacementProjections(k) += interval * projections(k)
            k += 1
          dz += interval * a
          time = end
          val previousA = a
          while index < breaks.size && breaks(index)._1 == end do
            val i = breaks(index)._2
            val d = direction(i)
            val z = (if d < 0.0 then l(i) else u(i)) - x(i)
            h -= g(i) * d
            dz -= d * z
            a = Math.max(0.0, a - d * d)
            k = 0
            while k < b.factors.size do
              projections(k) -= d * b.factors(k)._1(i)
              k += 1
            fixed(i) = true
            direction(i) = 0.0
            index += 1
          // Recondition when subtracting a large frozen contribution erased smaller free directions.
          if a <= 1e-12 * previousA then
            a = dot(direction, direction)
            val displacement = tabulate(n)(i =>
              if fixed(i) then (if g(i) > 0.0 then l(i) else u(i)) - x(i)
              else time * direction(i)
            )
            k = 0
            while k < b.factors.size do
              projections(k) = dot(b.factors(k)._1, direction)
              displacementProjections(k) = dot(b.factors(k)._1, displacement)
              k += 1
          // Still-free coordinates have d=-g and displacement=time*d along the projected path.
          h = -a
          dz = time * a
          if a == 0.0 then finished = true
    tabulate(n): i =>
      if fixed(i) then (if g(i) > 0.0 then l(i) else u(i))
      else Math.max(l(i), Math.min(u(i), x(i) + time * direction(i)))

  private def subspace(
      x: Array[Double],
      g: Array[Double],
      c: Array[Double],
      l: Array[Double],
      u: Array[Double],
      b: Model
  ): Array[Double] =
    val free = Array.tabulate(x.length)(i => c(i) > l(i) && c(i) < u(i))
    val indices = free.indices.filter(free(_)).toArray
    val z = tabulate(x.length)(i => c(i) - x(i))
    val bz = b(z)
    val rhs = tabulate(x.length)(i => if free(i) then -(g(i) + bz(i)) else 0.0)
    val v = new Array[Double](x.length)
    val residual = rhs.clone()
    val direction = residual.clone()
    var rr = dot(residual, residual)
    val target = 1e-20 * rr
    var iteration = 0
    val limit = Math.min(x.length, 2 * b.factors.size + 5)
    val product = new Array[Double](x.length)
    while rr > target && rr > 0.0 && iteration < limit do
      b.restricted(direction, indices, product)
      val denominator = dot(direction, product)
      if !denominator.isFinite || denominator <= 0.0 || !rr.isFinite then return c
      val alpha = rr / denominator
      var i = 0
      while i < x.length do
        v(i) += alpha * direction(i)
        residual(i) -= alpha * product(i)
        i += 1
      val next = dot(residual, residual)
      val beta = next / rr
      i = 0
      while i < x.length do
        direction(i) = residual(i) + beta * direction(i)
        i += 1
      rr = next
      iteration += 1
    var fraction = 1.0
    var i = 0
    while i < x.length do
      if v(i) > 0.0 then fraction = Math.min(fraction, (u(i) - c(i)) / v(i))
      else if v(i) < 0.0 then fraction = Math.min(fraction, (l(i) - c(i)) / v(i))
      i += 1
    val candidate = tabulate(x.length)(i => Math.max(l(i), Math.min(u(i), c(i) + fraction * v(i))))
    val nextZ = tabulate(x.length)(i => candidate(i) - x(i))
    val oldModel = dot(g, z) + 0.5 * dot(z, bz)
    val newModel = dot(g, nextZ) + 0.5 * dot(nextZ, b(nextZ))
    if allFinite(candidate) && newModel.isFinite && newModel <= oldModel then candidate else c

  def minimize(
      objective: DifferentiableObjective,
      bounds: MatrixBoxBounds,
      initial: DMat,
      config: LBFGSBConfig = LBFGSBConfig()
  ): Either[FirstOrderError, FirstOrderSolution] =
    minimizeWithScale(
      objective,
      bounds,
      initial,
      config,
      None,
      interpolate = false,
      approximateWolfe = false,
      onApproximateAcceptance = StrongWolfe.ignoreApproximateAcceptance
    )

  private[optim] def minimizeWithScale(
      objective: DifferentiableObjective,
      bounds: MatrixBoxBounds,
      initial: DMat,
      config: LBFGSBConfig,
      curvatureScale: Option[CurvatureScale],
      interpolate: Boolean,
      approximateWolfe: Boolean,
      onApproximateAcceptance: () => Unit
  ): Either[FirstOrderError, FirstOrderSolution] =
    if config.maxIterations < 0 || config.historySize <= 0 || config.maxLineSearch <= 0 ||
      !config.c1.isFinite || !config.c2.isFinite || config.c1 <= 0.0 || config.c1 >= config.c2 || config.c2 >= 1.0
    then Left(FirstOrderError.InvalidConfiguration("L-BFGS-B needs valid limits and 0 < c1 < c2 < 1"))
    else
      for
        _ <- validate(initial, objective.variableRows)
        _ <-
          if bounds.contains(initial) then Right(())
          else Left(FirstOrderError.InvalidConfiguration("initial point must satisfy the box bounds"))
        _ <- config.control.validate
        solution <- solve(
          objective,
          bounds,
          initial,
          config,
          curvatureScale,
          interpolate,
          approximateWolfe,
          onApproximateAcceptance
        )
      yield solution

  private def solve(
      objective: DifferentiableObjective,
      bounds: MatrixBoxBounds,
      initial: DMat,
      config: LBFGSBConfig,
      curvatureScale: Option[CurvatureScale],
      interpolate: Boolean,
      approximateWolfe: Boolean,
      onApproximateAcceptance: () => Unit
  ): Either[FirstOrderError, FirstOrderSolution] =
    val execution = new OptimizationExecution(config.control)
    val lower = array(bounds.lower)
    val upper = array(bounds.upper)
    var point = matrix(array(initial), initial.rows, initial.cols)
    var evaluation = execution.evaluate(objective, point) match
      case Left(error)  => return Left(error)
      case Right(value) => value
    var x = array(point)
    var gradient = array(evaluation.gradient)
    val reference = Math.max(1.0, infinity(projected(x, gradient, lower, upper)))
    if !config.tolerance.threshold(reference).isFinite then
      return Left(FirstOrderError.NumericalFailure("box tolerance overflow"))
    val settings = FirstOrderSettings(
      FirstOrderMethod.LBFGSB,
      config.maxIterations,
      config.tolerance,
      1.0,
      0.0,
      primalResidualScale = reference,
      algorithm =
        AlgorithmSettings.BoundedLimitedMemory(config.historySize, config.c1, config.c2, config.maxLineSearch),
      maxEvaluations = config.control.maxEvaluations
    )
    var history = Vector.empty[Pair]
    var iteration = 0
    var step = 0.0
    while true do
      val pg = projected(x, gradient, lower, upper)
      val residual = infinity(pg)
      val status =
        if residual <= config.tolerance.threshold(reference) then FirstOrderStoppingStatus.Converged
        else FirstOrderStoppingStatus.IterationLimit
      val latest = result(point, evaluation.value, residual, iteration, status, settings, execution, step)
      execution.record(latest)
      if !execution.continue then
        return execution.finish(Left(FirstOrderError.ExecutionStopped(FirstOrderStoppingStatus.Cancelled)))
      if status == FirstOrderStoppingStatus.Converged || iteration >= config.maxIterations then
        return execution.finish(Right(latest))
      val b = model(history, curvatureScale.fold(1.0)(s => 1.0 / s.value))
      val c = cauchy(x, gradient, lower, upper, b)
      val candidate = subspace(x, gradient, c, lower, upper, b)
      var direction = tabulate(x.length)(i => candidate(i) - x(i))
      var slope = dot(gradient, direction)
      if !allFinite(direction) || !slope.isFinite || slope >= 0.0 then
        history = Vector.empty
        curvatureScale.foreach(_.value = 1.0)
        direction = tabulate(pg.length)(i => -pg(i))
        slope = dot(gradient, direction)
      if infinity(direction) == 0.0 || !slope.isFinite || slope >= 0.0 then
        return execution.finish(Right(latest.copy(status = FirstOrderStoppingStatus.NumericalStagnation)))
      var feasibleStep = Double.PositiveInfinity
      var i = 0
      while i < x.length do
        if direction(i) > 0.0 then feasibleStep = Math.min(feasibleStep, (upper(i) - x(i)) / direction(i))
        else if direction(i) < 0.0 then feasibleStep = Math.min(feasibleStep, (lower(i) - x(i)) / direction(i))
        i += 1
      val maximum = Math.min(1e20, feasibleStep)
      if maximum <= 0.0 then
        return execution.finish(Right(latest.copy(status = FirstOrderStoppingStatus.NumericalStagnation)))
      StrongWolfe.search(
        objective,
        point,
        evaluation,
        matrix(direction, point.rows, point.cols),
        execution,
        config.c1,
        config.c2,
        config.maxLineSearch,
        maximum,
        feasibleBoundary = feasibleStep.isFinite && feasibleStep <= 1e20,
        bounds = Some(bounds),
        interpolate = interpolate,
        approximateWolfe = approximateWolfe,
        onApproximateAcceptance = onApproximateAcceptance
      ) match
        case Left(error)  => return execution.finish(Left(error))
        case Right(found) =>
          val s = difference(found.point, point)
          val y = difference(found.evaluation.gradient, evaluation.gradient)
          val sy = dot(s, y)
          val lengths = norm(s) * norm(y)
          if sy.isFinite && lengths.isFinite && sy > 1e-12 * lengths && sy > 0.0 then
            history = (history :+ Pair(s, y)).takeRight(config.historySize)
            curvatureScale.foreach(_.update(sy, dot(y, y)))
          point = found.point
          evaluation = found.evaluation
          x = array(point)
          gradient = array(evaluation.gradient)
          step = found.step
          iteration += 1
    execution.finish(Left(FirstOrderError.NumericalFailure("unreachable L-BFGS-B state")))
