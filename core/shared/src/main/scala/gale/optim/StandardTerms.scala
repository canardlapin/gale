package gale.optim

import gale.linalg.DMat

/** Small, explicitly scoped proximal terms and convex projection sets. */
object ProximalTerms:
  def zero(rows: Int): Either[FirstOrderError, ProximalTerm] =
    checkedRows(rows).map(r => SeparableTerm(r, _ => 0.0, (x, _) => x, (_, _) => 0.0))

  def l1(rows: Int, weight: Double): Either[FirstOrderError, ProximalTerm] =
    for _ <- checkedRows(rows); _ <- checkedWeight(weight)
    yield SeparableTerm(
      rows,
      x => weight * Math.abs(x),
      (x, step) => Math.signum(x) * Math.max(0.0, Math.abs(x) - step * weight),
      (x, _) => Math.max(-weight, Math.min(weight, x))
    )

  def weightedL1(weights: DMat): Either[FirstOrderError, ProximalTerm] =
    finiteMatrix("weighted L1 weights", weights).flatMap(_ =>
      if weights.rows <= 0 || weights.cols <= 0 then
        Left(FirstOrderError.InvalidConfiguration("weighted L1 weights must have positive dimensions"))
      else if hasNegative(weights) then
        Left(FirstOrderError.InvalidConfiguration("weighted L1 weights must be non-negative"))
      else Right(MatrixWeightedL1(weights))
    )

  def squaredL2(rows: Int, weight: Double): Either[FirstOrderError, ProximalTerm] =
    for _ <- checkedRows(rows); _ <- checkedWeight(weight)
    yield SeparableTerm(
      rows,
      x => 0.5 * weight * x * x,
      (x, step) => if (step * weight).isFinite then x / (1.0 + step * weight) else (x / weight) / step,
      (x, step) =>
        if weight == 0.0 then 0.0
        else if (step / weight).isFinite then x / (1.0 + step / weight)
        else (x * weight) / step
    )

  def conjugate(term: ProximalTerm): LinearCompositeFunctional = new LinearCompositeFunctional:
    def targetRows: Int = term.variableRows
    def value(at: DMat): Either[FirstOrderError, Double] = term.value(at)
    def proximalConjugate(at: DMat, step: Double): Either[FirstOrderError, DMat] =
      for
        _ <- checkedStep(step)
        _ <- input(at, targetRows)
        result <- term match
          case known: SeparableTerm    => map(at, "conjugate proximal")(x => known.dualProx(x, step))
          case known: MatrixWeightedL1 => known.dual(at)
          case _                       => moreau(term, at, step)
      yield result

  private def moreau(term: ProximalTerm, at: DMat, step: Double): Either[FirstOrderError, DMat] =
    for
      inverseStep <- finiteScalar("Moreau inverse step", 1.0 / step)
      scaled <- map(at, "Moreau scaling")(_ / step)
      primal <- term.proximal(scaled, inverseStep)
      _ <- sameShape("Moreau proximal result", primal, at)
      out <- zip(at, primal, "Moreau reconstruction")((z, p) => z - step * p)
      _ <-
        var unresolved = false
        var r = 0
        while r < at.rows do
          var c = 0
          while c < at.cols do
            val resolution = 8.0 * Math.ulp(Math.max(Math.abs(at(r, c)), Math.abs(step * primal(r, c))))
            if resolution > 1e-12 * Math.max(1.0, Math.abs(out(r, c))) then unresolved = true
            c += 1
          r += 1
        if unresolved then
          Left(
            FirstOrderError.NumericalFailure(
              "Moreau reconstruction loses numerical resolution; supply a direct conjugate proximal"
            )
          )
        else Right(())
    yield out

object ProjectionSets:
  def box(rows: Int, lower: Double, upper: Double): Either[FirstOrderError, ConvexProjectionSet] =
    for _ <- checkedRows(rows); _ <- checkedBounds(lower, upper)
    yield ScalarBox(rows, lower, upper)
  def box(lower: DMat, upper: DMat): Either[FirstOrderError, ConvexProjectionSet] =
    for
      _ <- sameShape("box bounds", lower, upper); _ <- finiteMatrix("box lower bound", lower);
      _ <- finiteMatrix("box upper bound", upper)
      _ <-
        if lower.rows > 0 && lower.cols > 0 && !anyLowerExceeds(lower, upper) then Right(())
        else Left(FirstOrderError.InvalidConfiguration("box bounds must have positive equal shape and lower <= upper"))
    yield MatrixBox(lower, upper)
  def nonnegative(rows: Int): Either[FirstOrderError, ConvexProjectionSet] = checkedRows(rows).map(Nonnegative(_))
  def simplex(rows: Int, mass: Double = 1.0): Either[FirstOrderError, ConvexProjectionSet] =
    for
      _ <- checkedRows(rows);
      _ <-
        if mass.isFinite && mass > 0.0 then Right(())
        else Left(FirstOrderError.InvalidConfiguration("simplex mass must be finite and positive"))
    yield Simplex(rows, mass)

/** A projection set that supplies convex-feasibility diagnostics. */
trait ConvexProjectionSet extends ProjectionSet:
  def violation(at: DMat): Either[FirstOrderError, Double]
  def indicator: ProximalTerm

private abstract class CheckedSet(val variableRows: Int) extends ConvexProjectionSet:
  final def indicator: ProximalTerm = new ProximalTerm:
    def variableRows: Int = CheckedSet.this.variableRows
    def value(at: DMat): Either[FirstOrderError, Double] =
      violation(at).flatMap(v =>
        if v <= feasibilityTolerance(at) then Right(0.0)
        else Left(FirstOrderError.OracleFailure("indicator", "point must be feasible"))
      )
    def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] = for _ <- checkedStep(step); out <- project(at)
    yield out

private final case class ScalarBox(override val variableRows: Int, lower: Double, upper: Double)
    extends CheckedSet(variableRows):
  def project(at: DMat) = for
    _ <- input(at, variableRows); out <- map(at, "box projection")(x => Math.max(lower, Math.min(upper, x)))
  yield out
  def violation(at: DMat) = for _ <- input(at, variableRows) yield maxViolation(at)(x => Math.max(lower - x, x - upper))
private final case class Nonnegative(override val variableRows: Int) extends CheckedSet(variableRows):
  def project(at: DMat) = for
    _ <- input(at, variableRows); out <- map(at, "nonnegative projection")(x => Math.max(0.0, x))
  yield out
  def violation(at: DMat) = for _ <- input(at, variableRows) yield maxViolation(at)(x => Math.max(0.0, -x))
private final case class MatrixBox(lower: DMat, upper: DMat) extends CheckedSet(lower.rows):
  def project(at: DMat) = for
    _ <- sameShape("box input", at, lower); _ <- finiteMatrix("box input", at);
    out <- zip3(at, lower, upper, "box projection")((x, l, u) => Math.max(l, Math.min(u, x)))
  yield out
  def violation(at: DMat) = for _ <- sameShape("box input", at, lower); _ <- finiteMatrix("box input", at)
  yield maxZipViolation(at, lower, upper)
private final case class Simplex(override val variableRows: Int, mass: Double) extends CheckedSet(variableRows):
  def project(at: DMat) = for
    _ <- input(at, variableRows)
    columns = Vector.tabulate(at.cols)(c => simplexColumn(at, c, mass))
    out <- map(DMat.tabulate(at.rows, at.cols)((r, c) => columns(c)(r)), "simplex projection")(identity)
    _ <-
      if simplexViolation(out, mass) <= feasibilityTolerance(out) then Right(())
      else Left(FirstOrderError.NumericalFailure("simplex projection cannot resolve the requested mass"))
  yield out
  def violation(at: DMat) = for _ <- input(at, variableRows) yield simplexViolation(at, mass)

private final case class SeparableTerm(
    variableRows: Int,
    f: Double => Double,
    prox: (Double, Double) => Double,
    dualProx: (Double, Double) => Double
) extends ProximalTerm:
  def value(at: DMat) = input(at, variableRows).flatMap(_ => finiteScalar("term value", sum(at, f)))
  def proximal(at: DMat, step: Double) = for
    _ <- checkedStep(step); _ <- input(at, variableRows); out <- map(at, "proximal result")(x => prox(x, step))
  yield out
private final case class MatrixWeightedL1(weights: DMat) extends ProximalTerm:
  def variableRows = weights.rows
  def dual(at: DMat): Either[FirstOrderError, DMat] = for
    _ <- sameShape("weighted L1 dual input", at, weights);
    out <- zip(at, weights, "weighted L1 conjugate")((x, w) => Math.max(-w, Math.min(w, x)))
  yield out
  def value(at: DMat) = for
    _ <- sameShape("weighted L1 input", at, weights); _ <- finiteMatrix("weighted L1 input", at);
    value <- finiteScalar("weighted L1 value", sumZip(at, weights)((x, w) => w * Math.abs(x)))
  yield value
  def proximal(at: DMat, step: Double) = for
    _ <- checkedStep(step); _ <- sameShape("weighted L1 input", at, weights);
    _ <- finiteMatrix("weighted L1 input", at);
    out <- zip(at, weights, "weighted L1 proximal")((x, w) => Math.signum(x) * Math.max(0.0, Math.abs(x) - step * w))
  yield out

private def checkedRows(rows: Int): Either[FirstOrderError, Int] =
  if rows > 0 then Right(rows) else Left(FirstOrderError.InvalidConfiguration("rows must be positive"))
private def checkedWeight(weight: Double) = if weight.isFinite && weight >= 0.0 then Right(())
else Left(FirstOrderError.InvalidConfiguration("weight must be finite and non-negative"))
private def checkedStep(step: Double) = if step.isFinite && step > 0.0 then Right(())
else Left(FirstOrderError.InvalidConfiguration("step must be finite and positive"))
private def input(at: DMat, rows: Int) =
  if at.rows != rows then Left(FirstOrderError.ShapeMismatch("input", rows, at.rows)) else finiteMatrix("input", at)
private def finiteMatrix(context: String, a: DMat): Either[FirstOrderError, Unit] =
  var i = 0
  while i < a.rows * a.cols do
    val value = a(i / a.cols, i % a.cols)
    if !value.isFinite then return Left(FirstOrderError.NonFiniteValue(context, i, value))
    i += 1
  Right(())
private def sameShape(context: String, a: DMat, b: DMat) = if (a.rows == b.rows && a.cols == b.cols) Right(())
else Left(FirstOrderError.InvalidConfiguration(s"$context requires exact ${b.rows}x${b.cols} shape"))
private def map(a: DMat, context: String)(f: Double => Double) = {
  val o = DMat.tabulate(a.rows, a.cols)((r, c) => f(a(r, c))); finiteMatrix(context, o).map(_ => o)
}
private def zip(a: DMat, b: DMat, context: String)(f: (Double, Double) => Double) = {
  val o = DMat.tabulate(a.rows, a.cols)((r, c) => f(a(r, c), b(r, c))); finiteMatrix(context, o).map(_ => o)
}
private def zip3(a: DMat, b: DMat, c: DMat, context: String)(f: (Double, Double, Double) => Double) = {
  val o = DMat.tabulate(a.rows, a.cols)((r, j) => f(a(r, j), b(r, j), c(r, j))); finiteMatrix(context, o).map(_ => o)
}
private def sum(a: DMat, f: Double => Double) = {
  var s = 0.0; var r = 0; while (r < a.rows) { var c = 0; while (c < a.cols) { s += f(a(r, c)); c += 1 }; r += 1 }; s
}
private def sumZip(a: DMat, b: DMat)(f: (Double, Double) => Double) = {
  var s = 0.0; var r = 0;
  while (r < a.rows) { var c = 0; while (c < a.cols) { s += f(a(r, c), b(r, c)); c += 1 }; r += 1 }; s
}
private def hasNegative(a: DMat) =
  a.rows * a.cols > 0 && (0 until a.rows).exists(r => (0 until a.cols).exists(c => a(r, c) < 0))
private def anyLowerExceeds(a: DMat, b: DMat) =
  (0 until a.rows).exists(r => (0 until a.cols).exists(c => a(r, c) > b(r, c)))
private def maxViolation(a: DMat)(f: Double => Double) = {
  var m = 0.0; var r = 0;
  while (r < a.rows) { var c = 0; while (c < a.cols) { m = Math.max(m, f(a(r, c))); c += 1 }; r += 1 }; m
}
private def maxZipViolation(a: DMat, l: DMat, u: DMat) = {
  var m = 0.0; var r = 0;
  while (r < a.rows) {
    var c = 0; while (c < a.cols) { m = Math.max(m, Math.max(l(r, c) - a(r, c), a(r, c) - u(r, c))); c += 1 }; r += 1
  }; m
}
private def feasibilityTolerance(a: DMat) = 8.0 * Math.ulp(Math.max(1.0, termMaxAbs(a))) * a.rows
private def termMaxAbs(a: DMat) = {
  var m = 0.0; var r = 0;
  while (r < a.rows) { var c = 0; while (c < a.cols) { m = Math.max(m, Math.abs(a(r, c))); c += 1 }; r += 1 }; m
}
private def finiteScalar(context: String, value: Double): Either[FirstOrderError, Double] =
  if (value.isFinite) Right(value) else Left(FirstOrderError.NumericalFailure(s"$context is non-finite"))
private def simplexColumn(a: DMat, col: Int, mass: Double): Array[Double] =
  val maximum = (0 until a.rows).map(a(_, col)).max
  val shifted = Array.tabulate(a.rows)(r => a(r, col) - maximum)
  val sorted = shifted.sorted(using Ordering[Double].reverse)
  var sum = 0.0
  var threshold = -mass
  var i = 0
  var running = true
  while i < sorted.length && running do
    // Stop before summing a very negative inactive tail (possibly -Infinity).
    if sorted(i) <= threshold then running = false
    else
      sum += sorted(i)
      threshold = (sum - mass) / (i + 1)
      i += 1
  Array.tabulate(a.rows)(r => Math.max(0.0, shifted(r) - threshold))

private def simplexViolation(a: DMat, mass: Double) = {
  var m = 0.0; var c = 0;
  while (c < a.cols) {
    var s = 0.0; var r = 0; while (r < a.rows) { m = Math.max(m, Math.max(0.0, -a(r, c))); s += a(r, c); r += 1 };
    m = Math.max(m, Math.abs(s - mass)); c += 1
  }; m
}

private def checkedBounds(lower: Double, upper: Double): Either[FirstOrderError, Unit] =
  if lower.isFinite && upper.isFinite && lower <= upper then Right(())
  else Left(FirstOrderError.InvalidConfiguration("box bounds must be finite and ordered"))
