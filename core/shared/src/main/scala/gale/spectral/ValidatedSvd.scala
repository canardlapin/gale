package gale.spectral

import gale.linalg.{Cols, DMat, DVec, LinAlgError, Rows, Shape}

/** An entrywise enclosure of a real dense matrix.  Endpoints are mathematical
  * bounds, rather than a tolerance convention: callers forming a matrix from
  * floating point arithmetic must account for every input operation before
  * constructing this value.
  */
final case class MatrixEnclosure private (lower: DMat, upper: DMat)

object MatrixEnclosure:
  def exact(value: DMat): Either[LinAlgError, MatrixEnclosure] =
    checked(value, value)

  def checked(lower: DMat, upper: DMat): Either[LinAlgError, MatrixEnclosure] =
    if lower.rows != upper.rows || lower.cols != upper.cols then
      Left(LinAlgError.DimensionMismatch(Shape(Rows(lower.rows), Cols(lower.cols)), Shape(Rows(upper.rows), Cols(upper.cols))))
    else
      var r = 0
      while r < lower.rows do
        var c = 0
        while c < lower.cols do
          if !lower(r, c).isFinite || !upper(r, c).isFinite || lower(r, c) > upper(r, c) then
            return Left(LinAlgError.InvalidArgument(s"invalid matrix enclosure endpoint at ($r,$c)"))
          c += 1
        r += 1
      Right(MatrixEnclosure(lower, upper))

/** Certified singular-value intervals obtained from finite candidate factors.
  * `factorError` bounds the Frobenius reconstruction error of every member of
  * the input enclosure; `leftDefect` and `rightDefect` bound the corresponding
  * Gram-matrix Frobenius defects.  The latter are deliberately Frobenius
  * bounds, hence also operator-norm bounds used by the polar-factor argument.
  */
final case class ValidatedSvdEnclosure(
    lower: DVec,
    upper: DVec,
    factorError: Double,
    leftDefect: Double,
    rightDefect: Double
)

/** Validates full/economy candidate SVD factors without trusting the solver.
  *
  * Let `W = U diag(s) Vt`.  Outward arithmetic bounds `e >= ||Z-W||_F`, and
  * `du`, `dv` bound the orthogonality defects.  For `du,dv < 1`, polar factors
  * give singular values of W between
  * `s_i sqrt(1-du) sqrt(1-dv)` and `s_i sqrt(1+du) sqrt(1+dv)`; Weyl's bound
  * then encloses every singular value of Z.  The supplied factors are merely
  * candidates.  A failed finite/arithmetic/defect check returns `Left` rather
  * than turning a heuristic SVD into a certificate.
  */
object ValidatedSvd:
  def enclosure(input: MatrixEnclosure, candidate: SVD): Either[LinAlgError, ValidatedSvdEnclosure] =
    val p = math.min(input.lower.rows, input.lower.cols)
    if candidate.size != p || candidate.u.rows != input.lower.rows || candidate.u.cols != p ||
        candidate.vt.rows != p || candidate.vt.cols != input.lower.cols then
      Left(LinAlgError.InvalidArgument("candidate factors do not have full/economy SVD shapes for the enclosure"))
    else
      val sigma = (0 until p).map(candidate.singularValues(_)).toVector
      if sigma.exists(s => !s.isFinite || s < 0.0) || sigma.sliding(2).exists(pair => pair.length == 2 && pair(0) < pair(1)) then
        Left(LinAlgError.InvalidArgument("candidate singular values must be finite, nonnegative and descending"))
      else
        for
          du <- gramDefect(candidate.u)
          dv <- gramDefect(candidate.vt.t)
          _ <- Either.cond(du < 1.0 && dv < 1.0, (), LinAlgError.InvalidArgument("candidate orthogonality defect is too large for a polar-factor enclosure"))
          e <- reconstructionError(input, candidate.u, sigma, candidate.vt)
          oneMinusDu <- downSub(1.0, du)
          onePlusDu <- upAdd(1.0, du)
          oneMinusDv <- downSub(1.0, dv)
          onePlusDv <- upAdd(1.0, dv)
          lowLeft <- sqrtLower(oneMinusDu)
          highLeft <- sqrtUpper(onePlusDu)
          lowRight <- sqrtLower(oneMinusDv)
          highRight <- sqrtUpper(onePlusDv)
          bounds <- singularBounds(sigma, lowLeft, highLeft, lowRight, highRight, e)
        yield ValidatedSvdEnclosure(
          DVec.tabulate(p)(i => bounds(i)._1),
          DVec.tabulate(p)(i => bounds(i)._2),
          e, du, dv
        )

  private def singularBounds(sigma: Vector[Double], lu: Double, uu: Double, lv: Double, uv: Double, e: Double)
      : Either[LinAlgError, Vector[(Double, Double)]] =
    val out = Vector.newBuilder[(Double, Double)]
    var i = 0
    while i < sigma.length do
      val bound = for
        firstLow <- downMul(sigma(i), lu)
        baseLow <- downMul(firstLow, lv)
        firstHigh <- upMul(sigma(i), uu)
        baseHigh <- upMul(firstHigh, uv)
        lower <- downSub(baseLow, e)
        upper <- upAdd(baseHigh, e)
      yield math.max(0.0, lower) -> upper
      bound match
        case Left(error) => return Left(error)
        case Right(value) => out += value
      i += 1
    Right(out.result())

  private def gramDefect(columns: DMat): Either[LinAlgError, Double] =
    val k = columns.cols
    var squares = 0.0
    var i = 0
    while i < k do
      var j = 0
      while j < k do
        var lo = 0.0
        var hi = 0.0
        var r = 0
        while r < columns.rows do
          intervalProduct(columns(r, i), columns(r, j)).flatMap { term =>
            for nextLo <- downAdd(lo, term._1); nextHi <- upAdd(hi, term._2) yield nextLo -> nextHi
          } match
            case Left(error) => return Left(error)
            case Right(value) =>
              lo = value._1
              hi = value._2
          r += 1
        val difference =
          if i == j then
            (downSub(lo, 1.0), upAdd(hi, -1.0)) match
              case (Right(deltaLo), Right(deltaHi)) => math.max(math.abs(deltaLo), math.abs(deltaHi))
              case (Left(error), _) => return Left(error)
              case (_, Left(error)) => return Left(error)
          else math.max(math.abs(lo), math.abs(hi))
        upMul(difference, difference).flatMap(term => upAdd(squares, term)) match
          case Left(error) => return Left(error)
          case Right(value) => squares = value
        j += 1
      i += 1
    sqrtUpper(squares)

  private def reconstructionError(input: MatrixEnclosure, u: DMat, sigma: Vector[Double], vt: DMat): Either[LinAlgError, Double] =
    var squares = 0.0
    var r = 0
    while r < input.lower.rows do
      var c = 0
      while c < input.lower.cols do
        var lo = 0.0
        var hi = 0.0
        var k = 0
        while k < sigma.length do
          val product = for
            first <- intervalProduct(u(r, k), sigma(k))
            term <- intervalProduct(first._1, first._2, vt(k, c), vt(k, c))
            nextLo <- downAdd(lo, term._1)
            nextHi <- upAdd(hi, term._2)
          yield nextLo -> nextHi
          product match
            case Left(error) => return Left(error)
            case Right(value) =>
              lo = value._1
              hi = value._2
          k += 1
        val distance = (downSub(input.lower(r, c), hi), upAdd(input.upper(r, c), -lo)) match
          case (Right(deltaLo), Right(deltaHi)) => math.max(math.abs(deltaLo), math.abs(deltaHi))
          case (Left(error), _) => return Left(error)
          case (_, Left(error)) => return Left(error)
        upMul(distance, distance).flatMap(term => upAdd(squares, term)) match
          case Left(error) => return Left(error)
          case Right(value) => squares = value
        c += 1
      r += 1
    sqrtUpper(squares)

  private def intervalProduct(a: Double, b: Double): Either[LinAlgError, (Double, Double)] =
    if !a.isFinite || !b.isFinite then Left(LinAlgError.InvalidArgument("candidate factor is nonfinite"))
    else
      for low <- downMul(a, b); high <- upMul(a, b) yield low -> high

  private def intervalProduct(aLo: Double, aHi: Double, bLo: Double, bHi: Double): Either[LinAlgError, (Double, Double)] =
    val products = Vector(aLo * bLo, aLo * bHi, aHi * bLo, aHi * bHi)
    if products.exists(value => !value.isFinite) then Left(LinAlgError.InvalidArgument("interval product overflowed or became nonfinite"))
    else for low <- down(products.min); high <- up(products.max) yield low -> high

  private def sqrtLower(value: Double): Either[LinAlgError, Double] =
    if !value.isFinite || value < 0.0 then Left(LinAlgError.InvalidArgument("square-root enclosure argument is invalid"))
    else if value == 0.0 then Right(0.0)
    else
      val candidate = math.sqrt(value)
      if !candidate.isFinite || candidate < 0.0 then Left(LinAlgError.InvalidArgument("square-root enclosure overflowed or was negative"))
      else
        // Do not assume the platform's sqrt was correctly rounded: decrease
        // until its square is a lower endpoint, then take one further endpoint.
        def seek(lower: Double, steps: Int): Either[LinAlgError, Double] =
          if !lower.isFinite || lower < 0.0 then Left(LinAlgError.InvalidArgument("could not validate square-root lower endpoint"))
          else upMul(lower, lower).flatMap { square =>
            if square <= value then down(lower)
            else if steps == 8 then Left(LinAlgError.InvalidArgument("could not validate square-root lower endpoint"))
            else seek(java.lang.Math.nextDown(lower), steps + 1)
          }
        seek(candidate, 0)

  private def sqrtUpper(value: Double): Either[LinAlgError, Double] =
    if !value.isFinite || value < 0.0 then Left(LinAlgError.InvalidArgument("square-root enclosure argument is invalid"))
    else if value == 0.0 then Right(0.0)
    else
      val candidate = math.sqrt(value)
      if !candidate.isFinite || candidate < 0.0 then Left(LinAlgError.InvalidArgument("square-root enclosure overflowed or was negative"))
      else
        // JS is permitted to approximate sqrt.  Validate by squaring and walk
        // outward, instead of asserting that one ulp is universally enough.
        def seek(upper: Double, steps: Int): Either[LinAlgError, Double] =
          if !upper.isFinite || upper < 0.0 then Left(LinAlgError.InvalidArgument("could not validate square-root upper endpoint"))
          else downMul(upper, upper).flatMap { square =>
            if square >= value then up(upper)
            else if steps == 8 then Left(LinAlgError.InvalidArgument("could not validate square-root upper endpoint"))
            else seek(java.lang.Math.nextUp(upper), steps + 1)
          }
        seek(candidate, 0)

  private def down(value: Double): Either[LinAlgError, Double] =
    if !value.isFinite then Left(LinAlgError.InvalidArgument("outward arithmetic overflowed or became nonfinite"))
    else
      val endpoint = java.lang.Math.nextDown(value)
      if endpoint.isFinite then Right(endpoint)
      else Left(LinAlgError.InvalidArgument("outward lower endpoint overflowed"))
  private def up(value: Double): Either[LinAlgError, Double] =
    if !value.isFinite then Left(LinAlgError.InvalidArgument("outward arithmetic overflowed or became nonfinite"))
    else
      val endpoint = java.lang.Math.nextUp(value)
      if endpoint.isFinite then Right(endpoint)
      else Left(LinAlgError.InvalidArgument("outward upper endpoint overflowed"))
  private def downAdd(a: Double, b: Double): Either[LinAlgError, Double] = down(a + b)
  private def upAdd(a: Double, b: Double): Either[LinAlgError, Double] = up(a + b)
  private def downSub(a: Double, b: Double): Either[LinAlgError, Double] = down(a - b)
  private def downMul(a: Double, b: Double): Either[LinAlgError, Double] = down(a * b)
  private def upMul(a: Double, b: Double): Either[LinAlgError, Double] = up(a * b)
