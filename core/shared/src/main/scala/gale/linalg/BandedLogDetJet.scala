package gale.linalg

import gale.platform.DoubleArray
import gale.platform.DoubleArray.*

/** Actual first and second derivatives of a symmetric banded matrix `A(θ)`,
  * `θ` of dimension 1 to 3, in the lower-band layout of
  * [[BandedCholesky.factorLower]]: `bands(i, d) = ∂A(i, i-d)`.
  *
  * `secondUpper(q*(q+1)/2 + p)` holds `∂²A/∂θp∂θq` for `p <= q`. These are
  * the second derivatives themselves, not Taylor half-coefficients. Every
  * matrix shares one `n`-by-`(b+1)` shape, which must equal the factor's
  * declared band even where a derivative entry is zero.
  */
final class BandedDerivatives private (val first: Vector[DMat], val secondUpper: Vector[DMat]):
  def dimension: Int = first.length

  def second(p: Int, q: Int): DMat =
    secondUpper(BandedDerivatives.pairIndex(math.min(p, q), math.max(p, q)))

object BandedDerivatives:
  val MaxDimension: Int = 3

  private[linalg] def pairIndex(p: Int, q: Int): Int = q * (q + 1) / 2 + p

  def apply(first: Vector[DMat], secondUpper: Vector[DMat]): Either[LinAlgError, BandedDerivatives] =
    val d = first.length
    if d < 1 || d > MaxDimension then
      Left(LinAlgError.InvalidArgument(s"banded derivative dimension must be 1..$MaxDimension, got $d"))
    else if secondUpper.length != d * (d + 1) / 2 then
      Left(LinAlgError.InvalidArgument(
        s"expected ${d * (d + 1) / 2} upper-packed second derivatives for dimension $d, got ${secondUpper.length}"))
    else
      val shape = first.head.shape
      (first ++ secondUpper).find(_.shape != shape) match
        case Some(other) => Left(LinAlgError.DimensionMismatch(shape, other.shape))
        case None => Right(new BandedDerivatives(first, secondUpper))

/** Gradient and full symmetric Hessian of `log det A(θ)` at an accepted
  * factor. `value` is the factor's own `logDet`; `factor` is the same instance
  * the jet was computed from. The recursion is floating-point evaluation of
  * exact identities, not a directed-rounding certificate.
  */
final class BandedLogDetJet private[linalg] (val factor: BandedCholesky, val gradient: DVec, val hessian: DMat):
  def value: Double = factor.logDet
  def dimension: Int = gradient.length

private[linalg] object BandedLogDetJet:
  /** Differentiates the accepted row-by-row band recursion
    * `L(i,j) L(j,j) = A(i,j) - Σ_{k<j} L(i,k) L(j,k)` once for every axis
    * and twice for each upper pair. All `d` first-derivative bands are
    * retained; each pair reuses one band of scratch, so the workspace is
    * `(d+1) n (b+1)` doubles. There are no RHS solves, inverses or pivots: the
    * factor's accepted pivot policy is inherited unchanged.
    */
  def compute(factor: BandedCholesky, derivatives: BandedDerivatives): Either[LinAlgError, BandedLogDetJet] =
    val n = factor.size
    val width = factor.bandwidth + 1
    val b = factor.bandwidth
    val d = derivatives.dimension
    val expected = Shape(Rows(n), Cols(width))
    val all = derivatives.first ++ derivatives.secondUpper
    all.find(_.shape != expected) match
      case Some(other) => return Left(LinAlgError.DimensionMismatch(expected, other.shape))
      case None => ()
    var index = 0
    while index < all.length do
      val m = all(index)
      var i = 0
      while i < n do
        var delta = 0
        while delta <= math.min(i, b) do
          if !entry(m, i, delta).isFinite then
            val role = if index < d then s"first derivative $index" else s"second derivative ${index - d}"
            return Left(LinAlgError.InvalidArgument(s"nonfinite active $role at ($i,$delta)"))
          delta += 1
        i += 1
      index += 1
    val cells = n.toLong * width
    val scratch = (d + 1).toLong * cells
    if scratch > Int.MaxValue.toLong - 8L then
      return Left(LinAlgError.InvalidArgument(s"banded log-determinant jet scratch of $scratch doubles exceeds array limits"))

    val bandCells = cells.toInt
    val l = factor.lowerBands.data
    val firstBands = DoubleArray.alloc(d * bandCells)
    val pairBand = DoubleArray.alloc(bandCells)
    val gradient = DoubleArray.alloc(d)
    val hessian = DoubleArray.alloc(d * d)

    var p = 0
    while p < d do
      val ap = derivatives.first(p)
      val lp = p * bandCells
      var sum = 0.0
      var i = 0
      while i < n do
        val start = math.max(0, i - b)
        var j = start
        while j < i do
          var value = entry(ap, i, i - j)
          var k = start
          while k < j do
            value -= firstBands(lp + i * width + i - k) * l(j * width + j - k) +
              l(i * width + i - k) * firstBands(lp + j * width + j - k)
            k += 1
          value -= l(i * width + i - j) * firstBands(lp + j * width)
          value /= l(j * width)
          if !value.isFinite then return Left(nonfinite(s"first derivative $p", i, j))
          firstBands(lp + i * width + i - j) = value
          j += 1
        var diagonal = entry(ap, i, 0)
        var k = start
        while k < i do
          diagonal -= 2.0 * l(i * width + i - k) * firstBands(lp + i * width + i - k)
          k += 1
        diagonal /= 2.0 * l(i * width)
        if !diagonal.isFinite then return Left(nonfinite(s"first derivative $p", i, i))
        firstBands(lp + i * width) = diagonal
        sum += diagonal / l(i * width)
        i += 1
      gradient(p) = 2.0 * sum
      if !gradient(p).isFinite then return Left(nonfinite(s"gradient $p", n - 1, n - 1))
      p += 1

    var q = 0
    while q < d do
      p = 0
      while p <= q do
        val apq = derivatives.secondUpper(BandedDerivatives.pairIndex(p, q))
        val lp = p * bandCells
        val lq = q * bandCells
        var sum = 0.0
        var i = 0
        while i < n do
          val start = math.max(0, i - b)
          var j = start
          while j < i do
            val ij = i * width + i - j
            val jj = j * width
            var value = entry(apq, i, i - j)
            var k = start
            while k < j do
              val ik = i * width + i - k
              val jk = j * width + j - k
              value -= pairBand(ik) * l(jk) + firstBands(lp + ik) * firstBands(lq + jk) +
                firstBands(lq + ik) * firstBands(lp + jk) + l(ik) * pairBand(jk)
              k += 1
            value -= firstBands(lp + ij) * firstBands(lq + jj) + firstBands(lq + ij) * firstBands(lp + jj) +
              l(ij) * pairBand(jj)
            value /= l(jj)
            if !value.isFinite then return Left(nonfinite(s"second derivative ($p,$q)", i, j))
            pairBand(ij) = value
            j += 1
          val ii = i * width
          var diagonal = entry(apq, i, 0)
          var k = start
          while k < i do
            val ik = i * width + i - k
            diagonal -= 2.0 * (firstBands(lp + ik) * firstBands(lq + ik) + l(ik) * pairBand(ik))
            k += 1
          diagonal -= 2.0 * firstBands(lp + ii) * firstBands(lq + ii)
          diagonal /= 2.0 * l(ii)
          if !diagonal.isFinite then return Left(nonfinite(s"second derivative ($p,$q)", i, i))
          pairBand(ii) = diagonal
          val pivot = l(ii)
          sum += diagonal / pivot - (firstBands(lp + ii) / pivot) * (firstBands(lq + ii) / pivot)
          i += 1
        val value = 2.0 * sum
        if !value.isFinite then return Left(nonfinite(s"Hessian ($p,$q)", n - 1, n - 1))
        hessian(p * d + q) = value
        hessian(q * d + p) = value
        p += 1
      q += 1

    Right(new BandedLogDetJet(factor, DVec.fromDoubleArrayOwned(gradient), DMat.fromDoubleArrayOwned(d, d, hessian)))

  /** Logical lower-band entry of a possibly strided view, without copying. */
  private inline def entry(m: DMat, row: Int, delta: Int): Double =
    m.data(m.offset.value + row * m.rowStride.value + delta * m.colStride.value)

  private def nonfinite(what: String, row: Int, col: Int): LinAlgError =
    LinAlgError.InvalidArgument(s"nonfinite banded log-determinant $what at ($row,$col)")
