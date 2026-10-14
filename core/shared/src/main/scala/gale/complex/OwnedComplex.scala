package gale.complex

import gale.linalg.DMat
import gale.spectral.Complex

enum ComplexFailure:
  case Invalid(detail: String)
  case Budget(required: BigInt, available: BigInt)
  case Numerical(detail: String)
  case Cancelled
  case NonConverged(sweeps: Int, offDiagonal: Double)
  case Indefinite(minimum: Double, allowedNegative: Double)
  case RankDeficient(rank: Int, dimension: Int)

/** Immutable split-complex row-major storage. Construction copies; no writable view escapes. */
final class OwnedComplex private[complex] (
    val rows: Int, val cols: Int, private[complex] val re: Array[Double], private[complex] val im: Array[Double]
):
  val storageBytes: BigInt = BigInt(rows) * cols * 16
  def apply(row: Int, column: Int): Complex =
    require(row >= 0 && row < rows && column >= 0 && column < cols, "complex index outside shape")
    Complex(re(row * cols + column), im(row * cols + column))
  def realCopy: Array[Double] = re.clone()
  def imaginaryCopy: Array[Double] = im.clone()
  def conjugateTranspose(budget: BigInt): Either[ComplexFailure, OwnedComplex] =
    if storageBytes > budget then Left(ComplexFailure.Budget(storageBytes, budget))
    else
      val rr = new Array[Double](re.length); val ii = new Array[Double](im.length)
      var r = 0
      while r < rows do
        var c = 0
        while c < cols do
          rr(c * rows + r) = re(r * cols + c); ii(c * rows + r) = -im(r * cols + c)
          c += 1
        r += 1
      Right(new OwnedComplex(cols, rows, rr, ii))
  /** Output payload budget; the immutable operands already belong to the caller. */
  def multiply(other: OwnedComplex, budget: BigInt): Either[ComplexFailure, OwnedComplex] =
    if cols != other.rows then Left(ComplexFailure.Invalid("complex product inner dimensions differ"))
    else OwnedComplex.inspect(rows, other.cols, budget).flatMap: _ =>
      val rr = new Array[Double](rows * other.cols); val ii = new Array[Double](rr.length)
      var r = 0
      var invalid = false
      while r < rows do
        var c = 0
        while c < other.cols do
          var real = 0.0; var imaginary = 0.0; var k = 0
          while k < cols do
            val a = r * cols + k; val b = k * other.cols + c
            real += re(a) * other.re(b) - im(a) * other.im(b)
            imaginary += re(a) * other.im(b) + im(a) * other.re(b)
            k += 1
          rr(r * other.cols + c) = real; ii(r * other.cols + c) = imaginary
          invalid ||= !real.isFinite || !imaginary.isFinite
          c += 1
        r += 1
      if invalid then Left(ComplexFailure.Numerical("nonfinite complex product"))
      else Right(new OwnedComplex(rows, other.cols, rr, ii))

object OwnedComplex:
  def inspect(rows: Int, cols: Int, budget: BigInt): Either[ComplexFailure, BigInt] =
    val cells = BigInt(rows) * cols
    if rows < 0 || cols < 0 || cells > Int.MaxValue then Left(ComplexFailure.Invalid("bounded nonnegative complex shape required"))
    else if cells * 16 > budget then Left(ComplexFailure.Budget(cells * 16, budget))
    else Right(cells * 16)
  def copying(rows: Int, cols: Int, real: Array[Double], imaginary: Array[Double], budget: BigInt): Either[ComplexFailure, OwnedComplex] =
    inspect(rows, cols, budget).flatMap: _ =>
      if real.length != rows * cols || imaginary.length != real.length then Left(ComplexFailure.Invalid("split-complex payload shape differs"))
      else if !real.forall(_.isFinite) || !imaginary.forall(_.isFinite) then Left(ComplexFailure.Numerical("nonfinite complex input"))
      else Right(new OwnedComplex(rows, cols, real.clone(), imaginary.clone()))
  def fromMatrices(real: DMat, imaginary: DMat, budget: BigInt): Either[ComplexFailure, OwnedComplex] =
    inspect(real.rows, real.cols, budget).flatMap: _ =>
      if real.rows != imaginary.rows || real.cols != imaginary.cols then Left(ComplexFailure.Invalid("split-complex matrix shapes differ"))
      else
        val rr = new Array[Double](real.rows * real.cols); val ii = new Array[Double](rr.length)
        var k = 0; var invalid = false
        while k < rr.length do
          rr(k) = real(k / real.cols, k % real.cols); ii(k) = imaginary(k / real.cols, k % real.cols)
          invalid ||= !rr(k).isFinite || !ii(k).isFinite
          k += 1
        if invalid then Left(ComplexFailure.Numerical("nonfinite complex input"))
        else Right(new OwnedComplex(real.rows, real.cols, rr, ii))
