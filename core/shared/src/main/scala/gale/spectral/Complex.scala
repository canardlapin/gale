package gale.spectral

/** A complex number, used only as a '''boundary value type''' at the spectral
  * API surface.
  *
  * gale's spectral inputs and storage are real throughout (fixed constraint 1 of
  * `docs/spectral-parity.md`): there is no complex storage tier and no complex
  * kernel tier. `Complex` exists solely so that nonsymmetric and generalized
  * eigenvalue accessors can hand back a conjugate-pair eigenvalue as a single
  * value. It carries only the arithmetic those accessors need — magnitude,
  * conjugation, and division by a real scalar (to form `α/β`) — deliberately not
  * a full complex-arithmetic tier.
  */
final case class Complex(re: Double, im: Double):
  /** The modulus `√(re² + im²)`, via [[scala.math.hypot]] to avoid intermediate
    * overflow.
    */
  def magnitude: Double =
    math.hypot(re, im)

  /** The complex conjugate `re − im·i`. */
  def conjugate: Complex =
    Complex(re, -im)

  /** This value divided by a real scalar, componentwise. Used to form a
    * generalized eigenvalue `α/β` from its projective `(α, β)` representation;
    * a zero divisor yields the IEEE infinities/NaNs that mark an infinite
    * generalized eigenvalue.
    */
  def /(scalar: Double): Complex =
    Complex(re / scalar, im / scalar)

object Complex:
  /** The additive identity `0 + 0i`. */
  val Zero: Complex = Complex(0.0, 0.0)

  /** A real number as a complex value with zero imaginary part. */
  def real(re: Double): Complex =
    Complex(re, 0.0)
