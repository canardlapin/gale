package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError

/** The result of a '''generalized nonsymmetric''' eigenproblem `A x = λ B x` (the
  * QZ / `ggev` family, § 1.3 of `docs/spectral-backend-boundary.md`). It mirrors
  * [[NonsymmetricEigenDecomposition]] — SoA real storage, real-Schur packed
  * vectors, typed accessors, `private[spectral]` constructor — but carries the
  * '''projective `(α, β)`''' spectrum, so an infinite eigenvalue from a singular /
  * rank-deficient `B` (`β = 0`) is representable, which a plain [[Complex]] cannot
  * express.
  *
  * Assembled only by the facade from a backend's raw QZ output; there is no
  * producer-less public path (every entry point requires a `given SpectralBackend`).
  *
  * '''Packing invariants''' (enforced at construction):
  *   - a '''finite''' complex eigenvalue `α = a ± b·i` (`b > 0`) owns two adjacent
  *     columns with the positive-imaginary member first, exact conjugate symmetry
  *     on `α`, '''and equal `β`''' — the `β` clause distinguishes two pairs that
  *     share an `α` but differ in `β`;
  *   - an '''infinite''' eigenvalue (`β = 0`) is always a real `1×1` block:
  *     `β(i) = 0 ⟹ alphaIm(i) = 0`. A complex-infinite eigenvalue is not
  *     representable and a backend that emits one is a conformance failure.
  */
final class GeneralizedEigenDecomposition private[spectral] (
    private[spectral] val alphaRe: DVec,
    private[spectral] val alphaIm: DVec,
    private[spectral] val beta: DVec,
    private[spectral] val rightVectorsPacked: DMat,
    private[spectral] val leftVectorsPacked: Option[DMat],
    val diagnostics: SpectralDiagnostics
):
  require(
    alphaRe.length == alphaIm.length && alphaIm.length == beta.length,
    "alphaRe, alphaIm, beta must have the same length"
  )
  require(
    rightVectorsPacked.cols == 0 || rightVectorsPacked.cols == alphaRe.length,
    s"right eigenvector columns (${rightVectorsPacked.cols}) must be 0 or size (${alphaRe.length})"
  )
  leftVectorsPacked.foreach: packed =>
    require(
      packed.cols == 0 || packed.cols == alphaRe.length,
      s"left eigenvector columns (${packed.cols}) must be 0 or size (${alphaRe.length})"
    )

  // Enforce the real-Schur packing + projective invariants: an infinite eigenvalue
  // is real (β=0 ⟹ alphaIm=0); a complex eigenvalue appears as an adjacent conjugate
  // pair (positive-imaginary first) with exact α symmetry AND equal β.
  private def validatePacking(): Unit =
    var i = 0
    while i < alphaRe.length do
      if beta(i) == 0.0 then
        require(alphaIm(i) == 0.0, s"infinite eigenvalue $i (β=0) must be real (alphaIm=0)")
      val imag = alphaIm(i)
      if imag > 0.0 then
        require(
          i + 1 < alphaRe.length && alphaIm(i + 1) == -imag && alphaRe(i + 1) == alphaRe(i) && beta(i + 1) == beta(i),
          s"complex eigenvalue $i must be followed by its conjugate (equal α real part, negated α imag, equal β)"
        )
      else if imag < 0.0 then
        require(
          i > 0 && alphaIm(i - 1) == -imag && alphaRe(i - 1) == alphaRe(i) && beta(i - 1) == beta(i),
          s"complex eigenvalue $i must be preceded by its positive-imaginary conjugate"
        )
      i += 1
  validatePacking()

  /** The number of generalized eigenvalues. */
  def size: Int =
    alphaRe.length

  /** The `i`-th eigenvalue as a projective [[GeneralizedEigenvalue]] `(α, β)`; use
    * [[isInfinite]] before reading `.value` (which is `α/β`).
    */
  def eigenvalue(i: Int): GeneralizedEigenvalue =
    checkIndex(i)
    GeneralizedEigenvalue(Complex(alphaRe(i), alphaIm(i)), beta(i))

  /** True when the `i`-th eigenvalue is infinite (`β = 0`) — a singular /
    * rank-deficient `B` direction.
    */
  def isInfinite(i: Int): Boolean =
    checkIndex(i)
    beta(i) == 0.0

  /** True when the `i`-th eigenvalue owns a single real column — a finite real `λ`
    * or an infinite eigenvalue (both have `alphaIm = 0`).
    */
  def isRealPair(i: Int): Boolean =
    checkIndex(i)
    alphaIm(i) == 0.0

  /** The `i`-th '''right''' eigenvector as `(realPart, imaginaryPart)`, decoded from
    * the real-Schur packing exactly like [[NonsymmetricEigenDecomposition]].
    */
  def eigenvector(i: Int): (DVec, DVec) =
    checkIndex(i)
    decode(rightVectorsPacked, "right eigenvectors", i)

  /** The `i`-th '''left''' eigenvector as `(realPart, imaginaryPart)`. */
  def leftEigenvector(i: Int): (DVec, DVec) =
    checkIndex(i)
    leftVectorsPacked match
      case Some(packed) => decode(packed, "left eigenvectors", i)
      case None         => throw LinAlgError.UnsupportedOperation("left eigenvectors")

  /** `Right(this)` when all requested pairs converged, else
    * `Left(`[[gale.linalg.LinAlgError.DidNotConverge]]`)`.
    */
  def requireConverged: Either[LinAlgError, GeneralizedEigenDecomposition] =
    diagnostics.requireConverged(this)

  /** `Right(this)` only when the requested global spectral extreme is
    * independently certified. See
    * [[SpectralDiagnostics.requireExtremeCertified]].
    */
  def requireExtremeCertified: Either[LinAlgError, GeneralizedEigenDecomposition] =
    diagnostics.requireExtremeCertified(this)

  private def decode(packed: DMat, what: String, i: Int): (DVec, DVec) =
    if packed.cols == 0 then throw LinAlgError.UnsupportedOperation(what)
    val imag = alphaIm(i)
    if imag == 0.0 then (packed.col(i).copy, DVec.zeros(packed.rows))
    else if imag > 0.0 then (packed.col(i).copy, packed.col(i + 1).copy)
    else (packed.col(i - 1).copy, packed.col(i) * -1.0)

  private def checkIndex(i: Int): Unit =
    if i < 0 || i >= size then throw LinAlgError.IndexOutOfBounds(i, size)
