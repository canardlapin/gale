package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError

/** The result of a '''symmetric''' eigenproblem `A V = V diag(λ)` — complex-free
  * by construction (fixed constraint 3 of `docs/spectral-parity.md`).
  *
  * Eigenvalues are returned '''ascending-algebraic always''', independent of the
  * selection order (§ 1): the selection chooses membership, this layout is fixed.
  * `eigenvectors` holds the orthonormal eigenvectors as columns, aligned with
  * `eigenvalues`; it is empty (zero columns) when only values were requested.
  */
final case class EigenDecomposition private[spectral] (
    eigenvalues: DVec,
    eigenvectors: DMat,
    diagnostics: SpectralDiagnostics
):
  /** The number of returned eigenvalues. */
  def size: Int =
    eigenvalues.length

  /** `Right(this)` when all requested pairs converged, else
    * `Left(`[[gale.linalg.LinAlgError.DidNotConverge]]`)`.
    */
  def requireConverged: Either[LinAlgError, EigenDecomposition] =
    diagnostics.requireConverged(this)

/** The result of a '''nonsymmetric''' eigenproblem, whose real input may have
  * complex eigenvalues in conjugate pairs (§ 2, § 7).
  *
  * Storage is real throughout: eigenvalues as structure-of-arrays real/imaginary
  * parts (`re`, `im`), right eigenvectors as a real `DMat` in LAPACK '''real-Schur
  * packed''' convention, and optional left eigenvectors in the same packing. That
  * packing is an implementation detail and is '''never''' exposed — the storage
  * fields are `private[gale]` and callers read the spectrum only through the
  * typed accessors below.
  *
  * '''Packing invariant.''' A real eigenvalue owns one column (its real
  * eigenvector, zero imaginary part). A complex-conjugate pair `λ = a ± b·i`
  * (`b > 0`) owns two '''adjacent''' columns `v_re, v_im`, with the
  * '''positive-imaginary member first''': the eigenvalue `a + b·i` at index `j`
  * has eigenvector `v_re + i·v_im` (columns `j`, `j+1`), and its conjugate
  * `a − b·i` at index `j+1` has eigenvector `v_re − i·v_im`. [[eigenvector]] and
  * [[leftEigenvector]] decode this, so no caller ever sees the shared columns.
  */
final class NonsymmetricEigenDecomposition private[spectral] (
    private[gale] val re: DVec,
    private[gale] val im: DVec,
    private[gale] val rightVectorsPacked: DMat,
    private[gale] val leftVectorsPacked: Option[DMat],
    val diagnostics: SpectralDiagnostics
):
  require(re.length == im.length, "re and im must have the same length")
  require(
    rightVectorsPacked.cols == 0 || rightVectorsPacked.cols == re.length,
    s"right eigenvector columns (${rightVectorsPacked.cols}) must be 0 or size (${re.length})"
  )
  leftVectorsPacked.foreach: packed =>
    require(
      packed.cols == 0 || packed.cols == re.length,
      s"left eigenvector columns (${packed.cols}) must be 0 or size (${re.length})"
    )
  // Enforce the real-Schur packing invariant at the SoA boundary: a complex
  // eigenvalue must appear as an adjacent conjugate pair with the positive-
  // imaginary member first, so the accessors can decode column ownership from a
  // single index. Both sides of each pair are checked, so a split pair or a
  // negative-first (or otherwise orphaned) member fails loudly here rather than
  // decoding to a silently wrong eigenvector. Short-circuit guards keep the
  // neighbour reads in bounds; one pass over `im`.
  private def validatePacking(): Unit =
    var i = 0
    while i < re.length do
      val imag = im(i)
      if imag > 0.0 then
        require(
          i + 1 < re.length && im(i + 1) == -imag && re(i + 1) == re(i),
          s"complex eigenvalue $i must be followed by its conjugate (positive-imaginary member first)"
        )
      else if imag < 0.0 then
        require(
          i > 0 && im(i - 1) == -imag && re(i - 1) == re(i),
          s"complex eigenvalue $i must be preceded by its positive-imaginary conjugate"
        )
      i += 1
  validatePacking()

  /** The number of eigenvalues. */
  def size: Int =
    re.length

  /** The `i`-th eigenvalue as a single [[Complex]]. An out-of-range `i` throws
    * [[gale.linalg.LinAlgError.IndexOutOfBounds]].
    */
  def eigenvalue(i: Int): Complex =
    checkIndex(i)
    Complex(re(i), im(i))

  /** True when the `i`-th eigenvalue is real (zero imaginary part), so it owns a
    * single real eigenvector column rather than a shared pair of columns. An
    * out-of-range `i` throws [[gale.linalg.LinAlgError.IndexOutOfBounds]].
    */
  def isRealPair(i: Int): Boolean =
    checkIndex(i)
    im(i) == 0.0

  /** The `i`-th '''right''' eigenvector as `(realPart, imaginaryPart)`, decoded
    * from the real-Schur packing. For a real eigenvalue the imaginary part is
    * zero; for the two members of a conjugate pair the imaginary parts are
    * negatives of each other. Throws
    * [[gale.linalg.LinAlgError.UnsupportedOperation]] if vectors were not
    * computed, or [[gale.linalg.LinAlgError.IndexOutOfBounds]] if `i` is out of
    * range.
    */
  def eigenvector(i: Int): (DVec, DVec) =
    checkIndex(i)
    decode(rightVectorsPacked, "right eigenvectors", i)

  /** The `i`-th '''left''' eigenvector as `(realPart, imaginaryPart)`, decoded
    * like [[eigenvector]]. Throws
    * [[gale.linalg.LinAlgError.UnsupportedOperation]] if left vectors were not
    * computed, or [[gale.linalg.LinAlgError.IndexOutOfBounds]] if `i` is out of
    * range.
    */
  def leftEigenvector(i: Int): (DVec, DVec) =
    checkIndex(i)
    leftVectorsPacked match
      case Some(packed) => decode(packed, "left eigenvectors", i)
      case None         => throw LinAlgError.UnsupportedOperation("left eigenvectors")

  /** `Right(this)` when all requested pairs converged, else
    * `Left(`[[gale.linalg.LinAlgError.DidNotConverge]]`)`.
    */
  def requireConverged: Either[LinAlgError, NonsymmetricEigenDecomposition] =
    diagnostics.requireConverged(this)

  /** Decode column(s) of a packed real-Schur eigenvector matrix at index `i` into
    * explicit `(realPart, imaginaryPart)` vectors. The imaginary part is a fresh
    * vector; the real part is a copy, so no view into the internal packing
    * escapes. Vector length is the ambient dimension `rightVectorsPacked.rows`,
    * which for a partial result exceeds [[size]].
    */
  private def decode(packed: DMat, what: String, i: Int): (DVec, DVec) =
    if packed.cols == 0 then
      throw LinAlgError.UnsupportedOperation(what)
    val imag = im(i)
    if imag == 0.0 then
      // Real eigenvalue: column i is the whole (real) eigenvector.
      (packed.col(i).copy, DVec.zeros(packed.rows))
    else if imag > 0.0 then
      // Positive-imaginary member: v_re + i·v_im over columns (i, i+1).
      (packed.col(i).copy, packed.col(i + 1).copy)
    else
      // Negative-imaginary member (conjugate): v_re − i·v_im, where v_re is the
      // pair's first column (i−1) and v_im is column i, negated.
      (packed.col(i - 1).copy, packed.col(i) * -1.0)

  private def checkIndex(i: Int): Unit =
    if i < 0 || i >= size then
      throw LinAlgError.IndexOutOfBounds(i, size)
