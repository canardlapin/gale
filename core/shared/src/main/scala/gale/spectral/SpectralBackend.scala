package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError

/** Raw numeric carriers a [[SpectralBackend]] returns — public, invariant-free
  * containers holding only what a native routine naturally produces. The facade
  * (in `gale.spectral`) canonicalizes each into the corresponding sealed result
  * type (`private[spectral]` constructor after migration step 0), imposing gale's
  * ordering/packing and re-deriving residuals; a backend never builds a result
  * type (§ 1.2 of `docs/spectral-backend-boundary.md`).
  */
final case class RawSymmetricEigen(values: DVec, vectors: DMat)

/** Raw nonsymmetric eigendecomposition, real-Schur SoA as `geev`/`dgeev` emit. */
final case class RawNonsymmetricEigen(re: DVec, im: DVec, rightPacked: DMat, leftPacked: Option[DMat])

/** Raw SVD triplets in any order — the facade fixes the descending layout.
  * Shape contract: for an `m×n` input with `p = min(m, n)`, `u` must carry AT
  * LEAST the leading `p` singular vectors as columns (`u.cols >= p`; a full
  * `m×m` LAPACK `jobz='A'` factor is fine) and `vt` at least the leading `p`
  * rows; the facade slices the leading block and discards any excess.
  */
final case class RawSvd(sigma: DVec, u: DMat, vt: DMat)

/** Raw generalized (QZ) eigendecomposition — the projective `(α, β)` SoA spectrum
  * (`beta == 0.0` marks an infinite eigenvalue), real-Schur packed vectors, and
  * the optional generalized-Schur factors (present iff the backend is
  * [[SpectralCapability.GeneralizedSchur]]-capable).
  *
  * '''Contract: `beta(i) ≥ 0` for every `i`''' (the LAPACK `ggev`/`gges`
  * convention). This is load-bearing: the facade's projective sorter compares
  * `α/β` by cross-multiplication, which stays a valid total order only for
  * non-negative `β`; a negative `β` is rejected by the facade with
  * `Left(InvalidArgument)` as a conformance violation, not silently mis-sorted.
  */
final case class RawGeneralizedEigen(
    alphaRe: DVec,
    alphaIm: DVec,
    beta: DVec,
    rightPacked: DMat,
    leftPacked: Option[DMat],
    schur: Option[QzSchur]
)

/** Raw GSVD factors — `c`, `s` as direct column norms (the facade classifies the
  * typed generalized singular values and imposes the descending order).
  */
final case class RawGsvd(u: DMat, v: DMat, x: DMat, c: DVec, s: DVec)

/** Optional generalized-Schur factors `A = Q·AA·Zᵀ`, `B = Q·BB·Zᵀ` — supplied only
  * by a [[SpectralCapability.GeneralizedSchur]]-capable backend; not required by the
  * boundary (§ 4, non-goal 6).
  */
final case class QzSchur(aa: DMat, bb: DMat, q: DMat, z: DMat)

/** The counts an iterative/rank-revealing backend alone knows, carried alongside a
  * raw carrier: iteration count, requested/converged pair counts, and numerical
  * rank where meaningful. Dense one-shot methods report `iterations = 0`,
  * `converged = requested`, and let the facade derive `rank` (§ 1.5).
  */
final case class BackendConvergence(requested: Int, converged: Int, iterations: Int, rank: Option[Int] = None)

/** The operation-level capabilities a [[SpectralBackend]] may advertise (§ 1.1).
  * Distinct from the general kernel-acceleration `Capability` enum: these are
  * "can this backend do QZ at all?", not "does it make `gemm` fast?".
  *
  * `ComplexShift` is deliberately absent — off-real-axis targeting needs a complex
  * tier, which is out (§ 4).
  */
enum SpectralCapability:
  case DenseSymmetricEigen, DenseNonsymmetricEigen, DenseSvd
  case GeneralizedNonsymmetricEigen
  case GeneralizedSchur
  case RankDeficientGsvd
  case ShiftInvertSolve
  case IterativeGeneralized

/** An optional provider of spectral computations the pure core does not ship
  * (§ 1.1 of `docs/spectral-backend-boundary.md`). Lives in gale-core as an
  * interface; every implementation is in a separate optional JVM module, brought
  * into scope by an explicit `given` import (§ 2). Methods take already-validated
  * dense inputs and return '''raw numeric carriers''' — never the sealed result
  * types, which the facade builds after canonicalization and residual re-derivation.
  *
  * '''Discovery vs invocation.''' [[capabilities]] answers the yes/no question; the
  * total `Either`-returning methods do the work. A backend that lists a
  * [[SpectralCapability]] '''MUST''' implement the corresponding method (return a
  * non-`UnsupportedOperation` for structurally valid input) and '''MUST NOT''' list
  * one it does not honour. This is a conformance obligation (verified by the laws
  * conformance suite), not a type-level guarantee.
  *
  * '''Thread-safety (a REQUIREMENT, G1).''' A given `SpectralBackend` is a shared
  * singleton resolved once and reused across all call sites and threads. An
  * implementation MUST be safe for concurrent invocation — either stateless, or
  * internally synchronized (a backend wrapping a non-reentrant native library MUST
  * serialize internally). Threading configuration does '''not''' flow through this
  * boundary (G2): a backend reads its thread policy from the shared `BackendConfig`
  * at construction, and the spectral facades take no per-call threading argument.
  */
trait SpectralBackend:
  def name: String
  def capabilities: Set[SpectralCapability]

  // Every operation defaults to Left(UnsupportedOperation); a backend overrides
  // only the ones it supports and lists exactly those in `capabilities`.
  def denseSymmetricEigen(a: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
    unsupported("dense symmetric eigen")

  def denseNonsymmetricEigen(a: DMat, vectors: EigenVectors): Either[LinAlgError, RawNonsymmetricEigen] =
    unsupported("dense nonsymmetric eigen")

  def denseSvd(a: DMat, wantVectors: Boolean): Either[LinAlgError, RawSvd] =
    unsupported("dense SVD")

  def generalizedNonsymmetricEigen(a: DMat, b: DMat, vectors: EigenVectors): Either[LinAlgError, RawGeneralizedEigen] =
    unsupported("generalized nonsymmetric eigen (QZ)")

  def rankDeficientGsvd(a: DMat, b: DMat, wantVectors: Boolean): Either[LinAlgError, RawGsvd] =
    unsupported("rank-deficient GSVD")

  def shiftInvertOperator(a: DMat, b: Option[DMat], sigma: Double): Either[LinAlgError, DoubleLinearOperator] =
    unsupported("shift-invert operator")

  protected final def unsupported(op: String): Either[LinAlgError, Nothing] =
    Left(LinAlgError.UnsupportedOperation(s"$op: $name backend does not provide it"))

object SpectralBackend:

  /** The always-in-scope, lowest-priority fallback: no capabilities, every
    * operation `Left(UnsupportedOperation)`. With no acceleration module imported
    * this is the only candidate for a `using SpectralBackend`, so every
    * backend-scoped call reproduces today's seam behaviour exactly (§ 2.1).
    */
  given none: SpectralBackend with
    def name: String = "none"
    def capabilities: Set[SpectralCapability] = Set.empty

  /** A composite whose [[SpectralBackend.capabilities]] is the '''union''' of the
    * parts, dispatching each operation to the '''first''' part whose `capabilities`
    * contains the needed capability (earlier parts win on overlap). This is how a
    * user combines several capability-disjoint engines: importing two `given`s is a
    * compile-time ambiguity, so instead they import the backend '''values''' and
    * declare one composite `given` (§ 2.2).
    *
    * Thread-safety follows from the parts' (G1): `compose` adds no mutable state.
    */
  def compose(primary: SpectralBackend, rest: SpectralBackend*): SpectralBackend =
    val parts = primary +: rest.toVector
    new SpectralBackend:
      def name: String = parts.map(_.name).mkString("compose(", ", ", ")")
      def capabilities: Set[SpectralCapability] = parts.foldLeft(Set.empty[SpectralCapability])(_ ++ _.capabilities)

      private def first(cap: SpectralCapability): Option[SpectralBackend] =
        parts.find(_.capabilities.contains(cap))

      override def denseSymmetricEigen(a: DMat, wantVectors: Boolean): Either[LinAlgError, RawSymmetricEigen] =
        first(SpectralCapability.DenseSymmetricEigen) match
          case Some(b) => b.denseSymmetricEigen(a, wantVectors)
          case None    => unsupported("dense symmetric eigen")

      override def denseNonsymmetricEigen(a: DMat, vectors: EigenVectors): Either[LinAlgError, RawNonsymmetricEigen] =
        first(SpectralCapability.DenseNonsymmetricEigen) match
          case Some(b) => b.denseNonsymmetricEigen(a, vectors)
          case None    => unsupported("dense nonsymmetric eigen")

      override def denseSvd(a: DMat, wantVectors: Boolean): Either[LinAlgError, RawSvd] =
        first(SpectralCapability.DenseSvd) match
          case Some(b) => b.denseSvd(a, wantVectors)
          case None    => unsupported("dense SVD")

      override def generalizedNonsymmetricEigen(
          a: DMat,
          b: DMat,
          vectors: EigenVectors
      ): Either[LinAlgError, RawGeneralizedEigen] =
        first(SpectralCapability.GeneralizedNonsymmetricEigen) match
          case Some(bk) => bk.generalizedNonsymmetricEigen(a, b, vectors)
          case None     => unsupported("generalized nonsymmetric eigen (QZ)")

      override def rankDeficientGsvd(a: DMat, b: DMat, wantVectors: Boolean): Either[LinAlgError, RawGsvd] =
        first(SpectralCapability.RankDeficientGsvd) match
          case Some(bk) => bk.rankDeficientGsvd(a, b, wantVectors)
          case None     => unsupported("rank-deficient GSVD")

      override def shiftInvertOperator(
          a: DMat,
          b: Option[DMat],
          sigma: Double
      ): Either[LinAlgError, DoubleLinearOperator] =
        first(SpectralCapability.ShiftInvertSolve) match
          case Some(bk) => bk.shiftInvertOperator(a, b, sigma)
          case None     => unsupported("shift-invert operator")
