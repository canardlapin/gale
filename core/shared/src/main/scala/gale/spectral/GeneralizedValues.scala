package gale.spectral

/** A generalized eigenvalue in projective `(α, β)` form, `λ = α/β`.
  *
  * The QZ / `ggev` family reports eigenvalues of the pencil `A x = λ B x` as a
  * ratio `α/β` rather than a single number, because a singular or rank-deficient
  * `B` yields '''infinite''' eigenvalues (`β = 0`) that a plain [[Complex]] cannot
  * represent (§ 5 of `docs/spectral-parity.md`). `α` is complex; `β` is real and,
  * by the QZ convention, is '''exactly zero''' for an infinite eigenvalue (a
  * producer snaps a negligible `β` to `0.0`).
  */
final case class GeneralizedEigenvalue(alpha: Complex, beta: Double):
  /** The eigenvalue `α/β`. When [[isInfinite]] this has non-finite components
    * (`±∞`/`NaN`); test [[isInfinite]] first.
    */
  def value: Complex =
    alpha / beta

  /** True when `β` is exactly zero — the pencil has an infinite eigenvalue at
    * this index (singular/rank-deficient `B`).
    */
  def isInfinite: Boolean =
    beta == 0.0

  /** True when the pair denotes a finite eigenvalue (`β ≠ 0`); the complement of
    * [[isInfinite]].
    */
  def isFinitePair: Boolean =
    beta != 0.0

/** The classification of a generalized singular value `c/s` from a GSVD
  * (§ 9 of `docs/spectral-parity.md`).
  *
  * The cosine–sine pair `(c, s)` with `c² + s² = 1` degenerates at the ends:
  * `s = 0` gives an '''infinite''' value and `c = 0` a '''zero''' value, both of
  * which the PRD requires be represented explicitly rather than as `±∞`/`0`
  * folded into a `Double`.
  */
enum GeneralizedSingularValue:
  /** A finite ratio `c/s`, with both `c` and `s` nonzero. */
  case Finite(ratio: Double)

  /** `s = 0` (`c = 1`): the generalized singular value is `+∞`. */
  case Infinite

  /** `c = 0` (`s = 1`): the generalized singular value is `0`. */
  case Zero

  /** The numeric value: the ratio for [[Finite]], `+∞` for [[Infinite]], `0` for
    * [[Zero]].
    */
  def value: Double =
    this match
      case Finite(ratio) => ratio
      case Infinite      => Double.PositiveInfinity
      case Zero          => 0.0
