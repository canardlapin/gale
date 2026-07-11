package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError

/** The result of a generalized singular value decomposition of a matrix pencil
  * `(A, B)` with `A` `m×n`, `B` `p×n` (§ 9 of `docs/spectral-parity.md`, and the
  * "Implications" item 12 struct):
  *
  * {{{
  *   A = U C Xᵀ ,   B = V S Xᵀ ,   CᵀC + SᵀS = I
  * }}}
  *
  * with generalized singular values `c_i / s_i`. `C = diag(c)` and `S = diag(s)`
  * are the cosine/sine diagonals (`c_i² + s_i² = 1`), so reconstruction is
  * `A ≈ u · diag(c) · x.t` and `B ≈ v · diag(s) · x.t` — note `x` holds '''X'''
  * (the reconstruction uses `x.t`, MATLAB's `[U,V,X,C,S]` orientation).
  *
  * '''Ordering.''' Entries are ordered '''descending by the generalized singular
  * value `c_i/s_i`''', with the [[GeneralizedSingularValue.Infinite]] values
  * (`s = 0`) '''first''' and the [[GeneralizedSingularValue.Zero]] values
  * (`c = 0`) '''last''' (§ 9 ordering table). `values` is the typed classification
  * aligned with `c`/`s`.
  *
  * '''Scope / vectors.''' v0.3.5 ships the '''pure full-column-rank''' pencil only
  * (`[A;B]` must have rank `n`); rank-deficient pencils are a `Left` (deferred to a
  * backend). `u` (`m×n`), `v` (`p×n`), and `x` (`n×n`) are empty when only values
  * were requested. `u` and `v` have '''orthonormal columns on the well-determined
  * entries only''' — a [[GeneralizedSingularValue.Zero]] leaves its `u` column
  * undetermined (stored as zero), an [[GeneralizedSingularValue.Infinite]] its `v`
  * column (this route cannot recover those from the Gram decomposition).
  */
final case class GeneralizedSVD private[gale] (
    u: DMat,
    v: DMat,
    x: DMat,
    c: DVec,
    s: DVec,
    values: IndexedSeq[GeneralizedSingularValue],
    diagnostics: SpectralDiagnostics
):
  /** The number of generalized singular values (`= n`). */
  def size: Int =
    c.length

  /** The `i`-th generalized singular value `c_i / s_i` as a `Double`
    * (`+∞` for [[GeneralizedSingularValue.Infinite]], `0` for
    * [[GeneralizedSingularValue.Zero]]).
    */
  def ratio(i: Int): Double =
    values(i).value

  /** `Right(this)` when the (dense one-shot) decomposition is complete. */
  def requireConverged: Either[LinAlgError, GeneralizedSVD] =
    diagnostics.requireConverged(this)
