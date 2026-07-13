package gale.syntax

import gale.backend.Backend
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix

/** Opt-in ergonomic syntax for `gale.linalg`, in focused import modules.
  *
  *   - [[all]] — the safe ASCII sugar: elementwise (`pointwise`) operators and
  *     `zipMap`, which the PRD's "Elementwise operations must be explicit" example
  *     assumes but core does not yet expose. Bring in with `import gale.syntax.all.*`.
  *   - [[unicode]] — the Unicode operator aliases `×` (matrix product) and `⋅`
  *     (dot), a '''separate''' opt-in (the PRD keeps symbol aliases out of the
  *     default surface). Bring in with `import gale.syntax.unicode.*`.
  *
  * Everything here is zero-cost (opaque wrappers / thin `*`-aliasing extensions),
  * adds no ambiguity with the core operators, the property wrappers, or the sized
  * layer under a combined import, and preserves the '''Either-first''' failure model
  * — the only throwing behavior is an elementwise shape mismatch, which matches the
  * core arithmetic (`DMat.+` / `DVec.+` already throw `DimensionMismatch` on a shape
  * mismatch; these are arithmetic primitives, not the Either-returning solve tier).
  *
  * `all` deliberately excludes the Unicode operators so a wildcard import stays
  * unsurprising; opt into those explicitly.
  */
object all:

  /** An elementwise (Hadamard) '''view''' of a matrix — `a.pointwise * b`,
    * `a.pointwise / b`, `a.pointwise.map(f)`. Opaque, so it never collides with the
    * core `DMat.*` (matrix product).
    */
  opaque type Pointwise = DMat

  extension (a: DMat)
    /** The elementwise view of `a`; combine with another matrix via `* ` / `/` or a
      * unary `map`.
      */
    def pointwise: Pointwise = a

    /** Elementwise combine of two '''same-shape''' matrices: `out(i,j) = f(a(i,j),
      * b(i,j))`. Throws `LinAlgError.DimensionMismatch` on a shape mismatch (as the
      * core `+`/`-` do).
      */
    def zipMap(b: DMat)(f: (Double, Double) => Double): DMat =
      zipShape(a, b, f)

  extension (p: Pointwise)
    /** Elementwise (Hadamard) product. */
    def *(b: DMat): DMat = zipShape(reveal(p), b, _ * _)

    /** Elementwise quotient. */
    def /(b: DMat): DMat = zipShape(reveal(p), b, _ / _)

    /** Apply `f` to every entry. */
    def map(f: Double => Double): DMat =
      val a = reveal(p)
      Matrix.tabulate(a.rows, a.cols)((i, j) => f(a(i, j)))

  private inline def reveal(p: Pointwise): DMat = p

  private def zipShape(a: DMat, b: DMat, f: (Double, Double) => Double): DMat =
    if a.rows != b.rows || a.cols != b.cols then throw LinAlgError.DimensionMismatch(a.shape, b.shape)
    Matrix.tabulate(a.rows, a.cols)((i, j) => f(a(i, j), b(i, j)))

/** Unicode operator aliases (PRD § "Symbol aliases … belong in opt-in syntax
  * modules"): `×` for the matrix product (matrix-vector and matrix-matrix) and `⋅`
  * for the inner product. Thin aliases of the core `*` / `dot`, so they are exact
  * synonyms with no new semantics. Kept out of [[all]] — import explicitly.
  */
object unicode:
  extension (a: DMat)
    /** Matrix-vector product (alias of `a * x`). */
    def ×(x: DVec): DVec = a * x

    /** Matrix-matrix product (alias of `a * b`). Forwards the ambient `given Backend` so the
      * alias stays a true synonym of `*` even when an accelerating backend is imported.
      */
    def ×(b: DMat)(using backend: Backend): DMat = a.*(b)(using backend)

  extension (x: DVec)
    /** Inner product (alias of `x dot y`). */
    def ⋅(y: DVec): Double = x.dot(y)
