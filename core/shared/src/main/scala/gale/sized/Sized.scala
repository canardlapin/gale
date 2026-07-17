package gale.sized

import gale.backend.Backend
import gale.linalg.Cols
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.linalg.Rows
import gale.linalg.Shape
import scala.annotation.targetName

/** An '''optional, zero-cost''' compile-time-dimension layer over the runtime-sized
  * `DVec`/`DMat` (PRD "Optional Sized Layer"). [[SVec]] and [[SMat]] are opaque
  * wrappers carrying the length / shape as Scala 3 literal-integer type parameters,
  * so the core operations (`+`, `-`, scale, `dot`, matrix-vector, matrix-matrix,
  * transpose) are '''shape-checked at compile time''' when the dimensions are
  * statically known, while lowering to the exact same `DVec`/`DMat` kernels.
  *
  * The layer allocates nothing (an opaque alias of `DVec`/`DMat`), is opt-in
  * (nothing in `gale.linalg` or the rest of core references it), and lowers
  * losslessly and '''zero-copy''' through [[SVec.toDVec]] / [[SMat.toDMat]] (the
  * identity on the wrapped value). Raising a runtime value is safe:
  * [[SVec.fromDVec]] / [[SMat.fromDMat]] return `Left` on a shape mismatch; the
  * trusted `unsafe` constructors are `private[sized]`.
  *
  * '''Naming.''' The PRD sketches `Matrix.Sized[A, R, C]` / `Vec.Sized[A, N]` inside
  * `object Matrix` / `object Vec`; those objects live in `gale.linalg`. This layer
  * ships in its own package `gale.sized` as [[SMat]]/[[SVec]] (Double-specialized,
  * matching gale's `DVec`/`DMat` op tier) so it stays fully additive and leaks
  * nothing into the core API.
  */
opaque type SVec[N <: Int] = DVec

/** A compile-time-sized matrix: `R` rows, `C` columns, over `DMat`. See [[SVec]]
  * for the layer's design (optional, zero-cost, opt-in).
  */
opaque type SMat[R <: Int, C <: Int] = DMat

object SVec:
  /** Wrap a `DVec` as an `N`-vector '''without''' a length check (trusted). */
  private[sized] inline def unsafe[N <: Int](v: DVec): SVec[N] = v

  /** Unwrap to the underlying `DVec` (for cross-type internals, e.g. matvec). */
  private[sized] inline def raw[N <: Int](v: SVec[N]): DVec = v

  /** Adopt a runtime `DVec` as a statically-sized `SVec[N]`, or `Left` if its length
    * is not `N` (zero-copy on success).
    */
  def fromDVec[N <: Int](v: DVec)(using n: ValueOf[N]): Either[LinAlgError, SVec[N]] =
    if v.length == n.value then Right(unsafe(v))
    else Left(LinAlgError.VectorLengthMismatch(n.value, v.length))

  /** Build an `SVec[N]` from exactly `N` values (`N` inferred from the expected
    * type). Throws `IllegalArgumentException` on an arity mismatch — a
    * programming error, not runtime data — so the size stays a compile-time fact.
    */
  def sized[N <: Int](values: Double*)(using n: ValueOf[N]): SVec[N] =
    require(values.length == n.value, s"expected ${n.value} values for SVec[${n.value}], got ${values.length}")
    unsafe(DVec.fromSeq(values.toIndexedSeq))

  /** The `N`-length zero vector. */
  def zeros[N <: Int](using n: ValueOf[N]): SVec[N] =
    unsafe(DVec.zeros(n.value))

  // Fixed-size factories (the parameter count IS the compile-time size check).
  def vec2(x: Double, y: Double): SVec[2] = unsafe(DVec.fromSeq(IndexedSeq(x, y)))
  def vec3(x: Double, y: Double, z: Double): SVec[3] = unsafe(DVec.fromSeq(IndexedSeq(x, y, z)))
  def vec4(x: Double, y: Double, z: Double, w: Double): SVec[4] = unsafe(DVec.fromSeq(IndexedSeq(x, y, z, w)))

  extension [N <: Int](v: SVec[N])
    /** Lower to the underlying `DVec` (zero-copy: the identity on the wrapper). */
    inline def toDVec: DVec = v

    /** The runtime length (`= N`). */
    inline def length: Int = v.toDVec.length

    inline def apply(i: Int): Double = v.toDVec(i)

    /** Elementwise sum with another `N`-vector (shape-checked: both are `SVec[N]`). */
    def +(that: SVec[N]): SVec[N] = unsafe(v.toDVec + that.toDVec)

    /** Elementwise difference with another `N`-vector. */
    def -(that: SVec[N]): SVec[N] = unsafe(v.toDVec - that.toDVec)

    /** Scale by a scalar. */
    def *(alpha: Double): SVec[N] = unsafe(v.toDVec * alpha)

    /** Inner product with another `N`-vector. */
    def dot(that: SVec[N]): Double = v.toDVec.dot(that.toDVec)

    /** Euclidean norm. */
    def norm2: Double = v.toDVec.norm2

object SMat:
  /** Wrap a `DMat` as an `R×C` matrix '''without''' a shape check (trusted). */
  private[sized] inline def unsafe[R <: Int, C <: Int](m: DMat): SMat[R, C] = m

  /** Adopt a runtime `DMat` as a statically-sized `SMat[R, C]`, or `Left` on a
    * shape mismatch (zero-copy on success).
    */
  def fromDMat[R <: Int, C <: Int](m: DMat)(using r: ValueOf[R], c: ValueOf[C]): Either[LinAlgError, SMat[R, C]] =
    if m.rows == r.value && m.cols == c.value then Right(unsafe(m))
    else Left(LinAlgError.DimensionMismatch(Shape(Rows(r.value), Cols(c.value)), m.shape))

  /** Build an `SMat[R, C]` from exactly `R·C` '''row-major''' values (`R`, `C`
    * inferred from the expected type). Throws on an arity mismatch.
    */
  def sized[R <: Int, C <: Int](values: Double*)(using r: ValueOf[R], c: ValueOf[C]): SMat[R, C] =
    require(
      values.length == r.value * c.value,
      s"expected ${r.value * c.value} values for SMat[${r.value}, ${c.value}], got ${values.length}"
    )
    unsafe(Matrix.dense(r.value, c.value, values.toIndexedSeq))

  /** The `R×C` zero matrix. */
  def zeros[R <: Int, C <: Int](using r: ValueOf[R], c: ValueOf[C]): SMat[R, C] =
    unsafe(Matrix.tabulate(r.value, c.value)((_, _) => 0.0))

  // Fixed-size row-major factories (the parameter count IS the compile-time size
  // check). These are the sized-layer counterpart to the field-based
  // `gale.linalg.Mat2`/`Mat3`/`Mat4`, adding `inverse` and general shape-checked ops.
  def mat2(a: Double, b: Double, c: Double, d: Double): SMat[2, 2] =
    unsafe(Matrix.dense(2, 2, IndexedSeq(a, b, c, d)))

  def mat3(
      a: Double, b: Double, c: Double,
      d: Double, e: Double, f: Double,
      g: Double, h: Double, i: Double
  ): SMat[3, 3] = unsafe(Matrix.dense(3, 3, IndexedSeq(a, b, c, d, e, f, g, h, i)))

  def mat4(
      a: Double, b: Double, c: Double, d: Double,
      e: Double, f: Double, g: Double, h: Double,
      i: Double, j: Double, k: Double, l: Double,
      m: Double, n: Double, o: Double, p: Double
  ): SMat[4, 4] = unsafe(Matrix.dense(4, 4, IndexedSeq(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p)))

  extension [R <: Int, C <: Int](m: SMat[R, C])
    /** Lower to the underlying `DMat` (zero-copy). */
    inline def toDMat: DMat = m

    inline def rows: Int = m.toDMat.rows
    inline def cols: Int = m.toDMat.cols
    inline def apply(i: Int, j: Int): Double = m.toDMat(i, j)

    /** Elementwise sum with another `R×C` matrix. */
    def +(that: SMat[R, C]): SMat[R, C] = unsafe(m.toDMat + that.toDMat)

    /** Elementwise difference. */
    def -(that: SMat[R, C]): SMat[R, C] = unsafe(m.toDMat - that.toDMat)

    /** Scale by a scalar. */
    def scale(alpha: Double): SMat[R, C] =
      val d = m.toDMat
      unsafe(Matrix.tabulate(d.rows, d.cols)((i, j) => d(i, j) * alpha))

    /** Transpose (`R×C ⇒ C×R`). */
    def t: SMat[C, R] = unsafe(m.toDMat.t)

    /** Matrix-vector product (shape-checked: the vector must be `SVec[C]`). */
    def *(x: SVec[C])(using backend: Backend): SVec[R] =
      SVec.unsafe(m.toDMat.*(SVec.raw(x))(using backend))

    /** Matrix-matrix product (shape-checked: the right factor's rows must be `C`). Forwards
      * the ambient `given Backend` to the underlying `DMat.*` so an imported accelerating
      * backend routes through the sized layer too (it would otherwise capture `Backend.pure`).
      */
    def *[C2 <: Int](that: SMat[C, C2])(using backend: Backend): SMat[R, C2] =
      unsafe(m.toDMat.*(that.toDMat)(using backend))

  // --- tiny fixed-size kernels: closed-form det / inverse (PRD § tiny kernels) ---

  extension (m: SMat[2, 2])
    /** Determinant of a `2×2` (unrolled). */
    @targetName("det2") def det: Double =
      m.toDMat(0, 0) * m.toDMat(1, 1) - m.toDMat(0, 1) * m.toDMat(1, 0)

    /** Inverse of a `2×2`, or `Left(SingularMatrix)` when the determinant is zero /
      * non-finite (unrolled adjugate).
      *
      * '''Conditioning caveat.''' The singularity test is '''exact''' (`det == 0`
      * or non-finite): a merely '''near'''-singular matrix returns `Right` with an
      * ill-conditioned, large-magnitude inverse in which roundoff is amplified by
      * `1/det`. This closed-form kernel does no pivoting or condition estimation —
      * for a numerically robust solve of an ill-conditioned system, lower to `DMat`
      * and use gale's LU/QR (`toDMat.solve` / `toDMat.leastSquares`).
      */
    @targetName("inverse2") def inverse: Either[LinAlgError, SMat[2, 2]] =
      val a = m.toDMat(0, 0)
      val b = m.toDMat(0, 1)
      val c = m.toDMat(1, 0)
      val dd = m.toDMat(1, 1)
      val det = a * dd - b * c
      if det == 0.0 || !det.isFinite then Left(LinAlgError.SingularMatrix(0))
      else
        val inv = 1.0 / det
        Right(unsafe(Matrix.dense(2, 2, IndexedSeq(dd * inv, -b * inv, -c * inv, a * inv))))

  extension (m: SMat[3, 3])
    /** Determinant of a `3×3` (unrolled cofactor expansion). */
    @targetName("det3") def det: Double =
      val a = m.toDMat(0, 0)
      val b = m.toDMat(0, 1)
      val c = m.toDMat(0, 2)
      val d = m.toDMat(1, 0)
      val e = m.toDMat(1, 1)
      val f = m.toDMat(1, 2)
      val g = m.toDMat(2, 0)
      val h = m.toDMat(2, 1)
      val i = m.toDMat(2, 2)
      a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

    /** Inverse of a `3×3`, or `Left(SingularMatrix)` when singular (unrolled
      * adjugate transpose divided by the determinant).
      *
      * '''Conditioning caveat.''' As with [[SMat.inverse2 the 2×2 case]], the
      * singularity test is '''exact''' (`det == 0` or non-finite); a near-singular
      * matrix yields a `Right` with a large-magnitude, roundoff-amplified inverse.
      * For ill-conditioned data, lower to `DMat` and solve via gale's LU/QR.
      */
    @targetName("inverse3") def inverse: Either[LinAlgError, SMat[3, 3]] =
      val a = m.toDMat(0, 0)
      val b = m.toDMat(0, 1)
      val c = m.toDMat(0, 2)
      val d = m.toDMat(1, 0)
      val e = m.toDMat(1, 1)
      val f = m.toDMat(1, 2)
      val g = m.toDMat(2, 0)
      val h = m.toDMat(2, 1)
      val i = m.toDMat(2, 2)
      val det3 = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
      if det3 == 0.0 || !det3.isFinite then Left(LinAlgError.SingularMatrix(0))
      else
        val inv = 1.0 / det3
        // Inverse = adjugate(m)ᵀ / det; cofactor Cᵢⱼ placed at (j, i).
        val values = IndexedSeq(
          (e * i - f * h) * inv, (c * h - b * i) * inv, (b * f - c * e) * inv,
          (f * g - d * i) * inv, (a * i - c * g) * inv, (c * d - a * f) * inv,
          (d * h - e * g) * inv, (b * g - a * h) * inv, (a * e - b * d) * inv
        )
        Right(unsafe(Matrix.dense(3, 3, values)))
