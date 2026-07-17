# Gale PRD: Scala 3 Linear Algebra Core

## Status

Draft product requirements for `gale`, a new Scala 3-first linear algebra
library. Gale is not Atlas, not a Breeze fork, and not a rename of the existing
`neuroim-linalg` fitting helper. It is a standalone, general-purpose algebraic
linear algebra core.

## Evidence Base

The design is grounded in the following external state as of July 2026:

- Breeze says it is "mostly retired" and lists 2.1.0 as its latest release,
  cross-built for Scala 3.1, 2.12, and 2.13:
  <https://github.com/scalanlp/breeze>
- Scala currently has two Scala 3 lines: Scala Next 3.8.4 and Scala LTS 3.3.8,
  with LTS recommended for published libraries:
  <https://www.scala-lang.org/download/>
- Scala.js cross-building is documented around `sbt-scalajs` 1.22.0 and
  `sbt-scalajs-crossproject` 1.3.2:
  <https://www.scala-js.org/doc/project/cross-build.html>
- Scala.js considers its WebAssembly backend stable since 1.22.0:
  <https://www.scala-js.org/doc/project/webassembly.html>
- JavaScript `Float64Array` is a standard 64-bit floating-point typed-array
  storage target:
  <https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Float64Array>
- JDK 26 delivered the Vector API as an eleventh incubator:
  <https://openjdk.org/jeps/529>
- JDK 22 finalized the Foreign Function & Memory API for native interop without
  JNI:
  <https://openjdk.org/jeps/454>
- Netlib describes BLAS as Level 1 scalar/vector/vector-vector operations,
  Level 2 matrix-vector operations, and Level 3 matrix-matrix operations:
  <https://www.netlib.org/blas/>
- Netlib's LAPACK guide states that LAPACK performs as much work as possible
  through BLAS calls, and that efficient machine-specific BLAS implementations
  are what let portable LAPACK routines achieve high performance:
  <https://www.netlib.org/lapack/lug/node11.html>
- OpenJDK describes JMH as the benchmark harness for building, running, and
  analyzing benchmarks on the JVM:
  <https://openjdk.org/projects/code-tools/jmh/>
- Oracle's Java 22 FFM guide documents foreign functions, foreign memory,
  memory segments, and arenas as the supported native interop model:
  <https://docs.oracle.com/en/java/javase/22/core/foreign-function-and-memory-api.html>
- MATLAB's spectral APIs cover ordinary and generalized `eig`, partial `eigs`,
  partial `svds`, and GSVD. Gale should target this capability space without
  copying MATLAB's syntax:
  <https://www.mathworks.com/help/matlab/ref/eig.html>
  <https://www.mathworks.com/help/matlab/ref/eigs.html>
  <https://www.mathworks.com/help/matlab/ref/svds.html>
  <https://www.mathworks.com/help/matlab/ref/gsvd.html>
- SciPy exposes dense generalized `eig`/`eigh`, sparse partial `eigs`/`eigsh`,
  and partial `svds`, including left/right vectors, selection modes,
  shift-invert-style targeting, and convergence failure semantics:
  <https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html>
  <https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html>
  <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html>
  <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html>
  <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html>

## Product Thesis

Gale should provide an elegant, functional, extremely well-typed linear algebra
library for Scala 3, with the same public API on the JVM, Scala.js, and
Scala.js Wasm, while keeping platform storage and hot-loop kernels private.

The public surface should feel mathematical:

```scala
import gale.linalg.*

val A = Matrix.dense(
  rows = 3,
  cols = 3,
  values = Seq(
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0,
    7.0, 8.0, 10.0
  )
)

val b = Vec(3.0, 3.0, 4.0)
val x = A.solve(b).orThrow
val r = A * x - b
val ok = r.norm2 < 1e-10
```

The implementation should remain storage-and-kernel shaped:

```scala
kernel.gemv(
  rows = A.rows.value,
  cols = A.cols.value,
  alpha = 1.0,
  a = A.data,
  aOffset = A.offset.value,
  rowStride = A.rowStride.value,
  colStride = A.colStride.value,
  x = x.data,
  xOffset = x.offset.value,
  xStride = x.stride.value,
  beta = -1.0,
  y = b.data,
  yOffset = b.offset.value,
  yStride = b.stride.value
)
```

The product bet is that Scala 3 can hide the second form behind the first
without turning the library into macro-heavy type art.

## Design Doctrine

Gale should be Oderskyan in the useful sense:

- Small core, large composability.
- Simple code remains simple.
- Power is available through typed escape hatches, not imposed everywhere.
- Public APIs model linear algebra.
- Internal APIs model layout, strides, ownership, and kernels.
- Generic interfaces describe capabilities; primitive concrete types do the
  high-throughput work.
- Type-level machinery protects real invariants, but runtime dimensions remain
  the default for normal numerical work.
- Functional values are the public default; local mutation is an explicit,
  disciplined implementation technique.

## Performance Doctrine

Gale's performance rule is:

> Generic at the API boundary; monomorphic and primitive inside every hot path.

The public interfaces may be generic:

```scala
trait Vec[A]
trait Matrix[A]
trait SparseMatrix[A]
```

The real dense and sparse implementations must be specialized by scalar family:

```scala
final class DVec extends Vec[Double]
final class FVec extends Vec[Float]
final class IVec extends Vec[Int]

final class DMat extends Matrix[Double]
final class FMat extends Matrix[Float]

final class CSRDouble extends SparseMatrix[Double]
final class CSRFloat extends SparseMatrix[Float]
```

Generic containers are not the fast path. This API should exist as a correct,
lawful fallback:

```scala
def dot[A](x: Vec[A], y: Vec[A])(using A: Field[A]): A
```

But primitive operations must dispatch to concrete kernels:

```scala
extension (x: DVec)
  def dot(y: DVec): Double =
    DoubleKernels.dot(x, y)

extension (A: DMat)
  def *(x: DVec): DVec =
    DoubleKernels.gemv(A, x)
```

Constructors must make the common numerical path fast by default:

```scala
Vec(1.0, 2.0, 3.0)           // DVec
Vec.float(1.0f, 2.0f, 3.0f)  // FVec
Vec.int(1, 2, 3)             // IVec
Vec.generic(BigDecimal(1))   // GVec[BigDecimal]
```

Do not rely on compiler specialization as the foundation. Gale should use
manual, explicit implementations for the primitive scalar families it supports.
That is less magical, but it is predictable across JVM, Scala.js, and Wasm.

Type classes describe capabilities and drive generic algorithms, exact
arithmetic, semiring sparse operations, construction, and algorithm selection.
They must not perform per-element arithmetic inside primitive dense kernels.

Acceptable operation-level dispatch:

```scala
given DenseOps[Double] with
  def dot(x: Vec[Double], y: Vec[Double]): Double =
    (x, y) match
      case (xd: DVec, yd: DVec) => DoubleKernels.dot(xd, yd)
      case _                   => GenericKernels.dot(x, y)
```

Unacceptable primitive kernel shape:

```scala
while i < n do
  acc = scalar.plus(acc, scalar.times(x(i), y(i)))
  i += 1
```

Acceptable primitive kernel shape:

```scala
while i < n do
  acc += xData(xi) * yData(yi)
  xi += xStride
  yi += yStride
  i += 1
```

Element access through `Vec[A]` is a convenience API, not a performance
contract. Hot loops must operate on validated primitive storage, offsets, and
strides.

`inline` is useful at the edges: constructors, tiny extension methods,
compile-time scalar dispatch, sized literals, and representation selection.
Macros are allowed only where they remove real runtime work: checked literals,
small fixed-size unrolled kernels, fused BLAS-shaped expressions, and generated
boilerplate for primitive kernel families. Macros must not become the foundation
of the algebra system or the solver stack.

Complex numbers need specialized primitive storage. A fast complex vector must
not be `Array[Complex]` when `Complex` is an object value. It should use packed
or structure-of-arrays primitive storage and materialize `Complex` only at the
API boundary.

Sparse matrices follow the same rule. Builder entries may be rich values, but
runtime CSR/CSC formats must store indices and primitive values in platform
arrays for primitive scalar families.

The project boxing policy is:

- Primitive dense and sparse kernels must not box per element.
- Primitive vector and matrix storage must use primitive platform arrays.
- Generic `Vec[A]` and `Matrix[A]` are allowed to box.
- Generic operations may dispatch once per operation.
- Element access through `Vec[A]` is not a performance contract.
- All high-throughput algorithms must have `Double` and `Float` fast paths.
- Complex numbers must have specialized primitive storage.
- Sparse primitive formats must store values in primitive arrays.
- No `Iterator`, `Seq`, `FunctionN`, collection combinators, or type-class
  arithmetic inside primitive kernels.

## Goals

1. Provide a Scala 3-first dense and sparse linear algebra core.
2. Cross-compile the core to JVM and Scala.js from day one.
3. Make Scala.js typed arrays and Wasm compatibility first-class, not fallback
   afterthoughts.
4. Keep the core pure Scala with no native, Breeze, Cats Effect, ZIO, Spire, or
   plotting dependency.
5. Offer specialized Double and Float fast paths, while still supporting a
   principled algebraic API for generic scalar work.
6. Make solver and factorization results typed, diagnostic, and reusable.
7. Expose allocation control through explicit destination APIs and limited,
   predictable fused operations.
8. Keep backend selection explicit through `given` imports.
9. Ship laws, conformance tests, and numerical regression tests as product
   artifacts, not as an afterthought.

## Non-Goals

- Do not modernize Breeze directly.
- Do not build a scientific-computing umbrella in v1.
- Do not include plotting, probability distributions, signal processing, or
  tensors in the core.
- Do not require native BLAS, LAPACK, SuiteSparse, Spire, Cats, ZIO, or Breeze.
- Do not cross-build for Scala 2 unless a concrete migration requirement appears.
- Do not make sparse `A.solve(b)` silently pick a major algorithm.
- Do not expose platform storage types as the public data model.
- Do not use macros, match types, or compile-time arithmetic as the foundation
  of ordinary matrix use.

## Target Users

- Scala developers who want a pleasant linear algebra surface without Breeze's
  maintenance risk.
- Library authors who need a pure, cross-platform numerical core.
- Browser and Node users who need meaningful numerical kernels through
  Scala.js, typed arrays, and Wasm.
- JVM numerical users who want optional BLAS, Vector API, or SuiteSparse
  acceleration without leaking those choices into ordinary code.
- Educators and researchers who benefit from algebraic clarity and readable
  examples, but still need predictable performance.

## Package Shape

The first public namespace is:

```scala
package gale.linalg
```

Recommended modules:

```text
gale-core                 JVM + JS
  Vec, Matrix, DenseMatrix, SparseMatrix, LinearOperator
  Shape, layout, slices, views, scalar type classes
  Pure algorithms and user-facing syntax

gale-platform             JVM + JS
  Private platform storage abstractions
  JVM: Array[Double], Array[Float], Array[Int]
  JS: Float64Array, Float32Array, Int32Array

gale-kernel               JVM + JS
  Low-level dot, axpy, scal, gemv, gemm, triangular solve
  Specialized Double and Float kernels

gale-factorization        JVM + JS
  LU, Cholesky, QR, triangular solve, least squares

gale-sparse               JVM + JS
  COO, CSR, CSC, diagonal, banded, identity, permutation
  SpMV, sparse addition, sparse scaling, sparse-dense multiply

gale-solvers              JVM + JS
  CG, MINRES, BiCGSTAB, restarted GMRES, LSQR, LSMR
  Preconditioner API and initial Jacobi / block-Jacobi support

gale-laws                 JVM + JS
  MUnit and ScalaCheck laws, backend conformance tests

gale-backend-jvm-vector   JVM only
  Optional JDK Vector API kernels

gale-backend-jvm-blas-ffm JVM only
  Optional BLAS/LAPACK backend through Java FFM

gale-backend-jvm-sparse   JVM only, post-v1
  Optional SuiteSparse backend boundary

gale-interop-breeze       JVM only
  Breeze conversion helpers and migration aids
```

## Core Domain Model

The public hierarchy should distinguish operators from matrices:

```scala
trait LinearOperator[A]:
  def rows: Rows
  def cols: Cols
  def applyTo(x: Vec[A], into: MutableVec[A]): Either[LinAlgError, Unit]

trait Matrix[A] extends LinearOperator[A]:
  def apply(row: RowIndex, col: ColIndex): A
  def row(i: RowIndex): Vec[A]
  def col(j: ColIndex): Vec[A]
  def t: Matrix[A]

trait DenseMatrix[A] extends Matrix[A]

trait SparseMatrix[A] extends Matrix[A]:
  def nnz: NonZeroCount
  def density: Double
```

Dense concrete types should be final and storage-aware:

```scala
final class DVec private[gale] (
  private[gale] val data: DoubleArray,
  val offset: Offset,
  val length: Length,
  val stride: Stride
)

final class DMat private[gale] (
  private[gale] val data: DoubleArray,
  val offset: Offset,
  val rows: Rows,
  val cols: Cols,
  val rowStride: Stride,
  val colStride: Stride,
  val layout: Layout
)
```

Sparse concrete types should make construction and canonicalization explicit:

```scala
final class COO[A] private[gale] (...)
final class CSR[A] private[gale] (...)
final class CSC[A] private[gale] (...)
```

Required sparse construction policy:

```scala
enum DuplicatePolicy:
  case Sum
  case Last
  case Error

val A =
  Sparse
    .coo[Double](rows = 5, cols = 5)
    .add(0, 0, 10.0)
    .add(0, 3, 2.0)
    .add(0, 3, 5.0)
    .toCSR(duplicates = DuplicatePolicy.Sum)
```

Sparse matrices must expose:

```scala
A.canonicalize
A.pruneZeros
A.prune(absBelow = 1e-12)
A.sortedIndices
A.hasCanonicalFormat
```

## Platform Storage Requirement

Gale must not expose `Array[Double]` as the public storage contract.

Shared code depends on private platform abstractions:

```scala
package gale.platform

opaque type DoubleArray
opaque type FloatArray
opaque type IndexArray
```

JVM implementation:

```scala
opaque type DoubleArray = Array[Double]
opaque type IndexArray = Array[Int]
```

Scala.js implementation:

```scala
opaque type DoubleArray = scala.scalajs.js.typedarray.Float64Array
opaque type IndexArray = scala.scalajs.js.typedarray.Int32Array
```

Public platform interop is copy-only so it does not become a storage contract:

```scala
// JVM only
Matrix.fromArrayCopy(rows, cols, values: Array[Double])

// JS only
Matrix.fromFloat64ArrayCopy(rows, cols, values: Float64Array)
```

## Type System Requirements

Gale should use Scala 3 types to protect real invariants.

### Opaque Domain Values

Required opaque value types:

```scala
opaque type Rows = Int
opaque type Cols = Int
opaque type Length = Int
opaque type Offset = Int
opaque type Stride = Int
opaque type RowIndex = Int
opaque type ColIndex = Int
opaque type NonZeroCount = Int
```

All constructors validate non-negativity, range, and overflow where applicable.
Kernel code may unwrap after validated boundaries.

### Optional Sized Layer

Runtime dimensions are the default. Statically sized matrices are an optional
zero-overhead wrapper:

```scala
object Matrix:
  opaque type Sized[A, R <: Int, C <: Int] = Matrix[A]

object Vec:
  opaque type Sized[A, N <: Int] = Vec[A]
```

Required behavior:

```scala
val A: Matrix.Sized[Double, 3, 4] = Matrix.sized[Double, 3, 4](...)
val x: Vec.Sized[Double, 4] = Vec.sized(1.0, 2.0, 3.0, 4.0)
val y: Vec.Sized[Double, 3] = A * x
```

The sized layer must not make ordinary `Matrix[Double]` usage noisy.

### Matrix Property Wrappers

Avoid boolean soup in solver APIs. Use explicit wrappers:

```scala
opaque type Symmetric[A] = A
opaque type PositiveDefinite[A] = A
opaque type LowerTriangular[A] = A
opaque type UpperTriangular[A] = A
opaque type CanonicalSparse[A] = A
```

Required API style:

```scala
val spd =
  A
    .assumeSymmetric
    .assumePositiveDefinite

val x = cg(spd, b, tolerance = 1e-10).orThrow
```

Checked alternatives must also exist:

```scala
A.verifySymmetric(tolerance = 1e-12)
A.verifyPositiveDefinite
A.verifyCanonicalSparse
```

`assume*` is explicit and cheap. `verify*` may be expensive and returns
`Either[LinAlgError, Wrapped[A]]`.

## Scalar API

Gale should not use `scala.math.Numeric` as its core numerical abstraction.

Required scalar hierarchy:

```scala
trait Semiring[A]:
  def zero: A
  def one: A
  def plus(x: A, y: A): A
  def times(x: A, y: A): A

trait Ring[A] extends Semiring[A]:
  def negate(x: A): A
  def minus(x: A, y: A): A = plus(x, negate(y))

trait Field[A] extends Ring[A]:
  def div(x: A, y: A): A

trait Real[A] extends Field[A]:
  def abs(x: A): A
  def sqrt(x: A): A
  def isNaN(x: A): Boolean
  def epsilon: A
```

Specialized kernels for `Double` and `Float` must bypass this abstraction in hot
paths. Type classes support generic algorithms, laws, exact or symbolic scalar
use, and optional interop.

Spire integration is allowed as an optional adapter, not a core dependency.

## Public API Requirements

Basic dense usage:

```scala
import gale.linalg.*

val x = Vec(1.0, 2.0, 3.0)
val y = Vec.fill(3)(1.0)
val z = x + 2.0 * y
val d = x dot z

val A = Matrix.dense(2, 3)(
  1.0, 2.0, 3.0,
  4.0, 5.0, 6.0
)

val b = Vec(1.0, 1.0)
val c = A.t * b
```

Required operators:

```scala
A * x          // matrix-vector
A * B          // matrix-matrix
x dot y
x + y
x - y
2.0 * x
A.t
A.solve(b)
A.lu
A.cholesky
A.qr
```

Elementwise operations must be explicit:

```scala
A.pointwise * B
A.pointwise / B
A.pointwise.map(math.exp)
A.zipMap(B)(_ + _)
```

Symbol aliases, including Unicode, belong in opt-in syntax modules:

```scala
import gale.linalg.syntax.unicode.*

val y = A × x
val d = x ⋅ y
```

## Functional Mutation Model

Gale is functional at the boundary and pragmatic in kernels.

Tier 1: pure, allocating:

```scala
val y = A * x + b
```

Tier 2: explicit destination:

```scala
val y = MutableVec.zeros[Double](A.rows)
A.mulInto(x, y)
y += b
```

Tier 3: limited fused operations:

```scala
y := alpha * (A * x) + beta * y
C := alpha * (A * B) + beta * C
```

The fused expression layer must be deliberately small and only recognize
BLAS-shaped patterns. It must not become a broad lazy expression graph.

Builders provide local mutation with immutable outputs:

```scala
val A =
  Matrix.build[Double](rows, cols): m =>
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        m(i, j) = f(i, j)
        j += 1
      i += 1
```

## Factorization Requirements

Factorizations are public values, not tuples:

```scala
final case class LU[A](
  packed: Matrix[A],
  pivots: IndexVec,
  parity: Int,
  diagnostics: FactorizationDiagnostics
):
  def solve(b: Vec[A]): Either[LinAlgError, Vec[A]]
  def det: Either[LinAlgError, A]

final case class Cholesky[A](
  lower: Matrix[A],
  diagnostics: FactorizationDiagnostics
):
  def solve(b: Vec[A]): Either[LinAlgError, Vec[A]]

final case class QR[A](
  q: Matrix[A],
  r: Matrix[A],
  diagnostics: FactorizationDiagnostics
):
  def solveLeastSquares(b: Vec[A]): Either[LinAlgError, Vec[A]]
```

Dense v1 must include:

- LU with partial pivoting.
- Cholesky for positive-definite matrices.
- QR via Householder reflections.
- Triangular solve.
- `solve` and `trySolve`.
- Least squares through QR.
- Determinant through LU.
- Rank and condition estimates.

Spectral methods should be a distinct phase rather than an incidental extension
of dense decompositions. They have different algorithmic and API risks from LU,
QR, and Cholesky.

Required spectral surface:

- Partial symmetric eigendecomposition for top and bottom eigenvalues.
- Partial nonsymmetric eigendecomposition for selected eigenvalues.
- Partial SVD for largest and smallest singular values.
- Generalized symmetric-definite eigenproblems, `A x = lambda B x`.
- Generalized nonsymmetric eigenproblems, where available through a backend.
- Generalized SVD for matrix pairs, with ordered generalized singular values.
- Residual and orthogonality diagnostics for every returned factorization.
- Explicit ordering schemes rather than ambiguous "top" APIs.

Ordering modes:

```scala
enum EigenOrder:
  case LargestMagnitude
  case SmallestMagnitude
  case LargestAlgebraic
  case SmallestAlgebraic
  case LargestRealPart
  case SmallestRealPart

enum SingularOrder:
  case Largest
  case Smallest
```

For symmetric real problems, `LargestAlgebraic` and `SmallestAlgebraic` should
mean the usual algebraic ordering on the real line. For nonsymmetric problems,
`LargestMagnitude` and real-part ordering should be explicit. For SVD and GSVD,
top means largest singular or generalized singular values; bottom means
smallest positive singular or generalized singular values, with rank-deficient
zeros reported explicitly.

The capability target is MATLAB/SciPy parity in the spectral space, expressed
through Gale's typed API rather than cloned syntax. This means:

- Dense ordinary eigenproblems with eigenvalues and optionally left and right
  eigenvectors.
- Dense generalized eigenproblems, including symmetric-definite problems and a
  backend boundary for QZ-style nonsymmetric generalized problems.
- Symmetric/Hermitian specializations that expose safer ordering and subset
  semantics than the general nonsymmetric API.
- Partial sparse and matrix-free eigenproblems equivalent in scope to
  `eigs`/`eigsh`: largest and smallest magnitude, largest and smallest
  algebraic values for symmetric problems, largest and smallest real parts for
  nonsymmetric problems, both-ends selection where meaningful, and target-near
  or shift-invert modes when the required linear solve is explicit.
- Partial SVD equivalent in scope to `svds`: largest and smallest singular
  values, caller-controlled subspace dimension, start vector, tolerance, maximum
  iterations, and explicit non-convergence results.
- GSVD for matrix pairs with the standard shape `A = U*C*X.t`,
  `B = V*S*X.t`, `C.t*C + S.t*S = I`, generalized singular values `c/s`, and
  explicit zero, infinite, and rank-deficient cases.

Stringly mode flags such as `LM`, `SM`, `LA`, `SA`, or `"largestabs"` should be
translated into typed options:

```scala
final case class SpectralSelection(
  count: Int,
  order: EigenOrder,
  target: Option[SpectralTarget] = None
)

enum SpectralTarget:
  case Around(value: Double)
  case ShiftInvert(sigma: Double, solver: LinearSolvePlan)

final case class SpectralOptions(
  tolerance: Double = 1e-10,
  maxIterations: Int = 1000,
  subspaceDimension: Option[Int] = None,
  startVector: Option[DVec] = None,
  returnVectors: EigenVectors = EigenVectors.Right
)

enum EigenVectors:
  case ValuesOnly
  case Right
  case Left
  case LeftAndRight
```

Suggested API:

```scala
val top = A.eigen(k = 5, order = EigenOrder.LargestAlgebraic)
val bot = A.eigen(k = 5, order = EigenOrder.SmallestAlgebraic)
val sv  = A.svd(k = 10, order = SingularOrder.Largest)

val gev = generalizedEigen(A, B, k = 6, order = EigenOrder.SmallestAlgebraic)
val gsv = generalizedSVD(A, B, k = 6, order = SingularOrder.Largest)
```

Result values should be first-class:

```scala
final case class EigenDecomposition(
  eigenvalues: DVec,
  eigenvectors: DMat,
  diagnostics: SpectralDiagnostics
)

final case class SVD(
  singularValues: DVec,
  u: DMat,
  vt: DMat,
  diagnostics: SpectralDiagnostics
)
```

Pure Scala algorithms should start with moderate-size, portable correctness:
Lanczos for symmetric partial eigenpairs, Arnoldi for nonsymmetric partial
eigenpairs, randomized or Lanczos bidiagonalization for partial SVD, and
inverse/shift-invert modes only when the required linear solves are explicit.
Production-scale spectral algorithms should be allowed to route to optional
native LAPACK/ARPACK/SLEPc-style backends later.

Later:

- LDLT.
- Full dense SVD.
- Full symmetric eigendecomposition.
- Full general eigendecomposition.
- Schur decomposition.

## Sparse Requirements

Sparse support is first-class but not overpromised.

Required v1 formats:

- COO for construction.
- CSR for row-oriented SpMV.
- CSC for transpose multiply and factorization-friendly column operations.
- Diagonal.
- Banded.
- Identity.
- Zero.
- Permutation.

Required v1 sparse operations:

- `A * x`.
- `A.t * x`.
- `A * B` for sparse-dense and sparse-sparse where practical.
- `A + B`, `A - B`, `alpha * A`.
- `A.mapValues(f)`.
- `A.zipValues(B)(f)`.
- `A.transpose`.
- `A.toCSR`, `A.toCSC`.
- `A.toDense` with size guard.
- `A.row(i)`, `A.col(j)`.
- `A.diagonal`, `A.trace`.
- `A.nnz`, `A.density`.
- `A.mulInto(x, y)` and `A.tMulInto(x, y)`.

Sparse direct solvers are not a v1 portability promise. The portable sparse
module should prioritize iterative methods and diagnostics.

Required v1 iterative solvers:

- CG for symmetric positive-definite operators.
- MINRES for symmetric indefinite operators, if implementation quality is high.
- BiCGSTAB for nonsymmetric systems.
- Restarted GMRES.
- LSQR / LSMR for least squares.
- Power iteration and Lanczos estimates if test infrastructure is ready.

All solvers return diagnostics:

```scala
sealed trait SolverResult[+A]:
  def iterations: Int
  def residualNorm: Double
  def value: A
  def orThrow: A

object SolverResult:
  final case class Converged[A](
    value: A,
    iterations: Int,
    residualNorm: Double
  ) extends SolverResult[A]

  final case class NotConverged[A](
    value: A,
    iterations: Int,
    residualNorm: Double
  ) extends SolverResult[A]
```

Preconditioners are first-class:

```scala
trait Preconditioner[A]:
  def solve(r: Vec[A], into: MutableVec[A]): Either[LinAlgError, Unit]
```

Initial preconditioners:

- Identity.
- Jacobi.
- Block Jacobi.
- Symmetric Gauss-Seidel.

Post-v1:

- ILU0.
- Incomplete Cholesky 0.
- ILUT.
- Algebraic multigrid.
- JVM-native SuiteSparse-backed sparse direct solvers.

## Backend Requirements

Backend selection is explicit:

```scala
import gale.backend.pure.given
// or
import gale.backend.jvm.blas.given
// or
import gale.backend.jvm.vector.given
```

Required backend model:

```scala
enum Capability:
  case NativeBlas
  case NativeSparse
  case Vectorized
  case Multithreaded
  case WasmFriendly
  case Deterministic

trait Backend:
  def name: String
  def capabilities: Set[Capability]
  def denseDouble: DenseDoubleKernel
```

Required kernel boundary:

```scala
trait DenseDoubleKernel:
  def dot(
    n: Int,
    x: DoubleArray,
    xOffset: Int,
    xStride: Int,
    y: DoubleArray,
    yOffset: Int,
    yStride: Int
  ): Double

  def axpy(...): Unit
  def scal(...): Unit
  def gemv(...): Unit
  def gemm(...): Unit
  def triangularSolve(...): Either[LinAlgError, Unit]
```

Core kernels must be pure Scala and cross-platform. JVM acceleration modules may
use the Vector API or FFM. Scala.js kernels should use typed arrays and avoid JS
interop inside hot loops.

## Backend Performance Strategy

Gale should treat pure Scala/JVM as a serious backend, not merely as a fallback.
It should also be honest that large dense BLAS and LAPACK kernels are exactly
where mature native implementations are hardest to beat.

The backend rule is:

> Pure Scala for portability, predictable memory-bound kernels, and tiny
> specialization; native BLAS/LAPACK for large dense throughput.

Expected backend strategy by operation class:

| Operation class | Preferred Gale strategy |
| --- | --- |
| `dot`, `axpy`, `scal`, `norm`, distance | Primitive pure kernels; optional Vector API backend |
| Pointwise dense operations | Fused primitive kernels |
| Small fixed matrices | Inline or generated unrolled kernels |
| Dense matrix-vector multiply | Layout-specialized primitive kernels |
| CSR sparse matrix-vector multiply | Primitive CSR kernels, row blocking where useful |
| Medium dense matrix multiplication | Pure blocked JVM kernel, later Vector API support |
| Large dense matrix multiplication | BLAS backend |
| Large dense LU, QR, Cholesky | LAPACK backend, pure portable fallback |
| SVD, eigen, Schur | Native backend for production-scale use |
| Generic `Matrix[A]` | Correct fallback, not the throughput path |

Pure Scala/JVM can be competitive when the operation is primitive,
monomorphic, allocation-free, simple enough for the JIT, and memory-bound. This
includes Level 1 BLAS-style vector operations, many elementwise fused kernels,
CSR SpMV, layout-specialized Level 2-style `gemv`, and tiny fixed-size kernels.

Native BLAS/LAPACK should be preferred when the operation is large dense
`gemm`, blocked dense factorization, SVD, eigendecomposition, or any kernel
where architecture-specific blocking, packing, SIMD microkernels, prefetching,
and threading dominate.

Backend selection must be thresholded and evidence-based, not hard-coded from
intuition:

```scala
trait BackendThresholds:
  def nativeDotMinLength: Int
  def nativeGemvMinWork: Long
  def nativeGemmMinFlops: Long
  def nativeFactorizationMinSize: Int
```

Thresholds must be benchmark-derived per backend family and conservative by
default. Calling native code for vectors of length three is a bug, not an
optimization.

Gale should support two dense storage modes:

```scala
final class DMat private[gale] (
  private[gale] val data: DoubleArray,
  val offset: Offset,
  val rows: Rows,
  val cols: Cols,
  val rowStride: Stride,
  val colStride: Stride
)

final class NativeDMat private[gale] (
  private[gale] val memory: java.lang.foreign.MemorySegment,
  val rows: Rows,
  val cols: Cols,
  val leadingDimension: Stride,
  val layout: Layout
)
```

`DMat` is the heap-backed portable representation. `NativeDMat` is the
off-heap FFM representation for BLAS/LAPACK-heavy code. Conversions must be
explicit or performed only above documented thresholds:

```scala
val nativeA = A.toNative
val x = nativeA.solve(b.toNative)
```

The native backend must account for copy cost. A native call that copies
heap-backed inputs into off-heap memory must include that cost in its threshold
decision.

Tiny fixed-size kernels should bypass general matrix code. Sized values such as
`Matrix.Sized[Double, 3, 3]`, `Vec.Sized[Double, 3]`, and dedicated `Mat2`,
`Mat3`, `Mat4` APIs may lower to unrolled code for operations such as
matrix-vector multiply, determinant, inverse, Cholesky, LU, and rigid-body
transforms.

Layout specialization is mandatory. Kernels should branch once:

```scala
if A.isContiguousRowMajor then
  DoubleKernels.gemvRowMajor(A, x, y)
else if A.isContiguousColMajor then
  DoubleKernels.gemvColMajor(A, x, y)
else
  DoubleKernels.gemvStrided(A, x, y)
```

Stride and layout checks do not belong inside innermost loops.

Pure JVM `gemm` should evolve through explicit stages:

```text
TinyUnrolledGemm
PureBlockedGemm
JvmVectorGemm
NativeBlasGemm
```

The v1 goal is a correct and respectable pure blocked kernel plus a clean
native backend path. Beating OpenBLAS, BLIS, MKL, or Accelerate on large GEMM is
not a v1 goal.

Decompositions need workspace-aware APIs:

```scala
final class Workspace private[gale] (
  val pivots: IndexArray,
  val panel: DoubleArray,
  val tau: DoubleArray,
  val work: DoubleArray
)

LU.factor(A, workspace = Workspace.forLU(A.rows))
```

Repeated solves, iterative algorithms, simulations, and optimization loops must
be able to reuse scratch memory.

Backend threading must be explicit. CPU-bound kernels should not use virtual
threads. Large JVM kernels may use a fixed compute pool. Native BLAS/LAPACK may
use its own native threads. Gale must not let both layers oversubscribe cores by
default:

```scala
final case class BackendConfig(
  jvmThreads: Int,
  nativeThreads: Int,
  allowNestedParallelism: Boolean = false
)
```

Expression fusion should stay BLAS-shaped:

```scala
y := alpha * x + y
y := alpha * (A * x) + beta * y
C := alpha * (A * B) + beta * C
z := x.pointwise * y + w
```

These should lower to `axpy`, `gemv`, `gemm`, or one fused pointwise pass.
Gale should not build a general lazy expression graph in v1.

Benchmarking must cover native crossover sizes, not only isolated throughput:

- `dot`, `axpy`, `norm`.
- fused pointwise operations.
- `gemv`.
- pure and native `gemm`.
- CSR SpMV.
- dense LU, Cholesky, QR.
- small fixed-size kernels.
- allocation rate.
- backend dispatch overhead.
- heap-to-native conversion cost.
- native threshold crossover by backend.
- warmup behavior and steady-state behavior.

JMH microbenchmarks are necessary but not sufficient. Gale should also keep a
small set of scenario benchmarks that resemble real solver and factorization
workloads, so JIT profiles are not only measured in artificial isolation.

## Error Model

Programmer shape errors are available as exceptions through convenience APIs,
but production APIs must have total alternatives.

Required error ADT:

```scala
sealed trait LinAlgError extends Exception:
  def message: String

object LinAlgError:
  final case class DimensionMismatch(expected: Shape, actual: Shape) extends LinAlgError
  final case class IndexOutOfBounds(index: Int, bound: Int) extends LinAlgError
  final case class NonSquareMatrix(shape: Shape) extends LinAlgError
  final case class SingularMatrix(index: Int) extends LinAlgError
  final case class NotPositiveDefinite(index: Int) extends LinAlgError
  final case class DidNotConverge(iterations: Int, residualNorm: Double) extends LinAlgError
  final case class UnsupportedOperation(operation: String, backend: String) extends LinAlgError
```

Required convenience pattern:

```scala
A.solve(b)       // Either[LinAlgError, Vec[Double]]
A.solveOrThrow(b)
A.lu.diagnostics
```

## Performance Requirements

Pure kernels must be serious, not merely fallback code.

Kernel rules:

- Use `while` loops.
- Use primitive storage.
- Avoid closures in hot loops.
- Avoid `Iterator`.
- Avoid `Range.foreach`.
- Avoid generic type classes in specialized dense kernels.
- Provide contiguous and strided paths.
- Provide transpose-view paths.
- Block `gemm`.
- Benchmark every primitive with JMH on JVM.
- Benchmark Scala.js fastLink and fullLink separately.
- Benchmark Scala.js Wasm on compute-heavy kernels.
- Add JVM allocation-rate tests for primitive kernels.
- Add JFR or equivalent allocation profiling for benchmark runs.
- Inspect bytecode or generated code for the smallest primitive kernel set.
- Fail performance regression checks when `DVec.dot(DVec)` or primitive SpMV
  allocates per element.
- Keep primitive dense and sparse kernels monomorphic after boundary validation.
- Measure native crossover thresholds instead of guessing them.
- Include heap-to-native conversion cost in native backend benchmarks.
- Include scenario benchmarks for solver/factorization workloads, not only
  isolated microbenchmarks.
- Record backend threading configuration during every benchmark run.

Shape checks happen at API boundaries. Unsafe kernels assume validated
dimensions and must be package-private.

## Cross-Platform Build Requirements

The initial build should use `sbt-crossproject`.

Recommended publish line:

```scala
ThisBuild / scalaVersion := "3.3.8"
```

CI should also test the current Scala Next line.

Required targets:

- JVM.
- Scala.js JavaScript backend.
- Scala.js fullLinkJS.
- Scala.js Wasm profile for selected kernel and solver tests.

Required artifact use:

```scala
libraryDependencies += "io.gale" %%% "gale-core" % version
```

JVM users can add:

```scala
libraryDependencies += "io.gale" %% "gale-backend-jvm-blas-ffm" % version
```

## Documentation Requirements

Documentation must teach Gale as a typed functional linear algebra library, not
as a list of LAPACK wrappers.

Required docs:

- Getting started.
- Dense vectors and matrices.
- Views, slices, and transpose.
- Functional updates and builders.
- Explicit destination APIs.
- Dense solving and factorization values.
- Partial spectral decompositions and ordering modes.
- Sparse construction and canonicalization.
- Iterative solvers and diagnostics.
- Matrix properties and `assume*` / `verify*`.
- Backend selection.
- Scala.js and Wasm notes.
- Breeze migration guide.
- Numerical accuracy and reproducibility guide.

## Verification Requirements

The test suite is part of the product.

Required law groups:

- Vector space laws with tolerances.
- Matrix-vector linearity.
- Matrix-matrix associativity where numerically stable.
- Transpose laws.
- Factorization reconstruction laws.
- Solve residual laws.
- Spectral residual, orthogonality, and ordering laws.
- Sparse canonicalization laws.
- Dense/sparse equivalence on small matrices.
- Backend conformance for pure, JS, Wasm, Vector API, and BLAS backends.

Required regression suites:

- Ill-conditioned dense systems.
- Singular matrices.
- Non-positive-definite Cholesky inputs.
- Clustered or repeated eigenvalues.
- Rank-deficient partial SVD inputs.
- Singular or indefinite generalized eigenproblem `B` matrices.
- Sparse duplicate entries.
- Explicit stored sparse zeros.
- Empty matrices and zero-length vectors.
- Strided views.
- Transpose views.
- JS typed-array indexing boundaries.

## Release Roadmap

### v0.1: Dense Core

- `DVec`, `DMat`, views, slices, transpose.
- Platform `DoubleArray`, `FloatArray`, `IndexArray`.
- Dot, norm, axpy, matrix-vector multiply, matrix-matrix multiply.
- Layout-specialized `gemv` paths.
- Initial tiny fixed-size kernel surface for `Mat2`, `Mat3`, `Mat4` or sized
  equivalents.
- Pure Scala kernels.
- MUnit and ScalaCheck laws.
- JVM and Scala.js tests.
- Initial JMH benchmark suite.

### v0.2: Dense Decompositions

- LU with partial pivoting.
- Cholesky.
- QR.
- Triangular solve.
- `solve`, `trySolve`, `leastSquares`.
- Determinant.
- Rank and condition estimates.
- Factorization diagnostics.

### v0.3: Sparse and Linear Operators

- `LinearOperator`.
- COO builder.
- CSR and CSC.
- Sparse transpose and SpMV.
- Sparse addition and scaling.
- Diagonal, identity, permutation.
- CG, BiCGSTAB, restarted GMRES.
- Jacobi preconditioner.
- Matrix Market read/write.

### v0.3.5: Spectral and Generalized Decompositions

- Spectral result types: `EigenDecomposition`, `SVD`, generalized spectral
  results, and `SpectralDiagnostics`.
- MATLAB/SciPy spectral parity matrix translated into Gale typed options rather
  than string flags.
- Explicit order enums for top/bottom and magnitude/algebraic/real-part
  selection.
- Partial symmetric eigendecomposition with largest and smallest algebraic
  eigenvalues.
- Partial nonsymmetric eigendecomposition with largest and smallest magnitude
  eigenvalues.
- Partial SVD with largest and smallest singular values.
- Generalized symmetric-definite eigenproblem support, `A x = lambda B x`.
- Generalized SVD API for matrix pairs, with largest and smallest generalized
  singular values.
- Residual, orthogonality, rank-deficiency, and convergence tests.
- Optional backend boundary for production spectral engines.

### v0.4: Type and API Polish

- Optional sized layer.
- Matrix property wrappers.
- Breeze interop.
- Syntax modules.
- Documentation and migration guide.
- Examples: regression, PCA, graph Laplacian, least squares.

### v0.5: Acceleration

- JVM Vector API backend.
- JVM BLAS/LAPACK FFM backend.
- Native memory-backed `NativeDMat` design.
- Backend threshold policy and measured native crossover defaults.
- Backend threading policy.
- Pure blocked `gemm` roadmap and benchmark comparison.
- Wasm kernel profile.
- Backend conformance suite.
- Benchmark dashboard.

### v1.0: Stability

- Binary compatibility policy.
- Complete dense real v1 API.
- Clear sparse v1 API.
- Documented numerical guarantees.
- Published migration guide.
- No exposed platform storage leakage.
- No residual `atlas` package, docs, or examples.

## Acceptance Criteria

Gale v1 is acceptable when:

- `gale-core`, `gale-kernel`, `gale-factorization`, `gale-sparse`,
  `gale-solvers`, and `gale-laws` cross-build for JVM and Scala.js.
- The public examples in this PRD compile with `import gale.linalg.*`.
- `gale-core` has no Breeze, Spire, Cats Effect, ZIO, native, plotting, or JVM-only
  dependency.
- Dense solve residual examples satisfy documented tolerance thresholds.
- Sparse iterative solver examples return convergence diagnostics.
- Partial spectral examples return ordered eigenvalue/singular-value results
  with residual diagnostics.
- Backend conformance tests pass for pure JVM and Scala.js.
- JVM acceleration backends are optional and selected only by explicit imports.
- JS storage uses typed arrays behind private platform abstractions.
- There is no public `Array[Double]` storage contract.
- `Vec(1.0, 2.0, 3.0)` returns the primitive `DVec` representation.
- Primitive dense and sparse kernels have no per-element boxing or collection
  allocation in benchmark allocation profiles.
- Generic `Vec[A]` fallback behavior is documented as correctness-oriented, not
  the primitive throughput path.
- Native BLAS/LAPACK routing is thresholded with benchmark-derived defaults.
- Native-backed dense storage is explicit and documented.
- Backend threading policy prevents JVM and native oversubscription by default.
- Tiny fixed-size kernels bypass the general strided matrix path.
- There are no `atlas.linalg`, `zephyr.linalg`, or placeholder package names in
  code or docs.

## Open Decisions

- Organization and artifact coordinates.
- Default owned dense layout: row-major, column-major, or explicit layout
  default with constructor value order separated from storage order.
- Whether `solve` should return `Either` by default or whether `trySolve` should
  be the total form with `solve` as throwing convenience.
- Whether `gale-laws` should be published as a normal module or testkit-only
  artifact.
- Minimum supported Node/browser versions for Scala.js Wasm examples.
- Whether Matrix Market IO belongs in `gale-sparse` or a small `gale-io` module.
- Initial native backend threshold defaults by operation and platform.
- Whether `NativeDMat` lives in the BLAS backend module or a small JVM-native
  storage module shared by BLAS and future sparse direct solvers.

## Product Warning

The central failure mode is overgeneralizing too early. Gale should not try to
be a tensor library, plotting system, probabilistic programming library, and
native numerical wrapper at once. The winning v1 is a small, typed, functional
linear algebra core with honest runtime backends and enough sparse support to be
useful without pretending to replace SuiteSparse.
