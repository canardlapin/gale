# gale — worked examples

A topic guide to `gale`'s public API, illustrated with runnable code. Every
snippet below reflects a real call pattern verified against the source; the
regression/PCA/graph-Laplacian examples referenced throughout have a
compiled, test-backed twin — see
["Where the numbers come from"](#where-the-numbers-come-from) at the end.

## Getting started

`gale` is not yet published; add it as a local `crossProject` dependency (e.g.
via `sbt publishLocal`, or as a source dependency inside a multi-project
build) under coordinates `"io.gale" %%% "gale-core" % "<version>"`. The two
imports you need almost everywhere are:

```scala
import gale.linalg.*      // DVec, DMat, Matrix, Vec, LinAlgError, ...
import gale.spectral.*    // Eigen, Svds, EigenSelection, EigenVectors, ...
```

A three-line hello world — build a small matrix and vector, multiply, print:

```scala
import gale.linalg.{Matrix, Vec}

val A = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
val x = Vec(1.0, 1.0)
println((A * x).toSeq) // Vector(3.0, 7.0)
```

## Dense vectors & matrices

Construction:

```scala
import gale.linalg.{DVec, Matrix, Vec}

val v  = Vec(1.0, 2.0, 3.0)               // DVec from literal values
val z  = Vec.zeros(4)                     // length-4 zero vector
val t  = Vec.tabulate(5)(i => i * i.toDouble)

val A  = Matrix.dense(2, 3)(              // row-major literal values
  1.0, 2.0, 3.0,
  4.0, 5.0, 6.0
)
val I  = Matrix.eye(3)                    // 3x3 identity
val M  = Matrix.tabulate(3, 3)((i, j) => if i == j then 2.0 else 0.0)
```

Arithmetic (`+`, `-`, scaling, matrix-vector, matrix-matrix, transpose) is
built into the core operators — no import beyond `gale.linalg` needed:

```scala
val sum   = A + A                 // elementwise, same shape required
val Av    = A * v                 // matrix-vector product -> DVec
val AtA   = A.t * A               // transpose is an O(1) strided view
val dot   = v.dot(Vec(1.0, 0.0, 0.0))
val norm  = v.norm2
```

`row`/`col`/`t` are `O(1)` aliasing views over the same backing storage;
`copy`, `updated`, and the `Seq`/`Array` exporters return independent data.
`DMat`/`DVec` arithmetic (`+`, `-`, `*`) throws `LinAlgError` on a shape
mismatch — these are primitive operators, not the `Either`-returning solve
tier described next.

**The Either-first convention.** Every operation that can fail structurally
(singular matrix, non-SPD input, rank deficiency, dimension mismatch on a
solve) returns `Either[LinAlgError, A]` rather than throwing:

```scala
A.solve(b) match
  case Right(x)    => println(s"solved: ${x.toSeq}")
  case Left(error) => println(s"failed: $error")

// or, compositionally:
val maybeDet: Either[LinAlgError, Double] = A.det
val scaled = maybeDet.map(_ * 2.0)
```

## Solvers & factorizations

`DMat` exposes the dense factorizations directly. `lu`, `cholesky`, `solve`,
`det`, and `leastSquares` are total but can fail structurally, so they return
`Either[LinAlgError, _]`; `qr` always succeeds structurally (every matrix has
a QR factorization), so it returns the factorization directly, not an
`Either`:

```scala
import gale.linalg.{LinAlgError, Matrix, Vec}

val A = Matrix.dense(3, 3)(4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0)
val b = Vec(1.0, 2.0, 3.0)

val solved = A.solve(b)          // Either[LinAlgError, DVec] via LU
val d      = A.det                // Either[LinAlgError, Double]
val chol   = A.cholesky           // Either[LinAlgError, Cholesky] (SPD only)

// QR — always total; least-squares is the Either-returning step.
val tall = Matrix.dense(4, 2)(1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0)
val qr   = tall.qr                            // QR (never fails)
val coeffs = qr.solveLeastSquares(Vec(1.0, 3.0, 5.0, 7.0))
// or, the DMat-level shortcut for the same thing:
val coeffs2 = tall.leastSquares(Vec(1.0, 3.0, 5.0, 7.0))
```

`qr.q` (the orthogonal factor) is rebuilt from the stored Householder
reflectors on first access and cached; `qr.r` is the upper-triangular factor;
`qr.diagnostics.rank` reports the numerical rank found during factorization.

Standalone triangular solves (for a matrix you already know is triangular,
without going through `lu`/`cholesky`):

```scala
import gale.linalg.TriangularSolve

val lowerSolve = TriangularSolve.lower(L, b)  // Either[LinAlgError, DVec]
val upperSolve = TriangularSolve.upper(U, b)
```

Iterative solvers live in `gale.solvers`: `cg` (SPD systems), `bicgstab` and
`gmres` (general nonsymmetric systems), and two least-squares routes, `cgnr`
(CG on the normal equations — simple, squares the condition number) and
`lsqr` (Golub–Kahan bidiagonalization — numerically better for
ill-conditioned `A`, never forms `AᵀA`):

```scala
import gale.solvers.*
import gale.linalg.{Matrix, Vec}

val A = Matrix.tabulate(100, 100)((i, j) => if i == j then 4.0 else 0.0) // a DMat is a DoubleLinearOperator
val b = Vec.tabulate(100)(_.toDouble)

val result = cg(A, b, config = SolverConfig(tolerance = 1e-10, maxIterations = 200))
result match
  case r if r.converged => println(s"converged in ${r.iterations} iters, residual ${r.residual}")
  case r                => println(s"did not converge: residual ${r.residual}")

// A preconditioned solve:
val jacobi = Preconditioner.Jacobi(A)
val precond = bicgstab(A, b, preconditioner = jacobi)
```

`SolverResult` is `Converged` or `NotConverged` (both carry `x`, `iterations`,
`residual`); `.converged: Boolean` distinguishes them, and `.orThrow` (a
throwing convenience, best kept to scripts/tests) unwraps to `x` or raises
`LinAlgError.DidNotConverge`.

## Spectral

The examples below share one symmetric positive-definite matrix — the
standard 1-D tridiagonal "stiffness" matrix, whose spectrum is known
analytically (`2 - 2*cos(kπ/(n+1))`):

```scala
import gale.linalg.Matrix

val n = 6
val A = Matrix.tabulate(n, n): (i, j) =>
  if i == j then 2.0
  else if math.abs(i - j) == 1 then -1.0
  else 0.0
```

**Symmetric eigendecomposition.** Eigenvalues are always returned
ascending-algebraic, independent of the selection's ordering criterion — the
selection chooses *membership*, never layout:

```scala
import gale.spectral.{Eigen, EigenSelection, EigenOrder, EigenVectors}

val full  = Eigen.eigSymmetric(A, EigenSelection.All, EigenVectors.Right)
val top3  = Eigen.eigSymmetric(A, EigenSelection.Count(3, EigenOrder.LargestAlgebraic), EigenVectors.Right)
val byIdx = Eigen.eigSymmetric(A, EigenSelection.IndexRange(0, 2))          // ranks 0..2, ascending
val byVal = Eigen.eigSymmetric(A, EigenSelection.ValueInterval(0.0, 3.0))   // every λ in (0, 3]

full.foreach { d =>
  println(d.eigenvalues.toSeq)      // ascending
  println(d.eigenvectors.cols)      // orthonormal columns, aligned with eigenvalues
  println((d.diagnostics.allConverged, d.diagnostics.worstResidual))
}
```

Only the matrix's lower triangle is read (the same convention as `cholesky`).
A partial, iterative solve (`eigsh`-style Lanczos) is available for a large or
matrix-free operator, restricted to `EigenSelection.Count` with `k < n`:

```scala
import gale.spectral.{Eigen, EigenSelection, EigenOrder, SpectralOptions}

val big = Matrix.tabulate(20, 20): (i, j) =>
  if i == j then 2.0 else if math.abs(i - j) == 1 then -1.0 else 0.0

val partial = Eigen.eigSymmetric(
  big,
  big.rows,
  EigenSelection.Count(3, EigenOrder.LargestMagnitude),
  SpectralOptions(tolerance = 1e-10, maxIterations = 50)
)
partial.foreach(d => println(d.requireConverged))
```

**Nonsymmetric eigendecomposition.** A real input can have complex
eigenvalues in conjugate pairs; the result never exposes gale's internal
packing directly — read it through the typed accessors. A 2x2 rotation
matrix is the simplest example with a genuinely complex spectrum (`±i`):

```scala
import gale.linalg.Matrix
import gale.spectral.{Eigen, EigenSelection, EigenVectors}

val rotation = Matrix.dense(2, 2)(0.0, -1.0, 1.0, 0.0)
val d = Eigen.eigNonsymmetric(rotation, EigenSelection.All, EigenVectors.LeftAndRight).toOption.get
for i <- 0 until d.size do
  val lambda = d.eigenvalue(i) // a Complex; here 0 + i and 0 - i
  if d.isRealPair(i) then
    val (vRe, _) = d.eigenvector(i) // real eigenvector
    ()
  else
    val (vRe, vIm) = d.eigenvector(i) // v = vRe + i*vIm; conjugate is the adjacent index
    ()
```

Conjugate pairs are stored as two adjacent columns (positive-imaginary member
first) and are always kept together — a `Count` selection whose boundary
would split a pair returns the whole pair instead. `EigenVectors.Left` /
`LeftAndRight` are available on the dense path only (left vectors are derived
from the full right-eigenvector matrix); the iterative (Arnoldi) path
supports right vectors only.

**Partial SVD** (`Svds`, v0.3.5 ships partial only — no full dense SVD yet;
`k` must satisfy `0 < k < min(rows, cols)`):

```scala
import gale.spectral.{Svds, SingularSelection, SingularOrder, EigenVectors}

val svd = Svds.svd(A, SingularSelection.Count(2, SingularOrder.Largest), EigenVectors.Right)
svd.foreach(s => println(s.singularValues.toSeq)) // always descending
```

**Generalized symmetric-definite** eigenproblem `A x = λ B x` (`B` SPD):

```scala
import gale.spectral.{Eigen, EigenSelection, EigenVectors}

val B = Matrix.tabulate(n, n)((i, j) => if i == j then (i + 1).toDouble else 0.0) // diag(1..n), SPD
val gen = Eigen.eigSymmetricGeneralized(A, B, EigenSelection.All, EigenVectors.Right)
// eigenvectors are B-orthonormal (XᵀBX = I), not Euclidean-orthonormal
```

The general nonsymmetric pencil (QZ) is a **backend-scoped seam**:
`Eigen.eigGeneralizedNonsymmetric(A, B)` needs a `given SpectralBackend` in
scope (found automatically via companion-object implicit search — no import
needed for the default), and with the always-available
`SpectralBackend.none` it returns `Left(LinAlgError.UnsupportedOperation(...))`
— the pure core ships no QZ engine:

```scala
import gale.spectral.Eigen

val qz = Eigen.eigGeneralizedNonsymmetric(A, B) // Left(UnsupportedOperation(...)) with no backend imported
```

## Sparse

Build via COO (coordinate) triples, then convert to a compressed format for
repeated use:

```scala
import gale.sparse.Sparse

val builder = Sparse.coo(rows = 3, cols = 3)
builder.add(0, 0, 2.0)
builder.add(0, 1, -1.0)
builder.add(1, 0, -1.0)
builder.add(1, 1, 2.0)
builder.add(2, 2, 5.0)

val csr = builder.toCSR()      // CSR, canonical (sorted, deduplicated) by default
val csc = builder.toCSC()
val coo = builder.toCOO()      // duplicate coordinates summed unless told otherwise

val y = csr * Vec.fill(3)(1.0) // sparse mat-vec, DVec result
```

`COOBuilder.toCSR`/`toCSC`/`toCOO` all take a `DuplicatePolicy` (`Sum`
default, `Last`, or `Error`). A `COO`/`CSR`/`CSC` value's `.canonicalize`
sorts, sums duplicates, and prunes exact zeros; `.hasCanonicalFormat` tells
you whether that work is already done. `CSR`/`CSC` support `+`, `-`, `*`
(scalar), `mapValues`, `zipValues`, `.t` (a zero-copy reinterpretation
between the two formats), `.toDense(maxEntries)`, `.diagonal`, `.trace`.
`Sparse.diagonal(...)`, `Sparse.identity(n)`, `Sparse.zero(rows, cols)`, and
`Sparse.permutation(...)` build the structural special cases directly.

Matrix Market I/O (`gale.sparse.MatrixMarket`, `coordinate real general`
only):

```scala
import gale.sparse.MatrixMarket

val text = MatrixMarket.writeCoordinate(csr)
val back = MatrixMarket.readCoordinate(text)
```

## Optional sized layer

`gale.sized` adds compile-time-checked shapes over the same runtime storage —
opaque wrappers (`SVec[N]`, `SMat[R, C]`) that lower to `DVec`/`DMat`
zero-copy and zero-allocation. Nothing in core depends on it; it is entirely
opt-in:

```scala
import gale.sized.{SVec, SMat}

val v: SVec[3] = SVec.sized[3](1.0, 2.0, 3.0)
val m: SMat[3, 3] = SMat.sized[3, 3](
  1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0
)
val mv: SVec[3] = m * v            // shape-checked at compile time

// Raising a runtime DVec/DMat is a checked, zero-copy operation:
val maybeSized: Either[gale.linalg.LinAlgError, SVec[3]] = SVec.fromDVec[3](v.toDVec)

// Lowering back to the runtime API is always safe and zero-copy:
val runtime: gale.linalg.DVec = v.toDVec
```

`SMat[2, 2]` and `SMat[3, 3]` additionally carry closed-form `det`/`inverse`
kernels (`Either[LinAlgError, SMat[N, N]]`, `Left(SingularMatrix)` on an
exactly-zero or non-finite determinant); for ill-conditioned systems, lower
to `DMat` and use gale's pivoted LU/QR instead.

## Syntax modules

Two focused, opt-in import modules add ergonomic sugar without touching the
default surface:

```scala
import gale.syntax.all.*   // elementwise (Hadamard) ops and zipMap

val hadamard = a.pointwise * b     // elementwise product, NOT a.* b (matrix product)
val quotient = a.pointwise / b
val mapped   = a.pointwise.map(x => x * x)
val combined = a.zipMap(b)((x, y) => x + y)
```

```scala
import gale.syntax.unicode.*  // × (matrix product) and ⋅ (dot product), ASCII synonyms

val Ax  = A × x     // same as A * x
val AB  = A × B     // same as A * B
val xy  = x ⋅ y      // same as x.dot(y)
```

`unicode` is a separate import from `all` so a wildcard `gale.syntax.all.*`
import stays unsurprising.

## Migrating from Breeze

`gale-interop-breeze` (JVM-only) is the single place gale meets
[Breeze](https://github.com/scalanlp/breeze); `gale-core` itself stays
Breeze-free.

**Copy vs view.** The direction is asymmetric, by construction: gale's
storage never escapes the package undoctored, so gale → Breeze is
copy-only; Breeze's layout is a strict subset of gale's `(offset, rowStride,
colStride)` layout, so Breeze → gale can be a zero-copy *view* that aliases
the same backing array.

```scala
import gale.interop.breeze.*
import gale.linalg.Matrix

val myGaleMatrix = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)

val bm: breeze.linalg.DenseMatrix[Double] = toBreezeCopy(myGaleMatrix)       // always copies
val fromCopy: gale.linalg.DMat = fromBreezeCopy(bm)                          // independent storage
val fromView: gale.linalg.DMat = fromBreezeView(bm)                         // aliases bm's storage!
```

Sparse conversions (`CSR`/`CSC` ↔ Breeze `CSCMatrix`) are always copies in
both directions — the two libraries differ in compressed axis and
canonical-ordering bookkeeping.

**Migration shims.** `BreezeMigration` mirrors common Breeze call sites
(Breeze types in and out, computed by gale underneath), throwing
`LinAlgError` on failure the way Breeze itself throws — useful for a
mechanical, one-call-site-at-a-time port, after which prefer gale's native
`Either`-returning API directly:

```scala
import gale.interop.breeze.BreezeMigration
import breeze.linalg.{DenseMatrix, DenseVector}

val a = DenseMatrix.tabulate(3, 3)((i, j) => if i == j then 4.0 else 0.1) // symmetric, diagonally dominant
val b = DenseVector(1.0, 2.0, 3.0)

val x      = BreezeMigration.solve(a, b)     // A \ b, square system
val d      = BreezeMigration.det(a)
val l      = BreezeMigration.cholesky(a)     // lower factor, throws if not SPD
val (w, v) = BreezeMigration.eigSym(a)       // (eigenvalues ascending, eigenvectors)

// leastSquares wants a tall (overdetermined) system:
val tallA  = DenseMatrix.tabulate(5, 2)((i, j) => if j == 0 then 1.0 else i.toDouble)
val tallB  = DenseVector.tabulate(5)(i => 2.0 + 3.0 * i)
val coeffs = BreezeMigration.leastSquares(tallA, tallB) // A \ b, overdetermined
```

Matrix right-hand sides are supported as well. These overloads factor `A` once
and reuse its LU or QR factors across all columns:

```scala
val manySolutions = BreezeMigration.solve(a, manyRightHandSides)
val manyCoeffs    = BreezeMigration.leastSquares(tallA, manyResponseColumns)
```

The exact supported replacement boundary—including rank/condition-number
differences, performance evidence, and deliberate exclusions—is documented in
[`breeze-equivalence.md`](breeze-equivalence.md).

## Accuracy & determinism

Every pure dense kernel is written against a shared strided algorithm, but the
platforms use different fused-multiply-add semantics: the JVM uses `Math.fma`,
while Scala.js uses plain `a*b + c`. Results from the pure single-threaded route
are deterministic for a fixed runtime and build, but JVM and browser results may
differ in the last ulps. Native vendor libraries and future multithreaded routes
may legally reassociate operations, so Gale does not promise bit identity from
an accelerated backend.

The failure model throughout is **Either-first**: structural or
precondition violations (`NonSquareMatrix`, `SingularMatrix`,
`NotPositiveDefinite`, `RankDeficient`, `DimensionMismatch`, ...) are
`Left(LinAlgError)`, never an exception, on every total entry point. The
`Either`-returning API is `gale`'s contract; treat any place that appears to
throw a `LinAlgError` (arithmetic operators like `DMat.+`, `DVec.dot`, sparse
`mulInto`) as a primitive whose precondition you are expected to have already
checked, not as a value-returning solve.

## Scala.js / WebAssembly

The exact same public API compiles and runs unchanged in the browser: `gale`
is a `crossProject` (`CrossType.Full`) with shared source in `core/shared`
and platform-specific storage/interop in `core/jvm` / `core/js`. A plain `sbt
coreJS/test` produces today's ordinary JavaScript build.

Scala.js can also target WebAssembly instead of JavaScript, gated behind an
opt-in environment variable so the default build is untouched:

```sh
GALE_WASM=1 sbt coreJS/test
GALE_WASM=1 sbt benchSmokeJSFull
```

The build supplies Node's required `--experimental-wasm-exnref` flag. Wasm is
correctness-tested but remains experimental and default-off: the current Node 24
profile is 23–43x slower than optimized JavaScript. See the README and backend
dashboard for the evidence; CI executes the tests/profile as allow-failure.

## Where the numbers come from

The four worked examples referenced throughout this guide — least-squares
regression, the normal-equations cross-check, PCA via the symmetric
eigendecomposition, and a graph Laplacian's Fiedler vector (plus a fifth,
sparse-matvec check on the same Laplacian) — are not just prose: they are
compiled and asserted in
[`gale.examples.WorkedExamplesSuite`](../core/shared/src/test/scala/gale/examples/WorkedExamplesSuite.scala),
which runs green on both `coreJVM/test` and `coreJS/test`. Read that file for
the exact, working code behind every number in this document.
