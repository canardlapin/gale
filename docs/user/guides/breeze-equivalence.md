# Moving a Breeze linear algebra workload to Gale

Gale can replace the linear algebra part of a Breeze program when the program
works mainly with real `Double` vectors and matrices. Gale provides dense
arithmetic and factorizations, compressed sparse matrices, iterative solvers,
and selected eigenvalue and singular-value routines. The shared API runs on the
JVM and Scala.js.

Gale is not source-compatible with Breeze, and it does not replace Breeze's
probability, optimization, statistics, signal-processing, plotting, machine
learning, or tensor modules. Use this guide to check the supported linear
algebra operations and to see how migrated code changes.

## Construct and inspect dense values

Gale uses zero-based indices. Slice endpoints are half-open: the starting index
is included and the ending index is not.

| Task | Breeze | Gale |
| --- | --- | --- |
| Construct a vector | `DenseVector(1.0, 2.0, 3.0)` | `Vec(1.0, 2.0, 3.0)` |
| Construct a matrix | `DenseMatrix((1.0, 2.0), (3.0, 4.0))` | `Matrix(2, 2)(1.0, 2.0, 3.0, 4.0)` |
| Read the shape | `a.rows`, `a.cols`, `x.length` | the same |
| Read an element | `a(i, j)`, `x(i)` | the same |
| Select a row | `a(i, ::).t` | `a.row(i)` |
| Select a column | `a(::, j)` | `a.col(j)` |
| Select a contiguous submatrix | `a(r0 until r1, c0 until c1)` | `a.slice(r0, r1, c0, c1)` |
| Select rows by index | Breeze indexing with `::` | `a.gatherRows(indices)` |
| Replace one value | mutate a copy with `copy(i, j) = value` | `a.updated(i, j, value)` |

`Matrix(rows, cols)(values*)` reads the values in row-major order.
`Matrix.dense(rows, cols)(values*)` is the more explicit spelling of the same
operation.

```scala mdoc
import gale.linalg.*

val a = Matrix(3, 3)(
  1.0, 2.0, 3.0,
  4.0, 5.0, 6.0,
  7.0, 8.0, 9.0
)

(a.rows, a.cols)
a(1, 2)
a.row(1).toSeq
a.col(0).toSeq
a.slice(0, 2, 1, 3).valuesRowMajor
```

`row`, `col`, `slice`, and `t` return views. A view does not copy the selected
values. `gatherRows` and `gatherColumns` return new matrices because an indexed
selection may reorder or repeat entries.

```scala mdoc
a.gatherRows(IndexedSeq(2, 0)).valuesRowMajor
a.gatherColumns(IndexedSeq(2, 2, 0)).valuesRowMajor

val changed = a.updated(0, 0, 10.0)
(a(0, 0), changed(0, 0))
```

Vectors follow the same rule. `x.slice(from, until)` returns a view,
`x.gather(indices)` returns a copy, and `x.updated(i, value)` returns a new
vector. Use a builder when filling many entries.

## Create generated and structured values

Use `tabulate` when each value can be computed from its row and column. Gale
also provides direct constructors for zero matrices, identity matrices, and
common sparse structures.

```scala mdoc
val zeroVector = Vec.zeros(4)
val identity3 = Matrix.eye(3)
val tridiagonal4 = Matrix.tabulate(4, 4): (row, col) =>
  if row == col then 2.0
  else if math.abs(row - col) == 1 then -1.0
  else 0.0

zeroVector.toSeq
(identity3 * Vec(1.0, 2.0, 3.0)).toSeq
(tridiagonal4 * Vec.fill(4)(1.0)).toSeq
```

`Matrix.tabulate` constructs a dense matrix. If most of the generated values
would be zero, build a sparse matrix or define an operator instead; both
choices are shown below.

## Multiply and transform values

Ordinary `*` performs an algebraic matrix product. Import
`gale.syntax.all.*` for pointwise matrix operations.

```scala mdoc
import gale.syntax.all.*

val x = Vec(1.0, 0.5, -1.0)
val ax = a * x
val gram = a.t * a
val scaled = a * 2.0
val squared = a.pointwise * a
val regularized = gram.addToDiagonal(1e-6)

ax.toSeq
gram.valuesRowMajor
scaled.valuesRowMajor
squared.valuesRowMajor
regularized(0, 0)
```

The main dense operations are:

| Task | Gale call |
| --- | --- |
| Vector addition, subtraction, scaling | `x + y`, `x - y`, `x * scalar` |
| Dot product and Euclidean norm | `x.dot(y)`, `x.norm2` |
| Matrix addition, subtraction, and scaling | `a + b`, `a - b`, `a * scalar`, `scalar * a` |
| Matrix-vector and matrix-matrix products | `a * x`, `a * b` |
| Transpose | `a.t` |
| Elementwise multiply, divide, or map | `a.pointwise * b`, `a.pointwise / b`, `a.pointwise.map(f)` |
| Add a scalar to the diagonal | `a.addToDiagonal(value)` |
| Average a matrix with its transpose | `a.symmetrizedAverage` |
| Kronecker product | `a.kron(b)` |

## Solve systems and least-squares problems

Call `solve` for a square system. A vector right-hand side returns one solution;
a matrix right-hand side returns one solution in each column. Gale factors the
coefficient matrix once for the call.

```scala mdoc
val system2 = Matrix(2, 2)(
  3.0, 1.0,
  1.0, 2.0
)
val rhs2 = Vec(9.0, 8.0)
val rhsMatrix2 = Matrix(2, 2)(
  9.0, 4.0,
  8.0, 5.0
)

system2.solve(rhs2).map(_.toSeq)
system2.solve(rhsMatrix2).map(_.valuesRowMajor)
```

Use `leastSquares` for a tall system. The following design has an intercept
column and a slope column. The observations follow `1 + 2x`, so the fitted
coefficients should be close to `(1, 2)`.

```scala mdoc
val design4x2 = Matrix(4, 2)(
  1.0, 0.0,
  1.0, 1.0,
  1.0, 2.0,
  1.0, 3.0
)
val observations4 = Vec(1.0, 3.0, 5.0, 7.0)

design4x2.leastSquares(observations4).map(_.toSeq)
```

Both methods return `Either[LinAlgError, A]`. A singular square matrix, a
rank-deficient design, or a right-hand side with the wrong shape returns
`Left(error)`.

Use `a.lu`, `a.cholesky`, or `a.qr` when you need the factor itself or when
separate calls reuse the same coefficient matrix. LU and Cholesky can fail, so
they return `Either`. QR factorization is defined for every dense matrix and
returns a `QR` directly.

```scala mdoc
val qr4x2 = design4x2.qr
(qr4x2.r.rows, qr4x2.r.cols, qr4x2.diagnostics.rank)
system2.cholesky.map(_.lower.valuesRowMajor)
```

See [Factorization capabilities](../reference/factorization-capabilities.md)
for the common interfaces implemented by LU, Cholesky, and QR factors.

## Compute eigenvalues and singular values

Use `Eigen.eigSymmetric` when a real square matrix is symmetric. The result
stores eigenvalues in ascending algebraic order. Ask for
`EigenVectors.ValuesOnly` when eigenvectors are not needed.

```scala mdoc
import gale.spectral.*

val symmetric3 = Matrix(3, 3)(
  2.0, -1.0, 0.0,
  -1.0, 2.0, -1.0,
  0.0, -1.0, 2.0
)

val symmetricSpectrum = Eigen.eigSymmetric(
  symmetric3,
  EigenSelection.All,
  EigenVectors.ValuesOnly
)

symmetricSpectrum.map(_.eigenvalues.toSeq)
```

`a.svd` computes a full economy-size singular value decomposition. If `a` is
`m × n`, Gale returns `U` with shape `m × k` and `Vᵀ` with shape `k × n`, where
`k = min(m, n)`.

```scala mdoc
symmetric3.svd.map: result =>
  (
    result.singularValues.toSeq,
    (result.u.rows, result.u.cols),
    (result.vt.rows, result.vt.cols)
  )
```

Use `Svds.svd` when only the largest or smallest singular values are needed.
Use the operator overloads of `Eigen` and `Svds` when the matrix is too large to
store. Iterative spectral results include residuals and convergence counts; call
`requireConverged` when the program requires every requested result.

Gale also supports:

- dense nonsymmetric eigendecomposition through `Eigen.eigNonsymmetric`;
- dense generalized symmetric-definite problems through
  `Eigen.eigSymmetricGeneralized`;
- matrix-free generalized symmetric-definite problems through LOBPCG or
  generalized block Lanczos; and
- pseudo-inverses through `a.pinv`.

Read [Matrix-free generalized symmetric eigensolving](../advanced/generalized-operator-eigen.md)
before choosing between LOBPCG and generalized block Lanczos.

## Build and use sparse matrices

Build a coordinate matrix by adding its stored entries, then convert it to CSR
or CSC. CSR is the usual choice for repeated row-oriented matrix-vector
products. CSC gives efficient access by column.

```scala mdoc
import gale.sparse.*

val sparseA = Sparse
  .coo(3, 3)
  .add(0, 0, 2.0)
  .add(0, 1, -1.0)
  .add(1, 0, -1.0)
  .add(1, 1, 2.0)
  .add(1, 2, -1.0)
  .add(2, 1, -1.0)
  .add(2, 2, 2.0)
  .toCSR()

(sparseA.nnz, (sparseA * Vec(1.0, 2.0, 3.0)).toSeq)
```

The builder sorts coordinates, handles duplicates according to a
`DuplicatePolicy`, and removes explicit zeros when it creates a canonical
compressed matrix. Use `Sparse.cooChecked` when malformed input or non-finite
values should be returned as `LinAlgError` rather than thrown.

`Sparse.diagonal`, `Sparse.identity`, `Sparse.zero`, and `Sparse.permutation`
construct common structures without passing through a coordinate builder.
CSR and CSC support element access, rows, columns, transpose, addition,
subtraction, scalar multiplication, and conversion to dense storage.

Gale does not provide sparse direct factorization. Use an iterative solver or
convert a small sparse matrix to dense storage when the conversion is known to
fit in memory. JVM and Scala.js both resolve the staged `SparseDirect` seam
with `SparseDirectProvider.none`. A capable provider is a later explicit
import, not a browser SuiteSparse; see
[Scala.js sparse-direct](../../sparse-direct-js.md).

## Solve with an iterative method

Dense matrices, sparse matrices, and custom operators all implement
`DoubleLinearOperator`. The iterative solvers accept that interface rather than
requiring a particular storage format.

Use conjugate gradient (`cg`) for a symmetric positive-definite operator.
`bicgstab` and `gmres` handle general nonsymmetric square systems. `lsqr` solves
least-squares problems without forming `AᵀA`; `cgnr` solves the normal equations
and is best reserved for well-conditioned problems.

```scala mdoc
import gale.solvers.*

val knownSolution = Vec(1.0, 2.0, 3.0)
val iterative = cg(
  sparseA,
  sparseA * knownSolution,
  SolverConfig(tolerance = 1e-12, maxIterations = 20)
)

(iterative.converged, iterative.iterations, iterative.x.toSeq)
```

An iterative solver returns `SolverResult.Converged` or
`SolverResult.NotConverged`. Both results contain the last iterate, iteration
count, and residual. Check `converged` before using the solution when the
requested tolerance is a requirement.

## Define a matrix-free operator

Use `LinearOperator.fromFunction` when code can compute `A x` without storing
the entries of `A`. The function receives an input vector and a mutable output
vector. It must write every output element.

This operator applies a one-dimensional second-difference stencil:

```scala mdoc
val stencilSize = 5
val stencil = LinearOperator.fromFunction(stencilSize, stencilSize):
  (input, output) =>
    var i = 0
    while i < stencilSize do
      val left = if i == 0 then 0.0 else input(i - 1)
      val right = if i + 1 == stencilSize then 0.0 else input(i + 1)
      output(i) = 2.0 * input(i) - left - right
      i += 1

val expectedStencilSolution = Vec.tabulate(stencilSize)(i => i.toDouble + 1.0)
val stencilResult = cg(
  stencil,
  stencil * expectedStencilSolution,
  SolverConfig(tolerance = 1e-12, maxIterations = 20)
)

(stencilResult.converged, stencilResult.x.toSeq)
```

Use `LinearOperator.fromFunctions` when an algorithm also needs `Aᵀx`, as LSQR
and partial SVD do. Operator composition, scaling, row restriction, column
restriction, adjoints, and matrix-free Kronecker products are available without
materializing the combined matrix.

## Differences that affect migrated code

### Rows, slices, and copies

Breeze has a general slicing language built around `::`. Gale uses separate
methods for a row, column, contiguous slice, and indexed gather. Rows, columns,
slices, and transposes are views. Gathers are copies because they may reorder or
repeat values.

### Updating a dense value

Gale matrices do not change after construction. `updated` returns a new matrix;
it does not alter the original. A builder is intended for bulk construction.
Views can share data safely because ordinary Gale values cannot modify it.

`unsafeFromBreezeView` is different. It shares data with a mutable Breeze value,
so changing the Breeze matrix also changes what Gale reads.

### Error handling

Breeze commonly throws when a numerical operation fails. Gale's factorization,
solve, and spectral entry points return `Either[LinAlgError, A]` for structural
failures. Library code can inspect the error. Tests and small applications can
call `.orThrow` when stopping immediately is appropriate.

Iterative methods use a different result because a non-converged iterate can
still be useful. They return a result with the current estimate and diagnostics;
the caller decides whether the reported convergence is sufficient.

### Numerical conventions

`conditionEstimate` is a 1-norm estimate. It is not Breeze's exact SVD-based
2-norm condition number. Gale uses a different name because the returned
quantity and its computational cost differ. LU reports `SingularMatrix` only
on an exact-zero (or NaN) pivot; a reconstructed near-singular plant is not a
portable `+∞` fixture. See the
[numerical contract](../advanced/numerical-contract.md#cross-platform-singularity-rank-and-backend-residuals).

Gale's full SVD returns economy-size factors rather than full square factors.
Symmetric eigenvalue and Cholesky routines read the lower triangle of the input
as the symmetric matrix. If a migrated program may receive asymmetric data,
validate or symmetrize it before the call.

Signs, pivot choices, and bases inside repeated eigenspaces may differ even when
two decompositions represent the same mathematical result. Compare
reconstruction errors, residuals, orthogonality, and subspaces rather than raw
factor entries.

## Coverage at a glance

Differential checks against Breeze 2.1 for these common operations live in the
[`parity/` module](../../../parity/README.md) (`sbt parityTest`). Operations
with no honest Breeze reference — generalized eigen / GSVD / QZ, sparse-direct
factorization, near-cutoff rank / `pinv` / `cond`, and Krylov algorithm
diagnostics — are checked against NumPy / SciPy fixtures in the same module.
The table below is the migration surface; the parity README tracks which rows
are cross-checked.

| Task | Gale API | Important difference |
| --- | --- | --- |
| Dense vectors and matrices | `DVec`, `DMat`, `Vec`, `Matrix` | Real `Double` values; ordinary results are immutable. |
| Dense arithmetic | `+`, `-`, `*`, `dot`, `norm2`, `pointwise` | Pointwise matrix operations require `gale.syntax.all.*`. |
| Square solve with vector or matrix RHS | `A.solve(b)`, `A.solve(B)` | Returns `Either`; matrix RHS is factored once. |
| Least squares | `A.leastSquares(b)` or `A.leastSquares(B)` | Tall, full-column-rank systems; rank deficiency is reported. |
| LU, Cholesky, QR, determinant | `A.lu`, `A.cholesky`, `A.qr`, `A.det` | Typed factors and failures; legal signs and pivots may differ. |
| Rank | `A.rankEstimate` | Numerical estimate with Gale's tolerance policy. |
| Condition number | `A.conditionEstimate` | Estimates the 1-norm condition number, not the exact SVD 2-norm. |
| Symmetric eigenproblem | `Eigen.eigSymmetric` | Eigenvalues are returned in ascending algebraic order. |
| Nonsymmetric eigenproblem | `Eigen.eigNonsymmetric` | Complex conjugate pairs use Gale's typed complex accessors. |
| Generalized symmetric-definite eigenproblem | `Eigen.eigSymmetricGeneralized` | Requires a symmetric `A` and positive-definite metric `B`. |
| Partial matrix-free eigenproblem | operator overloads of `Eigen` | Returns diagnostics and the converged pairs. |
| Full or partial SVD | `A.svd`, `Svds.svd` | Full factors use economy shapes; partial results include diagnostics. |
| Pseudo-inverse | `A.pinv` | Rank cutoff may differ near the numerical-rank boundary. |
| Kronecker product | `A.kron(B)` or operator Kronecker product | The operator form avoids materializing the product. |
| Sparse matrix-vector products | `COO`, `CSR`, `CSC`, `Banded`, `Diagonal` | Sparse formats and canonicalization rules differ from Breeze. |
| Iterative linear solve | `cg`, `bicgstab`, `gmres`, `lsqr`, `cgnr` | Non-convergence returns the last iterate and diagnostics. |
| Matrix-free computation | `DoubleLinearOperator` | The caller supplies `A x`, and sometimes `Aᵀx`, without exposing entries. |
| Dense and sparse conversion | `gale.interop.breeze.*` | Gale-to-Breeze conversion always copies. |

## Move one boundary at a time

Add `gale-interop-breeze` while a codebase still accepts or returns Breeze
values. Convert to Gale before the numerical work and convert back only where
an existing interface still requires Breeze:

```scala
import gale.interop.breeze.*

val galeA = fromBreezeCopy(breezeMatrix)
val galeB = fromBreezeCopy(breezeVector)

val result = galeA.solve(galeB)
val breezeResult = result.map(x => toBreezeCopy(x))
```

The module also contains `BreezeMigration` shims for old call sites whose
Breeze-typed signatures cannot yet change. Those shims throw on failure. Treat
them as temporary adapters, not as the Gale API for new code.

`unsafeFromBreezeView` creates a zero-copy view when the Breeze value has
positive strides. Use it only when aliasing is intentional and the Breeze owner
will not mutate the data during the Gale operation. Gale-to-Breeze conversion
always copies.

## Performance

Do not assume that Gale or Breeze is faster for every operation. Measure the
operations, sizes, JVM, and backend used by the application.

In Gale's recorded pure-JVM comparison against Breeze 2.1's Java fallback,
Gale was faster on 23 of 42 operation-and-size pairs, within the stated parity
band on 2, and slower on 17. The result does not describe Breeze with native
BLAS/LAPACK enabled, and it does not predict another processor or JVM.
See the [release-grade Breeze comparison](https://github.com/canardlapin/gale/blob/main/benchmarks/results/2026-07-11-breeze-release-grade.md).

Gale also has optional JVM Vector API and native FFM backends. They route only
supported operations, layouts, and measured problem sizes. A small or strided
operation may still use the portable implementation. See the
[JDK 22 Vector comparison](https://github.com/canardlapin/gale/blob/main/benchmarks/results/2026-07-17-breeze-jdk22-vector-enabled.md)
and [Numerical contract](../advanced/numerical-contract.md) before choosing a
backend.

## What Gale does not replace

Gale does not provide:

- Breeze collections, generic scalar operators, broadcasting, or compatible
  `::` slicing syntax;
- probability distributions, optimization, statistics, signal processing,
  plotting, machine learning, tensors, or other non-linear-algebra modules;
- general complex matrix storage and arithmetic;
- sparse direct factorization;
- negative-stride Gale views;
- source or binary compatibility with Breeze;
- identical exception classes, pivots, factor signs, eigenvector signs, or
  bit-identical floating-point results; or
- a guarantee that Gale is faster than a native-enabled Breeze installation.

Continue with [Worked examples](examples.md) for more complete programs. Read
[Advanced topics](../advanced/index.md) when allocation control, sparse
structure reuse, backend selection, or matrix-free generalized eigensolving
affects the application.
