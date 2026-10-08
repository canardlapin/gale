# Solve with sparse matrices and operators

Use a sparse matrix when stored positions matter. Use a
`DoubleLinearOperator` when code only needs to compute `A * x`. Both can feed
the same iterative solver APIs.

```scala mdoc:silent
import gale.linalg.*
import gale.solvers.*
import gale.sparse.*
```

## Build a compressed sparse matrix

Coordinate form is convenient for assembly. Convert to CSR for repeated
row-oriented matrix-vector products:

```scala mdoc
val stiffness = Sparse
  .coo(3, 3)
  .add(0, 0, 2.0)
  .add(0, 1, -1.0)
  .add(1, 0, -1.0)
  .add(1, 1, 2.0)
  .add(1, 2, -1.0)
  .add(2, 1, -1.0)
  .add(2, 2, 2.0)
  .toCSR()

(stiffness.nnz, (stiffness * Vec(1.0, 2.0, 3.0)).toSeq)
```

The default duplicate policy sums repeated coordinates. Checked ingestion via
`Sparse.cooChecked` keeps malformed coordinates, non-finite values, and a
requested `DuplicatePolicy.Error` in `Either[LinAlgError, A]` rather than the
throwing convenience path.

## Work with sparse vectors

`SparseVector` stores a fixed-length vector as strictly increasing indices and
their values. Build one from `(index, value)` pairs in any order; repeated
indices follow the same `DuplicatePolicy` as matrix assembly:

```scala mdoc
val features = SparseVector.fromEntries(8, Seq(5 -> 2.0, 1 -> -1.0, 5 -> 0.5, 3 -> 0.0))
val weights = SparseVector.fromDense(Vec(0.0, 3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 1.0))

(features.activeIndices, features.activeValues, features.dot(weights))
```

Explicit zeros stay stored, as in Breeze: index 3 above is active although its
value is `0.0`. `+`, `-`, scaling and `mapActive` also keep their patterns when
a value cancels. `compact` is the one operation that drops stored zeros, and
`fromDense` never stores them:

```scala mdoc
val cancelled = features - features
(cancelled.activeSize, cancelled.compact.activeSize)
```

Products (`dot`, `axpyInto`) visit active entries only. Reductions describe the
dense vector, so `max` and `min` include `0.0` whenever an implicit zero exists,
and a length-0 vector has no maximum:

```scala mdoc
val negative = SparseVector.fromEntries(4, Seq(0 -> -3.0, 2 -> -1.0))
(negative.max, negative.min, negative.norm2, SparseVector.zeros(0).sum)
```

`mapActive` applies its function to stored values only; it equals a dense map
only when the function sends `0.0` to `0.0`. Checked construction via
`SparseVector.tryFromEntries` returns `Either[LinAlgError, SparseVector]` for
out-of-range indices, a rejected duplicate, or a non-finite value under
`SparseValuePolicy.RequireFinite`.

CSR rows and CSC columns convert without a dense intermediate, and a CSR
matrix multiplies a sparse vector into a dense result:

```scala mdoc
val firstRow = stiffness.rowSparse(0)
val unitLoad = SparseVector.fromEntries(3, Seq(2 -> 1.0))
(firstRow.activeIndices, firstRow.activeValues, (stiffness * unitLoad).toSeq)
```

## Solve the sparse system iteratively

Conjugate gradient is appropriate when the operator is symmetric positive
definite:

```scala mdoc
val expected = Vec(1.0, 2.0, 3.0)
val result = cg(
  stiffness,
  stiffness * expected,
  SolverConfig(tolerance = 1e-12, maxIterations = 20)
)

(result.converged, result.iterations, result.x.toSeq)
```

`bicgstab` and restarted `gmres` handle general nonsymmetric square systems.
`lsqr` handles rectangular least squares without forming `A.t * A`; `cgnr`
solves the normal equations and therefore squares the condition number.

Iteration exhaustion is not a structural error. The result retains the last
iterate, residual, and iteration count with `converged = false`. Decide at the
application boundary whether that approximation is usable.

## Define a matrix-free operator

This second-difference stencil applies the same mathematical operator without
storing its entries:

```scala mdoc
val order = 5
val stencil = LinearOperator.fromFunction(order, order): (input, output) =>
  var row = 0
  while row < order do
    val left = if row == 0 then 0.0 else input(row - 1)
    val right = if row + 1 == order then 0.0 else input(row + 1)
    output(row) = 2.0 * input(row) - left - right
    row += 1

val known = Vec.tabulate(order)(index => index.toDouble + 1.0)
val matrixFreeResult = cg(
  stencil,
  stencil * known,
  SolverConfig(tolerance = 1e-12, maxIterations = 20)
)

(matrixFreeResult.converged, matrixFreeResult.x.toSeq)
```

The callback must write every destination element and must not retain mutable
destination storage. Use `LinearOperator.fromFunctions` when an algorithm also
needs a transpose action, as LSQR and partial SVD do.

## Reuse storage only after measurement

`CgWorkspace` keeps repeated same-size solve storage:

```scala mdoc
val cgWorkspace = CgWorkspace(order)
cgWith(
  stencil,
  stencil * known,
  cgWorkspace,
  SolverConfig(tolerance = 1e-12, maxIterations = 20)
)

val stableSolution = cgWorkspace.solution
val borrowedSolution = cgWorkspace.unsafeSolutionView

(stableSolution.toSeq, borrowedSolution.toSeq)
```

`solution` is an owned snapshot. `unsafeSolutionView` is allocation-free and
changes when the workspace is reused; its name is the lifetime warning.

Continue with [Compressed sparse patterns](../advanced/sparse-patterns.md) and
[Symbolic sparse plans](../advanced/sparse-plans.md) when structure stays fixed
across many numeric evaluations.
