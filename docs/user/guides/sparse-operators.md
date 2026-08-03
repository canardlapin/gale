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
