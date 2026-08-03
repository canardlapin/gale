# Solve dense systems and least squares

Use this guide when a coefficient matrix or design fits in dense memory. It
starts with allocating one-shot calls, then introduces retained factors,
pivoting, row scaling, and workspaces only when repeated work creates a reason
for them.

```scala mdoc:silent
import gale.linalg.*
```

## Solve a square system

`solve` factors the matrix once and accepts a vector or a matrix of right-hand
sides:

```scala mdoc
val system = Matrix(2, 2)(
  3.0, 1.0,
  1.0, 2.0
)
val rhs = Vec(9.0, 8.0)

system.solve(rhs).map(_.toSeq)
```

```scala mdoc
val rhsColumns = Matrix(2, 2)(
  9.0, 4.0,
  8.0, 5.0
)

system.solve(rhsColumns).map(_.valuesRowMajor)
```

Retain `system.lu.orThrow` when separate calls reuse the same general square
matrix. Retain `system.cholesky.orThrow` when the matrix is symmetric positive
definite and that premise is part of the application.

## Fit a full-rank design

The design below has an intercept and slope. `QRPivoting.Column` recomputes
candidate norms exactly, reports the selected permutation, and applies the
same deterministic rank policy on the JVM and Scala.js.

```scala mdoc
val design = Matrix(5, 2)(
  1.0, 0.0,
  1.0, 1.0,
  1.0, 2.0,
  1.0, 3.0,
  1.0, 4.0
)
val observations = Vec(1.0, 3.0, 5.0, 7.0, 9.0)
val options = QROptions(pivoting = QRPivoting.Column)
val qr = design.qr(options)

(
  qr.solveLeastSquares(observations).orThrow.toSeq,
  qr.diagnostics.rank,
  qr.columnPermutation.toIndexSeq
)
```

With pivoting, `qr.q * qr.r` reconstructs the input columns in
`qr.columnPermutation` order. Returned coefficients are unpermuted back to the
original design order. A rank-deficient or underdetermined least-squares solve
returns a typed `Left`; the factorization itself remains available for rank and
orthogonal-transform operations.

## Reuse a factor and scratch

`solveLeastSquaresWith` transforms the right-hand side in caller-owned scratch
and returns an independently owned coefficient vector or matrix:

```scala mdoc
val solveWorkspace = DenseWorkspace.forQRSolve(
  observations = design.rows,
  rightHandSides = 2
)
val responseMatrix = Matrix.tabulate(design.rows, 2): (row, col) =>
  observations(row) + col.toDouble

val firstFit = qr
  .solveLeastSquaresWith(responseMatrix, solveWorkspace)
  .orThrow
val secondFit = qr
  .solveLeastSquaresWith(responseMatrix, solveWorkspace)
  .orThrow

(firstFit.valuesRowMajor, secondFit.valuesRowMajor)
```

Reusing `solveWorkspace` cannot mutate `firstFit`. The workspace grows when a
later request needs more primitive cells and never shrinks. It is sequential
mutable state, so concurrent workers need separate instances.

## Apply row scales without a temporary matrix

An outer algorithm may need the algebraic system
`diag(scales) * design` and the matching response
`diag(scales) * observations`. Gale can construct and solve that QR without
materializing either scaled input:

```scala mdoc
val rowScales = Vec(1.0, 1.0, 0.5, 0.5, 0.25)
val pipelineWorkspace = DenseWorkspace.empty
val scaledQr = design
  .qrScaledRows(rowScales, options, pipelineWorkspace)
  .orThrow

val scaledFit = scaledQr
  .solveLeastSquaresScaledRowsWith(
    observations,
    rowScales,
    pipelineWorkspace
  )
  .orThrow

scaledFit.toSeq
```

These are **row multipliers**, not statistical weights. Zero and negative
finite scales are algebraically valid; non-finite or wrongly sized scales are
rejected. If an application starts from non-negative weights, it owns the
decision to pass their square roots.

## Consume a transient builder

When a design is filled once and immediately factored, a builder can transfer
its row-major storage directly into QR's working factor:

```scala mdoc
val transientDesign = design.toBuilder
val consumedQr = transientDesign.consumeQR(options, DenseWorkspace.empty)

consumedQr.diagnostics.rank
```

`consumeQR` permanently closes the builder. It does not leave a mutable alias
to the factor. Use `builder.result()` instead when the immutable matrix itself
must remain available.

## Choose the route deliberately

| Need | Route | Execution and ownership |
| --- | --- | --- |
| one general solve | `a.solve(rhs)` | allocating owned result; typed failure |
| several square solves | retain `a.lu.orThrow` or `a.cholesky.orThrow` | factor once, owned results |
| one least-squares fit | `design.leastSquares(y, options)` | portable explicit QR policy |
| several fits | retain `design.qr(options)` | reuse factor; owned results |
| reduce transformed-RHS allocation | `qr.solveLeastSquaresWith(rhs, workspace)` | borrowed scratch, owned result |
| avoid scaled design/RHS temporaries | `qrScaledRows` plus scaled solve | same row scales on both algebraic inputs |
| factor a transient fill buffer | `builder.consumeQR` | storage transfer; builder closes |

Read [Reusable workspaces](../advanced/workspaces.md) for requirement
composition and [Factorization capabilities](../reference/factorization-capabilities.md)
when a generic algorithm consumes retained factors.
