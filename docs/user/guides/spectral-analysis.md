# Compute and interpret a spectrum

Choose a dense route when the matrix fits in memory and all or many spectral
components are needed. Choose a partial operator route when only a few values
are needed and matrix application is cheaper than materialization.

```scala mdoc:silent
import gale.linalg.*
import gale.spectral.*
```

## Solve a dense symmetric eigenproblem

The dense symmetric facade reads the lower triangle and returns eigenvalues in
ascending algebraic order:

```scala mdoc
val symmetric = Matrix(3, 3)(
  2.0, -1.0, 0.0,
  -1.0, 2.0, -1.0,
  0.0, -1.0, 2.0
)

val eigen = Eigen
  .eigSymmetric(
    symmetric,
    EigenSelection.All,
    EigenVectors.Right
  )
  .orThrow

(
  eigen.eigenvalues.toSeq,
  eigen.diagnostics.converged,
  eigen.diagnostics.worstResidual
)
```

Signs of individual eigenvectors and bases inside repeated eigenspaces are not
identity guarantees. Check residuals, orthogonality, and invariant subspaces
rather than comparing raw vector columns.

## Compute an economy SVD

For an `m x n` dense matrix, the full facade returns economy factors with
`k = min(m, n)`: `U` is `m x k` and `V.t` is `k x n`.

```scala mdoc
val rectangular = Matrix(4, 2)(
  1.0, 0.0,
  1.0, 1.0,
  1.0, 2.0,
  1.0, 3.0
)
val svd = rectangular.svd.orThrow

(
  svd.singularValues.toSeq,
  (svd.u.rows, svd.u.cols),
  (svd.vt.rows, svd.vt.cols)
)
```

`Svds.svd` supplies partial singular values and an operator overload when a
transpose action is available. Singular values are returned in descending
order.

## Interpret partial results

Partial eigensolvers and SVD return diagnostics because a legal run can finish
before every requested component converges. `requireConverged` converts an
incomplete result to a typed failure when the application requires all pairs.

For a requested spectral extreme, residual convergence and global membership
are separate claims. `requireExtremeCertified` additionally requires Gale's
extremality certificate; it can reject a residual-converged Ritz pair that is
not proven to belong to the requested end.

## Choose the generalized route

For a dense symmetric-definite problem `A x = lambda B x`, use the dense
`Eigen.eigSymmetricGeneralized` facade when the matrices fit and a full dense
decomposition is acceptable.

For a few eigenpairs from operators, use LOBPCG when you can supply operator
applications and an optional preconditioner. Use generalized block Lanczos
when you can supply a genuine metric solve for `B`. The complete typed
selection, metric, convergence, and work-accounting contract is in
[Matrix-free generalized symmetric eigensolving](../advanced/generalized-operator-eigen.md).

Read the [Numerical contract](../advanced/numerical-contract.md) before relying
on ordering, tolerance, backend equivalence, or deterministic layout choices.
