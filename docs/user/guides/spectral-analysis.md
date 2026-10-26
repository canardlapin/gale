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

## Request a few eigenpairs and check the result

Passing an operator and its dimension selects the iterative overload. Here
only the largest eigenvalue of the earlier symmetric matrix is requested:

```scala mdoc
val partial = Eigen.eigSymmetric(
  symmetric, symmetric.rows,
  EigenSelection.Count(1, EigenOrder.LargestAlgebraic)
).orThrow

val complete = partial.requireConverged.orThrow
assert(complete.eigenvalues.length == 1)
assert(math.abs(complete.eigenvalues(0) - (2.0 + math.sqrt(2.0))) < 1e-10)
(partial.diagnostics.converged, partial.diagnostics.worstResidual)
```

An iterative run can return `Right` with fewer converged pairs than requested.
`requireConverged` checks completeness and returns a typed failure when pairs
are missing. Keep that `Either` in application code when partial output is
unacceptable; `.orThrow` above makes a failed tutorial expectation explicit.

Completeness does not always prove that the solver found the requested end of
the spectrum. A vector can satisfy `A * v = lambda * v` accurately while its
`lambda` is not the largest eigenvalue. An iterative method may explore a
subspace that misses that direction. Such an approximate eigenpair is called
a Ritz pair.

If the application also needs a certificate that the result belongs to the
requested spectral end, use the stricter check:

```scala mdoc
val globallyCertified = partial.requireExtremeCertified
globallyCertified.isRight
```

This check may reject residual-converged results. It does not rerun the solver
or improve the approximation; it checks the evidence already in the result.
Do not discard that failure if global extremality is required by your algorithm.

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
