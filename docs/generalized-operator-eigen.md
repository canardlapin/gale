# Matrix-free generalized symmetric eigensolving

Gale provides a partial solver for the real symmetric-definite pencil

```text
A x = λ B x
```

when `A` is symmetric and `B` is symmetric positive-definite. The portable
engine is LOBPCG. Both inputs are `DoubleLinearOperator`s: Gale applies them to
blocks of vectors and never constructs dense `A`, dense `B`, or an implicit
`B^-1`.

This is a numerical capability, not an alignment-specific abstraction. A
downstream alignment library can define its own operators and consume Gale's
eigenpairs without putting alignment concepts in Gale.

## Public contract

```scala
import gale.linalg.*
import gale.solvers.Preconditioner
import gale.spectral.*

val result = Eigen.eigSymmetricGeneralized(
  a.assumeSymmetricOperator,
  b.assumePositiveDefiniteOperator,
  n,
  EigenSelection.Count(k, EigenOrder.SmallestAlgebraic),
  GeneralizedSpectralOptions(
    tolerance = 1e-8,
    maxIterations = 200,
    returnVectors = EigenVectors.Right
  ),
  Preconditioner.Identity
)
```

The property wrappers are explicit evidence supplied by the caller. Gale cannot
inspect a matrix-free operator exhaustively, but it still detects encountered
non-positive `B` geometry, non-finite operator output, shape errors, and invalid
options through `LinAlgError`.

The operator overload accepts:

- `Count(k, SmallestAlgebraic)` or `Count(k, LargestAlgebraic)`;
- `1 <= k < n`;
- an optional owned `n × k` initial subspace;
- right eigenvectors or values only.

Magnitude, both-ends, interval, and shift-invert selections remain on other
spectral routes. The dense `DMat` overload supports the broader dense selection
surface but computes a dense generalized decomposition even for `Count`.

## Result geometry and convergence

Returned eigenvalues are ascending-algebraic. With vectors requested, the
columns of `X` satisfy the metric normalization

```text
Xᵀ B X = I
```

up to the reported `diagnostics.orthogonalityError`. Per-pair diagnostics are
the true ambient residual norms `||A x - λ B x||₂`, not projected residual
estimates.

Iteration exhaustion is an honest `Right` containing only pairs whose residual
passed the requested tolerance. Consequently `diagnostics.converged` may be
smaller than `diagnostics.requested`, including zero. `worstResidual` is zero
when no pair is returned; inspect `converged` before interpreting it.

Residual convergence alone does not prove that an arbitrary invariant initial
subspace contains the requested global spectral end. Gale performs one
deterministic complement exploration for an exactly invariant partial start,
but callers needing the stronger guarantee must inspect
`diagnostics.extremalityCertified` or call `requireExtremeCertified`.

## Preconditioner semantics

The preconditioner acts directly on each residual column. It is not a metric
solver and Gale never infers `B^-1` from it. `Preconditioner.Identity` is the
portable default. A useful preconditioner approximates the inverse action of the
operator near the requested spectral region while preserving the
symmetric-positive geometry assumed by LOBPCG; Jacobi and block-Jacobi are
common choices for suitable stiffness operators.

Application is columnwise and explicitly accounted. The kernel caches aligned
`A X` and `B X` blocks, transforms those images during Rayleigh-Ritz steps, and
does not reapply an operator merely to assemble ordinary LOBPCG diagnostics.
When an external iterative generalized backend returns raw Ritz pairs, the
facade deliberately reapplies `A` and `B` once per returned pair to establish
owned, independently derived final diagnostics.

## Backend capability

`SpectralCapability.IterativeGeneralized` is an explicit backend capability.
It is not inferred from dense generalized eigendecomposition, partial ordinary
eigen, or SVD support. A capable provider receives the same typed operators,
selection, options, and preconditioner. A provider `Left` declines and falls
back to portable LOBPCG. A malformed provider `Right` fails loudly after Gale
checks shape, finiteness, ordering, metric normalization, residuals, and
orthogonality.

## Performance evidence and limits

The fixed-seed JMH matrix covers `n = 128, 512, 2048`, `k = 4, 8, 16`,
clustered diagonal and finite-difference stiffness/mass pencils, and identity,
Jacobi, and block-Jacobi preconditioners. It records time and allocation with
JMH and exact iteration/operator work through a separate untimed receipt.
See the
[development baseline](../benchmarks/results/2026-07-25-generalized-lobpcg-baseline.md).

The current engine is LOBPCG. A generalized block-Lanczos engine requires an
explicit metric-solve contract so it can state whether and how `B^-1` is
applied. Production-scale native/provider performance remains a separate
backend evidence claim; the pure operator implementation does not imply one.
