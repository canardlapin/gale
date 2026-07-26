# Matrix-free generalized symmetric eigensolving

Gale provides a partial solver for the real symmetric-definite pencil

```text
A x = λ B x
```

when `A` is symmetric and `B` is symmetric positive-definite. The default
portable engine is LOBPCG. Gale also provides an explicitly named generalized
block-Lanczos engine for callers that can supply a metric solve. Both inputs are
`DoubleLinearOperator`s: Gale applies them to blocks of vectors and never
constructs dense `A`, dense `B`, or an implicit `B^-1`.

This is a numerical capability, not an alignment-specific abstraction. A
downstream alignment library can define its own operators and consume Gale's
eigenpairs without putting alignment concepts in Gale.

## Public contract

This diagonal pencil is intentionally small enough to inspect, but the call
uses only the operator interface. Replacing either dense value with a custom
`DoubleLinearOperator` does not change the eigensolver contract.

```scala mdoc
import gale.linalg.*
import gale.solvers.*
import gale.spectral.*

val n = 6
val k = 2
val expected = IndexedSeq(1.0, 2.0, 4.0, 7.0, 11.0, 16.0)
val metricDiagonal = IndexedSeq(1.0, 3.0, 0.5, 2.0, 4.0, 1.5)

val denseB = DMat.tabulate(n, n): (row, col) =>
  if row == col then metricDiagonal(row) else 0.0

val a = DMat.tabulate(n, n): (row, col) =>
  if row == col then expected(row) * metricDiagonal(row) else 0.0
val b = denseB

val lobpcg = Eigen.eigSymmetricGeneralized(
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
).orThrow

lobpcg.eigenvalues.toSeq
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

## Metric-solve contract

Generalized Lanczos needs the action of `B^-1`, which is a linear **solve**, not
an ordinary operator application and not a LOBPCG preconditioner. Gale represents
that distinction with `LinearSolveOperator`. Each call returns an owned solution
and `LinearSolveDiagnostics` containing convergence, inner iterations, measured
residual when available, and the exact number of system-operator applications.
Structural failures are `Left(LinAlgError)`; iterative exhaustion returns its
best iterate with `converged = false` so an outer eigensolver can identify the
failure as an inner solve rather than an outer convergence failure.

A direct solve adapts a factor the caller already created:

```scala mdoc:silent
val factor = denseB.cholesky.orThrow
val direct = LinearSolveOperator.direct(factor)
val metric = MetricSolveOperator
  .bind(denseB.assumePositiveDefinite, direct)
  .orThrow
```

No Gale spectral call factorizes `denseB` as a side effect. For a matrix-free
SPD metric, the caller can prepare a repeated CG solve:

```scala mdoc:silent
val iterativeSolve = LinearSolveOperator
  .conjugateGradient(
    b.assumePositiveDefiniteOperator,
    SolverConfig(tolerance = 1e-10, maxIterations = 500),
    toleranceMode = ToleranceMode.RelativeToRhs
  )
  .orThrow

val iterativeMetric = MetricSolveOperator
  .bind(b.assumePositiveDefiniteOperator, iterativeSolve)
  .orThrow
```

`MetricSolveOperator` binds the solve to the metric it claims to invert and
checks their dimensions. It never forms an inverse. If `A` is symmetric and the
solve implements `B^-1`, then `T = B^-1 A` is self-adjoint in the B inner
product:

```text
<x, T y>_B = xᵀ A y = (A x)ᵀ y = <T x, y>_B.
```

The same executable boundary is used by `SpectralTarget.ShiftInvert`.
`LinearSolvePlan.Use(solver)` supplies a prepared solve;
`LinearSolvePlan.Backend` explicitly asks a backend advertising
`ShiftInvertSolve`. There is no string method flag and no automatic
factorization.

## Generalized block Lanczos

Lanczos is a separate public entry point, so adding it does not change the
LOBPCG route or overload a generic `method = "..."` option:

```scala mdoc
val lanczos = Eigen.eigSymmetricGeneralizedLanczos(
  a.assumeSymmetricOperator,
  metric,
  n,
  EigenSelection.Count(k, EigenOrder.SmallestAlgebraic),
  GeneralizedLanczosOptions(
    tolerance = 1e-8,
    maxIterations = 100,
    subspaceDimension = Some(math.min(n, 4 * k))
  )
).orThrow

lanczos.eigenvalues.toSeq
```

The engine applies `A`, solves `B y = A x`, and then applies `B` to the new
directions for explicit B-orthogonalization. It uses two-pass
reorthogonalization, a multiplicity-safe block at least as wide as `k`,
deterministic complement replenishment after invariant breakdown, thick
restarts retaining the wanted Ritz block, and soft locking that keeps converged
vectors in the basis without expanding them again. Rayleigh-Ritz work is dense
only in the retained subspace.

`GeneralizedLanczosOptions` accepts the same algebraic ends, tolerance,
iteration, initial `n × k` block, and vector modes as LOBPCG, plus an explicit
Krylov `subspaceDimension`. `None` uses `min(n, max(4k, 20))`; an explicit value
must be in `[k + 1, n]`.

Outer iteration exhaustion follows Gale's normal partial-spectral contract:
`Right` contains exactly the true-residual-passing pairs. Inner work is separate
in `diagnostics.innerSolve`, aggregating solve calls, converged solves, inner
iterations, system-operator applications, and the worst measured residual.
An iterative inner solve that does not converge returns
`Left(InnerSolveDidNotConverge)` with both outer and inner work. A structural
inner failure returns `Left(InnerSolveFailed)` and retains its typed cause.

## Backend capability

`SpectralCapability.IterativeGeneralized` is an explicit backend capability.
It is not inferred from dense generalized eigendecomposition, partial ordinary
eigen, or SVD support. A capable provider receives the same typed operators,
selection, options, and preconditioner. A provider `Left` declines and falls
back to portable LOBPCG. A malformed provider `Right` fails loudly after Gale
checks shape, finiteness, ordering, metric normalization, residuals, and
orthogonality.

## Performance evidence and limits

Gale's fixed-seed benchmark matrix covers multiple problem sizes, requested
ranks, clustered and finite-difference pencils, and identity, Jacobi, and
block-Jacobi preconditioners. It records both time/allocation and exact
iteration/operator work. Treat those measurements as evidence for the measured
runtime and problem families, not as a universal engine ranking.

LOBPCG remains the default generalized engine because it requires only a
preconditioner and is generally more forgiving of inexact inverse information.
Generalized block Lanczos is available explicitly when the caller has a useful
metric solve. Gale's comparison found a narrow exact-solve Lanczos advantage on
some clustered requests, but no consistent lower-allocation crossover in
convergence-equivalent cells. LOBPCG was more robust across the full matrix,
while iterative metric solves could add substantial `B` work. Gale therefore
keeps the engines as separate named entry points and does not route
automatically from problem dimensions or operator types.

Production-scale native/provider performance remains a separate backend
evidence claim; the pure operator implementations do not imply one.
