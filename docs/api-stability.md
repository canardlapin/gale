# Gale 0.1 API stability boundary

This document identifies the public API that becomes compatible at
`v0.1.0-M1`. It is an inventory and review boundary, not a substitute for
Scaladoc or executable examples. The binary and TASTy baselines must be
generated from the published M1 artifacts after that immutable version exists.

## Admitted `gale-core` packages

| Package | Compatible contract from M1 |
| --- | --- |
| `gale.linalg` | dense values and views, builders, shapes, operators, factorizations, solve and least-squares entry points, typed errors, destinations, and workspaces |
| `gale.sparse` | COO/CSR/CSC and structured matrices, `SparseVector`, canonicalization, Matrix Market IO, compressed patterns, and symbolic replay plans |
| `gale.solvers` | iterative solver options, results, diagnostics, preconditioners, convergence semantics, and reusable workspaces |
| `gale.spectral` | dense and partial decompositions, typed selection, result ordering, convergence and extremality diagnostics, generalized operators, and explicit metric-solve contracts |
| `gale.optim` | optimization contracts, L-BFGS/L-BFGS-B, dense nonlinear least squares/LM, accelerated proximal gradient, standard terms/sets, work controls, diagnostics, and constrained-Rayleigh helpers |
| `gale.sized` | optional compile-time sized wrappers and their checked conversion boundary |
| `gale.backend` | caller-visible capabilities, backend selection, fallback behavior, configuration, thresholds, and factorization-provider contracts |
| `gale.syntax` | the `all` and `unicode` opt-in extension modules, including `zipMapExact`, pointwise operations, matrix product, and dot aliases |

The `gale.kernel` implementation and platform-storage operations are
`private[gale]` and are not public extension points. JVM Vector, native, FFM,
and Breeze modules are separate provisional artifacts; their APIs acquire no
0.1 compatibility promise until an admission record names a baseline.

## Observable semantic contracts

Compatibility covers more than method descriptors:

- `DVec` and `DMat` results own their returned mutable storage unless a method
  is explicitly documented as a view.
- Builders, destinations, and workspaces are single-owner mutable resources.
  Their methods do not expose Gale backing arrays.
- Total factorization and solve entry points return `Either[LinAlgError, A]`.
  Throwing conveniences remain explicitly named.
- Iterative and spectral results distinguish a returned approximation from
  residual convergence and from certification of a requested global extreme.
- Sparse canonicalization, duplicate handling, stored-zero behavior, and
  structure-preserving replay retain their documented meanings.
- Selection chooses spectral membership; result types retain their documented
  canonical ordering independently of the selected backend.
- Backends may change floating-point association but must satisfy the same
  shape, validation, ownership, residual, and conformance contracts.
- Work and allocation counters retain their units and inclusion rules. A field
  must not silently change from measured work to an estimate or proxy.

An incompatible change to one of these contracts after M1 requires `0.2.0`,
even when a binary compatibility checker cannot detect it.

## Published `gale-laws` API

The M1 law module admits these reusable entry points:

- `VecLaws` and `MatrixLaws`;
- `SparseLaws` and `SolverLaws`;
- `SpectralLaws` and `GeneralizedOperatorLaws`;
- `SpectralBackendLaws`; and
- `BackendConformanceSuite`.

Their public method signatures and assertion meanings become part of the 0.1
contract. Test fixtures and private helpers do not.

## Pre-M1 API-diff receipt

Known source migrations that must be settled before M1, not restored as
permanent aliases:

1. `DMatBuilder.updateRowMajor` became `writeLinear`. Callers in multivar and
   ScalaFIM still need to move.
2. `gale.syntax.all.zipMap` became `zipMapExact`. No remaining direct consumer
   of the old Gale name is known.

Before M1, compare the exact candidate commit against the last shared consumer
pin, review every exported addition and removal, and retain the complete diff
in the release evidence. Do not describe the candidate as binary-verified until
M1 artifacts exist and a post-M1 compatibility gate runs against them.

## Optimization additions before M1

Existing four first-order entry points and callback traits remain available.
Primal-dual calls add an optional dual start; configurations add optional work
controls. New stopping/error/method cases require consumers with exhaustive
matches to handle the additions. Primal-dual extrapolation is restricted to 1.
Summary-only certificates no longer bind matrices, and supplied certificates
must bind the exact result and matching objective. Rayleigh matrices must be
symmetric, with a positive-definite denominator. These are deliberate numerical
contract repairs, not binary-compatibility claims against a published baseline.

The LM/L-BFGS-B extension adds `LeastSquaresObjective`, `LeastSquaresSolution`,
`LevenbergMarquardtConfig`, `LevenbergMarquardt`, `MatrixBoxBounds`, `LBFGSBConfig`, and
`LBFGSB`. It adds two `FirstOrderMethod` cases, two `AlgorithmSettings` cases,
and residual/Jacobian callback counters to `EvaluationCounts`. Consumers that
pattern-match these enums or case-class products must review the additions.
The LM result uses raw gradient stationarity; L-BFGS-B uses a unit projected
gradient. Neither changes the meaning of existing fixed-step residuals.

CMA-ES adds separate `CMAES`, `CMAESSession`, `CMAESBatch`, configuration, control,
distribution, receipt, and result types. It reuses `MatrixBoxBounds` for vector
bounds but does not extend first-order method or stopping enums. Its stopping
reasons are sampling diagnostics, not stationarity certificates. Ask/tell has
one outstanding owned batch and rejects invalid receipts without consuming it.
Seeded reproducibility is runtime/version specific, not a cross-platform bitwise
compatibility promise.

Augmented Lagrangian adds `NonlinearConstraints`, `ConstraintEvaluation`,
`AugmentedLagrangian`, its configuration/result/status types, and
`ConstrainedDiagnostics`. It accepts a single-column `DMat` with analytic
derivatives, reuses `MatrixBoxBounds` and `SolverControl`, and does not extend
existing method/stopping enums. Diagnostics require raw and scaled feasibility,
box-projected Lagrangian stationarity, and complementarity. Fused constraint
callbacks use the existing `jacobians` work counter. These are first-order
conditions, not a certificate of a minimum or of infeasibility.

## Checked sparse builder finalization

`COOBuilder.tryToCSR` and `tryToCSC` now keep explicit stored zeros, including
duplicates that sum to zero, exactly as the total `toCSR` and `toCSC` always
have. Earlier releases pruned them on the checked path only, so the two paths
could disagree on `nnz` and on `hasCanonicalFormat`. Callers that relied on
the pruning should call `pruneZeros` on the result.

## Sparse vector additions before M1

`gale.sparse` adds `SparseVector`, its `SparseVectorEntryConsumer` callback,
`CSR.rowSparse`, `CSC.colSparse`, and `CSR * SparseVector`. The existing dense
`row`/`col` accessors and matrix products are unchanged. `gale.linalg` adds
`LinAlgError.EmptyInput`, raised by `SparseVector.max`/`min` on a length-0
vector; consumers that match `LinAlgError` exhaustively must handle it.

The stored-zero rule is part of the contract: construction from entries,
scaling, `mapActive`, `+` and `-` keep explicit zeros and cancellations, and
only `compact` (or `fromDense`) drops them. Products visit active entries
only, so a non-finite value facing an implicit zero is not propagated, while
`max`/`min` include the implicit zero whenever one exists.

Vector names follow the matrix vocabulary where the concept is shared:
`activeSize` is the vector counterpart of `nnz`, `compact` of `pruneZeros`,
`mapActive` of `mapValues`, `foreachActive` of `foreachStoredEntry`, and
`fromEntries`/`tryFromEntries` reuse `DuplicatePolicy` and `SparseValuePolicy`
from COO assembly. Unlike Breeze's in-place `compact()`, `compact` returns a
new immutable vector. `CSR * SparseVector` follows the `dot` rule: implicit
zeros never multiply stored values, so it can differ from
`A * x.toDense` when `A` stores a non-finite value.

## Everyday-ops additions before M1

`DVec` adds `sum`, `sumExact`, `mean`, `max`, `min`, `argmax`, `argmin`,
`norm1`, and `normInf`. `DMat` adds `sum`, `sumExact`, `mean`, `max`, `min`,
`argmax`/`argmin` (returning `(row, col)`), per-axis `sum`/`mean`/`max`/`min`
taking the new `gale.linalg.Axis` enum, and `norm1`, `normInf`, and
`normFrobenius`. `gale.linalg.Numerics` adds `exp`, `log`, `log1p`, `expm1`,
`sigmoid`, `logSumExp`, `softmax`, and `logSoftmax` for vectors and matrices,
with per-axis matrix forms of the last three. `Axis.Rows` means one result per
row and `Axis.Cols` one result per column.

`LinAlgError` gains the `EmptyInput(operation)` case. Consumers with exhaustive
matches on `LinAlgError` must handle it. Order statistics and means of empty
inputs throw it; sums and norms of empty inputs return `0.0`.

`DVec.norm2` now returns `+Inf`, rather than NaN, for a vector holding several
infinite entries and no NaN. The private dense 1-norm used by
`conditionEstimate` and the Frobenius norm used by the nonsymmetric
left-eigenvector residual guard now share the public implementations; the
Frobenius guard is now overflow-safe. For a row-major input the shared 1-norm
sums each column in the same order as before; for column-major or other strided
inputs it uses the unrolled kernel, so `conditionEstimate` can change in the
last ulp. Compact pivoted QR (width 8 or less) applies the same infinite-entry
repair as `norm2`, so its pivot choice matches the wider screened path. These are numerical-contract repairs, not
binary-compatibility claims against a published baseline.

## Deferred dense spectral diagnostics

`SpectralDiagnostics` is no longer a case class. It is a final class with the
same field names, the same eager `SpectralDiagnostics(...)` constructor and
defaults, and a `copy` with the same parameters. The synthesized `unapply`,
`Product` members and `canEqual` are gone; `equals`, `hashCode` and `toString`
keep the former case-class meaning (`residuals` still compares as a `DVec`).
The dense one-shot facades — `Eigen.eigSymmetric`, `Eigen.eigSymmetricWith`
and the full dense `Svds.svd` path (which `pinv` and minimum-norm least squares
share) — now measure `residuals` and `orthogonalityError` on first access and
cache them, instead of forming `A·V` and `VᵀV` (and, for SVD, `AᵀU`) on every
call. No convergence or certification decision reads either value on these
paths, so `converged`, `requireConverged` and `requireExtremeCertified` are
unchanged, and the measured values are bit-identical to the former eager ones.
The measurement uses snapshots taken when the result is built (the mirrored
lower triangle of the input, or a copy of the SVD input), so it is unaffected
by later mutation of a borrowed input view. The first read pays the two matrix
products; iterative, generalized and nonsymmetric solvers still report values
computed during the solve.

## Divide-and-conquer dense symmetric eigenvectors

From order 48, `Eigen.eigSymmetric` and `Eigen.eigSymmetricWith` with
eigenvectors solve the tridiagonal problem by divide and conquer (Cuppen with
`dlaed`-style deflation, a safeguarded secular solver and Gu–Eisenstat
vectors), then form `V = Q Z`. Below order 48, for values only, and for the
Lanczos and block-Krylov projected problems, the implicit QL solver is
unchanged. This changes result bits at those orders. The tested bounds are:
eigenvalues within `n·ε·max|λ|` of QL; residual `‖AV − VΛ‖_F ≤ n·ε·‖A‖_F` and
at most `max(4 × QL residual, 2√n·ε·‖A‖_F)`; orthogonality
`‖VᵀV − I‖_F ≤ 4n·ε`. The probes include clusters, repeated zeros, grading,
glued Wilkinson matrices and 1e±300 scaling. Both routes still run one kernel and agree
exactly, and the JVM and Scala.js agree because the merge products use
unfused arithmetic.

`Eigen.symmetricScratchRequirement(n, EigenVectors.Right)` therefore grows
for `n ≥ 48` from `n` doubles to `n + 2n² + 105n` doubles and `7n` indices
(the transposed tridiagonal eigenvectors, a packing region and per-merge
vectors and panels — the same order as LAPACK `dsyevd`). Callers that size a
`DenseWorkspace` from this requirement need no change; callers that
hard-coded the former `n` must re-query it. Values-only requirements are
unchanged. Above the largest order whose divide-and-conquer scratch is
addressable (`n + 2n² + 105n ≤ Int.MaxValue`, i.e. `n ≤ 32,741`), both
routes fall back to QL and the requirement is again `n` doubles, so vector
solves keep working up to QL's own `n² ≤ Int.MaxValue` limit. Non-finite input keeps the QL path and its `DidNotConverge`
behaviour.
