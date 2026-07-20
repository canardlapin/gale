# JVM sparse-direct provider boundary

Status: **contract only; no sparse-direct implementation is shipped or
advertised.** The public boundary exists only in the JVM build under
`gale.sparse.direct`. Scala.js has no corresponding API or provider.

`SparseDirectProvider.none` is the capability-less default. Its capability set
is empty and workspace creation returns `Left(UnsupportedOperation)`. Merely
having these types on the classpath does not mean sparse LU, Cholesky, or QR is
available.

## Staged contract

The boundary separates five lifetimes:

1. `SparseDirectProvider` is an explicitly selected, thread-safe singleton. Its
   capabilities state which factorization and optional solve features exist.
2. `SparseDirectWorkspace` owns mutable provider scratch and native temporary
   resources. It is single-user, provider-specific, grow/reuse policy is owned by
   the provider, and it is never safe for concurrent calls.
3. `SparseDirectSymbolicAnalysis` owns ordering and symbolic factor structure for
   one exact immutable `CSRPattern`.
4. `SparseDirectNumericFactor` binds changing values to that analysis and owns
   the resulting numeric factor handle.
5. `SparseDirect.solve` and `solveInto` apply a numeric factor to vector or matrix
   right-hand sides and return per-solve diagnostics.

The intended flow is:

```scala
import gale.sparse.direct.*

// A provider module supplies an explicit given SparseDirectProvider.
val workspace = SparseDirect.newWorkspace().toOption.get
val symbolic = SparseDirect
  .analyze(a.pattern, SparseDirectFactorization.LU, workspace)
  .toOption
  .get

val factor1 = SparseDirect.factor(symbolic, a, workspace).toOption.get
val factor2 = SparseDirect.factor(symbolic, aWithNewValues, workspace).toOption.get

val solved = SparseDirect.solve(factor2, rhs, workspace).toOption.get
```

`aWithNewValues` must share the exact compressed storage analyzed. Build it with
`symbolic.pattern.bind(values)`, `CSR.rebind`, or another pattern-preserving
numeric operation. A structurally similar but independently allocated CSR is
rejected before provider code runs, preventing stale native mappings from being
applied to reordered indices.

## Capability and selection rules

Capabilities are operation-specific:

- `LU`, `Cholesky`, and `QR` authorize the corresponding symbolic analysis;
- `UserOrdering` authorizes `SparseDirectOrdering.User`;
- `TransposeSolve` and `MultipleRhs` authorize those solve forms.

Provider selection is an explicit `given SparseDirectProvider`. It is deliberately
separate from the portable `Backend` and from `Capability.NativeSparse`: no Gale
sparse operation routes here implicitly, and composing dense backends cannot
silently activate sparse direct code. An eventual optional module may expose a
provider value and a `given`, but users must import it intentionally.

LU and Cholesky analysis is square-only at the facade. QR may be rectangular;
its factor reports the accepted right-hand-side and solution dimensions for
normal and transpose operations. Providers remain responsible for family-specific
preconditions such as structural/numeric symmetry and positive definiteness.

## Ordering and permutation ownership

`Natural` is the identity column ordering. `ProviderDefault` lets an
implementation choose AMD, COLAMD, nested dissection, or another named internal
policy. `User(columnPermutation)` maps new column position to original column.

Symbolic analysis reports the column permutation actually used. Numeric factors
report both the symbolic column permutation and any value-dependent row
permutation (for example LU pivoting). These are Gale `Permutation` values. A
provider must copy native index buffers before returning them; neither native
storage nor a mutable aliased array may escape.

## Workspace, lifecycle, and threads

Provider, workspace, analysis, and factor ownership is validated by object
identity. A workspace created by one provider cannot be passed to another.
`close()` is idempotent; every facade call rejects a closed resource with a typed
error before entering provider code.

A shared provider must be safe for concurrent invocation and fixes its
`BackendConfig` at construction. Symbolic analyses and numeric factors may be
used concurrently only when each invocation has a distinct workspace. A single
workspace is sequential mutable state. A numeric factor has an independent
lifetime: closing its symbolic analysis after factorization does not invalidate
the factor. Native implementations must use independent handles or reference
counting to uphold that rule.

The default configuration is single-threaded. An implementation that uses native
threads must declare them in `BackendConfig`, prevent accidental nested
oversubscription, and document whether the linked solver also reads process-wide
thread environment variables.

## Failure and diagnostics

Structural and lifecycle failures use existing `LinAlgError` values:

- unavailable family/solve feature: `UnsupportedOperation`;
- changed pattern, wrong provider/workspace, malformed returned handle, or
  invalid ordering: `InvalidArgument`;
- nonsquare LU/Cholesky: `NonSquareMatrix`;
- numeric failures: the specific existing error when available
  (`SingularMatrix`, `NotPositiveDefinite`, `RankDeficient`), otherwise a
  provider-labelled `InvalidArgument` rather than an exception code leak.

Symbolic, numeric, and solve diagnostics are distinct immutable records. They
carry provider/family identity, fill/rank/pivot/condition information when the
implementation can support it, and solve residual/refinement information. No
mutable “last solve” diagnostics live on a shared factor.

## Implementation go/no-go gate

The current decision is **no-go** for enabling any provider: no candidate has yet
supplied the evidence below. A SuiteSparse wrapper, a faer-inspired implementation,
or a new pure provider becomes eligible only when all applicable gates pass.

| Gate | Required evidence before capability is advertised |
| --- | --- |
| Packaging | Separate optional JVM module; no native/JNI/FFM dependency in shared core or Scala.js; supported OS/architecture matrix and licence audit recorded. |
| Capability honesty | Advertise only independently implemented families and optional features. Missing symbols or unsupported matrix classes produce empty/missing capabilities, not late linkage failures. |
| Contract conformance | Deterministic provider suite plus independent dense/reference oracles for ordering, pivots, normal/transpose and multiple-RHS solves, singular/rank-deficient/indefinite failures, zero-sized matrices, and malformed inputs. |
| Symbolic reuse | Changing values on one exact pattern demonstrably reuses symbolic state; changed patterns fail before native execution; factor lifetime remains independent from analysis lifetime. |
| Resource safety | Repeated create/factor/solve/close stress shows no native-memory, file-descriptor, or thread leak; idempotent close and wrong-workspace rejection are tested. |
| Threading | Concurrent solves with distinct workspaces pass race/stress tools; configured thread counts avoid JVM/native oversubscription. |
| Performance | Two-fork scenario sweeps across representative size, density, fill, and RHS counts include analysis and conversion costs. A measured reuse crossover and conservative routing guidance are documented; one favourable matrix is insufficient. |
| Distribution | CI loads the actual packaged artifacts on every supported OS, verifies missing-library diagnostics, and exercises an installed consumer rather than only the source checkout. |

For SuiteSparse specifically, symbol discovery, index width, ABI, library
version, and ordering availability must be explicit. “Faer-inspired” means
adopting sound symbolic/numeric separation and ownership ideas; calling Rust
through an unmeasured bridge or transliterating algorithms without independent
oracles does not pass the gate. A pure JVM provider must meet the same correctness
and crossover evidence, not a lower bar because it avoids native packaging.
