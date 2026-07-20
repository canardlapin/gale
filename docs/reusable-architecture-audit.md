# Reusable structure and allocation-control acceptance audit

Date: 2026-07-20

This audit closes the epic that adopted the most useful structural lessons from
faer and nalgebra while retaining Gale's `Double`-specialized JVM/Scala.js
model. The result is a small opt-in reuse tier beneath the existing allocating
APIs, not a Rust-shaped ownership, scalar, or storage hierarchy.

## Outcome audit

| Required outcome | Delivered contract | Primary evidence |
| --- | --- | --- |
| Compressed sparse structure independent from values | Immutable checked `CSRPattern` and `CSCPattern`; owned rebinding; single-owner value builders; structural transpose | `CompressedPatternSuite`, [compressed-pattern contract](sparse-patterns.md) |
| Safe dense destination tier | `gemmInto` and `linearCombinationInto` on open `DMatBuilder` destinations, with shape, lifetime, and alias checks | `DenseDestinationSuite`, [dense-destination contract](dense-destinations.md) |
| General reusable workspaces | Checked `ScratchRequirement` composition and grow-only `DenseWorkspace`; QR, symmetric eigen, and CSR canonicalization reuse | `WorkspaceSuite`, `DenseSymmetricWorkspaceSuite`, `SparseWorkspaceSuite`, [workspace contract](workspaces.md) |
| Exact-solve and least-squares capabilities | Narrow `ExactSolveFactor` and `LeastSquaresFactor` traits without leaking concrete-only operations | `FactorizationCapabilitySuite`, [capability contract](factorization-capabilities.md) |
| Symbolic sparse plans | Checked `CSRUnionPlan` and `CSRProductPlan`, exact-pattern validation, reusable numeric destinations, explicit-zero preservation | `SparsePlanSuite`, [symbolic-plan contract](sparse-plans.md) |
| Optional JVM sparse-direct seam | Explicit provider, symbolic, numeric, workspace, solve, capability, diagnostic, and lifecycle contracts; capability-less default | `SparseDirectProviderSuite`, [provider boundary and go/no-go gates](sparse-direct-provider.md) |

All ordinary allocating APIs remain available. The reuse tier is explicit and
does not silently alter backend selection, numerical ordering, diagnostics, or
failure behavior.

## Cross-platform and safety matrix

| Acceptance area | Evidence |
| --- | --- |
| Row-major, transposed, and strided dense inputs | Dense destination tests cover all three layouts and compare with allocating dense oracles on JVM and Scala.js. |
| Ownership and alias safety | Pattern constructors copy caller arrays; numeric rebinding owns one copy; closed builders, closed provider workspaces, and dense aliases fail before execution; plan snapshots remain unchanged after destination reuse. |
| Invalid inputs | Tests cover negative/overflowing scratch requirements, malformed compressed structure, incompatible dimensions, noncanonical plans, changed pattern identity, wrong destinations/workspaces/providers, and unsupported direct capabilities. |
| Zero-sized inputs | Portable dense destinations, workspace requirements, compressed patterns, and symbolic plans cover legal empty shapes and rejected incompatible empties on both JVM and Scala.js. |
| Repeated reuse | Tests reuse dense destinations, workspace regions, symbolic plans, changing numeric values, and sparse-direct symbolic analyses without mutating prior owned results. |
| Determinism and diagnostics | Pure paths preserve existing deterministic behavior and typed `LinAlgError` failures. Factorization and sparse-direct capability layers retain immutable operation-specific diagnostics. |
| Backend routing | Dense destination GEMM uses the same pure, Vector, and native routing policy as allocating GEMM. Caller-workspace symmetric eigen deliberately selects the pure route because providers do not promise caller scratch reuse. |

The sparse-direct seam is intentionally JVM-only. Scala.js tests cover every
portable outcome; they do not emulate or advertise a JVM provider boundary.

## Storage-encapsulation audit

- Public pattern constructors accept ordinary arrays only as copied input; no
  caller array is adopted.
- Compressed offsets, indices, numeric arrays, dense destination storage, and
  workspace regions remain private or package-private.
- Public destinations expose checked cell operations and traversal, never a
  mutable backing array or aliased mutable matrix view.
- Builder `result()` transfers only storage allocated and exclusively owned by
  Gale, then permanently closes the builder.
- Public pivot/permutation array accessors return fresh copies. A future native
  sparse provider is contractually required to copy native index buffers before
  returning Gale values.
- JVM dense array interop remains explicitly copy-based (`fromArrayCopy`,
  `toArray`, and `toArrayRowMajor`).

No new public API exposes or adopts aliased backing storage.

## Verification receipt

The final source state passed:

```text
sbt 'coreJVM/test' 'coreJS/test'
  500 JVM tests, 491 Scala.js tests

sbt testAllFull compileAll
  core tests above, 30 JVM laws, 30 Scala.js laws,
  full optimized Scala.js linking, Scala 3.4 consumer compilation

sbt vectorBackendTest nativeBackendTest blasFfmBackendTest parityTest interopBreezeTest benchCompile
  27 Vector, 3 native-boundary, 21 FFM BLAS/LAPACK,
  35 Breeze parity, and 24 Breeze interop tests; benchmark compilation passed

sbt 'coreJVM/doc' 'coreJS/doc'
  both public API documentation builds passed
```

Focused JVM JMH and optimized Scala.js results are recorded in the
[allocation receipt](../benchmarks/results/2026-07-20-faer-nalgebra-allocation-baseline.md).
At `n=128`, reusable dense GEMM measures 5.984 B/op, fused dense linear
combination 0.259 B/op, symbolic union replay 0.682 B/op, and symbolic product
replay 1.550 B/op. Reused values-only symmetric eigen scratch reduces allocation
from 141,264.218 to 9,201.132 B/op (93.5%). Union analysis amortizes by roughly
the sixth replay in the focused JVM run, inside the tenth-replay acceptance
budget. Optimized Scala.js counters report the expected destination, workspace,
and plan reuse with no owned result on destination replay paths.

These are focused development receipts, not universal release-grade crossover
claims. Automatic dispatch remains limited to existing measured backend seams;
the new sparse-direct boundary advertises no implementation or capability.
