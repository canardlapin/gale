# Reusable-architecture allocation baseline

Date: 2026-07-20

Machine: Apple ARM64

Runtime: OpenJDK 25.0.1, pure Gale backend; optimized Scala.js on Node

Harness: JMH 1.37, one fork, one 250 ms warmup, two 250 ms measurements;
Scala.js uses five warmup batches and the median of nine measured batches

This is a focused development receipt, not a release-grade performance claim.
It establishes the allocation baseline and acceptance budgets for the
destination, workspace, compressed-pattern, and symbolic-plan work. Inputs are
created in setup and excluded from timed invocations.

## Reproduction

```bash
sbt 'benchmarksJVM/Jmh/run -i 2 -wi 1 -f 1 -r 250ms -w 250ms -prof gc -p n=64 gale.bench.AllocationArchitectureJmh.*'
sbt 'benchmarksJVM/Jmh/run -i 2 -wi 1 -f 1 -r 250ms -w 250ms -prof gc -p n=128 gale.bench.AllocationArchitectureJmh.*'
sbt benchSmokeJSFull
```

## JVM allocation baseline

| Operation | n=64 time | n=64 bytes/op | n=128 time | n=128 bytes/op |
| --- | ---: | ---: | ---: | ---: |
| Dense add then scale | 3.886 us | 65,680 | 15.434 us | 262,288 |
| Dense GEMM, allocating | 26.656 us | 32,825 | 214.037 us | 131,134 |
| Symmetric eigenvalues | 186.423 us | 35,654 | 1,137.414 us | 141,268 |
| QR, fresh workspace | 54.141 us | 67,416 | 427.955 us | 338,388 |
| QR, reused workspace | 59.393 us | 66,869 | 429.894 us | 264,612 |
| Sparse add, allocating | 0.994 us | 6,040 | 1.378 us | 11,928 |
| Sparse value map, allocating | 0.328 us | 3,736 | 0.595 us | 7,320 |
| Sparse matvec, allocating | 0.188 us | 592 | 0.333 us | 1,104 |
| Sparse matvec, reused destination | 0.160 us | 0.004 | 0.304 us | 0.008 |

The dense GEMM allocation is effectively its one owned `n x n` Double result:
131,072 payload bytes at n=128 versus 131,134 measured. The add/scale route owns
two such matrices and measures 262,288 bytes, only 144 bytes above the two
payloads. A fused destination operation can therefore remove one full result,
while a GEMM destination can remove the other when the caller already owns
storage.

Reusing the existing QR workspace saves 73,776 bytes/op at n=128 (21.8%) with
neutral time in this short run. The values-only symmetric eigen call allocates
141,268 bytes at n=128 despite returning only 128 scalar eigenvalues, making its
temporary storage a higher-leverage workspace target. Sparse add and value-map
rebuild structure on every call. By contrast, the existing sparse matvec
destination contract reduces 1,104 bytes/op to measurement noise, demonstrating
that explicit reuse works across the current storage boundary.

## Scala.js construction and timing baseline

| Scenario | Median ns/op | Owned results/op | Destination reuses/op | Workspace reuses/op | Pattern reuses/op |
| --- | ---: | ---: | ---: | ---: | ---: |
| dot, 65,536 | 32,000.0 | 0 | 0 | 0 | 0 |
| GEMV, reused destination | 130,000.0 | 0 | 1 | 0 | 0 |
| GEMM, allocating | 475,000.0 | 1 | 0 | 0 | 0 |
| add then scale | 66,666.7 | 2 | 0 | 0 | 0 |
| QR, fresh workspace | 125,000.0 | 1 | 0 | 0 | 0 |
| QR, reused workspace | 125,000.0 | 1 | 0 | 1 | 0 |
| symmetric eigenvalues | 200,000.0 | 1 | 0 | 0 | 0 |
| sparse add | 12,500.0 | 1 | 0 | 0 | 0 |
| sparse value map | 4,600.0 | 1 | 0 | 0 | 0 |
| sparse matvec, allocating | 1,200.0 | 1 | 0 | 0 | 0 |
| sparse matvec, reused destination | 1,000.0 | 0 | 1 | 0 | 0 |

The baseline predates the compressed-pattern implementation, so its sparse map
row reports zero pattern reuse. Current runs report one pattern reuse/op for that
scenario. Counters describe Gale's public construction contract; they do not
inspect or estimate hidden JavaScript-engine allocation. Millisecond timer
granularity also makes these timings directional rather than release-grade.

## Ranked opportunities

1. Dense destination operations remove result-sized allocation with a direct,
   already-proven ownership model. GEMM and fused add/scale are the first gates.
2. Immutable compressed patterns and symbolic sparse plans avoid reconstructing
   row pointers and column indices in repeated map/add/factorization pipelines.
3. Spectral scratch planning has the largest observed mismatch between a small
   values-only result and temporary allocation.
4. General workspace reuse should retain the existing QR gain and extend it to
   at least one independently measured spectral or sparse-plan family.

## Acceptance budgets for dependent work

These are development gates at the measured sizes. A later release/default-route
claim still requires the repository's normal two-fork sweep.

| Path | Allocation/construction budget | Timing/crossover gate |
| --- | --- | --- |
| Dense GEMM or fused AXPBY into caller storage | JVM `gc.alloc.rate.norm <= 16 B/op` after destination setup; Scala.js reports zero owned results and one destination reuse/op | No more than 5% slower than the equivalent allocating kernel in a two-fork sweep |
| Sparse matvec into caller storage | Retain `<= 16 B/op` on JVM and zero owned results on Scala.js | No material regression from this baseline |
| Pattern-preserving sparse value map | High-level owned result allocates at most the numeric values plus 512 bytes at n=128; a planned destination route is `<= 32 B/op` | Planned/reused route no more than 5% slower than rebuilding structure |
| Symbolically planned sparse add/numeric replay | Replayed destination is `<= 32 B/op`; Scala.js reports one plan/workspace reuse and zero owned structural results per replay | Planning cost must amortize by the tenth replay at n=128 |
| QR with reusable workspace | Preserve at least a 15% allocation reduction at n=128 | No more than 5% slower than fresh workspace in a two-fork sweep |
| Values-only symmetric eigen with reusable scratch | Reduce bytes/op by at least 40% at n=128 while preserving the owned result contract | No more than 5% slower after workspace setup |

Budgets apply only after destination, workspace, or symbolic-plan construction.
Allocating high-level APIs remain the simple default surface and are measured
separately from their explicit reuse tier.

## Compressed-pattern follow-up

After `CSR.mapValues` was moved onto shared immutable pattern storage, the same
n=128 focused command measured 2,104 B/op and 0.101 us/op. The numeric payload is
2,048 bytes, so the result is 56 bytes above its one required value array and
comfortably inside the payload-plus-512-byte budget. This is a 71.3% allocation
reduction from the 7,320 B/op baseline; the short-run timing is directional only.
The optimized Scala.js smoke likewise moved from 4,600 to 600 ns/op while
retaining one owned numeric result and reusing the immutable pattern internally.

```bash
sbt 'benchmarksJVM/Jmh/run -i 2 -wi 1 -f 1 -r 250ms -w 250ms -prof gc -p n=128 gale.bench.AllocationArchitectureJmh.sparseMapAllocating'
```

## Dense-destination follow-up

The reusable GEMM and fused linear-combination paths both meet the `<= 16 B/op`
JVM budget after destination setup.

| n=128 path | Time | Allocation |
| --- | ---: | ---: |
| GEMM, allocating | 238.324 us/op | 131,134.596 B/op |
| GEMM, reused destination | 217.420 us/op | 5.984 B/op |
| Add then scale, allocating | 15.988 us/op | 262,288.443 B/op |
| Fused linear combination, reused destination | 9.524 us/op | 0.259 B/op |

The optimized Scala.js run reports zero owned results and one destination
reuse/op for both new paths. GEMM measured 450,000 ns/op for both allocating and
reused forms; fused add/scale measured 14,000 ns/op versus the allocating
pipeline's 66,666.7 ns/op. These short runs support the absence of a timing
regression, but remain development evidence rather than release-grade crossover
claims.

```bash
sbt 'benchmarksJVM/Jmh/run -i 2 -wi 1 -f 1 -r 250ms -w 250ms -prof gc -p n=128 gale.bench.AllocationArchitectureJmh.dense(AddScalePipeline|GemmAllocating|GemmReusedDestination|LinearCombinationReusedDestination)'
sbt benchSmokeJSFull
```

## Reusable-workspace follow-up

The generalized `ScratchRequirement`/`DenseWorkspace` tier now reuses both
Double and index scratch. QR reports its existing requirement through the same
checked model; dense symmetric eigenvalues and CSR canonicalization are the two
new independently measured families.

| n=128 path | Time | Allocation | Change from fresh |
| --- | ---: | ---: | ---: |
| Symmetric eigenvalues, fresh | 1,126.623 us/op | 141,264.218 B/op | — |
| Symmetric eigenvalues, reused workspace | 1,119.055 us/op | 9,201.132 B/op | -93.5% |
| CSR canonicalize, fresh workspace | 2.719 us/op | 18,136.075 B/op | — |
| CSR canonicalize, reused workspace | 2.591 us/op | 11,928.072 B/op | -34.2% |

The eigen route exceeds its 40% reduction gate and is slightly faster in this
short run. Its remaining allocation is owned result/facade storage, not the
`n*n + n` primitive scratch. Reused CSR canonicalization removes the full
Double/index sorting buffers and is likewise slightly faster; exact-sized output
structure and values remain owned by each result.

The optimized Scala.js smoke records one workspace reuse/op for both new reuse
paths while preserving one owned result/op. Eigen measured 200,000 ns/op for
both forms; CSR canonicalization measured 9,000 ns/op reused versus 10,000 ns/op
fresh. These counters validate the public ownership contract rather than hidden
engine allocation.

```bash
sbt 'benchmarksJVM/Jmh/run -i 2 -wi 1 -f 1 -r 250ms -w 250ms -prof gc -p n=128 gale.bench.AllocationArchitectureJmh.(denseSymmetricEigenvalues|denseSymmetricEigenvaluesReusedWorkspace|sparseCanonicalizeFreshWorkspace|sparseCanonicalizeReusedWorkspace)'
sbt benchSmokeJSFull
```

## Symbolic sparse-plan follow-up

`CSRUnionPlan` and `CSRProductPlan` now separate checked symbolic analysis from
numeric replay. Both retain canonical result patterns and exact input mappings;
replay writes into a caller-owned `CSRValuesDestination` and preserves explicit
numeric zeros until an explicit prune.

| n=128 path | Time | Allocation |
| --- | ---: | ---: |
| Sparse add, allocating/rebuilding structure | 1.481 us/op | 11,928.041 B/op |
| Union symbolic analysis | 6.526 us/op | 27,056.186 B/op |
| Union planned destination replay | 0.284 us/op | 0.682 B/op |
| Product symbolic analysis | 8.492 us/op | 24,072.239 B/op |
| Product planned destination replay | 0.775 us/op | 1.550 B/op |

Both replays are comfortably below the 32 B/op gate and allocate no owned
numeric or structural result. For union, one analysis plus ten replays costs
9.366 us versus 14.810 us for ten ordinary sparse additions. The measured
break-even is approximately six replays, satisfying the by-the-tenth-replay
amortization gate. Product replay has no former CSR-by-CSR public baseline; its
measurement is an allocation and liveness contract, not a crossover claim.

The optimized Scala.js run reports zero owned results, one destination reuse,
and one plan reuse per replay. At n=256, union analysis/replay measured
38,000/1,200 ns/op and product analysis/replay measured 55,000/3,500 ns/op. The
ordinary sparse add measured 12,500 ns/op, putting the directional JS union
break-even at four replays.

```bash
sbt 'benchmarksJVM/Jmh/run -i 2 -wi 1 -f 1 -r 250ms -w 250ms -prof gc -p n=128 gale.bench.AllocationArchitectureJmh.(sparseAddAllocating|sparseUnionPlanAnalysis|sparseUnionPlannedReplay|sparseProductPlanAnalysis|sparseProductPlannedReplay)'
sbt benchSmokeJSFull
```
