# W2.1 SIMD level-1 and reduction spike

This spike tests whether ADR A-2b ("the Vector-API backend does not seam L1/L2",
`docs/backend-architecture.md`, §A.2.1) should be reopened. It is item W2.1 of
`docs/breeze-competitiveness-plan.md`. It covers measurement only: no facade or
public API routes to these kernels yet. That wiring is W2.2, and it happens only if
the gate below passes.

## What exists

| Piece | Location |
|---|---|
| Kernels | `backend-jvm-vector/src/main/scala/gale/backend/jvm/vector/VectorL1Kernels.scala` |
| Correctness tests | `backend-jvm-vector/src/test/scala/gale/backend/jvm/vector/VectorL1KernelsSuite.scala` (`sbt vectorBackendTest`) |
| JMH | `benchmarks/jvm/src/main/scala/gale/bench/SimdL1SpikeJmh.scala`: `SimdL1SpikeJmh` (ops) and `BackendDispatchJmh` (A-R1) |

Each kernel has the strided signature of its `DoubleKernels` twin. The unit-stride
path uses `DoubleVector.SPECIES_PREFERRED` vectors plus a scalar tail. A strided
input, or a runtime narrower than two doubles, forwards to the pure kernel.

| Kernel | Semantics relative to `DoubleKernels` |
|---|---|
| `ddot` | Four vector fma accumulators. Reassociates: within `γ_n·Σ|x_i y_i|` |
| `daxpy` | One fma per element. **Bit-identical** |
| `dsum` | Four vector accumulators. Reassociates: within `γ_n·Σ|x_i|` |
| `dnrm2` | Optimistic fma sum of squares. If that is non-finite or below `1e-280`, it rescans as a SIMD max-abs pass plus a scaled pass (the `dnrmFrobenius` strategy). Keeps the contract: NaN iff an entry is NaN, else `+Inf` iff an entry is infinite, `0` for an empty input |
| `dmaxIndex` | Blocked: a SIMD `MAX` per 512-element block (NaN-propagating), strict `>` across blocks, then a scan of the winning block for the first equal element. **Exact**: first maximum, first NaN wins, `±0` tie to the first, all `-Inf` gives 0 |
| `dexpInto` | `lanewise(EXP)` plus a `Math.exp` tail. Neither bit-identical nor tier-stable (see below) |

## Correctness evidence (JDK 25.0.1, Apple M3 Max, 2 lanes)

`VectorL1KernelsSuite` has 10 tests and all pass. Lengths run from 0 to 65537,
including ragged tails and offsets.

- `daxpy` is bitwise equal to the pure kernel, both contiguous and strided.
- `ddot` and `dsum` stay within `2·(n+2)·u·Σ|terms|` of the pure kernel.
- `dnrm2` agrees to `(n+10)·u` relative at scales of 1, `1e300`, `1e-300` and
  `1e-170`. `4096 × 1e300` gives `64e300`, and `4096 × 1e-300` gives `64e-300`.
- `dmaxIndex` is equal to the pure kernel on random data, coarse ties, all `-Inf`,
  alternating `±0`, one or two NaNs, tied `+Inf`, and ties or NaNs at the
  block boundaries.
- `dexpInto` was checked over 2,000,003 points in `[-745, 710]` against
  `StrictMath.exp`. The **maximum error is 1.000 ulp**. On this JDK `Math.exp`
  is 0 ulp from `StrictMath.exp`. The specials (`±0`, `±1`, NaN, `±Inf`, overflow
  and underflow edges, subnormals, `±MaxValue`) are exact where the result is
  NaN, `±Inf`, 0 or 1, and otherwise within 1 ulp.
- **`exp` is not tier-stable.** Across five reruns of the same sweep, about 1.1–1.2%
  of results (115,649 and 122,502 of 10,000,015 in two runs) differ by 1 ulp from the first, partly
  interpreted, run. They are still within 1 ulp of `StrictMath.exp`. A W2.2 `exp`
  contract can promise an accuracy bound. It cannot promise run-to-run
  determinism within one JVM.

## Gate (from the plan, §W2.1)

The gate is evaluated **per op**:

- SIMD ≥ **1.25×** the pure scalar kernel at every `n ≥ 4K`, **and**
- SIMD ≥ **1.0×** Breeze lane B: VectorBLAS for dot, axpy and nrm2; `breeze.linalg`
  `sum` and `argmax`; `breeze.numerics.exp`.

A ratio counts only when the 99.9% CIs do not overlap. `exp` is expected to fail,
and that is an acceptable result. If every op fails, close the ADR as "keep A-2b".

Requirement A-R1 for dispatch: at n = 16..256, the `(using Backend)` call through
`DenseDoubleKernel` must be allocation-free (`-prof gc`: `gc.alloc.rate.norm ≈ 0`)
and within JMH noise of the direct static `DoubleKernels.ddot` call. Check this both
with a monomorphic site (`implementations=1`) and with a site polluted by three
backends (`implementations=3`).

## Quiet-machine run (not yet executed)

Run this on an otherwise idle machine. Check `uptime` first: the load average should
be below 2. Use JDK 25 and lane B (`--add-modules` is inherited from
`Jmh / javaOptions`).

```sh
export JAVA_HOME=/opt/homebrew/opt/openjdk@25
OUT="$PWD/docs/verification/w21-simd-spike"
sbt -java-home "$JAVA_HOME" \
  "benchmarksJVM/breezeNetlibReset" \
  "benchmarksJVM/Jmh/run -jvmArgsAppend -Dgale.bench.lane=B -f 3 -wi 5 -i 10 -w 1s -r 1s -rf json -rff $OUT/ops-laneB.json gale.bench.SimdL1SpikeJmh" \
  "benchmarksJVM/Jmh/run -jvmArgsAppend -Dgale.bench.lane=B -f 3 -wi 5 -i 10 -w 1s -r 1s -prof gc -rf json -rff $OUT/dispatch.json gale.bench.BackendDispatchJmh"
cp benchmarks/jvm/target/breeze-netlib.jsonl "$OUT/netlib.jsonl"
```

`netlib.jsonl` must show `blas=dev.ludovic.netlib.blas.VectorBLAS`. Record the
commit SHA, JDK, machine and load average with the results.

## Results (quiet run — to fill)

Times are ns/op, lower is better. Speedup is `pure / simd` and `breeze / simd`.

| op | n | simd | pure | Breeze | ×pure | ×Breeze | gate |
|---|---:|---:|---:|---:|---:|---:|---|
| dot | 4096 | | | | | | |
| dot | 65536 | | | | | | |
| dot | 1048576 | | | | | | |
| axpy | 4096 | | | | | | |
| axpy | 65536 | | | | | | |
| axpy | 1048576 | | | | | | |
| nrm2 | 4096 | | | | | | |
| nrm2 | 65536 | | | | | | |
| nrm2 | 1048576 | | | | | | |
| sum | 4096 | | | | | | |
| sum | 65536 | | | | | | |
| sum | 1048576 | | | | | | |
| argmax | 4096 | | | | | | |
| argmax | 65536 | | | | | | |
| argmax | 1048576 | | | | | | |
| exp | 4096 | | | | | | |
| exp | 65536 | | | | | | |
| exp | 1048576 | | | | | | |

(n = 1024 is also measured. It is report-only.)

| A-R1 dispatch | n | impls | direct (ns) | via `using Backend` (ns) | alloc B/op | within noise? |
|---|---:|---:|---:|---:|---:|---|
| | 16 | 1 | | | | |
| | 16 | 3 | | | | |
| | 64 | 1 | | | | |
| | 64 | 3 | | | | |
| | 256 | 1 | | | | |
| | 256 | 3 | | | | |

## Smoke run (harness check only — NOT evidence)

This run used 1 fork, 1–2 warmup iterations and 2 iterations of 300 ms. The machine
was shared, with a load average of about 27. These numbers prove that the harness
runs. They must not be used for the gate.

| op | n | simd/pure | simd/Breeze |
|---|---:|---:|---:|
| dot | 4096 / 64K / 1M | 2.4× / 1.2× / 1.3× | 3.4× / 2.2× / 2.5× |
| axpy | 4096 / 64K / 1M | 2.0× / 1.3× / 1.5× | 1.1× / 0.9× / 1.0× |
| nrm2 | 4096 / 64K / 1M | 7× / 14× / 29× | 3.8× / 4.0× / 2.2× |
| sum | 4096 / 64K / 1M | 2.1× / 1.8× / 0.7× | 6.7× / 5.9× / 2.6× |
| argmax (blocked) | 4096 / 64K / 1M | 3.7× / 3.6× / 3.6× | not re-run |
| exp | 4096 / 64K / 1M | 1.4× / 2.2× / 2.1× | 1.4× / 4.5× / 3.7× |

The smoke dispatch run gave `implementations=1` at about +0–5% over direct, and
`implementations=3` at about +16–21% over direct at n = 16..128 and +5% at n = 256. Allocation was not measured in the smoke run. A megamorphic site is a real A-R1 risk if several backends are ever live in one process.

Several things in this run are noise or need explaining:

- `nrm2`'s large ratios reflect the pure kernel's per-element LAPACK
  scale/divide recurrence. The ratio is not SIMD width.
- `sum` at 1M below 1× is unexplained in this noisy run.
- The first argmax design (lane-wise compare/blend with a loop-carried mask) was
  about 10× *slower* than scalar on 2-lane NEON. The two-pass max-then-find design
  was about 2× slower at 64K. Both were replaced by the blocked design above.
