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
| `ddot` | Four vector fma accumulators, then a fixed-order pairwise lane combine. Reassociates: within `γ_n·Σ|x_i y_i|` |
| `daxpy` | One fma per element. **Bit-identical** |
| `dsum` | Four vector accumulators and the same fixed lane combine. Reassociates: within `γ_n·Σ|x_i|` |
| `dnrm2` | Optimistic fma sum of squares. If that is non-finite or below `1e-280`, it rescans as a SIMD max-abs pass plus a scaled pass (the `dnrmFrobenius` strategy). Keeps the contract: NaN iff an entry is NaN, else `+Inf` iff an entry is infinite, `0` for an empty input |
| `dmaxIndex` | Blocked: a SIMD `MAX` per 512-element block (NaN-propagating), strict `>` across blocks, then a scan of the winning block for the first equal element. **Exact**: first maximum, first NaN wins, `±0` tie to the first, all `-Inf` gives 0 |
| `dexpInto` | `lanewise(EXP)` plus a `Math.exp` tail. Neither bit-identical nor tier-stable (see below) |

### Reproducibility and species width

The sums do not use `reduceLanes(ADD)`: the Vector API leaves its lane order
unspecified, so its rounding could change between the interpreter, C1 and the C2
intrinsic on 4- or 8-lane hardware. Each horizontal sum extracts the lanes and adds
them in a fixed pairwise order instead. `ddot`, `dsum` and `dnrm2` are therefore
bit-reproducible within one JVM for a given species width, and a test checks this
across JIT tiers.

The result still depends on the species width (`SPECIES_PREFERRED.length()`): a
2-lane NEON, a 4-lane AVX2 and an 8-lane AVX-512 machine reassociate differently.
This is permitted by the numerical contract's clause that "vector and vendor
BLAS/LAPACK backends may reassociate operations"
(`docs/user/advanced/numerical-contract.md`). It is not covered by the
pure-path determinism promise.

The `MAX` reductions in `dnrm2` and `dmaxIndex` are order-independent, so they are
exact regardless of order.

## Correctness evidence (JDK 25.0.1, Apple M3 Max, 2 lanes)

`VectorL1KernelsSuite` has 12 tests and all pass. Lengths run from 0 to 65537,
including ragged tails and offsets.

- `daxpy` is bitwise equal to the pure kernel, both contiguous and strided.
- `ddot` and `dsum` stay within `2·(n+2)·u·Σ|terms|` of the pure kernel. They also
  stay within `(n+2)·u·Σ|terms| + u·|exact|` of an independent correctly rounded
  reference: `BigDecimal` products for the dot, and `gale.numeric.ExactSum` for the sum.
- `ddot`, `dsum` and `dnrm2` (both the optimistic and the rescan path) give
  bit-identical results across 20,000 reruns starting from the first, interpreted
  call (tier-stability test).
- `dnrm2` agrees with the pure kernel to `(n+10)·u` relative at scales of 1,
  `1e300`, `1e-300` and `1e-170`. `4096 × 1e300` gives `64e300`, and
  `4096 × 1e-300` gives `64e-300`.
- `dmaxIndex` is equal to the pure kernel on random data, coarse ties, all `-Inf`,
  alternating `±0`, one or two NaNs, tied `+Inf`, and ties or NaNs at the
  block boundaries.
- `dexpInto` was checked over 2,000,003 points in `[-745, 710]`. The measured error is **within
  1 ulp of `StrictMath.exp`**. The test asserts ≤ 2 ulp, to allow other platforms. Because fdlibm's
  `StrictMath.exp` is itself within 1 ulp, the error against the true value is
  **≤ 2 ulp**. On this JDK, `Math.exp` is 0 ulp from `StrictMath.exp`. The
  specials (`±0`, `±1`, NaN, `±Inf`, overflow and underflow edges, subnormals,
  `±MaxValue`) are exact where the result is NaN, `±Inf`, 0 or 1, and otherwise
  within 1 ulp of `StrictMath.exp`.
- **`exp` is not tier-stable.** Across five reruns of the same sweep, about 1–2.5%
  of results (the rate varies by run) differ by 1 ulp from the first, partly
  interpreted, run. In three local runs this was 115,649, 122,502 and 141,389 of
  10,000,015 results. They stay within 1 ulp of `StrictMath.exp`. The test asserts
  ≤ 2 ulp against `StrictMath.exp` (x86 SVML may sit 2 ulp from fdlibm) and prints
  the measured maximum.

### Recommendation: exclude `exp` from W2.2 default routing

SIMD `exp` breaks the determinism promise in
`docs/user/advanced/numerical-contract.md:17`, which says the pure single-threaded
implementation is deterministic for a fixed build and runtime. A selected backend
should not weaken that promise for an elementwise map. It breaks it in two ways:

- **JIT tier flips.** The same input can round differently before and after C2
  compilation.
- **Position dependence.** Elements in the scalar tail go through `Math.exp`, and
  elements in the vector body go through the `EXP` stub. The same value can
  therefore give different bits depending on its index and on `n mod lanes`.

W2.2 should not route `expInto` to this kernel by default. It could be offered later
as an explicit opt-in fast-math capability, with its accuracy bound documented
(≤ 2 ulp) and no reproducibility promise. The spike still measures `exp`, for
information only.

## Gate (from the plan, §W2.1)

The gate is evaluated **per op**:

- SIMD ≥ **1.25×** the scalar baseline at every `n ≥ 4K`, **and**
- SIMD ≥ **1.0×** Breeze lane B: VectorBLAS for dot, axpy and nrm2; `breeze.linalg`
  `sum` and `argmax`; `breeze.numerics.exp`.

A ratio counts only when the 99.9% CIs do not overlap. If every op fails, close the
ADR as "keep A-2b". `exp` is information-only (see the recommendation above).

The scalar baseline is the pure kernel with one exception. For nrm2 it is
`pureOptimisticNrm2`: the scalar single-pass optimistic algorithm
(`DoubleKernels.dnrmFrobenius` on a `1×n` block), so the ratio isolates the SIMD
gain from the algorithm change. `pureNrm2`, the LAPACK-style recurrence, is
report-only. Its large ratio is mostly algorithmic, and W1.0 will bring the
optimistic algorithm to the pure `dnrm2`. The bench data is in `[-1, 1)`, so its
sum of squares is about `n/3`, and **neither nrm2 kernel ever takes the rescan
path** in these benchmarks. The rescan path is covered by tests only.

**Scope.** The ADR decision covers only the species width and CPU that were
measured. With only the aarch64 run, the decision is scoped to aarch64 (2-lane
NEON). Covering x86 needs the AVX2 run below. Each fork logs
`[w21-simd] … lanes=… arch=… jdk=…` to stderr; copy those values into the
results table.

Requirement A-R1 for dispatch: at n = 16..256, the `(using Backend)` call through
`DenseDoubleKernel` must be allocation-free (`-prof gc`: `gc.alloc.rate.norm ≈ 0`)
and within JMH noise of the direct static `DoubleKernels.ddot` call. Check this both
with a monomorphic site (`implementations=1`) and with a site polluted by three
backends (`implementations=3`).

## Quiet-machine run, aarch64 (not yet executed)

Run this on an otherwise idle machine. Check `uptime` first: the load average should
be below 2. Use JDK 25 and lane B (`--add-modules` is inherited from
`Jmh / javaOptions`).

```sh
export JAVA_HOME=/opt/homebrew/opt/openjdk@25
OUT="$PWD/docs/verification/w21-simd-spike/aarch64"
mkdir -p "$OUT"
{ uptime; sysctl -n machdep.cpu.brand_string; git rev-parse HEAD; } > "$OUT/env.txt"
sbt -java-home "$JAVA_HOME" \
  "benchmarksJVM/breezeNetlibReset" \
  "benchmarksJVM/Jmh/run -jvmArgsAppend -Dgale.bench.lane=B -f 3 -wi 5 -i 10 -w 1s -r 1s -rf json -rff $OUT/ops-laneB.json gale.bench.SimdL1SpikeJmh" \
  "benchmarksJVM/Jmh/run -jvmArgsAppend -Dgale.bench.lane=B -f 3 -wi 5 -i 10 -w 1s -r 1s -prof gc -rf json -rff $OUT/dispatch.json gale.bench.BackendDispatchJmh" \
  2>&1 | tee "$OUT/run.log"
grep -E '\[w21-simd\]|\[breeze-netlib\]' "$OUT/run.log" > "$OUT/forks.log"
cp benchmarks/jvm/target/breeze-netlib.jsonl "$OUT/netlib.jsonl"
```

## Quiet-machine run, x86-64 AVX2 (cluster node; not yet executed)

Use an exclusive node (for example `srun --exclusive -N1 -c <all cores>`) with
JDK 25 and sbt on the `PATH`. On an AVX-512 host, `-XX:UseAVX=2` caps the preferred
species at 4 lanes (256 bit), so the result describes AVX2 hardware. It is passed
through `JAVA_TOOL_OPTIONS` so that every JMH fork inherits it. Add a second run
without the flag for an 8-lane AVX-512 data point.

```sh
export JAVA_HOME=/path/to/jdk-25            # e.g. module load openjdk/25
OUT="$PWD/docs/verification/w21-simd-spike/x86-avx2"
mkdir -p "$OUT"
export JAVA_TOOL_OPTIONS=-XX:UseAVX=2
{ uptime; lscpu; git rev-parse HEAD; } > "$OUT/env.txt"
sbt -java-home "$JAVA_HOME" \
  "benchmarksJVM/breezeNetlibReset" \
  "benchmarksJVM/Jmh/run -jvmArgsAppend -Dgale.bench.lane=B -f 3 -wi 5 -i 10 -w 1s -r 1s -rf json -rff $OUT/ops-laneB.json gale.bench.SimdL1SpikeJmh" \
  "benchmarksJVM/Jmh/run -jvmArgsAppend -Dgale.bench.lane=B -f 3 -wi 5 -i 10 -w 1s -r 1s -prof gc -rf json -rff $OUT/dispatch.json gale.bench.BackendDispatchJmh" \
  2>&1 | tee "$OUT/run.log"
grep -E '\[w21-simd\]|\[breeze-netlib\]' "$OUT/run.log" > "$OUT/forks.log"
cp benchmarks/jvm/target/breeze-netlib.jsonl "$OUT/netlib.jsonl"
grep -m1 'lanes=4' "$OUT/forks.log"     # must report 4 lanes under UseAVX=2
```

On Linux, netlib tries native JNI BLAS first. If `netlib.jsonl` shows a native class
instead of `VectorBLAS`, the Breeze column measures native BLAS. Record which class
was used. To compare against VectorBLAS, re-run on a node without a system
`libblas.so.3`, or with it hidden from the loader.

Run `vectorBackendTest` on the same node as well, so that tier stability and the
`exp` bound are checked at that species width.

## Results (quiet runs — to fill)

Times are ns/op, lower is better. `×scalar` is `scalar baseline / simd`, and
`×Breeze` is `breeze / simd`. Fill in one table per platform.

Platform: CPU `______` · arch `______` · lanes `__` · JDK `______` · SHA `______` · load `____`

| op | n | lanes | simd | scalar baseline | Breeze | ×scalar | ×Breeze | gate |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| dot | 4096 | | | | | | | |
| dot | 65536 | | | | | | | |
| dot | 1048576 | | | | | | | |
| axpy | 4096 | | | | | | | |
| axpy | 65536 | | | | | | | |
| axpy | 1048576 | | | | | | | |
| nrm2 (vs `pureOptimisticNrm2`) | 4096 | | | | | | | |
| nrm2 (vs `pureOptimisticNrm2`) | 65536 | | | | | | | |
| nrm2 (vs `pureOptimisticNrm2`) | 1048576 | | | | | | | |
| sum | 4096 | | | | | | | |
| sum | 65536 | | | | | | | |
| sum | 1048576 | | | | | | | |
| argmax | 4096 | | | | | | | |
| argmax | 65536 | | | | | | | |
| argmax | 1048576 | | | | | | | |
| exp (info only) | 4096 | | | | | | | n/a |
| exp (info only) | 65536 | | | | | | | n/a |
| exp (info only) | 1048576 | | | | | | | n/a |

(n = 1024 is also measured and is report-only, as is `pureNrm2`.)

| A-R1 dispatch | n | impls | lanes | direct (ns) | via `using Backend` (ns) | alloc B/op | within noise? |
|---|---:|---:|---:|---:|---:|---:|---|
| | 16 | 1 | | | | | |
| | 16 | 3 | | | | | |
| | 64 | 1 | | | | | |
| | 64 | 3 | | | | | |
| | 256 | 1 | | | | | |
| | 256 | 3 | | | | | |

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
