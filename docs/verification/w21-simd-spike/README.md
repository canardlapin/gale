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

## Results (quiet runs)

Times are ns/op, lower is better. `×scalar` is `scalar baseline / simd`, and
`×Breeze` is `breeze / simd`. Fill in one table per platform.

### aarch64 (to fill)

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

### x86-64 AVX-512 (SciNet trillium, executed 2026-10-09)

Platform: CPU `AMD EPYC 9655 (Zen 5), 2 × 96 cores, NPS4` · arch `amd64` · lanes `8`
(AVX-512, no `UseAVX` cap) · JDK `Temurin 25+36` · SHA `17813d0d761467187ba0764f1081fd30747ad18e`
· load at start `18.8` (1-min, decaying from the previous job, falling to 3.5 two minutes
into the run; the only busy process afterwards was this JVM).

How it was run: Slurm job 2524302, partition `compute`, `--exclusive --nodes=1`, node
`tri0379`. No sbt: the exported `benchmarksJVM/Jmh/fullClasspath` was copied to the
node and run as `java --add-modules=jdk.incubator.vector -cp … org.openjdk.jmh.Main`
with the recipe's options (`-jvmArgsAppend -Dgale.bench.lane=B -f 3 -wi 5 -i 10 -w 1s
-r 1s`, plus `-prof gc` for dispatch). The JVM was pinned to one CCD and its memory node with
`numactl --physcpubind=0-7 --membind=0`. The site's `JAVA_TOOL_OPTIONS=-Xmx2g` was
unset. Every fork logged `lanes=8 arch=amd64 jdk=25`. Breeze resolved
`dev.ludovic.netlib.blas.VectorBLAS` (no native BLAS is visible to the loader on these
nodes). The runner script, environment, fork logs, sidecar and both JMH JSON files are in
[`x86-avx512/`](x86-avx512/).

Correctness at this species width: the `vectorBackendTest` suites
(`VectorL1KernelsSuite`, `VectorGemmSuite`, `VectorBackendConformanceSuite`), run
with JUnitCore on the same node, gave **42/42 OK** (`x86-avx512/vtest.log`).
`dexpInto` at 8 lanes: max 1.000 ulp against `StrictMath.exp` over 2,000,003
points. 223,717 of 10,000,015 rerun results differed by 1 ulp from the first run, which
confirms the tier instability recorded above.

Values are mean ± 99.9% CI half-width, in ns/op. A `(CI overlap)` ratio does not count.

| op | n | lanes | simd | scalar baseline | Breeze | ×scalar | ×Breeze | gate |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| dot | 1024 | 8 | 61.1 ± 0.1 | 264.8 ± 0.5 | 64.2 ± 0.1 | 4.33× | 1.05× | report-only |
| dot | 4096 | 8 | 325.9 ± 4.8 | 1,031.5 ± 2.9 | 415.3 ± 0.1 | 3.17× | 1.27× | pass |
| dot | 65536 | 8 | 6,954.4 ± 417.8 | 16,620 ± 198.0 | 7,411.7 ± 4.0 | 2.39× | 1.07× | pass |
| dot | 1048576 | 8 | 141,307 ± 7,808.9 | 242,626 ± 3,076.1 | 133,957 ± 1,954.6 | 1.72× | 0.95× (CI overlap) | FAIL (×Breeze) |
| axpy | 1024 | 8 | 66.3 ± 5.7 | 285.1 ± 0.5 | 61.9 ± 11.7 | 4.30× | 0.93× (CI overlap) | report-only |
| axpy | 4096 | 8 | 361.2 ± 8.2 | 1,137.6 ± 4.5 | 368.8 ± 10.5 | 3.15× | 1.02× (CI overlap) | FAIL (×Breeze) |
| axpy | 65536 | 8 | 6,732.9 ± 321.9 | 18,131 ± 55.2 | 6,699.3 ± 276.8 | 2.69× | 1.00× (CI overlap) | FAIL (×Breeze) |
| axpy | 1048576 | 8 | 150,071 ± 7,825.6 | 299,092 ± 3,807.9 | 147,483 ± 6,785.1 | 1.99× | 0.98× (CI overlap) | FAIL (×Breeze) |
| nrm2 (vs `pureOptimisticNrm2`) | 1024 | 8 | 62.8 ± 0.1 | 267.1 ± 0.2 | 60.6 ± 0.0 | 4.25× | 0.97× (CI overlap) | report-only |
| nrm2 (vs `pureOptimisticNrm2`) | 4096 | 8 | 235.6 ± 2.1 | 1,034.2 ± 2.7 | 401.8 ± 0.2 | 4.39× | 1.71× | pass |
| nrm2 (vs `pureOptimisticNrm2`) | 65536 | 8 | 4,492.9 ± 0.8 | 16,272 ± 49.7 | 7,235.7 ± 3.6 | 3.62× | 1.61× | pass |
| nrm2 (vs `pureOptimisticNrm2`) | 1048576 | 8 | 79,444 ± 88.3 | 239,450 ± 2,787.7 | 116,686 ± 35.5 | 3.01× | 1.47× | pass |
| sum | 1024 | 8 | 30.0 ± 0.0 | 120.2 ± 0.2 | 426.3 ± 0.0 | 4.00× | 14.20× | report-only |
| sum | 4096 | 8 | 98.1 ± 16.9 | 485.6 ± 0.7 | 1,789.9 ± 0.2 | 4.95× | 18.26× | pass |
| sum | 65536 | 8 | 3,685.9 ± 18.6 | 7,747.2 ± 27.3 | 29,067 ± 6.3 | 2.10× | 7.89× | pass |
| sum | 1048576 | 8 | 64,581 ± 100.5 | 135,265 ± 4,084.7 | 465,765 ± 58.5 | 2.09× | 7.21× | pass |
| argmax | 1024 | 8 | 258.3 ± 0.4 | 295.9 ± 2.9 | 788.8 ± 0.8 | 1.15× | 3.05× | report-only |
| argmax | 4096 | 8 | 604.2 ± 0.5 | 1,315.0 ± 140.8 | 3,106.6 ± 7.1 | 2.18× | 5.14× | pass |
| argmax | 65536 | 8 | 8,102.5 ± 4.1 | 20,135 ± 1,348.4 | 29,375 ± 28.3 | 2.49× | 3.63× | pass |
| argmax | 1048576 | 8 | 165,876 ± 17,407 | 296,080 ± 114.3 | 471,985 ± 10,058 | 1.78× | 2.85× | pass |
| exp (info only) | 1024 | 8 | 308.7 ± 2.5 | 2,811.4 ± 2.5 | 2,911.7 ± 5.7 | 9.11× | 9.43× | report-only |
| exp (info only) | 4096 | 8 | 1,325.2 ± 41.8 | 11,828 ± 4.7 | 11,929 ± 29.8 | 8.93× | 9.00× | n/a |
| exp (info only) | 65536 | 8 | 19,237 ± 184.0 | 183,154 ± 2,147.0 | 198,890 ± 1,658.3 | 9.52× | 10.34× | n/a |
| exp (info only) | 1048576 | 8 | 319,885 ± 4,127.4 | 2,930,642 ± 2,045.2 | 4,577,037 ± 5,664.7 | 9.16× | 14.31× | n/a |

Per-op gate verdict at n ≥ 4K on this platform:

- dot: FAIL (2/3 sizes)
- axpy: FAIL (0/3 sizes)
- nrm2 (vs `pureOptimisticNrm2`): PASS (3/3 sizes)
- sum: PASS (3/3 sizes)
- argmax: PASS (3/3 sizes)

So nrm2, sum and argmax pass. dot fails only at n = 1M, where it is
statistically tied with VectorBLAS (0.95×, CIs overlap); it passes at 4K and 64K. axpy is
2.0–3.2× scalar, but at every size it is statistically tied with VectorBLAS
(0.98–1.02×, CIs overlap) rather than ahead of it. On this platform the ADR question (reopen
A-2b) therefore has three passing ops. The axpy criterion of "≥ 1.0× Breeze" is
not met in the strict CI sense, but neither is gale behind.

Report-only: SIMD nrm2 against the LAPACK-style `pureNrm2` (mostly algorithmic, see above).

| report-only | n | simd nrm2 | pureNrm2 (LAPACK recurrence) | ratio |
|---|---:|---:|---:|---:|
| nrm2 vs pureNrm2 | 1024 | 62.8 | 921.0 | 14.66× |
| nrm2 vs pureNrm2 | 4096 | 235.6 | 3,660.7 | 15.54× |
| nrm2 vs pureNrm2 | 65536 | 4,492.9 | 59,245 | 13.19× |
| nrm2 vs pureNrm2 | 1048576 | 79,444 | 956,800 | 12.04× |

A-R1 dispatch (`-prof gc`, 3 forks × 10 × 1 s):

| A-R1 dispatch | n | impls | lanes | direct (ns) | via `using Backend` (ns) | overhead | alloc B/op (direct / via) | within noise? |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| | 16 | 1 | 8 | 7.83 ± 0.06 | 8.21 ± 0.01 | +4.8% | 0.000 / 0.000 | no |
| | 16 | 3 | 8 | 7.79 ± 0.07 | 10.51 ± 0.06 | +34.8% | 0.000 / 0.000 | no |
| | 32 | 1 | 8 | 11.91 ± 0.07 | 12.52 ± 0.01 | +5.1% | 0.000 / 0.000 | no |
| | 32 | 3 | 8 | 11.89 ± 0.03 | 14.20 ± 0.01 | +19.5% | 0.000 / 0.000 | no |
| | 64 | 1 | 8 | 19.92 ± 0.02 | 20.48 ± 0.02 | +2.8% | 0.000 / 0.000 | no |
| | 64 | 3 | 8 | 19.92 ± 0.02 | 22.14 ± 0.02 | +11.2% | 0.000 / 0.000 | no |
| | 128 | 1 | 8 | 35.94 ± 0.02 | 36.44 ± 0.04 | +1.4% | 0.000 / 0.000 | no |
| | 128 | 3 | 8 | 35.94 ± 0.02 | 38.01 ± 0.10 | +5.8% | 0.000 / 0.000 | no |
| | 256 | 1 | 8 | 67.99 ± 0.03 | 68.50 ± 0.08 | +0.7% | 0.000 / 0.000 | no |
| | 256 | 3 | 8 | 67.94 ± 0.03 | 70.13 ± 0.16 | +3.2% | 0.000 / 0.000 | no |

A-R1 reading: the call through `(using Backend)` is **allocation-free** at every n and
site shape (`gc.alloc.rate.norm` ≈ 0 B/op). The monomorphic site (`impls=1`) costs
+0.7% to +5%, about 0.4–0.6 ns. That is small, but with 3 × 10 one-second iterations the
CIs are tight enough that it is *not* strictly within noise. The megamorphic site
(`impls=3`) costs +35% at n = 16, +20% at 32, +11% at 64, +6% at 128 and +3% at
256, which is about 2.1–2.7 ns per call. That is the same megamorphic risk the smoke
run flagged.

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
