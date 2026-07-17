# Gale Vector versus Breeze VectorBLAS — ARM64 / JDK 22 (2026-07-17)

## Scope and environment

- Apple ARM64, macOS 14.3
- Oracle OpenJDK 22, `--add-modules=jdk.incubator.vector`
- `DoubleVector.SPECIES_PREFERRED.length == 2`
- Breeze 2.1; native JNI BLAS unavailable, so dev.ludovic netlib selected its
  pure-Java VectorBLAS implementation
- identical deterministic values, separate native layouts, result allocation and
  Gale GEMM packing included in the timed public-facade methods
- average-time mode, 5 warmups, 8 measurements, 2 forks, one thread

Commands:

```sh
sbt "benchmarksJVM/Jmh/run -f 2 -wi 5 -i 8 -w 500ms -r 500ms -p n=128,256,512 .*VectorGemmJmh.*"
sbt "benchmarksJVM/Jmh/run -f 2 -wi 5 -i 8 -w 500ms -r 500ms -p n=128,256,1024 .*VectorGemvJmh.*"
```

## GEMM

Lower is better. Ratios compare operation time, using the reported JMH means.

| n | Breeze VectorBLAS | pure Gale | Vector Gale | Vector / pure speedup | Vector vs Breeze |
|---:|---:|---:|---:|---:|---:|
| 128 | 175.578 us | 206.098 us | 144.176 us | 1.43x | 1.22x faster |
| 256 | 1,228.730 us | 2,106.915 us | 1,319.746 us | 1.60x | 7.4% slower |
| 512 | 10,146.650 us | 17,655.212 us | 10,549.609 us | 1.67x | 4.0% slower |

The packed-column 3x3 dot-product kernel fixes the original two-lane regression.
The public Vector route is faster than pure Gale at every dispatched size and is
competitive with Breeze's mature VectorBLAS implementation.

## GEMV

| n | Breeze VectorBLAS | pure Gale | Vector Gale | Vector / pure speedup | Vector vs Breeze |
|---:|---:|---:|---:|---:|---:|
| 128 | 1.432 us | 3.546 us | 1.787 us | 1.98x | 24.8% slower |
| 256 | 6.821 us | 14.554 us | 8.205 us | 1.77x | 20.3% slower |
| 1024 | 114.012 us | 278.990 us | 146.992 us | 1.90x | 28.9% slower |

The four-output row kernel materially improves Gale, but Breeze retains the GEMV
lead on this runtime. This is a documented optimization opportunity, not hidden by
the broader equivalence claim.

## Dispatch decision

For runtimes with at least two preferred double lanes:

- contiguous GEMM dispatches at `rows * cols * shared >= 128^3`;
- contiguous GEMV dispatches at `rows * cols >= 128^2`;
- strided/transposed layouts fall back to the pure kernel;
- factorization remains pure.

These thresholds are evidence for this CPU/JVM family. JDK 21/22 CI establishes
compatibility and conformance, not cross-machine performance; new platform data
may justify platform-aware policy later.
