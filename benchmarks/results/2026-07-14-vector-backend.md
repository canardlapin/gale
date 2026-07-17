# Vector backend crossover — ARM64 / JDK 22 (2026-07-14)

> Historical rejected-kernel record. The packed-column GEMM and four-output GEMV
> replacements are measured in
> [the 2026-07-17 receipt](2026-07-17-breeze-jdk22-vector-enabled.md), which
> supersedes the dispatch conclusion below.

## Environment

- Apple ARM64, macOS 14.3
- Oracle OpenJDK 22
- `DoubleVector.SPECIES_PREFERRED.length == 2`
- Committed JMH protocol: average-time mode, 5 warmups, 8 measurements, 2 forks

## Finding

The original explicit-SIMD GEMM was slower than Gale's tuned pure GEMM throughout
the tested range (roughly 0.53–0.76x). Two replacement kernels were measured:

1. a four-row by preferred-species register tile retaining C across the K loop;
2. a packed-A kernel vectorized across rows.

At only two preferred double lanes, neither beat the pure blocked kernel. The
register tile reached roughly 0.60–0.73x and packed A roughly 0.45–0.60x. Gale
therefore does not dispatch explicit SIMD on runtimes narrower than four doubles:
`VectorThresholds.nativeGemmMinFlops` is `Long.MaxValue` there.

With that policy, warmed public-facade probes were performance-neutral
(approximately 0.99–1.03x pure) rather than a regression. A focused JMH run gave:

| n | pure | Vector backend | interpretation |
|---:|---:|---:|---|
| 128 | 184.664 us | 185.016 us | neutral; both use pure kernel |
| 256 | stable samples about 2.19–2.22 ms | 2.190 ms | neutral; both use pure kernel |

The committed `VectorGemmJmh` sweep covers 64, 128, 256, and 512.
Four-or-more-lane runtimes retain the register-tiled SIMD path with a conservative
`128^3` element-product threshold. CI on JDK 21 and 22 verifies compatibility and
conformance; each platform's JMH evidence remains authoritative for moving that
threshold.
