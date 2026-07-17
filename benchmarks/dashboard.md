# Backend conformance and performance dashboard

Last refreshed: 2026-07-17. Performance numbers are evidence for the named
machine, JDK, and native library only. CI proves compatibility and conformance;
it does not establish cross-machine speed.

## Conformance and availability

| Route | Runtime | Correctness gate | CI status | Default policy |
| --- | --- | --- | --- | --- |
| Pure Gale | JVM 21/22 | core + laws | required | active without an import |
| Pure Gale | Scala.js / Node 22 | core + laws | required | active without an import |
| Vector API | JVM 21/22 | shared backend conformance + backend suites | required on both JDKs | opt-in import; measured coarse ops only |
| FFM BLAS/LAPACK | JVM 22 / OpenBLAS | shared backend conformance + native storage + loader/LAPACK suites | required | opt-in import; family-specific thresholds |
| Scala.js Wasm | Node 22+ / experimental linker | 373 core tests + kernel profile | allow-failure | explicit `GALE_WASM=1` only |

The FFM capability is conditional. `NativeBlas` requires a conforming CBLAS
candidate; `NativeLapack` additionally requires all bound factorization and
symmetric-eigen symbols. The shared conformance suite uses independent
`BigDecimal` oracles for dense kernels and typed reconstruction/error laws for
factorizations. Breeze parity is a separate differential gate and is not counted
as the backend's only oracle.

## Measured automatic routes

| Backend and platform | Operation | Shipped threshold | Evidence at routed sizes | Decision |
| --- | --- | ---: | --- | --- |
| Vector, Apple ARM64, JDK 22 | contiguous GEMM | `128^3` work | 1.43–1.67x over pure Gale | enabled |
| Vector, Apple ARM64, JDK 22 | contiguous GEMV | `128^2` work | 1.77–1.98x over pure Gale | enabled |
| FFM Accelerate, Apple ARM64, JDK 22 | square GEMM | `512^3` work | 2.84x at n=512 | enabled above discontinuity |
| FFM Accelerate, Apple ARM64, JDK 22 | LU / solve | n=128 | 1.53–1.60x solve at n=128–256 | enabled |

## Available but deliberately default-off

| Route | Evidence | Reason |
| --- | --- | --- |
| FFM Accelerate GEMV | 0.30–0.85x at n=64–2048 | heap/native copies never amortized |
| FFM Accelerate Cholesky | 0.28–3.88x at n=64–512 | non-monotone; one minimum cannot represent it |
| FFM Accelerate QR / least squares | 0.13–1.34x | severe regressions at n>=128 |
| FFM symmetric eigen | 1.18–2.77x | real raw provider; no public backend-routed eigen facade yet |
| OpenBLAS/MKL automatic routes | not measured here | loadable and explicit, but no borrowed Accelerate policy |
| FFM/reference BLAS automatic routes | not measured | correctness does not imply performance |
| Scala.js Wasm on Node 24 | 0.02–0.04x optimized JS | 23–43x slower; experimental/default-off |

The declined-native control is neutral after adequate warmup: n=32 LU measured
3.952 us/op with the FFM backend selected versus 3.905 pure, and QR measured
20.865 versus 20.944 us/op. Selecting an accelerator therefore does not impose a
material warmed penalty below its threshold.

## Evidence index

- [Vector versus Breeze, JDK 22](results/2026-07-17-breeze-jdk22-vector-enabled.md)
- [FFM GEMM crossover](results/2026-07-14-ffm-blas-crossover.md)
- [FFM GEMV crossover](results/2026-07-17-ffm-gemv-crossover.md)
- [FFM LAPACK and solver scenarios](results/2026-07-17-ffm-lapack-crossover.md)
- [Scala.js Wasm profile](results/2026-07-17-wasm-profile.md)
- [Breeze equivalence guide](../docs/breeze-equivalence.md)

## Refresh commands

```bash
sbt testAll parityTest interopBreezeTest vectorBackendTest benchCompile
sbt nativeBackendTest blasFfmBackendTest benchFfmCompile  # JDK 22+
sbt "benchmarksJVM/Jmh/run .*Vector.*Jmh.*"
sbt "benchmarksFfm/Jmh/run .*FfmGemvJmh.*"
sbt "benchmarksFfm/Jmh/run .*FfmLapackJmh.*"
GALE_WASM=1 sbt coreJS/test benchSmokeJSFull
```
