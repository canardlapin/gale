# Scala.js WebAssembly kernel profile (2026-07-17)

## Environment and protocol

- Apple ARM64, macOS 14.3
- Node.js 24.1.0 / V8, launched with `--experimental-wasm-exnref`
- Scala.js 1.22 experimental WebAssembly backend, ECMAScript 2022 modules
- identical Gale source and deterministic inputs under `fullLinkJS`
- five complete warmup batches followed by the median of nine timed batches;
  checksums are consumed and printed
- optimized JavaScript is Gale's normal default; Wasm requires `GALE_WASM=1`

This is a controlled wall-clock profile, not JMH. Its purpose is to determine
whether the experimental backend is ready to become a supported performance
route, not to generalize across browsers or V8 releases.

## Result

| Kernel | Optimized JS | Wasm | JS / Wasm speed |
| --- | ---: | ---: | ---: |
| dot, n=65,536 | 32,000 ns/op | 1,374,000 ns/op | 42.9x faster JS |
| GEMV, 512 x 512 | 135,000 ns/op | 4,965,000 ns/op | 36.8x faster JS |
| GEMM, 128 x 128 | 475,000 ns/op | 10,725,000 ns/op | 22.6x faster JS |

A second full Wasm execution after the v1 storage-contract audit measured
1,372,000 / 4,840,000 / 9,775,000 ns/op respectively. The 0.1–8.9% run-to-run
movement does not change the decision: optimized JavaScript remained 20.6–42.9x
faster in the repeat.

The Wasm-linked core suite passed all 373 tests, so this is not a broken-link or
wrong-result artifact. It is a clear performance rejection on the measured
runtime. Wasm remains experimental, explicit, and default-off; optimized
JavaScript remains the supported Scala.js performance route. CI now executes
both the Wasm core tests and this profile as an allow-failure visibility job.
