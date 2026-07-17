# FFM BLAS copy-inclusive crossover — Accelerate / JDK 22 (2026-07-14)

> The focused one-fork dispatch conclusion below is superseded by the full
> two-fork sweep added on 2026-07-17. That sweep exposed a severe n=256
> discontinuity and moved the Accelerate default to `512^3`; see the addendum.

## Environment and protocol

- Apple ARM64, macOS 14.3
- Oracle OpenJDK 22
- runtime-discovered macOS Accelerate CBLAS
- public `DMat.*` facade on both sides; native timing includes heap-to-native
  copies, native GEMM, native-to-heap copy, and result allocation
- focused JMH validation: 3 warmups, 5 measurements, 1 fork
- committed authoritative harness: 5 warmups, 8 measurements, 2 forks, sizes
  64/128/256/512/1024

## Focused result

| n | Pure Gale | FFM Accelerate | Native speedup |
|---:|---:|---:|---:|
| 128 | 182.091 us | 98.510 us | 1.85x |
| 256 | 2,131.592 us | 400.962 us | 5.32x |
| 512 | 14,605.279 us | 2,881.599 us | 5.07x |

The 512 native measurements were noisier than the smaller sizes but remained far
below pure in every sample. This focused run originally suggested `256^3`; the
authoritative follow-up below rejected that policy. Unknown/reference `libblas`
candidates stay at `Long.MaxValue` until measured.

The kernel also checks arithmetic intensity against the actual heap arrays copied.
This closes the product-only threshold hole where a `1 x k` by `k x 1` operation
could exceed the cubic threshold without enough computation to amortize copying.

## Copy-free NativeDMat control

The committed harness also keeps A/B/C in long-lived native storage. A same-fork
control at n=128 measured 74.131 us/op for copy-free `NativeDMat` GEMM, versus
147.485 us/op for the heap-copy FFM route and 181.737 us/op pure. An earlier
isolated copy-free-only sweep was erratic and implausibly slow; it was rejected
rather than averaged into the result. The full two-fork harness is the required
follow-up before assigning a separate copy-free automatic threshold.

## Authoritative two-fork follow-up (2026-07-17)

Protocol: 5 x 500 ms warmups, 8 x 500 ms measurements, 2 forks, one thread.

| n | Pure Gale | Heap-copy FFM | NativeDMat control | Heap-copy speedup |
| ---: | ---: | ---: | ---: | ---: |
| 128 | 193.061 us | 126.083 us | 96.462 us | 1.53x |
| 256 | 2,246.015 us | 9,083.703 us | 2,869.916 us | 0.25x |
| 512 | 18,691.067 us | 6,582.583 us | 7,093.006 us | 2.84x |

Accelerate is sharply non-monotone at n=256 in this run, including in the
copy-free control. A lower-bound threshold therefore cannot safely capture the
n=128 win without also routing the n=256 loss. The default is `512^3`, the first
size after the discontinuity with a complete two-fork win. OpenBLAS and MKL no
longer inherit this threshold: their automatic routes stay disabled until each
library family has its own copy-inclusive sweep.
