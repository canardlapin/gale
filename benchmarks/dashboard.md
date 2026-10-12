# Backend conformance and performance dashboard

Last refreshed: 2026-10-11 (Breeze scoreboard, Milestone 3 receipts); the
backend route tables below were last refreshed 2026-07-25. Performance numbers
are evidence for the named machine, JDK, and native library only. CI proves
compatibility and conformance; it does not establish cross-machine speed.

## gale vs Breeze 2.1: two-lane scoreboard (Milestone 3)

Receipts: source `96c95e577f6b437c3de995a081b8ba075a6054e4` (`96c95e5`);
SciNet trillium nodes `tri0463` (lane A) and `tri0501` (lane B), 2 × AMD EPYC
9655 (Zen 5, AVX-512), pinned to one CCD; Temurin 25+36; Breeze netlib
`dev.ludovic.netlib.blas.Java11BLAS` (lane A) and
`dev.ludovic.netlib.blas.VectorBLAS` (lane B), with
`dev.ludovic.netlib.lapack.F2jLAPACK` in both lanes.

This section is assembled by hand from the output of
`tools/bench/breeze_scoreboard.py`, which has no dashboard mode. Each table is
the tool's output, filtered to the rows named in its heading; nothing is
re-measured or re-scored here. Every scoreboard ran with `--strict` and
`--validation` against the receipt's `.validation.jsonl` sidecar (all records
passed and are bound to `96c95e5`). To reproduce, from `benchmarks/results/`:

```sh
P=2026-10-11-breeze-milestone3-trillium
python3 -I ../../tools/bench/breeze_scoreboard.py --lane A --strict \
  --commit 96c95e577f6b437c3de995a081b8ba075a6054e4 \
  --netlib $P-laneA.netlib.jsonl --validation $P-laneA.validation.jsonl \
  $P-laneA.json
```

Use `--lane B` and the `laneB` files for lane B, and the `-eig16` or
`-hugemethods` files for the supplements. Ratios are gale speed / Breeze speed
(> 1 means gale is faster); verdicts need non-overlapping 99.9% confidence
intervals. A `withheld` pair is not like-for-like and carries no verdict.
Speed claims hold for this machine, JDK and netlib build only.

| Lane | Comparison | Ahead | Tie | Behind | Withheld | Pairs |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| A (scalar) | gale pure vs Breeze `Java11BLAS` | 124 | 13 | 4 | 10 | 151 |
| B (SIMD), headline | gale vector + backend-insensitive vs Breeze `VectorBLAS` | 118 | 14 | 9 | 10 | 151 |
| B (SIMD), secondary | gale pure vs Breeze `VectorBLAS` | 46 | 0 | 11 | 1 | 58 |

Full receipts:
[lane A](results/2026-10-11-breeze-milestone3-trillium-laneA.md) and
[lane B](results/2026-10-11-breeze-milestone3-trillium-laneB.md), each with
its delta against Milestone 2 and its provenance.

### Still behind

| Gap | Lane A | Lane B headline | Cause or owner |
| --- | --- | --- | --- |
| `sigmoid`, n = 1K and 64K | 0.96×, 0.97× | 0.96×, 0.95× | Accepted: the per-element range test that keeps subnormal results below −36; see [Documented gaps](../docs/breeze-competitiveness-plan.md#documented-gaps) |
| `max`, n = 1M | 0.76× | 0.90× (tie) | In progress: HC-3 (`bd-01M4GY1QTAE17DV1A1990V88B8`) |
| vector `gemvT`, n = 1024 | — | 0.94–0.95× | Documented gap: follow-up `breeze/hc-gemvt-vec` abandoned; see [Documented gaps](../docs/breeze-competitiveness-plan.md#documented-gaps); [receipt](results/2026-10-11-hc-gemv-vec-trillium.md) |
| vector `gemv`/`gemvT` kernels, 1000×50 (L-BFGS design) | — | 0.91× / 0.81× (kernel harness, vs VectorBLAS) | gemvT: documented gap; gemv: open (HC-4); [receipt](results/2026-10-11-hc-gemv-vec-trillium.md) |
| `dot`, n = 64K | 0.97× | 0.97× | Paired rerun pending before kernel work (`bd-01M4PAPKQTBWWQSVC5PJ5R9CFT`) |
| sparse `dot` with a dense vector, length = 100000, nnz = 10000 | 1.00× (tie) | 0.96× | No tracker item yet |

Since Milestone 3, HC-4 (`7c460eb`, [receipt](results/2026-10-11-hc-gemv-vec-trillium.md))
closed the lane-B square vector `gemv` rows (n = 2048 1.08–1.10× ahead;
n = 256 and 1024 ties) and the gemv-bound L-BFGS logistic objective at a fixed
budget (0.95× to a 0.99–1.00× tie). The rows above are the gaps that remain.

No lane-A L2, L3, factorization or eigen pair at n ≥ 64 is behind. The
secondary lane-B rows (pure gale against SIMD Breeze) are report-only: they
measure what a caller gets without `import gale.backend.jvm.vector.given`.

### Lane A: rows not ahead

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 178469 ± 6.959e+04 | 187164 ± 3.602e+04 | ops/s | 0.95x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 29410.6 ± 3057 | 29185.5 ± 2766 | ops/s | 1.01x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7172.44 ± 784.6 | 7609.52 ± 763.4 | ops/s | 0.94x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 65460.4 ± 784.2 | 67497.1 ± 76.53 | ops/s | 0.97x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16209.9 ± 979.7 | 16950.7 ± 149.4 | ops/s | 0.96x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4202.93 ± 14.31 | 4243.21 ± 41.85 | ops/s | 0.99x | tie |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1894.91 ± 18.61 | 1766.24 ± 32.49 | us/op | 0.93x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 1238.45 ± 4.886 | 1871.51 ± 18.25 | us/op | 1.51x | withheld | not like-for-like: each library stops on its own convergence test |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 3574.05 ± 60.5 | 3083.48 ± 435.2 | ops/s | 1.16x | tie | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 340926 ± 808.1 | 341208 ± 1259 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5021.92 ± 55.77 | 5012.11 ± 56.02 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 300307 ± 930 | 301260 ± 1141 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 5.62406e+06 ± 5.047e+05 | 5.22969e+06 ± 1.505e+05 | ops/s | 1.08x | tie |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 78077.1 ± 1.686e+04 | 72557.5 ± 7569 | ops/s | 1.08x | tie |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4225.91 ± 3.417 | 5543.38 ± 332 | ops/s | 0.76x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 278961 ± 1639 | 290778 ± 1031 | ops/s | 0.96x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 4164.28 ± 35.82 | 4289.44 ± 15.44 | ops/s | 0.97x | behind |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140691 ± 0.0003203 | 1.31426 ± 0.01737 | ms/op | 9.34x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.6044 ± 0.09684 | 133.883 ± 0.1513 | ms/op | 9.17x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19677 ± 0.007083 | 12.0672 ± 0.03205 | ms/op | 10.08x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 145.603 ± 5.146 | 1274.15 ± 82.58 | ms/op | 8.75x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00497075 ± 8.995e-06 | 0.0407935 ± 5.192e-05 | ms/op | 8.21x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.543703 ± 0.002865 | 3.79305 ± 0.004615 | ms/op | 6.98x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0532194 ± 0.0001792 | 0.369588 ± 0.001226 | ms/op | 6.94x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.62393 ± 0.01572 | 37.6054 ± 0.1011 | ms/op | 6.69x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 40373.4 ± 2.278e+04 | 35138.2 ± 53.5 | ops/s | 1.15x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 222646 ± 3242 | 223760 ± 140.9 | ops/s | 1.00x | tie |  |

<details>
<summary>Lane A: all 151 pairs</summary>

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 178469 ± 6.959e+04 | 187164 ± 3.602e+04 | ops/s | 0.95x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 29410.6 ± 3057 | 29185.5 ± 2766 | ops/s | 1.01x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7172.44 ± 784.6 | 7609.52 ± 763.4 | ops/s | 0.94x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 65460.4 ± 784.2 | 67497.1 ± 76.53 | ops/s | 0.97x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16209.9 ± 979.7 | 16950.7 ± 149.4 | ops/s | 0.96x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4202.93 ± 14.31 | 4243.21 ± 41.85 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Norm | n=65536 | backend-insensitive | thrpt | 66484.5 ± 1730 | 4647.43 ± 13.88 | ops/s | 14.31x | ahead |  |
| BlasL1BreezeJmh | Norm | n=262144 | backend-insensitive | thrpt | 16828.8 ± 265.8 | 1239.4 ± 258 | ops/s | 13.58x | ahead |  |
| BlasL1BreezeJmh | Norm | n=1048576 | backend-insensitive | thrpt | 4199.86 ± 24.93 | 357.421 ± 0.1816 | ops/s | 11.75x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 98684.9 ± 477.3 | 79078.4 ± 1944 | ops/s | 1.25x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6512.89 ± 12.99 | 5097.73 ± 13.94 | ops/s | 1.28x | ahead |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1612.59 ± 19.38 | 1220.66 ± 10.01 | ops/s | 1.32x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 82874.9 ± 195.7 | 68880.8 ± 45.89 | ops/s | 1.20x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5727.07 ± 8.637 | 4219.04 ± 10.26 | ops/s | 1.36x | ahead |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1420.81 ± 5.403 | 1046.77 ± 1.367 | ops/s | 1.36x | ahead |  |
| BlasL3BreezeJmh | AtA | n=16 | pure | thrpt | 774071 ± 1.855e+04 | 397292 ± 9261 | ops/s | 1.95x | ahead |  |
| BlasL3BreezeJmh | AtA | n=64 | pure | thrpt | 15278.1 ± 78.97 | 6258.28 ± 43.47 | ops/s | 2.44x | ahead |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 215.715 ± 5.578 | 113.371 ± 0.0688 | ops/s | 1.90x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=16 | pure | thrpt | 1.78606e+06 ± 8011 | 1.50336e+06 ± 1.622e+04 | ops/s | 1.19x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=64 | pure | thrpt | 31726.5 ± 321.8 | 28332.1 ± 184.8 | ops/s | 1.12x | ahead |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 517.783 ± 5.022 | 417.692 ± 3.203 | ops/s | 1.24x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=16 | pure | thrpt | 436170 ± 1043 | 397773 ± 1460 | ops/s | 1.10x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7857.56 ± 61.99 | 6099.58 ± 36.11 | ops/s | 1.29x | ahead |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 128.777 ± 1.685 | 105.108 ± 0.9533 | ops/s | 1.23x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=16 | pure | thrpt | 1.3677e+06 ± 4070 | 944059 ± 1.844e+04 | ops/s | 1.45x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=64 | pure | thrpt | 47651.8 ± 264.4 | 31534.7 ± 211.1 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Det | n=256 | pure | thrpt | 870.107 ± 79.62 | 475.344 ± 2.519 | ops/s | 1.83x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33197.7 ± 299.5 | 12878.2 ± 61.74 | ops/s | 2.58x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1084.49 ± 0.941 | 720.489 ± 6.567 | ops/s | 1.51x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 18.9002 ± 0.05056 | 1.31477 ± 0.02058 | ops/s | 14.38x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=16 | pure | thrpt | 465548 ± 1.24e+04 | 340246 ± 8001 | ops/s | 1.37x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=64 | pure | thrpt | 17386.4 ± 108.2 | 9789.1 ± 24.46 | ops/s | 1.78x | ahead |  |
| DenseDecompositionBreezeJmh | Inv | n=256 | pure | thrpt | 362.227 ± 4.136 | 164.439 ± 0.04758 | ops/s | 2.20x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=16 | backend-insensitive | thrpt | 45538.7 ± 1173 | 15310.1 ± 213.6 | ops/s | 2.97x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=64 | backend-insensitive | thrpt | 1374.92 ± 11.86 | 1076.96 ± 12.82 | ops/s | 1.28x | ahead |  |
| DenseDecompositionBreezeJmh | Pinv | n=256 | backend-insensitive | thrpt | 27.6434 ± 0.04689 | 24.2384 ± 0.03626 | ops/s | 1.14x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=16 | backend-insensitive | thrpt | 51389.5 ± 836.1 | 18562.9 ± 179.1 | ops/s | 2.77x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=64 | backend-insensitive | thrpt | 1686.35 ± 11.75 | 1292.03 ± 11.96 | ops/s | 1.31x | ahead |  |
| DenseDecompositionBreezeJmh | Svd | n=256 | backend-insensitive | thrpt | 32.9521 ± 0.0447 | 28.1554 ± 0.07188 | ops/s | 1.17x | ahead |  |
| ElementwiseBreezeJmh | Add | n=256 | backend-insensitive | thrpt | 44260.3 ± 438.7 | 37666.6 ± 199.2 | ops/s | 1.18x | ahead |  |
| ElementwiseBreezeJmh | Add | n=1024 | backend-insensitive | thrpt | 2063.31 ± 17.33 | 1608.02 ± 22.59 | ops/s | 1.28x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=256 | backend-insensitive | thrpt | 44195.7 ± 379.1 | 38052.9 ± 416.5 | ops/s | 1.16x | ahead |  |
| ElementwiseBreezeJmh | Hadamard | n=1024 | backend-insensitive | thrpt | 2044.28 ± 8.556 | 1615.05 ± 24.92 | ops/s | 1.27x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=256 | backend-insensitive | thrpt | 43610.4 ± 419.4 | 38206.4 ± 95.82 | ops/s | 1.14x | ahead |  |
| ElementwiseBreezeJmh | Sub | n=1024 | backend-insensitive | thrpt | 2034.87 ± 9.786 | 1621.36 ± 10 | ops/s | 1.26x | ahead |  |
| FactorizationBreezeJmh | Chol | n=16 | pure | thrpt | 1.46053e+06 ± 3355 | 405552 ± 6558 | ops/s | 3.60x | ahead |  |
| FactorizationBreezeJmh | Chol | n=64 | pure | thrpt | 55876.9 ± 353.9 | 20205 ± 77.21 | ops/s | 2.77x | ahead |  |
| FactorizationBreezeJmh | Chol | n=256 | pure | thrpt | 1072.65 ± 12.28 | 653.943 ± 6.526 | ops/s | 1.64x | ahead |  |
| FactorizationBreezeJmh | Lu | n=16 | pure | thrpt | 1.40703e+06 ± 1.495e+04 | 1.04176e+06 ± 5787 | ops/s | 1.35x | ahead |  |
| FactorizationBreezeJmh | Lu | n=64 | pure | thrpt | 47588.5 ± 276.3 | 32276.7 ± 333 | ops/s | 1.47x | ahead |  |
| FactorizationBreezeJmh | Lu | n=256 | pure | thrpt | 924.305 ± 10.29 | 475.784 ± 2.433 | ops/s | 1.94x | ahead |  |
| FactorizationBreezeJmh | Qr | n=16 | pure | thrpt | 773808 ± 1.205e+04 | 368128 ± 3770 | ops/s | 2.10x | ahead |  |
| FactorizationBreezeJmh | Qr | n=64 | pure | thrpt | 16777.1 ± 171.3 | 15381.9 ± 248.7 | ops/s | 1.09x | ahead |  |
| FactorizationBreezeJmh | Qr | n=256 | pure | thrpt | 340.76 ± 0.812 | 230.765 ± 1.982 | ops/s | 1.48x | ahead |  |
| FactorizationBreezeJmh | Solve | n=16 | pure | thrpt | 1.07547e+06 ± 2.136e+04 | 787011 ± 2.143e+04 | ops/s | 1.37x | ahead |  |
| FactorizationBreezeJmh | Solve | n=64 | pure | thrpt | 41265 ± 250.3 | 29943.2 ± 202.8 | ops/s | 1.38x | ahead |  |
| FactorizationBreezeJmh | Solve | n=256 | pure | thrpt | 787.147 ± 5.251 | 469.665 ± 2.42 | ops/s | 1.68x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=512 | pure | avgt | 8.2104 ± 0.01123 | 10.4887 ± 0.00906 | ms/op | 1.28x | ahead |  |
| FactorizationLargeBreezeJmh | Chol | n=1024 | pure | avgt | 71.6164 ± 0.1905 | 74.2817 ± 0.2921 | ms/op | 1.04x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=512 | backend-insensitive | avgt | 104.319 ± 1.854 | 144.215 ± 0.2098 | ms/op | 1.38x | ahead |  |
| FactorizationLargeBreezeJmh | EigSym | n=1024 | backend-insensitive | avgt | 772.359 ± 5.799 | 1078.17 ± 3.141 | ms/op | 1.40x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=512 | pure | avgt | 77.5877 ± 0.7946 | 93.1493 ± 2.513 | ms/op | 1.20x | ahead |  |
| FactorizationLargeBreezeJmh | Lstsq | n=1024 | pure | avgt | 535.852 ± 5.149 | 736.376 ± 23.07 | ms/op | 1.37x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=512 | pure | avgt | 8.84362 ± 0.03005 | 16.4023 ± 0.01208 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Lu | n=1024 | pure | avgt | 69.6041 ± 0.1938 | 128.798 ± 0.1417 | ms/op | 1.85x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=512 | pure | avgt | 25.9314 ± 0.0754 | 36.7898 ± 0.1385 | ms/op | 1.42x | ahead |  |
| FactorizationLargeBreezeJmh | Qr | n=1024 | pure | avgt | 194.495 ± 3.625 | 289.078 ± 7.962 | ms/op | 1.49x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=512 | pure | avgt | 9.13491 ± 0.0408 | 16.5078 ± 0.02352 | ms/op | 1.81x | ahead |  |
| FactorizationLargeBreezeJmh | Solve | n=1024 | pure | avgt | 70.5868 ± 0.3249 | 129.271 ± 0.1598 | ms/op | 1.83x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 815.697 ± 7.908 | 932.429 ± 15.33 | us/op | 1.14x | ahead |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1894.91 ± 18.61 | 1766.24 ± 32.49 | us/op | 0.93x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=fixed | backend-insensitive | avgt | 45.1146 ± 0.2034 | 77.7914 ± 0.4456 | us/op | 1.72x | ahead |  |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 1238.45 ± 4.886 | 1871.51 ± 18.25 | us/op | 1.51x | withheld | not like-for-like: each library stops on its own convergence test |
| LeastSquaresBreezeJmh | Lstsq | n=16 | pure | thrpt | 134782 ± 2233 | 105840 ± 767.8 | ops/s | 1.27x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=64 | pure | thrpt | 2940.3 ± 15.1 | 2292.61 ± 4.808 | ops/s | 1.28x | ahead |  |
| LeastSquaresBreezeJmh | Lstsq | n=256 | pure | thrpt | 45.6575 ± 0.2991 | 38.0479 ± 0.8664 | ops/s | 1.20x | ahead |  |
| MatrixReductionBreezeJmh | LogSumExpRows | n=1024 | backend-insensitive | thrpt | 291.365 ± 0.07458 | 103.76 ± 0.9512 | ops/s | 2.81x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxCols | n=1024 | backend-insensitive | thrpt | 3574.05 ± 60.5 | 3083.48 ± 435.2 | ops/s | 1.16x | tie | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | MaxRows | n=1024 | backend-insensitive | thrpt | 4621.02 ± 18.26 | 187.912 ± 1.132 | ops/s | 24.59x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | NormFrobenius | n=1024 | backend-insensitive | thrpt | 4244.48 ± 19.06 | 357.162 ± 0.09516 | ops/s | 11.88x | ahead |  |
| MatrixReductionBreezeJmh | SoftmaxRows | n=1024 | backend-insensitive | thrpt | 253.436 ± 0.68 | 64.956 ± 2.263 | ops/s | 3.90x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass; gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MatrixReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7428.67 ± 4.382 | 2145.95 ± 0.7686 | ops/s | 3.46x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| MatrixReductionBreezeJmh | SumCols | n=1024 | backend-insensitive | thrpt | 3540.26 ± 70.11 | 2269.14 ± 8.674 | ops/s | 1.56x | ahead | gale row-major: columns strided for gale, contiguous for Breeze (column-major) |
| MatrixReductionBreezeJmh | SumRows | n=1024 | backend-insensitive | thrpt | 7187.75 ± 3.312 | 390.895 ± 20.04 | ops/s | 18.39x | ahead | gale row-major: rows contiguous for gale, strided for Breeze (column-major) |
| MultiRhsBreezeJmh | CholSolve | k=64, n=256 | pure | thrpt | 658.509 ± 8.711 | 319.728 ± 0.8647 | ops/s | 2.06x | ahead |  |
| MultiRhsBreezeJmh | LuSolve | k=64, n=256 | pure | thrpt | 619.664 ± 3.346 | 245.285 ± 6.426 | ops/s | 2.53x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1024 | backend-insensitive | thrpt | 3.40713e+06 ± 5.067e+04 | 1.28014e+06 ± 1630 | ops/s | 2.66x | ahead |  |
| ReductionBreezeJmh | Argmax | n=65536 | backend-insensitive | thrpt | 53793.2 ± 44.57 | 32555.8 ± 4522 | ops/s | 1.65x | ahead |  |
| ReductionBreezeJmh | Argmax | n=1048576 | backend-insensitive | thrpt | 3387.62 ± 3.073 | 2140.71 ± 3.67 | ops/s | 1.58x | ahead |  |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 340926 ± 808.1 | 341208 ± 1259 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5021.92 ± 55.77 | 5012.11 ± 56.02 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=1048576 | backend-insensitive | thrpt | 323.608 ± 0.4378 | 216.999 ± 0.3323 | ops/s | 1.49x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 300307 ± 930 | 301260 ± 1141 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | LogSumExp | n=65536 | backend-insensitive | thrpt | 4828.42 ± 13.91 | 4608.12 ± 11.94 | ops/s | 1.05x | ahead |  |
| ReductionBreezeJmh | LogSumExp | n=1048576 | backend-insensitive | thrpt | 299.096 ± 4.256 | 283.547 ± 0.8746 | ops/s | 1.05x | ahead |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 5.62406e+06 ± 5.047e+05 | 5.22969e+06 ± 1.505e+05 | ops/s | 1.08x | tie |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 78077.1 ± 1.686e+04 | 72557.5 ± 7569 | ops/s | 1.08x | tie |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4225.91 ± 3.417 | 5543.38 ± 332 | ops/s | 0.76x | behind |  |
| ReductionBreezeJmh | Mean | n=1024 | backend-insensitive | thrpt | 7.49719e+06 ± 1069 | 245562 ± 246.6 | ops/s | 30.53x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=65536 | backend-insensitive | thrpt | 118161 ± 21.64 | 3789 ± 0.7185 | ops/s | 31.19x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Mean | n=1048576 | backend-insensitive | thrpt | 7352.53 ± 31.08 | 237.264 ± 0.07126 | ops/s | 30.99x | ahead | Breeze `stats.mean` is a running mean (a division per element); gale is sum/n |
| ReductionBreezeJmh | Norm1 | n=1024 | backend-insensitive | thrpt | 4.36118e+06 ± 1297 | 107642 ± 49.98 | ops/s | 40.52x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=65536 | backend-insensitive | thrpt | 68566 ± 10.98 | 1674.33 ± 3.147 | ops/s | 40.95x | ahead |  |
| ReductionBreezeJmh | Norm1 | n=1048576 | backend-insensitive | thrpt | 4273.67 ± 1.5 | 104.361 ± 0.04779 | ops/s | 40.95x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1024 | backend-insensitive | thrpt | 2.86964e+06 ± 9796 | 2.31916e+06 ± 4993 | ops/s | 1.24x | ahead |  |
| ReductionBreezeJmh | NormInf | n=65536 | backend-insensitive | thrpt | 55597 ± 17.58 | 34140.7 ± 30.39 | ops/s | 1.63x | ahead |  |
| ReductionBreezeJmh | NormInf | n=1048576 | backend-insensitive | thrpt | 3498.56 ± 1.982 | 2129.72 ± 0.7619 | ops/s | 1.64x | ahead |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 278961 ± 1639 | 290778 ± 1031 | ops/s | 0.96x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 4164.28 ± 35.82 | 4289.44 ± 15.44 | ops/s | 0.97x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=1048576 | backend-insensitive | thrpt | 264.022 ± 0.6472 | 186.397 ± 0.2735 | ops/s | 1.42x | ahead |  |
| ReductionBreezeJmh | Softmax | n=1024 | backend-insensitive | thrpt | 269738 ± 4456 | 147852 ± 2273 | ops/s | 1.82x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=65536 | backend-insensitive | thrpt | 3980.49 ± 24.22 | 2170.46 ± 255.4 | ops/s | 1.83x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Softmax | n=1048576 | backend-insensitive | thrpt | 256.029 ± 0.7024 | 117.21 ± 0.4665 | ops/s | 2.18x | ahead | Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass |
| ReductionBreezeJmh | Sum | n=1024 | backend-insensitive | thrpt | 7.66707e+06 ± 1.741e+04 | 2.3442e+06 ± 4433 | ops/s | 3.27x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=65536 | backend-insensitive | thrpt | 122890 ± 816.5 | 34389.5 ± 30.77 | ops/s | 3.57x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| ReductionBreezeJmh | Sum | n=1048576 | backend-insensitive | thrpt | 7241.14 ± 25.36 | 2146.41 ± 0.9092 | ops/s | 3.37x | ahead | gale multi-accumulator sum vs Breeze's single-accumulator loop |
| SmallDenseBreezeJmh | Det | n=3 | pure | thrpt | 2.8239e+07 ± 7.097e+04 | 7.65217e+06 ± 1.305e+05 | ops/s | 3.69x | ahead |  |
| SmallDenseBreezeJmh | Det | n=4 | pure | thrpt | 1.88645e+07 ± 1.958e+06 | 6.63636e+06 ± 1.722e+04 | ops/s | 2.84x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=3 | pure | thrpt | 3.12498e+07 ± 2.051e+05 | 2.48481e+07 ± 7.515e+05 | ops/s | 1.26x | ahead |  |
| SmallDenseBreezeJmh | Gemm | n=4 | pure | thrpt | 3.81162e+07 ± 1.601e+06 | 1.77847e+07 ± 3.469e+05 | ops/s | 2.14x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=3 | pure | thrpt | 9.87179e+06 ± 1.177e+04 | 3.32575e+06 ± 7.207e+04 | ops/s | 2.97x | ahead |  |
| SmallDenseBreezeJmh | Inv | n=4 | pure | thrpt | 8.28318e+06 ± 2.625e+05 | 2.78387e+06 ± 2.909e+04 | ops/s | 2.98x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=3 | pure | thrpt | 1.583e+07 ± 2.535e+05 | 5.46582e+06 ± 5.176e+04 | ops/s | 2.90x | ahead |  |
| SmallDenseBreezeJmh | Solve | n=4 | pure | thrpt | 1.12167e+07 ± 8.052e+04 | 4.75671e+06 ± 1.586e+05 | ops/s | 2.36x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.230103 ± 0.0006735 | 1.32992 ± 0.008756 | ms/op | 5.78x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 24.3302 ± 0.09494 | 133.904 ± 0.412 | ms/op | 5.50x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.85066 ± 0.005471 | 12.0577 ± 0.03633 | ms/op | 6.52x | ahead |  |
| SparseMatrixBreezeJmh | CscMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 222.827 ± 0.9248 | 1276.58 ± 74.46 | ms/op | 5.73x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00578915 ± 2.347e-05 | 0.0406048 ± 0.0004838 | ms/op | 7.01x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.595366 ± 0.002092 | 3.86725 ± 0.1502 | ms/op | 6.50x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0465821 ± 0.0002329 | 0.370606 ± 0.002 | ms/op | 7.96x | ahead |  |
| SparseMatrixBreezeJmh | CscMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 6.1913 ± 0.001807 | 37.5547 ± 0.02451 | ms/op | 6.07x | ahead |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140691 ± 0.0003203 | 1.31426 ± 0.01737 | ms/op | 9.34x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.6044 ± 0.09684 | 133.883 ± 0.1513 | ms/op | 9.17x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19677 ± 0.007083 | 12.0672 ± 0.03205 | ms/op | 10.08x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 145.603 ± 5.146 | 1274.15 ± 82.58 | ms/op | 8.75x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00497075 ± 8.995e-06 | 0.0407935 ± 5.192e-05 | ms/op | 8.21x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.543703 ± 0.002865 | 3.79305 ± 0.004615 | ms/op | 6.98x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0532194 ± 0.0001792 | 0.369588 ± 0.001226 | ms/op | 6.94x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.62393 ± 0.01572 | 37.6054 ± 0.1011 | ms/op | 6.69x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Add | length=100000, nnz=1000 | backend-insensitive | thrpt | 211949 ± 1532 | 191978 ± 9497 | ops/s | 1.10x | ahead |  |
| SparseVectorBreezeJmh | Add | length=100000, nnz=10000 | backend-insensitive | thrpt | 17717.5 ± 51.76 | 15977.3 ± 275.6 | ops/s | 1.11x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.0777e+06 ± 2.306e+04 | 314889 ± 823.1 | ops/s | 6.60x | ahead |  |
| SparseVectorBreezeJmh | Axpy | length=100000, nnz=10000 | backend-insensitive | thrpt | 207477 ± 6220 | 32537.3 ± 1371 | ops/s | 6.38x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=1000 | backend-insensitive | thrpt | 1.13297e+06 ± 6.088e+04 | 324106 ± 557 | ops/s | 3.50x | ahead |  |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 40373.4 ± 2.278e+04 | 35138.2 ± 53.5 | ops/s | 1.15x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | backend-insensitive | thrpt | 2.26188e+06 ± 3923 | 2.25174e+06 ± 804 | ops/s | 1.00x | ahead |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 222646 ± 3242 | 223760 ± 140.9 | ops/s | 1.00x | tie |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97212.8 ± 835.8 | 34065.2 ± 517.5 | ops/s | 2.85x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=64 | backend-insensitive | thrpt | 3706.14 ± 16.28 | 1766.47 ± 8.402 | ops/s | 2.10x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=128 | backend-insensitive | thrpt | 584.76 ± 16.31 | 314.132 ± 2.936 | ops/s | 1.86x | ahead |  |

</details>

### Lane B headline: rows not ahead (gale vector and backend-insensitive)

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL1BreezeJmh | Axpy | n=65536 | backend-insensitive | thrpt | 186230 ± 3.68e+04 | 188835 ± 3.623e+04 | ops/s | 0.99x | tie |  |
| BlasL1BreezeJmh | Axpy | n=262144 | backend-insensitive | thrpt | 29574.1 ± 2365 | 29616.6 ± 3327 | ops/s | 1.00x | tie |  |
| BlasL1BreezeJmh | Axpy | n=1048576 | backend-insensitive | thrpt | 7688.42 ± 974.2 | 7824.41 ± 691.1 | ops/s | 0.98x | tie |  |
| BlasL1BreezeJmh | Dot | n=65536 | backend-insensitive | thrpt | 65853.8 ± 647.3 | 68035.1 ± 560.3 | ops/s | 0.97x | behind |  |
| BlasL1BreezeJmh | Dot | n=262144 | backend-insensitive | thrpt | 16949.7 ± 111.8 | 16835.2 ± 455.5 | ops/s | 1.01x | tie |  |
| BlasL1BreezeJmh | Dot | n=1048576 | backend-insensitive | thrpt | 4199.25 ± 45.13 | 4238.06 ± 42.41 | ops/s | 0.99x | tie |  |
| BlasL2BreezeJmh | Gemv | n=256 | vector | thrpt | 225633 ± 4394 | 234667 ± 946.5 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | vector | thrpt | 12607.6 ± 41 | 13287.1 ± 43.12 | ops/s | 0.95x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | vector | thrpt | 2724.76 ± 14.22 | 2843.01 ± 22.01 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | vector | thrpt | 253835 ± 1.673e+04 | 258440 ± 1.802e+04 | ops/s | 0.98x | tie |  |
| BlasL2BreezeJmh | GemvT | n=1024 | vector | thrpt | 12645.8 ± 35.44 | 13224.9 ± 126.6 | ops/s | 0.96x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | vector | thrpt | 2956.65 ± 14.89 | 2973.64 ± 9.758 | ops/s | 0.99x | tie |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | vector | avgt | 579.51 ± 7.848 | 556.114 ± 2.054 | us/op | 0.96x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | vector | avgt | 1358.84 ± 30.91 | 1047.9 ± 5.737 | us/op | 0.77x | withheld | not like-for-like: each library stops on its own convergence test |
| LbfgsBreezeJmh | Rosenbrock | budget=tolerance | backend-insensitive | avgt | 1237.37 ± 8.69 | 1855.32 ± 7.69 | us/op | 1.50x | withheld | not like-for-like: each library stops on its own convergence test |
| ReductionBreezeJmh | Exp | n=1024 | backend-insensitive | thrpt | 342333 ± 1594 | 341468 ± 1264 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | Exp | n=65536 | backend-insensitive | thrpt | 5061.57 ± 118.2 | 5041.31 ± 102.4 | ops/s | 1.00x | tie |  |
| ReductionBreezeJmh | LogSumExp | n=1024 | backend-insensitive | thrpt | 298826 ± 3593 | 302869 ± 656.9 | ops/s | 0.99x | tie |  |
| ReductionBreezeJmh | Max | n=1024 | backend-insensitive | thrpt | 5.11082e+06 ± 1.766e+04 | 5.18955e+06 ± 1.34e+05 | ops/s | 0.98x | tie |  |
| ReductionBreezeJmh | Max | n=65536 | backend-insensitive | thrpt | 89352.6 ± 56.75 | 87719.5 ± 9906 | ops/s | 1.02x | tie |  |
| ReductionBreezeJmh | Max | n=1048576 | backend-insensitive | thrpt | 4921.85 ± 1082 | 5457.77 ± 531 | ops/s | 0.90x | tie |  |
| ReductionBreezeJmh | Sigmoid | n=1024 | backend-insensitive | thrpt | 280270 ± 666.6 | 291829 ± 2122 | ops/s | 0.96x | behind |  |
| ReductionBreezeJmh | Sigmoid | n=65536 | backend-insensitive | thrpt | 4097.21 ± 58.95 | 4303.39 ± 18.45 | ops/s | 0.95x | behind |  |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=1000 | backend-insensitive | avgt | 0.140046 ± 0.0003984 | 1.27823 ± 0.008975 | ms/op | 9.13x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.01, n=10000 | backend-insensitive | avgt | 14.5458 ± 0.09704 | 121.97 ± 0.4268 | ms/op | 8.39x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=1000 | backend-insensitive | avgt | 1.19333 ± 0.007258 | 12.2149 ± 0.007124 | ms/op | 10.24x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatmul | density=0.1, n=10000 | backend-insensitive | avgt | 145.81 ± 2.87 | 1228.01 ± 50.63 | ms/op | 8.42x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=1000 | backend-insensitive | avgt | 0.00494833 ± 3.096e-06 | 0.0403577 ± 0.000279 | ms/op | 8.16x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.01, n=10000 | backend-insensitive | avgt | 0.543913 ± 0.001658 | 3.80498 ± 0.02886 | ms/op | 7.00x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=1000 | backend-insensitive | avgt | 0.0529125 ± 0.000249 | 0.370746 ± 0.001161 | ms/op | 7.01x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseMatrixBreezeJmh | CsrMatvec | density=0.1, n=10000 | backend-insensitive | avgt | 5.64577 ± 0.002769 | 37.7333 ± 0.06931 | ms/op | 6.68x | withheld | not like-for-like: Breeze has no CSR, its twin runs the CSC product |
| SparseVectorBreezeJmh | Dot | length=100000, nnz=10000 | backend-insensitive | thrpt | 41171.1 ± 2.501e+04 | 35306.5 ± 25.54 | ops/s | 1.17x | tie |  |
| SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | backend-insensitive | thrpt | 214957 ± 1693 | 224836 ± 38.42 | ops/s | 0.96x | behind |  |

### Lane B secondary: rows not ahead (gale pure)

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| BlasL2BreezeJmh | Gemv | n=256 | pure | thrpt | 98494.1 ± 476.7 | 234667 ± 946.5 | ops/s | 0.42x | behind |  |
| BlasL2BreezeJmh | Gemv | n=1024 | pure | thrpt | 6518.39 ± 12.64 | 13287.1 ± 43.12 | ops/s | 0.49x | behind |  |
| BlasL2BreezeJmh | Gemv | n=2048 | pure | thrpt | 1609.73 ± 15.19 | 2843.01 ± 22.01 | ops/s | 0.57x | behind |  |
| BlasL2BreezeJmh | GemvT | n=256 | pure | thrpt | 83452.6 ± 502.9 | 258440 ± 1.802e+04 | ops/s | 0.32x | behind |  |
| BlasL2BreezeJmh | GemvT | n=1024 | pure | thrpt | 5790.73 ± 32.04 | 13224.9 ± 126.6 | ops/s | 0.44x | behind |  |
| BlasL2BreezeJmh | GemvT | n=2048 | pure | thrpt | 1425.41 ± 10.74 | 2973.64 ± 9.758 | ops/s | 0.48x | behind |  |
| BlasL3BreezeJmh | AtA | n=256 | pure | thrpt | 214.672 ± 5.491 | 356.925 ± 9.868 | ops/s | 0.60x | behind |  |
| BlasL3BreezeJmh | Gemm | n=256 | pure | thrpt | 518.193 ± 3.501 | 1143.61 ± 20.08 | ops/s | 0.45x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=64 | pure | thrpt | 7893.88 ± 54.32 | 12469.2 ± 337.9 | ops/s | 0.63x | behind |  |
| BlasL3BreezeJmh | GemmTall | n=256 | pure | thrpt | 129.443 ± 1.533 | 294.196 ± 8.022 | ops/s | 0.44x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=fixed | pure | avgt | 814.5 ± 10.11 | 556.114 ± 2.054 | us/op | 0.68x | behind |  |
| LbfgsBreezeJmh | Logistic | budget=tolerance | pure | avgt | 1900.06 ± 16.79 | 1047.9 ± 5.737 | us/op | 0.55x | withheld | not like-for-like: each library stops on its own convergence test |

### JIT-fair nonsymmetric eig (lane A, `-XX:-DontCompileHugeMethods`)

Breeze's F2J `dgeev` has methods above HotSpot's huge-method limit that stay
partly interpreted by default, which inflates gale's lane-A eig margin
(2.58× at n = 16 and 14.38× at n = 256 in the sweep). With the flag passed to
every fork on both sides, gale stays ahead. These are the ratios to quote for
nonsymmetric eig.

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33186.5 ± 297.1 | 21857.1 ± 141.9 | ops/s | 1.52x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=64 | backend-insensitive | thrpt | 1081.71 ± 0.8313 | 1023.36 ± 8.694 | ops/s | 1.06x | ahead |  |
| DenseDecompositionBreezeJmh | Eig | n=256 | backend-insensitive | thrpt | 18.8832 ± 0.04402 | 14.2852 ± 0.7325 | ops/s | 1.32x | ahead |  |

### Controlled n = 16 eigen recheck (4 forks × 10 × 10 iterations)

Lane A:

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 34172.8 ± 519.3 | 12913.1 ± 32.47 | ops/s | 2.65x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 97869.5 ± 243.1 | 34417.5 ± 79.68 | ops/s | 2.84x | ahead |  |

Lane B:

| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |
|---|---|---|---|---|---:|---:|---|---:|---|---|
| DenseDecompositionBreezeJmh | Eig | n=16 | backend-insensitive | thrpt | 33388.1 ± 98.88 | 12974.5 ± 32.78 | ops/s | 2.57x | ahead |  |
| SymEigenBreezeJmh | EigSym | n=16 | backend-insensitive | thrpt | 98066.4 ± 64.86 | 34707.5 ± 32.65 | ops/s | 2.83x | ahead |  |

## Conformance and availability

| Route | Runtime | Correctness gate | CI status | Default policy |
| --- | --- | --- | --- | --- |
| Pure Gale | JVM 21/22 | core + laws | required | active without an import |
| Pure Gale | Scala.js / Node 22 | core + laws | required | active without an import |
| Vector API | JVM 21/22 | shared backend conformance + backend suites | required on both JDKs | opt-in import; measured coarse ops only |
| FFM BLAS/LAPACK | JVM 22 / OpenBLAS | shared backend conformance + native storage + loader/LAPACK suites | required | opt-in import; family-specific thresholds |
| Scala.js Wasm | Node 25+ / Wasm 3.0 linker | 373 core tests + kernel profile | allow-failure | explicit `GALE_WASM=1` only |

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
| Vector, Apple ARM64, JDK 22 | contiguous GEMV | `128^2` work | 1.77–1.98x over pure Gale (measured before HC-4 `7c460eb`; not re-measured since, see the [ARM directional check](results/2026-10-11-hc-gemv-vec-trillium.md#arm-directional-check)) | enabled |
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

- [Breeze two-lane scoreboard, Milestone 3, lane A](results/2026-10-11-breeze-milestone3-trillium-laneA.md)
- [Breeze two-lane scoreboard, Milestone 3, lane B](results/2026-10-11-breeze-milestone3-trillium-laneB.md)
- [Sigmoid path attribution (documented gap)](results/2026-10-11-hc-sigmoid2-path-trillium.md)
- [Matrix-free generalized LOBPCG development baseline](results/2026-07-25-generalized-lobpcg-baseline.md)
- [Reusable-architecture allocation baseline](results/2026-07-20-faer-nalgebra-allocation-baseline.md)
- [Scalafim migration workloads](results/2026-07-19-scalafim-migration.md)
- [Vector versus Breeze, JDK 22](results/2026-07-17-breeze-jdk22-vector-enabled.md)
- [FFM GEMM crossover](results/2026-07-14-ffm-blas-crossover.md)
- [FFM GEMV crossover](results/2026-07-17-ffm-gemv-crossover.md)
- [FFM LAPACK and solver scenarios](results/2026-07-17-ffm-lapack-crossover.md)
- [Scala.js Wasm profile](results/2026-07-17-wasm-profile.md)
- [Breeze migration guide](../docs/user/guides/breeze-equivalence.md)

## Refresh commands

```bash
sbt testAll parityTest interopBreezeTest vectorBackendTest benchCompile
sbt nativeBackendTest blasFfmBackendTest benchFfmCompile  # JDK 22+
sbt "benchmarksJVM/Jmh/run .*Vector.*Jmh.*"
sbt "benchmarksJVM/Jmh/run -prof gc gale.bench.GeneralizedLobpcgJmh.solve"
sbt "benchmarksJVM/Jmh/runMain gale.bench.GeneralizedLobpcgWorkReceipt"
sbt "benchmarksFfm/Jmh/run .*FfmGemvJmh.*"
sbt "benchmarksFfm/Jmh/run .*FfmLapackJmh.*"
sbt "benchmarksJVM/Jmh/run -prof gc gale.bench.ScalafimMigrationJmh.*"
GALE_WASM=1 sbt coreJS/test benchSmokeJSFull
```
