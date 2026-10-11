# hc-lstsq16 paired run on trillium (Zen 5, lanes A and B)

`7ffd558` (breeze/hc-lstsq16: fused least squares, narrow-panel QR shared by pivoted and unpivoted QR) against the Milestone 2 baseline. Slurm job 2537701, one exclusive compute node (`tri0319`, 2 × AMD EPYC 9655, Zen 5, AVX-512).

## Protocol

- **Bundles:** base is the Milestone 2 bundle compiled from `0ff6440`, which is source-identical to `5d3950e`, the branch base. New is `7ffd558`, compiled locally in a clean worktree (local branch `bench/hcls-7ffd558`) with sbt 1.11.7, JDK 25.0.1 and Scala 3.7.4. No sbt ran on the cluster.
- **Harness:** `docs/verification/w21-simd-spike/x86-avx512/ls.sh`, which is the reviewed `run.sh` lane paths (lane A `-jvmArgs`/`backend=pure`; lane B inherited `--add-modules`/`backend=pure,vector`), with `-Xms4g -Xmx4g`, a bundle-bound `gale.bench.sourceSha`, `-foe true`, fresh netlib and validation sidecars per run, and the annotation defaults (2 forks). Temurin 25+36, `JAVA_TOOL_OPTIONS` unset, `numactl --physcpubind=0-7 --membind=0`. It has two modes:
  - `ls*`: `LeastSquaresBreezeJmh` (n = 16, 64, 256), `FactorizationBreezeJmh.(gale|breeze)Qr` (n = 16, 64, 256) and `FactorizationLargeBreezeJmh.(gale|breeze)Qr` (n = 512, 1024).
  - `pq*`: `TallPivotedQrJmh`, the copy-inclusive pivoted factorization over n ∈ {512 … 10000} and p ∈ {3 … 24}. It is gale-only, since Breeze has no pivoted QR.
- **Order:** interleaved per round and lane: base then new for each mode. Full order: ls-laneA-base-r1, ls-laneA-new-r1, pq-laneA-base-r1, pq-laneA-new-r1, ls-laneB-base-r1, ls-laneB-new-r1, pq-laneB-base-r1, pq-laneB-new-r1, ls-laneA-base-r2, ls-laneA-new-r2, pq-laneA-base-r2, pq-laneA-new-r2, ls-laneB-base-r2, ls-laneB-new-r2, pq-laneB-base-r2, pq-laneB-new-r2.
- **Netlib** (blas, lapack, vectorModule): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False), ('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- **Validation:** every `ls*` run wrote a `BreezeBenchmarkValidation` sidecar bound to its bundle's source SHA: base `0ff6440`, new `7ffd558`. Every record passed, with 32 records per lane-A run and 48 per lane-B run and no failures. Each `ls*` run's scoreboard (`breeze_scoreboard.py --lane {A,B} --strict --validation <its sidecar>`) exited 0. It validated 16 result rows per lane-A run and 24 per lane-B run, with 0 unpaired rows. `TallPivotedQrJmh` is gale-only and not a Breeze pair, so it is neither preflight-covered nor scoreboard-able (it has no netlib records). Its rows are reported as new/base only.

## gale / Breeze (speed ratio, >1 means gale is faster; verdict from 99.9% CIs)

| lane | class | op | params | backend | base r1 | base r2 | new r1 | new r2 |
|---|---|---|---|---|---|---|---|---|
| A | FactorizationBreezeJmh | Qr | n=16 | pure | 1.42x ahead | 1.39x ahead | 2.08x ahead | 2.12x ahead |
| A | FactorizationBreezeJmh | Qr | n=64 | pure | 1.10x ahead | 1.09x ahead | 1.09x ahead | 1.09x ahead |
| A | FactorizationBreezeJmh | Qr | n=256 | pure | 1.48x ahead | 1.48x ahead | 1.48x ahead | 1.48x ahead |
| A | FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.42x ahead | 1.40x ahead | 1.42x ahead | 1.42x ahead |
| A | FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.54x ahead | 1.53x ahead | 1.54x ahead | 1.51x ahead |
| A | LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.87x behind | 0.88x behind | 1.27x ahead | 1.27x ahead |
| A | LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.23x ahead | 1.23x ahead | 1.28x ahead | 1.29x ahead |
| A | LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.22x ahead | 1.22x ahead | 1.22x ahead | 1.20x ahead |
| B | FactorizationBreezeJmh | Qr | n=16 | pure | 1.41x ahead | 1.43x ahead | 2.09x ahead | 2.12x ahead |
| B | FactorizationBreezeJmh | Qr | n=16 | vector | 1.38x ahead | 1.43x ahead | 2.09x ahead | 2.16x ahead |
| B | FactorizationBreezeJmh | Qr | n=64 | pure | 1.09x ahead | 1.09x ahead | 1.09x ahead | 1.08x ahead |
| B | FactorizationBreezeJmh | Qr | n=64 | vector | 1.09x ahead | 1.09x ahead | 1.09x ahead | 1.08x ahead |
| B | FactorizationBreezeJmh | Qr | n=256 | pure | 1.48x ahead | 1.48x ahead | 1.48x ahead | 1.48x ahead |
| B | FactorizationBreezeJmh | Qr | n=256 | vector | 1.49x ahead | 1.49x ahead | 1.49x ahead | 1.49x ahead |
| B | FactorizationLargeBreezeJmh | Qr | n=512 | pure | 1.40x ahead | 1.42x ahead | 1.41x ahead | 1.42x ahead |
| B | FactorizationLargeBreezeJmh | Qr | n=512 | vector | 1.43x ahead | 1.48x ahead | 1.48x ahead | 1.47x ahead |
| B | FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 1.53x ahead | 1.55x ahead | 1.55x ahead | 1.54x ahead |
| B | FactorizationLargeBreezeJmh | Qr | n=1024 | vector | 2.11x ahead | 2.30x ahead | 2.15x ahead | 2.20x ahead |
| B | LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 0.88x behind | 0.87x behind | 1.27x ahead | 1.27x ahead |
| B | LeastSquaresBreezeJmh | Lstsq | n=16 | vector | 0.87x behind | 0.88x behind | 1.27x ahead | 1.27x ahead |
| B | LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 1.22x ahead | 1.23x ahead | 1.29x ahead | 1.29x ahead |
| B | LeastSquaresBreezeJmh | Lstsq | n=64 | vector | 1.23x ahead | 1.23x ahead | 1.28x ahead | 1.29x ahead |
| B | LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 1.22x ahead | 1.21x ahead | 1.22x ahead | 1.22x ahead |
| B | LeastSquaresBreezeJmh | Lstsq | n=256 | vector | 1.27x ahead | 1.28x ahead | 1.26x ahead | 1.27x ahead |

## gale new / base (speed ratio; `*` = disjoint 99.9% CIs)

| lane | class | op | params | backend | base r1 (raw) | new r1 (raw) | r1 | r2 | call |
|---|---|---|---|---|---:|---:|---|---|---|
| A | FactorizationBreezeJmh | Qr | n=16 | pure | 5.157e+05 ± 7.5e+03 ops/s | 7.729e+05 ± 2.2e+04 ops/s | 1.50x* | 1.49x* | improved |
| A | FactorizationBreezeJmh | Qr | n=64 | pure | 1.688e+04 ± 2.1e+02 ops/s | 1.68e+04 ± 2e+02 ops/s | 1.00x | 0.99x |  |
| A | FactorizationBreezeJmh | Qr | n=256 | pure | 341.8 ± 0.45 ops/s | 341.7 ± 0.93 ops/s | 1.00x | 1.00x |  |
| A | FactorizationLargeBreezeJmh | Qr | n=512 | pure | 25.89 ± 0.14 ms/op | 25.85 ± 0.066 ms/op | 1.00x | 1.00x |  |
| A | FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 190.1 ± 3.4 ms/op | 190.3 ± 2.6 ms/op | 1.00x | 1.01x |  |
| A | LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 9.268e+04 ± 1.3e+03 ops/s | 1.348e+05 ± 2.1e+03 ops/s | 1.45x* | 1.46x* | improved |
| A | LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 2,830 ± 33 ops/s | 2,947 ± 9.6 ops/s | 1.04x* | 1.04x* | improved |
| A | LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 45.63 ± 0.27 ops/s | 45.97 ± 0.11 ops/s | 1.01x | 1.00x |  |
| B | FactorizationBreezeJmh | Qr | n=16 | pure | 5.257e+05 ± 1.4e+04 ops/s | 7.626e+05 ± 2.8e+03 ops/s | 1.45x* | 1.49x* | improved |
| B | FactorizationBreezeJmh | Qr | n=16 | vector | 5.14e+05 ± 3.1e+03 ops/s | 7.632e+05 ± 9e+03 ops/s | 1.48x* | 1.52x* | improved |
| B | FactorizationBreezeJmh | Qr | n=64 | pure | 1.688e+04 ± 1.9e+02 ops/s | 1.68e+04 ± 2e+02 ops/s | 1.00x | 1.00x |  |
| B | FactorizationBreezeJmh | Qr | n=64 | vector | 1.691e+04 ± 1.8e+02 ops/s | 1.682e+04 ± 1.8e+02 ops/s | 0.99x | 1.00x |  |
| B | FactorizationBreezeJmh | Qr | n=256 | pure | 341.1 ± 1.6 ops/s | 341.2 ± 0.95 ops/s | 1.00x | 1.00x |  |
| B | FactorizationBreezeJmh | Qr | n=256 | vector | 344.8 ± 1.6 ops/s | 345.3 ± 1.3 ops/s | 1.00x | 1.00x |  |
| B | FactorizationLargeBreezeJmh | Qr | n=512 | pure | 25.87 ± 0.1 ms/op | 25.97 ± 0.25 ms/op | 1.00x | 1.00x |  |
| B | FactorizationLargeBreezeJmh | Qr | n=512 | vector | 25.2 ± 0.19 ms/op | 24.78 ± 0.15 ms/op | 1.02x* | 0.99x |  |
| B | FactorizationLargeBreezeJmh | Qr | n=1024 | pure | 192.2 ± 5.1 ms/op | 189.9 ± 3.9 ms/op | 1.01x | 1.00x |  |
| B | FactorizationLargeBreezeJmh | Qr | n=1024 | vector | 139 ± 3.2 ms/op | 136.6 ± 0.84 ms/op | 1.02x | 0.96x |  |
| B | LeastSquaresBreezeJmh | Lstsq | n=16 | pure | 9.322e+04 ± 3.1e+02 ops/s | 1.349e+05 ± 1.5e+03 ops/s | 1.45x* | 1.45x* | improved |
| B | LeastSquaresBreezeJmh | Lstsq | n=16 | vector | 9.28e+04 ± 5.7e+02 ops/s | 1.352e+05 ± 2.1e+03 ops/s | 1.46x* | 1.44x* | improved |
| B | LeastSquaresBreezeJmh | Lstsq | n=64 | pure | 2,805 ± 33 ops/s | 2,949 ± 13 ops/s | 1.05x* | 1.05x* | improved |
| B | LeastSquaresBreezeJmh | Lstsq | n=64 | vector | 2,817 ± 23 ops/s | 2,946 ± 17 ops/s | 1.05x* | 1.05x* | improved |
| B | LeastSquaresBreezeJmh | Lstsq | n=256 | pure | 45.75 ± 0.31 ops/s | 45.79 ± 0.31 ops/s | 1.00x | 1.00x |  |
| B | LeastSquaresBreezeJmh | Lstsq | n=256 | vector | 47.58 ± 1.1 ops/s | 47.47 ± 1 ops/s | 1.00x | 0.99x |  |
| A | TallPivotedQrJmh | factorPivotedQr | n=512,p=3 | - | 12.2 ± 0.83 us/op | 10.08 ± 0.24 us/op | 1.21x* | 1.16x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=512,p=5 | - | 22.24 ± 0.61 us/op | 19.91 ± 0.48 us/op | 1.12x* | 1.11x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=512,p=6 | - | 29.02 ± 0.48 us/op | 26.46 ± 0.84 us/op | 1.10x* | 1.09x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=512,p=8 | - | 46.26 ± 1.6 us/op | 42.18 ± 0.8 us/op | 1.10x* | 1.10x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=512,p=16 | - | 96.42 ± 1.6 us/op | 82.78 ± 1.8 us/op | 1.16x* | 1.16x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=512,p=24 | - | 168.9 ± 2.6 us/op | 153.7 ± 2.5 us/op | 1.10x* | 1.10x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=1024,p=3 | - | 22.77 ± 0.58 us/op | 19.66 ± 0.62 us/op | 1.16x* | 1.16x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=1024,p=5 | - | 45.12 ± 0.94 us/op | 41.2 ± 0.83 us/op | 1.10x* | 1.09x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=1024,p=6 | - | 59.5 ± 1.2 us/op | 55.54 ± 1.1 us/op | 1.07x* | 1.07x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=1024,p=8 | - | 90.27 ± 1.5 us/op | 85.25 ± 1.8 us/op | 1.06x* | 1.06x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=1024,p=16 | - | 192 ± 3.1 us/op | 164.8 ± 3.4 us/op | 1.17x* | 1.16x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=1024,p=24 | - | 336.6 ± 4.1 us/op | 307.1 ± 4 us/op | 1.10x* | 1.10x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=2048,p=3 | - | 47.58 ± 1.2 us/op | 41.27 ± 1.1 us/op | 1.15x* | 1.14x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=2048,p=5 | - | 91.79 ± 2.2 us/op | 83.97 ± 1.9 us/op | 1.09x* | 1.09x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=2048,p=6 | - | 118.6 ± 2.3 us/op | 110.4 ± 2.5 us/op | 1.07x* | 1.07x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=2048,p=8 | - | 187 ± 3.4 us/op | 172.2 ± 4.1 us/op | 1.09x* | 1.07x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=2048,p=16 | - | 388.4 ± 7.8 us/op | 327.2 ± 7.6 us/op | 1.19x* | 1.18x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=2048,p=24 | - | 670.5 ± 11 us/op | 600.7 ± 11 us/op | 1.12x* | 1.12x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=4096,p=3 | - | 95.02 ± 2.3 us/op | 82.91 ± 2.1 us/op | 1.15x* | 1.15x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=4096,p=5 | - | 180.4 ± 3.4 us/op | 164.1 ± 3.9 us/op | 1.10x* | 1.10x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=4096,p=6 | - | 234.6 ± 4 us/op | 217.1 ± 4.7 us/op | 1.08x* | 1.08x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=4096,p=8 | - | 363.4 ± 10 us/op | 338.7 ± 7.3 us/op | 1.07x* | 1.08x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=4096,p=16 | - | 775 ± 23 us/op | 662.5 ± 23 us/op | 1.17x* | 1.17x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=4096,p=24 | - | 1,365 ± 43 us/op | 1,246 ± 42 us/op | 1.10x* | 1.09x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=10000,p=3 | - | 231.4 ± 7.2 us/op | 202 ± 6.6 us/op | 1.15x* | 1.15x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=10000,p=5 | - | 440.9 ± 14 us/op | 401.9 ± 13 us/op | 1.10x* | 1.05x |  |
| A | TallPivotedQrJmh | factorPivotedQr | n=10000,p=6 | - | 571.4 ± 12 us/op | 530.8 ± 12 us/op | 1.08x* | 1.08x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=10000,p=8 | - | 918.8 ± 22 us/op | 854.9 ± 22 us/op | 1.07x* | 1.07x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=10000,p=16 | - | 1,892 ± 27 us/op | 1,835 ± 14 us/op | 1.03x* | 1.06x* | improved |
| A | TallPivotedQrJmh | factorPivotedQr | n=10000,p=24 | - | 3,418 ± 72 us/op | 3,361 ± 73 us/op | 1.02x | 1.02x |  |
| B | TallPivotedQrJmh | factorPivotedQr | n=512,p=3 | - | 11.75 ± 0.31 us/op | 10.11 ± 0.31 us/op | 1.16x* | 1.17x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=512,p=5 | - | 22.35 ± 0.57 us/op | 20.02 ± 0.69 us/op | 1.12x* | 1.11x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=512,p=6 | - | 29.07 ± 0.44 us/op | 26.6 ± 0.75 us/op | 1.09x* | 1.10x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=512,p=8 | - | 46.25 ± 1.5 us/op | 42.34 ± 0.87 us/op | 1.09x* | 1.10x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=512,p=16 | - | 96.41 ± 1.6 us/op | 83.05 ± 1.9 us/op | 1.16x* | 1.16x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=512,p=24 | - | 168.9 ± 2.9 us/op | 153.5 ± 2.3 us/op | 1.10x* | 1.10x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=1024,p=3 | - | 22.76 ± 0.74 us/op | 19.62 ± 0.63 us/op | 1.16x* | 1.16x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=1024,p=5 | - | 45.36 ± 1 us/op | 41.41 ± 1 us/op | 1.10x* | 1.09x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=1024,p=6 | - | 59.92 ± 1.2 us/op | 55.79 ± 1.4 us/op | 1.07x* | 1.07x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=1024,p=8 | - | 90.59 ± 1.5 us/op | 85.71 ± 1.8 us/op | 1.06x* | 1.06x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=1024,p=16 | - | 192.1 ± 3.2 us/op | 165.3 ± 3.4 us/op | 1.16x* | 1.16x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=1024,p=24 | - | 336.5 ± 4.8 us/op | 306.9 ± 3.9 us/op | 1.10x* | 1.10x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=2048,p=3 | - | 47.58 ± 1.1 us/op | 41.27 ± 1.3 us/op | 1.15x* | 1.15x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=2048,p=5 | - | 92.02 ± 2.2 us/op | 84.26 ± 2.1 us/op | 1.09x* | 1.10x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=2048,p=6 | - | 119.1 ± 2.4 us/op | 110.8 ± 2.7 us/op | 1.07x* | 1.07x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=2048,p=8 | - | 182.6 ± 3.1 us/op | 171.9 ± 3.3 us/op | 1.06x* | 1.06x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=2048,p=16 | - | 388.2 ± 7.6 us/op | 328.8 ± 7.2 us/op | 1.18x* | 1.18x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=2048,p=24 | - | 670.9 ± 12 us/op | 601.3 ± 11 us/op | 1.12x* | 1.12x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=4096,p=3 | - | 95.16 ± 2.4 us/op | 82.96 ± 2.7 us/op | 1.15x* | 1.15x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=4096,p=5 | - | 181.1 ± 3.8 us/op | 164.9 ± 4.3 us/op | 1.10x* | 1.10x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=4096,p=6 | - | 235.8 ± 3.8 us/op | 217.8 ± 4.2 us/op | 1.08x* | 1.08x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=4096,p=8 | - | 360.5 ± 7.1 us/op | 339.7 ± 7.6 us/op | 1.06x* | 1.07x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=4096,p=16 | - | 777.3 ± 23 us/op | 666.5 ± 23 us/op | 1.17x* | 1.17x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=4096,p=24 | - | 1,373 ± 44 us/op | 1,265 ± 46 us/op | 1.08x* | 1.09x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=10000,p=3 | - | 231.6 ± 6.7 us/op | 201.9 ± 7.7 us/op | 1.15x* | 1.14x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=10000,p=5 | - | 441.7 ± 13 us/op | 403.1 ± 12 us/op | 1.10x* | 1.08x |  |
| B | TallPivotedQrJmh | factorPivotedQr | n=10000,p=6 | - | 571.6 ± 13 us/op | 534.7 ± 12 us/op | 1.07x* | 1.08x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=10000,p=8 | - | 915.2 ± 20 us/op | 863 ± 22 us/op | 1.06x* | 1.07x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=10000,p=16 | - | 1,892 ± 23 us/op | 1,807 ± 57 us/op | 1.05x* | 1.05x* | improved |
| B | TallPivotedQrJmh | factorPivotedQr | n=10000,p=24 | - | 3,422 ± 65 us/op | 3,314 ± 1e+02 us/op | 1.03x | 1.01x |  |

## Summary

**There is no QR regression in either lane, on either backend, at any size, pivoted or unpivoted.** No new-vs-base gale change is CI-separated and slower in both rounds, and none is CI-separated and slower even in a single round.

**Improvements, CI-separated in both rounds:**

| row | new / base | gale, base r1 → new r1 | gale/Breeze, base → new |
|---|---|---|---|
| lstsq n = 16 | 1.44–1.46× (A and B, both backends) | 92,680 ± 1,300 → 134,800 ± 2,100 ops/s (lane A) | **0.87–0.88× behind → 1.27× ahead** |
| lstsq n = 64 | 1.04–1.05× | 2,830 ± 33 → 2,947 ± 10 ops/s (lane A) | 1.22–1.23× → 1.28–1.29× |
| QR n = 16 | 1.45–1.52× | 515,700 ± 7,500 → 772,900 ± 22,000 ops/s (lane A) | 1.38–1.43× → 2.08–2.16× |
| pivoted QR (`TallPivotedQrJmh`) | 1.03–1.21× at 28 of 30 (n, p) shapes, in both lanes | e.g. n = 4096, p = 6: 1.08× | gale-only |

**Unchanged** (no CI-separated difference in both rounds): QR at n = 64, 256, 512 and 1024, lstsq n = 256, and pivoted QR at n = 10000 with p = 5 or 24. Those last two were faster with separated CIs in round 1 (1.10×) and overlapping in round 2. QR n = 1024 vector in lane B measured 1.02× and 0.96× with overlapping CIs, which is noise at ±3 ms/op.

**gale vs Breeze:** every least-squares and QR pair is now ahead in both lanes. The lane-A scoreboard goes from 7 ahead / 1 behind to 8 / 0, and lane B from 14 / 2 to 16 / 0. The two former losses were lstsq n = 16 on pure and vector.

Raw data: `2026-10-11-hc-lstsq16-trillium.json`. It holds the per-run JMH JSON, netlib and validation records, source SHA and load, keyed by mode, lane, variant and round.
