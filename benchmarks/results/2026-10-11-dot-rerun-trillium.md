# dot 64K / DotDense rerun on trillium (Zen 5, lanes A and B)

Rerun item `bd-01M4PAPKQTBWWQSVC5PJ5R9CFT`: are the Milestone 3 `dot` n = 64K and `DotDense` nnz = 10K losses kernel gaps or noise? Source: main `88dde85`. Slurm job 2540116, one exclusive compute node (`tri1081`, 2 × AMD EPYC 9655, Zen 5, AVX-512).

## Protocol

- **Build:** `88dde85`, compiled locally in a clean bench worktree (branch `bench/dot-88dde85`) with `benchmarksJVM/clean` first. sbt 1.11.7, JDK 25.0.1, Scala 3.7.4. No sbt ran on the cluster.
- **Runner:** `docs/verification/w21-simd-spike/x86-avx512/dot.sh`, which uses the Milestone 3 `run.sh` lane flags. Lane A: `-jvmArgs "-Xms4g -Xmx4g …lane=A…"` (no Vector module) with `-p backend=pure`. Lane B: inherited `--add-modules`, `-jvmArgsAppend "-Xms4g -Xmx4g …lane=B…"` and `-p backend=pure,vector`. Both are bundle SHA-bound, use `-foe true`, Temurin 25+36, `JAVA_TOOL_OPTIONS` unset and `numactl --physcpubind=0-7 --membind=0`.
- **Benchmarks:** `BlasL1BreezeJmh.(galeDot|breezeDot)$` at n = 65536, 262144 and 1048576. These are the class's only sizes; there is no n = 1K parameter. Also `SparseVectorBreezeJmh.(galeDotDense|breezeDotDense)$` at length 100000 with nnz = 1000 and 10000. All runs use `-f 8 -wi 5 -i 10` with the annotation iteration times.
- **Order:** laneA-r1, laneB-r1, laneA-r2, laneB-r2.
- **Netlib** (blas, lapack, vectorModule): [('dev.ludovic.netlib.blas.Java11BLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', False), ('dev.ludovic.netlib.blas.VectorBLAS', 'dev.ludovic.netlib.lapack.F2jLAPACK', True)].
- **Validation:** neither `BlasL1BreezeJmh` nor `SparseVectorBreezeJmh` is preflight-covered (they are not in `PREFLIGHT_CLASSES`), and no validation records were written. Each run's scoreboard (`breeze_scoreboard.py --strict --validation <empty>`) exited 0 with 0 unpaired rows. Its counts per run are: `laneA-r1`: 2 ahead, 1 tie, 2 behind, 0 n/a, 0 withheld; `laneA-r2`: 1 ahead, 2 tie, 2 behind, 0 n/a, 0 withheld; `laneB-r1`: Headline, gale vector vs Breeze VectorBLAS: 1 ahead, 2 tie, 2 behind, 0 n/a, 0 withheld; `laneB-r2`: Headline, gale vector vs Breeze VectorBLAS: 1 ahead, 2 tie, 2 behind, 0 n/a, 0 withheld.
- **Bimodality rule:** a series is called bimodal when the largest gap between consecutive sorted per-fork means exceeds 4% of the median.

## gale / Breeze (ops/s; ratio is the speed ratio with conservative CI bounds, `*` = disjoint 99.9% CIs)

| lane | round | class | op | params | gale | Breeze | gale/Breeze |
|---|---:|---|---|---|---:|---:|---|
| A | 1 | BlasL1BreezeJmh | Dot | n=65536 | 65,659 ± 362 | 67,887 ± 79 | 0.967× [0.961, 0.974]* behind |
| A | 1 | BlasL1BreezeJmh | Dot | n=262144 | 16,717 ± 84 | 16,380 ± 438 | 1.021× [0.989, 1.054] tie |
| A | 1 | BlasL1BreezeJmh | Dot | n=1048576 | 4,133 ± 16 | 4,238 ± 21 | 0.975× [0.967, 0.984]* behind |
| A | 1 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | 2,272,802 ± 1,735 | 2,258,289 ± 2,061 | 1.006× [1.005, 1.008]* ahead |
| A | 1 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | 224,724 ± 105 | 223,449 ± 644 | 1.006× [1.002, 1.009]* ahead |
| A | 2 | BlasL1BreezeJmh | Dot | n=65536 | 65,077 ± 954 | 67,861 ± 81 | 0.959× [0.944, 0.974]* behind |
| A | 2 | BlasL1BreezeJmh | Dot | n=262144 | 16,750 ± 112 | 17,008 ± 93 | 0.985× [0.973, 0.997]* behind |
| A | 2 | BlasL1BreezeJmh | Dot | n=1048576 | 4,169 ± 22 | 4,198 ± 15 | 0.993× [0.984, 1.002] tie |
| A | 2 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | 2,274,400 ± 1,953 | 2,260,228 ± 1,810 | 1.006× [1.005, 1.008]* ahead |
| A | 2 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | 224,808 ± 235 | 224,412 ± 375 | 1.002× [0.999, 1.004] tie |
| B | 1 | BlasL1BreezeJmh | Dot | n=65536 | 65,341 ± 402 | 67,745 ± 27 | 0.965× [0.958, 0.971]* behind |
| B | 1 | BlasL1BreezeJmh | Dot | n=262144 | 16,791 ± 82 | 17,041 ± 82 | 0.985× [0.976, 0.995]* behind |
| B | 1 | BlasL1BreezeJmh | Dot | n=1048576 | 4,190 ± 23 | 4,210 ± 18 | 0.995× [0.986, 1.005] tie |
| B | 1 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | 2,277,065 ± 1,101 | 2,259,463 ± 2,152 | 1.008× [1.006, 1.009]* ahead |
| B | 1 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | 224,705 ± 79 | 224,217 ± 430 | 1.002× [1.000, 1.004] tie |
| B | 2 | BlasL1BreezeJmh | Dot | n=65536 | 65,404 ± 479 | 67,790 ± 18 | 0.965× [0.957, 0.972]* behind |
| B | 2 | BlasL1BreezeJmh | Dot | n=262144 | 16,701 ± 83 | 16,842 ± 190 | 0.992× [0.976, 1.008] tie |
| B | 2 | BlasL1BreezeJmh | Dot | n=1048576 | 4,123 ± 44 | 4,190 ± 20 | 0.984× [0.969, 0.999]* behind |
| B | 2 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=1000 | 2,273,068 ± 1,053 | 2,256,845 ± 3,424 | 1.007× [1.005, 1.009]* ahead |
| B | 2 | SparseVectorBreezeJmh | DotDense | length=100000, nnz=10000 | 224,678 ± 101 | 224,428 ± 205 | 1.001× [1.000, 1.002] tie |

## Per-fork means (ops/s, 10 iterations each, fork order)

| lane | round | op | params | side | forks 1…8 | max/min | modes |
|---|---:|---|---|---|---|---:|---|
| A | 1 | Dot | n=65536 | gale | 65,526, 65,856, 65,329, 65,397, 66,009, 65,055, 66,162, 65,936 | 1.017 | unimodal |
| A | 1 | Dot | n=65536 | breeze | 67,746, 67,809, 67,734, 67,861, 67,741, 68,162, 68,288, 67,757 | 1.008 | unimodal |
| A | 1 | Dot | n=262144 | gale | 16,384, 16,911, 16,479, 16,977, 16,752, 16,904, 16,624, 16,700 | 1.036 | unimodal |
| A | 1 | Dot | n=262144 | breeze | 16,142, 17,125, 17,126, 17,135, 13,552, 16,168, 17,023, 16,772 | 1.264 | bimodal (1 low / 7 high, gap 15%) |
| A | 1 | Dot | n=1048576 | gale | 4,204, 4,146, 4,124, 4,137, 4,131, 4,060, 4,153, 4,112 | 1.035 | unimodal |
| A | 1 | Dot | n=1048576 | breeze | 4,221, 4,108, 4,289, 4,229, 4,257, 4,275, 4,287, 4,239 | 1.044 | unimodal |
| A | 1 | DotDense | length=100000, nnz=1000 | gale | 2,272,635, 2,270,495, 2,269,654, 2,280,347, 2,270,841, 2,270,141, 2,268,052, 2,280,252 | 1.005 | unimodal |
| A | 1 | DotDense | length=100000, nnz=1000 | breeze | 2,250,629, 2,260,465, 2,258,877, 2,261,938, 2,261,580, 2,248,710, 2,260,557, 2,263,557 | 1.007 | unimodal |
| A | 1 | DotDense | length=100000, nnz=10000 | gale | 224,733, 224,831, 224,494, 224,459, 224,958, 224,467, 225,154, 224,700 | 1.003 | unimodal |
| A | 1 | DotDense | length=100000, nnz=10000 | breeze | 221,273, 224,883, 224,743, 224,748, 225,098, 221,019, 224,073, 221,752 | 1.018 | unimodal |
| A | 2 | Dot | n=65536 | gale | 65,702, 58,718, 65,083, 65,858, 65,976, 66,215, 66,245, 66,820 | 1.138 | bimodal (1 low / 7 high, gap 10%) |
| A | 2 | Dot | n=65536 | breeze | 67,725, 68,394, 67,735, 67,807, 67,918, 67,787, 67,766, 67,755 | 1.010 | unimodal |
| A | 2 | Dot | n=262144 | gale | 16,074, 16,632, 16,903, 16,922, 16,944, 16,892, 16,888, 16,743 | 1.054 | unimodal |
| A | 2 | Dot | n=262144 | breeze | 17,140, 16,778, 17,153, 17,138, 16,458, 17,132, 17,130, 17,134 | 1.042 | unimodal |
| A | 2 | Dot | n=1048576 | gale | 4,222, 4,123, 4,158, 4,222, 4,073, 4,135, 4,194, 4,224 | 1.037 | unimodal |
| A | 2 | Dot | n=1048576 | breeze | 4,244, 4,132, 4,207, 4,180, 4,233, 4,148, 4,212, 4,229 | 1.027 | unimodal |
| A | 2 | DotDense | length=100000, nnz=1000 | gale | 2,272,451, 2,270,689, 2,271,771, 2,278,837, 2,282,383, 2,267,178, 2,271,565, 2,280,325 | 1.007 | unimodal |
| A | 2 | DotDense | length=100000, nnz=1000 | breeze | 2,262,982, 2,260,578, 2,256,799, 2,264,437, 2,264,705, 2,256,582, 2,264,098, 2,251,639 | 1.006 | unimodal |
| A | 2 | DotDense | length=100000, nnz=10000 | gale | 224,920, 225,153, 224,952, 224,947, 224,642, 225,041, 224,837, 223,969 | 1.005 | unimodal |
| A | 2 | DotDense | length=100000, nnz=10000 | breeze | 224,806, 225,009, 224,155, 221,948, 224,826, 224,859, 225,027, 224,663 | 1.014 | unimodal |
| B | 1 | Dot | n=65536 | gale | 65,353, 66,194, 65,904, 65,559, 65,342, 65,123, 66,024, 63,234 | 1.047 | unimodal |
| B | 1 | Dot | n=65536 | breeze | 67,686, 67,798, 67,777, 67,783, 67,712, 67,725, 67,754, 67,724 | 1.002 | unimodal |
| B | 1 | Dot | n=262144 | gale | 16,903, 16,927, 16,927, 16,489, 16,607, 16,661, 16,913, 16,898 | 1.027 | unimodal |
| B | 1 | Dot | n=262144 | breeze | 17,135, 16,482, 17,051, 17,142, 17,139, 17,123, 17,120, 17,135 | 1.040 | unimodal |
| B | 1 | Dot | n=1048576 | gale | 4,216, 4,197, 4,229, 4,121, 4,191, 4,233, 4,245, 4,090 | 1.038 | unimodal |
| B | 1 | Dot | n=1048576 | breeze | 4,260, 4,253, 4,144, 4,250, 4,246, 4,148, 4,185, 4,191 | 1.028 | unimodal |
| B | 1 | DotDense | length=100000, nnz=1000 | gale | 2,270,612, 2,278,509, 2,276,323, 2,279,749, 2,277,483, 2,277,929, 2,279,911, 2,276,000 | 1.004 | unimodal |
| B | 1 | DotDense | length=100000, nnz=1000 | breeze | 2,264,047, 2,265,444, 2,260,961, 2,264,014, 2,250,133, 2,262,591, 2,251,395, 2,257,122 | 1.007 | unimodal |
| B | 1 | DotDense | length=100000, nnz=10000 | gale | 224,514, 224,618, 224,914, 224,567, 224,670, 224,749, 224,611, 224,998 | 1.002 | unimodal |
| B | 1 | DotDense | length=100000, nnz=10000 | breeze | 224,745, 224,815, 224,701, 224,599, 221,327, 224,850, 224,648, 224,049 | 1.016 | unimodal |
| B | 2 | Dot | n=65536 | gale | 62,417, 65,620, 65,075, 65,625, 66,003, 66,468, 66,262, 65,762 | 1.065 | bimodal (1 low / 7 high, gap 4%) |
| B | 2 | Dot | n=65536 | breeze | 67,845, 67,846, 67,733, 67,771, 67,771, 67,834, 67,730, 67,792 | 1.002 | unimodal |
| B | 2 | Dot | n=262144 | gale | 16,700, 16,909, 16,644, 16,873, 16,461, 16,712, 16,926, 16,384 | 1.033 | unimodal |
| B | 2 | Dot | n=262144 | breeze | 17,145, 17,052, 17,127, 16,107, 15,881, 17,131, 17,150, 17,144 | 1.080 | bimodal (2 low / 6 high, gap 6%) |
| B | 2 | Dot | n=1048576 | gale | 4,207, 4,151, 3,883, 4,016, 4,209, 4,118, 4,190, 4,208 | 1.084 | unimodal |
| B | 2 | Dot | n=1048576 | breeze | 4,172, 4,094, 4,192, 4,199, 4,154, 4,205, 4,291, 4,213 | 1.048 | unimodal |
| B | 2 | DotDense | length=100000, nnz=1000 | gale | 2,276,326, 2,276,622, 2,273,916, 2,270,343, 2,269,719, 2,272,212, 2,269,937, 2,275,470 | 1.003 | unimodal |
| B | 2 | DotDense | length=100000, nnz=1000 | breeze | 2,266,752, 2,238,501, 2,260,827, 2,262,102, 2,262,742, 2,249,093, 2,262,767, 2,251,977 | 1.013 | unimodal |
| B | 2 | DotDense | length=100000, nnz=10000 | gale | 224,328, 224,859, 224,745, 224,764, 224,193, 224,698, 225,029, 224,811 | 1.004 | unimodal |
| B | 2 | DotDense | length=100000, nnz=10000 | breeze | 224,686, 224,766, 223,281, 225,018, 224,149, 224,472, 224,783, 224,268 | 1.008 | unimodal |

## Summary: dot 64K is a real kernel gap; DotDense nnz = 10K is noise

**`dot` n = 65536: a small, real kernel gap.** gale is **0.959–0.967× of Breeze in all four runs**, both lanes, with disjoint 99.9% CIs every time: lane A 0.967× and 0.959×, lane B 0.965× and 0.965×. Per fork, **every one of the 32 gale forks is below every one of the 32 Breeze forks**. The best gale fork reaches 0.977–0.987× of the worst Breeze fork, and median-to-median the gap is 0.966–0.973×.

Breeze is extremely stable (max/min ≤ 1.010). gale has one slow fork in two of the four runs: 58,718 ops/s in lane A r2 and 62,417 in lane B r2, against about 65–66K for the rest. Those outliers widen gale's CI but are not the cause of the gap, since gale's typical forks are still about 3% slower. The deficit is the same in lane A (Java11BLAS) and lane B (VectorBLAS); gale's dot is backend-insensitive.

**`dot` n = 262144 and 1048576 (context):** within ±2% of Breeze, with verdicts that flip between runs (0.985–1.021×). These are noise-level. The Breeze side is bimodal at 256K in two runs, with one or two slow forks.

**`DotDense` nnz = 10000: not a gap.** gale/Breeze is **1.001–1.006×** in all four runs: one CI-separated "ahead" and three ties. Both sides are unimodal, and gale is the steadier of the two (max/min ≤ 1.005 against Breeze's ≤ 1.018). The Milestone 3 lane-B loss (0.96×) did not reproduce, so it was run-to-run noise. **`DotDense` nnz = 1000** is 1.006–1.008×, ahead with separated CIs in all four runs.

Raw data: `2026-10-11-dot-rerun-trillium.json` (per-run JMH JSON with per-fork rawData, netlib records, source SHA, load).
