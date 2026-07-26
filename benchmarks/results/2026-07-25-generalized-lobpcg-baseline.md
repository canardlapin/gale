# Matrix-free generalized LOBPCG development baseline

Date: 2026-07-25

Machine: Apple ARM64

Runtime: OpenJDK 25.0.1, pure Gale backend

Harness: JMH 1.37. The committed benchmark defaults to one fork, one 200 ms
warmup, and two 200 ms measurements. The numbers below are a fast development
receipt run with no fork, no warmup, and one 100 ms measurement. They prove the
scenario and allocation instrumentation is live; they are not release-grade
performance or crossover claims.

## Fixed fixtures

`GeneralizedLobpcgJmh` covers the full product:

- `n` in `{128, 512, 2048}`;
- `k` in `{4, 8, 16}`;
- clustered diagonal and finite-difference stiffness/mass pencils;
- identity, Jacobi, and pairwise block-Jacobi preconditioners.

Inputs are deterministic. The clustered pencil uses an analytic diagonal metric
and a repeated, tightly clustered low end. The stiffness/mass pencil uses
Dirichlet-like tridiagonal operators with a variable SPD stiffness diagonal and
a consistent-mass metric. All construction is outside timed invocations.

## Reproduction

```bash
sbt benchCompile
sbt "benchmarksJVM/Jmh/runMain gale.bench.GeneralizedLobpcgWorkReceipt"
sbt "benchmarksJVM/Jmh/run -f 0 -wi 0 -i 1 -r 100ms -prof gc gale.bench.GeneralizedLobpcgJmh.solve"
```

For trustworthy performance numbers, omit the short-run overrides and retain
the committed fork/warmup/measurement settings:

```bash
sbt "benchmarksJVM/Jmh/run -prof gc gale.bench.GeneralizedLobpcgJmh.solve"
```

## Directional time and allocation

The tables compare the useful extremes. Diagonal block-Jacobi reduces exactly
to Jacobi, so its small timing differences are noise. Plain Jacobi is only a
scalar rescaling on the stiffness fixture, so block-Jacobi is the meaningful
comparison there. Allocation is MB/op using decimal MB for readability.

### Clustered diagonal: identity versus Jacobi

| n | k | Identity ms/op | Jacobi ms/op | Identity MB/op | Jacobi MB/op |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 128 | 4 | 36.017 | 2.403 | 9.963 | 2.125 |
| 128 | 8 | 7.019 | 3.009 | 7.315 | 3.521 |
| 128 | 16 | 71.409 | 79.163 | 50.827 | 52.411 |
| 512 | 4 | 12.130 | 5.393 | 14.352 | 7.161 |
| 512 | 8 | 28.151 | 14.703 | 28.899 | 14.421 |
| 512 | 16 | 99.854 | 49.175 | 58.033 | 29.375 |
| 2048 | 4 | 53.654 | 19.932 | 62.422 | 27.764 |
| 2048 | 8 | 128.031 | 57.913 | 119.768 | 55.797 |
| 2048 | 16 | 469.387 | 213.068 | 238.062 | 112.627 |

The first identity cell includes cold non-forked startup and must not be read as
a stable ratio. In the converged `n >= 512` cases, Jacobi halves iterations and
roughly halves allocation. At `n=128, k=16`, both routes exhaust the shared
40-iteration budget on the deliberately ambiguous repeated cluster; time is not
a convergence-equivalent comparison.

### Stiffness/mass: identity versus block-Jacobi

| n | k | Identity ms/op | Block-Jacobi ms/op | Identity MB/op | Block-Jacobi MB/op |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 128 | 4 | 13.401 | 11.763 | 12.545 | 11.081 |
| 128 | 8 | 21.540 | 17.708 | 22.774 | 20.294 |
| 128 | 16 | 66.900 | 58.732 | 46.647 | 42.003 |
| 512 | 4 | 28.642 | 27.106 | 40.690 | 40.644 |
| 512 | 8 | 87.333 | 86.288 | 82.533 | 82.518 |
| 512 | 16 | 316.144 | 304.797 | 171.249 | 169.477 |
| 2048 | 4 | 161.574 | 94.849 | 158.727 | 133.952 |
| 2048 | 8 | 377.217 | 311.172 | 318.413 | 268.579 |
| 2048 | 16 | 1400.277 | 1354.251 | 642.612 | 622.800 |

The allocation scale is a baseline opportunity, not an allocation budget.
LOBPCG currently owns block temporaries at each iteration; fewer iterations and
a later explicit workspace tier are the two distinct ways to reduce this cost.

## Exact work accounting

Each tuple is `converged/requested @ iterations; A/B/P`, where `P` is the number
of residual columns sent through the preconditioner. Zero converged pairs are
honest exhaustion results, not failures; their `worstResidual` is structurally
zero because no pair was returned.

### Clustered diagonal: identity / Jacobi

| n | k | Identity | Jacobi |
| ---: | ---: | --- | --- |
| 128 | 4 | `4/4 @ 13; 103/55/51` | `4/4 @ 7; 53/29/25` |
| 128 | 8 | `8/8 @ 12; 191/103/95` | `8/8 @ 6; 96/56/48` |
| 128 | 16 | `0/16 @ 40; 1204/603/587` | `2/16 @ 40; 1248/649/612` |
| 512 | 4 | `4/4 @ 14; 112/60/56` | `4/4 @ 7; 56/32/28` |
| 512 | 8 | `8/8 @ 14; 224/120/112` | `8/8 @ 7; 112/64/56` |
| 512 | 16 | `16/16 @ 14; 433/225/209` | `16/16 @ 7; 220/124/108` |
| 2048 | 4 | `4/4 @ 15; 120/64/60` | `4/4 @ 7; 56/32/28` |
| 2048 | 8 | `8/8 @ 15; 240/128/120` | `8/8 @ 7; 112/64/56` |
| 2048 | 16 | `16/16 @ 15; 473/249/233` | `16/16 @ 7; 224/128/112` |

### Stiffness/mass: identity / block-Jacobi

| n | k | Identity | Block-Jacobi |
| ---: | ---: | --- | --- |
| 128 | 4 | `0/4 @ 40; 320/164/160` | `3/4 @ 40; 286/141/137` |
| 128 | 8 | `4/8 @ 40; 611/300/292` | `6/8 @ 40; 522/241/233` |
| 128 | 16 | `12/16 @ 40; 1080/493/477` | `15/16 @ 40; 936/400/384` |
| 512 | 4 | `0/4 @ 40; 320/164/160` | `0/4 @ 40; 320/164/160` |
| 512 | 8 | `0/8 @ 40; 640/328/320` | `0/8 @ 40; 640/328/320` |
| 512 | 16 | `0/16 @ 40; 1280/656/640` | `6/16 @ 40; 1262/651/629` |
| 2048 | 4 | `0/4 @ 40; 320/164/160` | `4/4 @ 34; 270/138/134` |
| 2048 | 8 | `0/8 @ 40; 640/328/320` | `8/8 @ 34; 539/275/267` |
| 2048 | 16 | `0/16 @ 40; 1280/656/640` | `0/16 @ 40; 1226/627/611` |

The counters verify that `A`, `B`, and preconditioning are separate work
dimensions. Cached `A X` and `B X` images avoid a hidden final-diagnostic pass
on the pure engine. The backend facade has a separate conformance test proving
that externally supplied raw pairs cost exactly one final `A` and `B`
application per returned pair.

## Interpretation

This receipt supports three bounded conclusions:

1. The requested size/block/preconditioner scenario matrix is executable and
   compile-checked.
2. Preconditioning can materially change iterations, allocation, and time, but
   a weak preconditioner does not guarantee convergence within a fixed budget.
3. Partial or zero convergence must remain visible in every benchmark
   comparison; timing rows are not comparable when the completed work differs.

It does not establish a production-scale backend crossover, a native provider
claim, or a release threshold.
