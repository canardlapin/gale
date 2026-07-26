# Generalized operator eigensolver comparison

Date: 2026-07-25

Machine: Apple ARM64, macOS 14.3

Runtime: OpenJDK 25.0.1, pure Gale backend

Harness: JMH 1.37. The runtime/allocation numbers below are a short development
receipt with no fork, no warmup, and one 100 ms measurement. They verify the
regression matrix and provide directional evidence; they are not
release-quality crossover thresholds.

## Question and decision

This comparison asks whether generalized block Lanczos should replace LOBPCG as
Gale's default matrix-free solver, or whether it provides a useful explicit
alternative when a caller already has a metric solve.

The evidence supports the latter:

- `Eigen.eigSymmetricGeneralized` remains the LOBPCG route;
- `Eigen.eigSymmetricGeneralizedLanczos` remains an explicit named route;
- Gale does not automatically select an engine from dimensions or operator
  types.

Preconditioned LOBPCG was the most robust route across the matrix. Exact-solve
Lanczos was useful on some clustered or small-wide problems, but showed no
consistent allocation win in convergence-equivalent cells and was less robust
on the stiffness/mass pencil. Iterative metric solves add a material inner-work
tax. Those results do not justify implicit routing.

## Fixed scenario matrix

`GeneralizedEigenComparisonJmh` covers:

- `n` in `{128, 512, 2048}`;
- `k` in `{4, 16}`;
- clustered diagonal, ill-conditioned diagonal, and sparse-structured
  matrix-free finite-difference stiffness/mass pencils;
- LOBPCG with identity or problem-aware preconditioning;
- Lanczos with a reusable exact tridiagonal/diagonal metric solve or
  Jacobi-preconditioned CG.

The exact solve is constructed outside the timed invocation and uses O(n)
storage and work. The iterative solve is also constructed once. Inputs and
initial blocks are deterministic.

## Reproduction

```bash
sbt benchCompile
sbt "benchmarksJVM/Jmh/runMain gale.bench.GeneralizedEigenComparisonWorkReceipt"
sbt "benchmarksJVM/Jmh/run -f 0 -wi 0 -i 1 -r 100ms -prof gc gale.bench.GeneralizedEigenComparisonJmh.solve"
```

For a stable performance study, omit the short-run overrides and retain the
committed fork/warmup/measurement settings:

```bash
sbt "benchmarksJVM/Jmh/run -prof gc gale.bench.GeneralizedEigenComparisonJmh.solve"
```

## Convergence receipt

Each cell reports pairs converged within the configured budget. The columns are
LOBPCG identity (`LI`), LOBPCG problem-aware preconditioned (`LP`), Lanczos
exact solve (`LE`), and Lanczos iterative solve (`LCG`). A zero is an honest
exhaustion result, not a failure.

| n | k | Pencil | LI | LP | LE | LCG |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 128 | 4 | clustered | 4 | 4 | 4 | 4 |
| 128 | 4 | ill-conditioned | 0 | 0 | 0 | 0 |
| 128 | 4 | stiffness/mass | 0 | 3 | 0 | 0 |
| 128 | 16 | clustered | 0 | 2 | 16 | 4 |
| 128 | 16 | ill-conditioned | 0 | 4 | 14 | 16 |
| 128 | 16 | stiffness/mass | 12 | 15 | 0 | 0 |
| 512 | 4 | clustered | 4 | 4 | 4 | 4 |
| 512 | 4 | ill-conditioned | 0 | 0 | 0 | 0 |
| 512 | 4 | stiffness/mass | 0 | 0 | 0 | 0 |
| 512 | 16 | clustered | 16 | 16 | 0 | 0 |
| 512 | 16 | ill-conditioned | 0 | 0 | 0 | 0 |
| 512 | 16 | stiffness/mass | 0 | 6 | 0 | 0 |
| 2048 | 4 | clustered | 4 | 4 | 4 | 4 |
| 2048 | 4 | ill-conditioned | 0 | 0 | 0 | 0 |
| 2048 | 4 | stiffness/mass | 0 | 4 | 0 | 0 |
| 2048 | 16 | clustered | 16 | 16 | 16 | 1 |
| 2048 | 16 | ill-conditioned | 0 | 0 | 0 | 0 |
| 2048 | 16 | stiffness/mass | 0 | 0 | 0 | 0 |

The matrix deliberately includes difficult cells. Runtime from differently
converged cells is not treated as an engine comparison.

## Exact work accounting

Each tuple is
`converged/requested @ outer iterations; A/B/P; solves/inner iterations/inner applications`.
`P` is outer-preconditioner applications. LOBPCG has no inner solve; the exact
metric solve reports solve calls but zero iterative work.

| n | k | Pencil | Engine | Exact work |
| ---: | ---: | --- | --- | --- |
| 128 | 4 | clustered | LI | `4/4 @ 13; 103/55/51; 0/0/0` |
| 128 | 4 | clustered | LP | `4/4 @ 7; 53/29/25; 0/0/0` |
| 128 | 4 | clustered | LE | `4/4 @ 3; 112/112/0; 48/0/0` |
| 128 | 4 | clustered | LCG | `4/4 @ 3; 112/208/0; 48/48/96` |
| 128 | 16 | clustered | LP | `2/16 @ 40; 1248/649/612; 0/0/0` |
| 128 | 16 | clustered | LE | `16/16 @ 7; 800/807/0; 337/0/0` |
| 512 | 16 | clustered | LP | `16/16 @ 7; 220/124/108; 0/0/0` |
| 512 | 16 | clustered | LE | `0/16 @ 20; 2256/2256/0; 960/0/0` |
| 2048 | 4 | clustered | LP | `4/4 @ 7; 56/32/28; 0/0/0` |
| 2048 | 4 | clustered | LE | `4/4 @ 3; 112/112/0; 48/0/0` |
| 2048 | 4 | clustered | LCG | `4/4 @ 3; 112/208/0; 48/48/96` |
| 2048 | 4 | stiffness/mass | LP | `4/4 @ 34; 270/138/134; 0/0/0` |
| 2048 | 4 | stiffness/mass | LE | `0/4 @ 20; 724/724/0; 320/0/0` |
| 2048 | 4 | stiffness/mass | LCG | `0/4 @ 20; 724/6491/0; 320/5447/5767` |
| 2048 | 16 | clustered | LP | `16/16 @ 7; 224/128/112; 0/0/0` |
| 2048 | 16 | clustered | LE | `16/16 @ 18; 2032/2110/0; 894/0/0` |

Lanczos's fewer outer iterations for `k = 4` do not mean fewer operator actions:
building and stabilizing its retained basis costs more `A` and `B` work than
preconditioned LOBPCG. On the `n = 2048`, `k = 4` stiffness cell, iterative
metric solves increase metric applications from 724 to 6491 without producing a
converged pair.

## Directional runtime and allocation

The following clustered cells completed the requested work for every listed
engine. Allocation is decimal MB/op.

| n | k | Engine | ms/op | MB/op |
| ---: | ---: | --- | ---: | ---: |
| 128 | 4 | LP | 1.219 | 1.856 |
| 128 | 4 | LE | 1.981 | 2.240 |
| 128 | 4 | LCG | 1.831 | 2.561 |
| 512 | 4 | LP | 4.536 | 7.101 |
| 512 | 4 | LE | 5.802 | 7.810 |
| 512 | 4 | LCG | 5.539 | 8.834 |
| 2048 | 4 | LP | 17.887 | 27.765 |
| 2048 | 4 | LE | 25.759 | 30.380 |
| 2048 | 4 | LCG | 27.920 | 35.537 |
| 2048 | 16 | LP | 210.694 | 112.720 |
| 2048 | 16 | LE | 1184.729 | 715.431 |

In these convergence-equivalent cells, Lanczos did not establish a lower-memory
crossover. Its explicit route is still useful: at `n = 128`, `k = 16`, exact
Lanczos completed all 16 clustered pairs where the two LOBPCG configurations
completed zero and two. Conversely, at `n = 512`, `k = 16`, preconditioned
LOBPCG completed all pairs and both Lanczos routes exhausted. These
complementary outcomes are the reason to expose both engines without silently
choosing between them.

## Regression findings and evidence boundary

The comparison found two basis-management defects that now have direct
regression coverage: cancellation in near-dependent B-orthogonalization could
be misclassified as an indefinite metric, and a final block could exceed a
non-multiple subspace cap. The Lanczos engine now refreshes suspicious metric
images, stabilizes the complete basis before Rayleigh-Ritz, and respects the
exact cap.

This receipt establishes portable correctness/work behavior for the committed
fixtures. It does not establish a production-scale native/provider crossover,
an automatic-routing threshold, or a release performance promise.
