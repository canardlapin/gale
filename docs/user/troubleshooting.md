# Troubleshoot a Gale workflow

Start from the observable symptom. Keep structural failure, numerical
non-convergence, and runtime configuration separate; they require different
responses.

| Symptom | Likely cause | What to inspect or change |
| --- | --- | --- |
| dependency does not resolve | Gale has no public release | use a pinned source `ProjectRef` or publish the exact checkout locally |
| `DimensionMismatch` or `VectorLengthMismatch` | incompatible shapes | print matrix rows/cols and RHS length or columns at the call boundary |
| `SingularMatrix` | square system has no unique solve under the pivot policy | inspect rank/conditioning; reformulate rather than loosening an unrelated tolerance |
| `NotPositiveDefinite` | Cholesky or metric premise failed | validate the lower triangle and the application's SPD assumption |
| `RankDeficient` from least squares | design lacks full numerical column rank | use pivoted QR diagnostics, remove dependent columns, or define a different model |
| iterative result has `converged = false` | legal iteration budget ended first | inspect residual and iterations; improve scaling/preconditioning or change the budget deliberately |
| partial spectral result contains fewer pairs | only a subset passed residual tolerance | inspect `diagnostics.converged`; call `requireConverged` when partial output is unacceptable |
| `SpectralExtremeNotCertified` | residual passed but global-end membership was not proved | change the start/subspace/engine or accept only residual convergence if scientifically appropriate |
| backend import does not accelerate a call | shape, layout, routine, or library family is outside measured dispatch | inspect the backend dashboard; do not infer a threshold from another machine or BLAS family |
| FFM load or native-access failure | JDK, library discovery, symbols, or native access is missing | use JDK 22+, enable native access, and verify the selected library explicitly |
| a retained CG solution changes | code retained `unsafeSolutionView` | use `workspace.solution` for an owned snapshot |
| builder operation fails after `result()` or `consumeQR` | ownership was already transferred | allocate a new builder for a new lifetime |

## Distinguish `Either` from thrown preconditions

Solves, factorizations that can fail, and spectral facades report structural
failure as `Either[LinAlgError, A]`. Keep the `Either` in reusable code and use
`.orThrow` only where aborting is the intended boundary.

Primitive arithmetic methods such as matrix addition, multiplication, dot
product, and destination writes validate programmer preconditions and may
throw `LinAlgError`. Check shapes before calling them when input dimensions are
not already guaranteed by your program.

## Report a numerical problem usefully

Include the Gale revision, Scala version, runtime, backend import, input shapes,
selection/options, the exact typed error or diagnostics, and a deterministic
small reproducer. For floating-point disagreement, report a residual or
reconstruction error rather than only raw factor entries.

Continue with [Core concepts](core-concepts.md) for the general model or the
[Numerical contract](advanced/numerical-contract.md) for precise guarantees.
