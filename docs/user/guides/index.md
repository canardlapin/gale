# Choose a task guide

Start from the problem you need to solve:

- [Dense systems and least squares](dense-systems.md) — solve square systems,
  fit one or many responses, inspect pivoted QR, apply row scales, and reuse
  scratch.
- [Sparse matrices and operator solves](sparse-operators.md) — construct CSR,
  run iterative methods, and define a matrix-free operator.
- [Spectral analysis](spectral-analysis.md) — compute dense eigenvalues or SVD
  and interpret partial convergence.
- [Minimum-norm solves and subspaces](subspaces-and-minimum-norm.md) — choose
  a singular-value cutoff, reuse projections, and report null equations.
- [First-order composite optimization](first-order-optimization.md) — choose
  proximal, projected, primal-dual, or exact-reduction capabilities.
- [Nonlinear fitting and bounded optimization](nonlinear-and-bounded-optimization.md) —
  fit residual models with LM or minimize smooth objectives with L-BFGS-B.
- [Derivative-free optimization with CMA-ES](cma-es.md) — evaluate continuous
  black-box objectives with seeded populations, work limits, and optional bounds.
- [General constrained optimization](augmented-lagrangian.md) — minimize smooth
  objectives with nonlinear equalities, inequalities, and optional box bounds.
- [Moving a Breeze workload to Gale](breeze-equivalence.md) — migrate the
  supported real-`Double` linear algebra slice without assuming source
  compatibility.

The [worked-example map](examples.md) summarizes the learning path and routes
to advanced ownership, sparse-structure, and generalized-eigensolver material.

If you have not completed one Gale call yet, begin with
[Get your first result](../getting-started.md). If you know the call but not the
failure, use [Troubleshooting](../troubleshooting.md).
