# Gale guide

Gale provides real-`Double` linear algebra for Scala 3 on the JVM and
Scala.js. Its ordinary API combines immutable-facing dense and sparse values,
typed numerical failures, matrix-free operators, and diagnostics that qualify
iterative results. Optional allocation-control and JVM acceleration stay
explicit.

This site describes Gale @VERSION@. The first intended public checkpoint is
`v0.1.0-M1`, so installation currently uses a pinned source revision or locally
published artifacts.

## Reach a first result

Start with [Getting started](getting-started.md). It shows the current source
dependency, fits a small regression, and explains when to retain
`Either[LinAlgError, A]` rather than call `.orThrow`.

Then read [Core concepts](core-concepts.md) for the five rules that explain
later behavior: immutable-facing values, explicit shapes, typed failure,
operator contracts, and opt-in execution policy.

Every `mdoc` block in this guide is compiled or executed by `docsCheck` against
the live `gale-core` API.

## Continue with your task

- [Dense systems and least squares](guides/dense-systems.md) — square solves,
  pivoted QR, matrix right-hand sides, row scaling, and reusable scratch.
- [Sparse matrices and operator solves](guides/sparse-operators.md) — CSR,
  iterative methods, and matrix-free application.
- [Spectral analysis](guides/spectral-analysis.md) — dense eigenproblems, SVD,
  partial results, and diagnostics.
- [First-order composite optimization](guides/first-order-optimization.md) —
  proximal, projected, primal-dual, and exact-reduction methods.
- [Moving from Breeze](guides/breeze-equivalence.md) — the supported migration
  boundary, conversions, and deliberate differences.

The [worked-example map](guides/examples.md) routes by problem rather than by
package name.

## Make an informed advanced choice

Use [Advanced topics](advanced/index.md) when you need precise numerical
guarantees, matrix-free generalized eigensolving, reusable destinations or
workspaces, sparse structure reuse, or ownership details.

Use [Reference](reference/index.md) for modules, platforms, factorization
capabilities, and local Scaladoc. If a call fails or does not converge, start
with [Troubleshooting](troubleshooting.md).

Internal implementation plans, acceptance audits, and release receipts remain
versioned in the repository but are intentionally outside this public site.
