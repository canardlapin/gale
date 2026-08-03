# Find a worked example

Choose the example whose inputs and result match your problem. Each linked
guide uses executable `mdoc` blocks against the public API; this page does not
duplicate their setup.

| I have... | I need... | Follow... |
| --- | --- | --- |
| a small square dense matrix | one or many direct solutions | [Dense systems and least squares](dense-systems.md#solve-a-square-system) |
| a tall design matrix | regression coefficients and rank evidence | [Dense systems and least squares](dense-systems.md#fit-a-full-rank-design) |
| repeated response matrices | one retained QR and reusable scratch | [Dense systems and least squares](dense-systems.md#reuse-a-factor-and-scratch) |
| row multipliers from an outer algorithm | factor and solve `diag(scales) * A` without materializing it | [Dense systems and least squares](dense-systems.md#apply-row-scales-without-a-temporary-matrix) |
| coordinate triples | a canonical CSR matrix and matrix-vector product | [Sparse matrices and operator solves](sparse-operators.md#build-a-compressed-sparse-matrix) |
| only an `A * x` procedure | an iterative solve without dense storage | [Sparse matrices and operator solves](sparse-operators.md#define-a-matrix-free-operator) |
| a symmetric dense matrix | eigenvalues and residual diagnostics | [Spectral analysis](spectral-analysis.md#solve-a-dense-symmetric-eigenproblem) |
| a rectangular dense matrix | economy SVD factors | [Spectral analysis](spectral-analysis.md#compute-an-economy-svd) |
| symmetric operators `A` and SPD `B` | a few generalized eigenpairs | [Matrix-free generalized symmetric eigensolving](../advanced/generalized-operator-eigen.md) |
| a composite differentiable/proximal objective | a typed first-order method and certificate | [First-order composite optimization](first-order-optimization.md) |
| Breeze values at an old boundary | explicit copy/view conversion and a migration plan | [Moving a Breeze workload to Gale](breeze-equivalence.md) |

For storage reuse after profiling, continue with [Dense destinations](../advanced/dense-destinations.md),
[Reusable workspaces](../advanced/workspaces.md), or
[Sparse patterns and plans](../advanced/sparse-patterns.md). For copy, view,
builder, and borrowed-lifetime rules, read [Ownership](../advanced/ownership.md).

The browser PCA demo under `demo/` is a runnable integration example rather
than the first tutorial. Build it with `sbt demoBuild`, then open
`demo/index.html`; no hosted demo is currently claimed.
