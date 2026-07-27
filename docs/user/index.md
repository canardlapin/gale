# Gale

Gale is a cross-platform linear algebra library for Scala 3. Its shared dense,
sparse, solver, and spectral APIs run on the JVM and Scala.js, with optional
JVM acceleration kept behind explicit backend modules.

This site describes Gale @VERSION@. Gale is still preparing its first public
release, so the guides currently assume a source checkout or locally published
artifacts rather than promising coordinates that are not yet available.

## Your first ten minutes

1. Follow [Getting started](getting-started.md) to build a matrix, solve a
   system, and compute a symmetric eigendecomposition.
2. Read [Core concepts](core-concepts.md) for Gale's value, error, operator,
   and backend model.
3. Continue with the [task-oriented guides](guides/index.md) when you have a concrete
   problem to solve.

Every code block on the getting-started path is compiled and executed during
the documentation build.

## Choose a path

- **Learn by example:** [Worked examples](guides/examples.md) covers dense
  operations, solvers, factorizations, spectral routines, sparse matrices, and
  the optional sized layer.
- **Move from Breeze:** [Migrating a focused Breeze workload](guides/breeze-equivalence.md)
  explains what Gale replaces, what differs, and where copying or aliasing
  occurs.
- **Solve large operator problems:** [Matrix-free generalized symmetric
  eigensolving](advanced/generalized-operator-eigen.md) introduces LOBPCG and
  generalized block Lanczos without materializing dense operators.
- **Solve composite optimization problems:** [First-order composite
  optimization](guides/first-order-optimization.md) covers proximal, projected,
  primal-dual, and exact-reduction methods with typed certificates.
- **Optimize a measured hot path:** [Advanced topics](advanced/index.md) covers
  numerical guarantees, reusable destinations and workspaces, sparse structure,
  and ownership boundaries.

## Look something up

Use the [reference section](reference/index.md) for stable capability summaries and the
generated Scaladoc for symbol-level signatures. Internal implementation plans,
acceptance audits, and release evidence remain in the repository but are not
part of this user guide.
