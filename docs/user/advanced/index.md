# Advanced topics

These pages assume you have completed [Getting started](../getting-started.md)
and understand [Core concepts](../core-concepts.md). They are for users who
need to make a numerical, performance, storage, or ownership decision:

- [Numerical, sparse, and backend contract](numerical-contract.md) states the
  guarantees shared by portable and accelerated routes.
- [Matrix-free generalized symmetric eigensolving](generalized-operator-eigen.md)
  covers LOBPCG, generalized block Lanczos, metric solves, and convergence
  evidence.
- [Dense destinations](dense-destinations.md) and
  [reusable workspaces](workspaces.md) reduce allocation in measured repeated
  pipelines.
- [Compressed sparse patterns](sparse-patterns.md) and
  [symbolic sparse plans](sparse-plans.md) separate reusable structure from
  changing values.
- [Ownership and mutable boundaries](ownership.md) explains when results copy,
  share immutable storage, or deliberately expose an unsafe borrowed view.

Implementation plans, backend design records, acceptance audits, and raw
release evidence are intentionally not part of this section.
