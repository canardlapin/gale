# Ownership and mutable boundaries

Gale's ordinary dense and sparse values are immutable-facing. Once a result
reaches caller code, later work through another public handle does not change
the elements it observes.

## Copies, views, and builders

- Export and ordinary interop functions return independent copies unless their
  names explicitly say otherwise.
- Matrix transpose, row, column, and slice views may share storage with another
  immutable Gale value. Immutable sharing does not permit mutation.
- `DVecBuilder`, `DMatBuilder`, sparse value destinations, and workspaces are
  single-owner mutable resources. They are not safe for concurrent use.
- A builder's `result()` transfers its storage to an immutable value and closes
  the builder. Later reads, writes, or a second result request fail.
- `DMatBuilder.consumeQR(...)` closes the builder and transfers its owned
  row-major storage into QR construction. The returned factor retains no mutable
  alias to the builder.

## Deliberately unsafe borrowed views

An API that permits later mutation to be observed carries an `unsafe` name.
Examples include `CgWorkspace.unsafeSolutionView` and
`unsafeFromBreezeView`. Keep such a view within the mutable owner's lifetime and
do not retain it across workspace reuse or external mutation.

The safe defaults are `CgWorkspace.solution` and `fromBreezeCopy`, which return
stable snapshots.

Dense solve workspaces follow the same rule: `QR.solveLeastSquaresWith` and
`QR.solveLeastSquaresScaledRowsWith` use caller-owned scratch for transformed
right-hand sides, but the returned coefficient vector or matrix owns its
storage. Reusing the workspace cannot change an earlier result.

## Choosing the boundary

Prefer immutable results until profiling demonstrates that allocation matters.
Then keep mutation local:

1. allocate one destination or workspace;
2. reuse it sequentially inside a clearly owned scope;
3. publish only an immutable `result()` or `snapshot()`.

This preserves Gale's ordinary reasoning model while still allowing explicit
allocation control in hot numerical pipelines.
