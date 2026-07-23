# Immutable vector ownership

`DVec` is an immutable value. Once a `DVec` reaches public caller code through
an unqualified API, later mutation through another public handle cannot change
its observed elements. Zero-copy slices, matrix rows, and matrix columns may
share storage with other immutable values; sharing alone does not weaken this
contract. Deliberate borrowed exceptions must carry an `unsafe` name.

## Mutable conversion boundary

`MutableDVec.toVec` is the only public conversion from a mutable vector to
`DVec`. It returns an independently owned snapshot. The allocation-free
`MutableDVec.asVec` doorway is `private[gale]` and exists only for borrowed
internal reads while Gale retains the mutable owner.

`CgWorkspace.solution` likewise returns a stable snapshot. The intentionally
live, allocation-free exception is named `unsafeSolutionView`; its documentation
states that reusing the workspace changes values observed through an earlier
view.

## Public-return audit

| Return path | Mutable origin | Ownership result |
| --- | --- | --- |
| `DMat`, `CSR`, `CSC`, and `Banded` allocating multiplication | fresh local `MutableDVec`; no destination callback | safe exclusive transfer through Gale's private borrowed view |
| allocating iterative solvers | fresh local mutable recurrence state | safe after return because no mutable owner escapes |
| `DoubleLinearOperator.apply` | destination passed to an arbitrary operator implementation | copied before return so a retained destination cannot mutate the result |
| `Preconditioner.apply` | destination passed to an arbitrary preconditioner implementation | copied before return so a retained destination cannot mutate the result |
| `SparseDirect.solve` for vector right-hand sides | destination passed to an arbitrary provider | copied before return so a retained provider destination cannot mutate the result |
| `CgWorkspace.solution` | persistent caller-owned workspace | copied on every access; stable across workspace reuse |
| `CgWorkspace.unsafeSolutionView` | persistent caller-owned workspace | explicitly borrowed and unstable across workspace reuse |
| `unsafeFromBreezeView` in `gale-interop-breeze` | externally mutable Breeze storage | explicitly unsafe borrowed matrix/vector view; `fromBreezeCopy` is the stable default |
| `DVecBuilder.result` / `DMatBuilder.result` | single-owner builder storage | zero-copy ownership transfer; builder closes permanently |

Regression tests exercise consumer-level inaccessibility of `asVec`, retained
operator/preconditioner/provider destinations, builder closure, and workspace
reuse. New APIs returning `DVec` from mutable storage must either prove exclusive
ownership transfer or take an independent snapshot before returning.
