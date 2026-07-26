# Dense matrix destinations

Gale's ordinary matrix operators remain allocating and immutable. Callers with a
measured repeated-operation hot path can instead allocate a `DMatBuilder` once
and reuse it as a single-owner row-major destination.

## Matrix multiplication

```scala
import gale.linalg.*

val destination = DMatBuilder.zeros(a.rows, b.cols)
a.gemmInto(b, destination, alpha = 1.0, beta = 0.0)

// Reuse the same open destination for another product.
a.gemmInto(b, destination, alpha = 0.5, beta = 1.0)

val result = destination.result() // transfers storage and closes the builder
```

`gemmInto` implements

```text
destination := alpha * A * B + beta * destination
```

- `beta == 0.0` is replacement semantics. The old destination is not read, so
  existing NaN or infinity cannot contaminate the result through `0 * NaN`.
- Any nonzero `beta` is accumulation semantics and reads the previous
  destination.
- Row-major, transposed, and general positive-stride inputs are accepted.
- The same pure, Vector, or native GEMM routing policy used by `A * B` is used by
  `gemmInto`. The `A.t * A` assign-only case retains the SYRK route.

## Fused scale and add

```scala
val destination = DMatBuilder.zeros(a.rows, a.cols)
a.linearCombinationInto(
  b,
  destination,
  alpha = 0.5,
  beta = 0.5
)
```

This is replacement semantics:

```text
destination := alpha * A + beta * B
```

It never reads the old destination. A zero coefficient suppresses reading that
input as well, giving an explicit and testable NaN policy. The operation fuses
the add and scale into one pass and creates no intermediate matrix.

## Safety and lifetime

Destination shape is checked before computation. A closed builder and any
internal alias between a destination and an input are rejected deterministically
with `LinAlgError`; no kernel runs after rejection. Public Gale APIs do not expose
the builder's storage or provide a mutable matrix view, so aliasing is normally
impossible without internal/test-only access.

`DMatBuilder.result()` transfers the destination storage to an immutable `DMat`
without copying and permanently closes the builder. Allocate a new builder for a
new ownership lifetime.

Use this tier only when destination reuse is part of the caller's design. For
one-off operations, `A * B`, `A + B`, and the other high-level methods remain the
clear default.
