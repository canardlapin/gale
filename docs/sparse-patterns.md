# Compressed sparse patterns

`CSRPattern` and `CSCPattern` represent compressed structure independently from
numeric values. They are immutable values available on the JVM and Scala.js.
Offsets and indices remain private; callers use dimensions, `nnz`, equality,
lane queries, and allocation-free traversal.

## Checked construction and ownership

```scala
import gale.sparse.CSRPattern

val pattern = CSRPattern.checked(
  rows = 3,
  cols = 4,
  rowOffsets = Array(0, 2, 3, 5),
  columnIndices = Array(0, 2, 1, 0, 3)
).toOption.get
```

Construction checks non-negative dimensions, exact offset length, a zero first
offset, monotone in-range lanes, a final offset equal to `nnz`, and every minor
index. Input arrays are defensively copied. Mutating them after construction
cannot affect the pattern.

`pattern.hasCanonicalFormat` means minor indices are strictly increasing within
each lane. Checked construction accepts safe non-canonical ordering and repeated
indices so imported compressed data can be represented honestly; algorithms that
require a canonical pattern must test this property.

## Numeric rebinding

Binding a `DVec` or `Array[Double]` validates `nnz` and makes exactly one owned
copy of the logical values:

```scala
val first = pattern.bind(Array(1.0, 2.0, 3.0, 4.0, 5.0))
val second = pattern.bind(Array(5.0, 4.0, 3.0, 2.0, 1.0))
```

Both matrices share immutable compressed structure while owning independent
numeric storage. For direct construction, a single-owner builder avoids the
copy:

```scala
val values = pattern.valuesBuilder()
values(0) = 1.0
values(1) = 2.0
// fill remaining stored positions
val matrix = values.result()
```

`result()` transfers the builder's storage and permanently closes the builder;
later reads, writes, or a second `result()` fail with a typed `LinAlgError`.

Every `CSR` and `CSC` exposes its immutable `.pattern` and a checked `.rebind`.
`CSRPattern.t` produces a `CSCPattern` for the transposed shape (and vice versa)
without copying structure. Matrix transpose shares numeric storage as before.

## Explicit zeros and structural change

Numeric-only operations preserve the pattern. In particular, `mapValues`,
scalar multiplication, `rebind`, and values-builder construction retain stored
positions even when a value becomes exactly zero. This makes a fixed pattern
safe to reuse for graph weights, iterative algorithms, and later symbolic
plans.

`pruneZeros`, `prune(absBelow)`, and `canonicalize` are the explicit
structure-changing operations. A matrix containing explicit zeros is not in the
numeric `hasCanonicalFormat`, even when its `.pattern.hasCanonicalFormat` is
true.

## Traversal without storage exposure

```scala
pattern.foreachStoredPosition { (row, col) =>
  // deterministic CSR storage order
}

pattern.foreachRow(0) { (column, storedIndex) =>
  // storedIndex addresses the matching numeric builder cell
}
```

CSC provides the dual `foreachColumn`, `columnNnz`, and `rowIndexAt` surface.
No public API adopts caller arrays, exposes backing arrays, or returns a mutable
view of immutable sparse storage.
