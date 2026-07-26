# Symbolic sparse operation plans

When compressed structure stays fixed while numeric values change, Gale can
analyze that structure once and replay the numeric operation through precomputed
index mappings. The ordinary sparse operators remain the simple choice for
one-off work.

## Union plans

`CSRUnionPlan` merges two canonical patterns into a deterministic canonical
union. Its mappings identify the left and right stored position corresponding
to every output position.

```scala
import gale.sparse.*

val plan = CSRUnionPlan.analyze(a.pattern, b.pattern).toOption.get
val destination = plan.newDestination()

plan.evaluateInto(a, b, destination).toOption.get

// Rebind values without changing either analyzed pattern.
val a2 = plan.leftPattern.bind(nextAValues).toOption.get
val b2 = plan.rightPattern.bind(nextBValues).toOption.get
plan.evaluateInto(a2, b2, destination, leftScale = 0.5, rightScale = 2.0)

val owned = destination.snapshot()
```

Numeric replay computes `leftScale * A + rightScale * B`. A zero scale suppresses
reading that input, including NaN propagation from it. `evaluate` is the
allocating convenience form: it creates one owned numeric result while sharing
the analyzed result pattern.

## Product plans

`CSRProductPlan` computes the Boolean matrix-product pattern. For every output
position it records the exact pairs of left/right stored indices that contribute
to the numeric dot product.

```scala
val product = CSRProductPlan.analyze(a.pattern, b.pattern).toOption.get
val destination = product.newDestination()

product.evaluateInto(a, b, destination).toOption.get
```

`contributionCount` reports the number of stored input pairs retained by the
plan. This mapping can be materially larger than the output pattern; product
plans are intended for repeated replay where avoiding sparse lookup repays that
one-time memory and analysis cost.

## Exact-pattern contract

Analysis is checked and returns `Either[LinAlgError, Plan]`. Union requires equal
shapes; product requires compatible inner dimensions. Both inputs must have
canonical patterns. Result and contribution counts are checked against Gale's
`Int`-addressable storage limit.

Numeric replay accepts values bound to the exact immutable compressed storage
used during analysis. `pattern.bind`, `CSR.rebind`, and pattern-preserving numeric
updates retain that storage. A same-shaped matrix with changed structure—or even
an independently allocated lookalike pattern—is rejected with a typed error in
O(1), before the destination is written. A destination from another plan is
rejected the same way.

This exact identity rule prevents a stale mapping from being applied to newly
ordered indices. To use independently loaded but structurally equal data, bind
its numeric values through the plan's `leftPattern` or `rightPattern` first.

## Destinations, ownership, and explicit zeros

`CSRValuesDestination` owns one mutable numeric array and exposes stored-index
and coordinate access plus allocation-free entry traversal; it never exposes or
adopts backing arrays. It is reusable single-owner execution state and is not
thread-safe. `snapshot()` makes an independent immutable CSR that later replay
cannot mutate.

The symbolic result pattern never changes during numeric evaluation. A
cancellation therefore remains an explicit stored zero, and a product whose dot
sum vanishes retains its Boolean-product position. Call `snapshot().pruneZeros`
or `prune(absBelow)` explicitly when numeric compaction is desired. This is the
same separation between numeric update and structural change used by compressed
patterns elsewhere in Gale.

No numeric workspace is needed during replay: union uses direct output mappings,
and product uses grouped contribution mappings. These APIs do not implement a
sparse direct factorization; provider integration is a separate boundary.
