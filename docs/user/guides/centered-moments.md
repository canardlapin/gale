# Centered moments with frequency weights

`gale.numeric.CenteredMoments` is a worker-local owned builder for exact finite
binary64 moments. It accepts an observation vector and a nonnegative integral
frequency: literal replication of that observation. It does not imply independent
observations or choose masks, sensor identities, shrinkage or a noise model.

```scala mdoc
import gale.numeric.CenteredMoments
val moments = CenteredMoments.make(2).toOption.get
moments.add(Array(1e16, 0.0))
moments.add(Array(1e16 + 2.0, 2.0))
moments.add(Array(1e16 + 4.0, 4.0))
val sample = moments.finish(degreesRemoved = 1).toOption.get
sample.mean.toSeq
sample.covariance.valuesRowMajor
```

Every finite double is an integer multiple of 2^-1074. Sums use those exact
integer units; product sums use 2^-2148. For frequency mass W, sample/population
centered covariance is `(W*sum(w*x*y)-sum(w*x)*sum(w*y))/(W*(W-ddof))`.
The subtraction operates on exact integers and loses no bits. Means and covariance
entries round once to nearest binary64, ties to even. Constant maximum finite
values therefore have finite means and zero covariance even when unnormalized
floating products would overflow. Genuinely unrepresentable normalized entries
refuse. Subnormal underflow to zero is a rounded representation, not an added
noise term or a promise of preserving every positive eigenvalue.

`finish(0)` uses population normalization and `finish(1)` uses sample normalization.
The frequency mass must exceed the removed degree. Other divisors refuse. The
maximum total frequency is 2^60; zero-frequency inputs contribute no moments but
must still have valid shape and finite values. Inputs are not retained. Validation
and capacity failures leave the prefix unchanged. `copy` has independent mutable
array slots; immutable integer limbs may safely share. `merge` preserves its source,
refuses self-merge and accumulates exact state across arbitrary partitions/trees.
Owned final means/matrices do not alias future state.

Exact real covariance is PSD. Final entrywise Double rounding can produce tiny
negative eigenvalues in a rank-deficient matrix, so downstream code must declare
its symmetry/PSD/rank tolerance and cannot invent a noise floor from this API.

`resources(dimension)` inspects capacities without allocating the builder. Dense
output and triangular shapes must fit Int indexing. Source-level estimates bound
sum/product/readout integer limbs by finite binary64 range and total frequency,
plus staged numerical arrays, temporary operations and owned output. They are
conservative numerical capacities, excluding object/reference/host/allocator/GC
and RSS behavior. Time is O(dimension²) per positive observation with arbitrary-
precision products, and there is no cheap-loop or whole-workload performance claim.

The implementation shares ties-to-even exact ratio readout with `ExactSum`; its
existing behavior and resource model are unchanged. An ephys4s consumer adds common
observation masks, actual weight/count/window receipts, scientific space/units,
insufficiency rules and explicit covariance-transform preconditions.
