# Compare two linear fits

Use this comparison when two fitted models have nearly equal residual energies
and subtracting their rounded totals could obscure which fit is better.

`gale.numeric.PairedResidualWorkspace` bounds the difference in squared-residual
energy between two stored row-major linear predictions, using outward rounding
instead of subtracting independently rounded residual energies. The arrays
store each design in row-major order. In this example the previous prediction
is `(1, 0)` and the candidate prediction is `(1, 1)`, exactly matching the
response. Both supplied coefficients minimize their respective residual sums.
The squared-residual energy falls from 1 to 0:

```scala mdoc
import gale.numeric.*

val workspace = PairedResidualWorkspace(2, 1).toOption.get
val comparison = workspace.compare(
  response = Array(1.0, 1.0),
  previousDesign = Array(1.0, 0.0),
  candidateDesign = Array(1.0, 1.0),
  previousCoefficients = Array(1.0),
  candidateCoefficients = Array(1.0)
).toOption.get
assert(comparison.certifiesProfileDecrease)
```

## Interpret the certificate

A profile comparison asks whether the candidate improves on the previous
design even after allowing its coefficients to be refitted.
The default profile check bounds the previous coefficients' possible
suboptimality using the normal-equation residual and a positive Gershgorin
lower bound for the previous Gram matrix. It certifies only a strictly negative
profile-difference upper bound. A nonpositive sufficient Gram bound or an
inconclusive sign returns explicit `Unresolved`; shape, nonfinite input and
unbounded arithmetic return typed errors. Candidate rank need not be assumed.

## Reuse the workspace

A workspace owns `2 * rows + 2 * columns` primitive scratch cells;
`requiredArrayCells` checks that shape before allocation.

Inputs must stay unchanged during the synchronous call; workspaces are not
thread-safe and results retain no input/scratch aliases. Common additive energy
offsets are absent from this API. Priors, optimizer acceptance policy, callback
identity and evaluation quotas remain the caller's responsibility. A data-only
certificate cannot establish decrease of an augmented objective.
