# First-order composite optimization

`gale.optim` provides portable first-order kernels for objectives assembled by
downstream libraries. The API owns numerical iteration, stopping, typed
operator failures, and result certificates. It does not define statistical
models, penalties, coordinate systems, or application-specific solver plans.

The same implementation runs on the JVM and Scala.js. Variables are `DMat`
values: rows are optimization coordinates and columns are independent right
hand sides or parameter blocks advanced together.

## Choose the method from the objective structure

| Objective structure | Gale method |
| --- | --- |
| smooth term plus a directly proximable term | `proximalGradient` |
| smooth term plus a projection set | `projectedGradient` |
| proximable primal term plus `g(Kx)` | `linearCompositePrimalDual` |
| smooth term plus direct proximal term plus `g(Kx)` | `smoothCompositePrimalDual` |
| exact linear null-space parameterization | `ExactLinearReduction.verify` |

`FirstOrderCapabilities.select` lets a host library separate model lowering
from the methods available in a runtime. Automatic selection is deterministic;
a required but unavailable method is a typed `MissingCapability`.

## Proximal-gradient example

This example minimizes

```text
0.5 ||x - c||² + 0.25 ||x||₁
```

whose independent solution is soft-thresholding `c` by `0.25`.

```scala mdoc
import gale.linalg.DMat
import gale.optim.*

val center = DMat.dense(2, 1, Seq(1.0, -2.0))

val smooth = new SmoothObjective:
  val variableRows = 2
  val lipschitz = 1.0

  def value(at: DMat): Either[FirstOrderError, Double] =
    var squared = 0.0
    var row = 0
    while row < at.rows do
      val difference = at(row, 0) - center(row, 0)
      squared += difference * difference
      row += 1
    Right(0.5 * squared)

  def gradient(at: DMat): Either[FirstOrderError, DMat] =
    Right(DMat.tabulate(at.rows, at.cols): (row, column) =>
      at(row, column) - center(row, column)
    )

val l1 = new ProximalTerm:
  val variableRows = 2

  def value(at: DMat): Either[FirstOrderError, Double] =
    var total = 0.0
    var row = 0
    while row < at.rows do
      total += math.abs(at(row, 0))
      row += 1
    Right(0.25 * total)

  def proximal(at: DMat, step: Double): Either[FirstOrderError, DMat] =
    Right(DMat.tabulate(at.rows, at.cols): (row, column) =>
      val value = at(row, column)
      math.signum(value) * math.max(0.0, math.abs(value) - 0.25 * step)
    )

val fit = FirstOrderSolvers
  .proximalGradient(smooth, l1, DMat.zeros(2, 1))
  .fold(error => throw new IllegalStateException(error.message), identity)

fit.primal(0, 0)
fit.primal(1, 0)
fit.status
```

Objective and proximal callbacks return `Either[FirstOrderError, _]`. A
non-finite callback value, mismatched result shape, invalid configuration, or
failed Gale operator application remains a typed failure.

## Linear operators and norm bounds

Primal-dual methods accept a `BoundedLinearOperator`: a
`DoubleLinearOperator` paired with a finite non-negative upper bound on its
induced norm. Dense `DMat` values already implement `DoubleLinearOperator`;
matrix-free, sparse, block, and composed operators use the same contract.

The bound is evidence supplied by the caller, not estimated by the solver. It
controls the safe primal and dual steps, so an invalid underestimate voids the
method's convergence premise. Construction validates the representation of the
claim, while the model layer remains responsible for how the bound was derived.

Operator application returns `FirstOrderError.OperatorFailure` with the
original `LinAlgError` retained as `cause`. Gale never flattens that failure
into an untyped exception string.

## Read the result and certificate

A successful numerical execution returns `FirstOrderSolution`, including:

- the owned primal value and optional dual value;
- the final objective;
- `Converged` or `IterationLimit`;
- fixed-point residuals, objective change, iteration count, and exact settings;
- value summaries that bind the certificate to the returned matrices.

Iteration exhaustion is a successful execution with
`FirstOrderStoppingStatus.IterationLimit`, not a false convergence claim.
`FirstOrderCertificate.binds` detects substitution of either returned value.
Residuals certify the implemented fixed-point equations at the reported
settings; application-specific optimality or statistical validity remains the
host library's responsibility.

`ExactLinearReduction.verify` separately checks that a proposed basis lies in
the declared constraint null space. Its `LinearReductionCertificate` records
the basis image, constraint image, measured maximum residual, and tolerance
threshold.
