# Owned complex Hermitian factors

`gale.complex` is a bounded portable complex tier alongside Gale's existing
real-storage spectral API. It provides owned split-complex matrices, conjugate
transpose, products, and certified PSD Hermitian eigensystems. It is not a
general nonsymmetric complex decomposition API or an admitted native backend.

```scala mdoc
import gale.complex.*

val budget = BigInt(1024 * 1024)
val matrix = OwnedComplex.copying(2, 2,
  Array(2.0, 0.0, 0.0, 2.0), Array(0.0, 1.0, -1.0, 0.0), budget).toOption.get
val policy = HermitianPolicy(HermitianCutoff.Relative(1e-10))
val factor = Hermitian.factor(matrix, policy, budget).toOption.get
assert(factor.record.rank == 2)
assert(math.abs(factor.record.eigenvalues.head - 1.0) < 1e-12)
assert(math.abs(factor.record.eigenvalues.last - 3.0) < 1e-12)
val rhs = OwnedComplex.copying(2, 1, Array(1.0, 0.0), Array(0.0, 1.0), budget).toOption.get
val solved = factor.solve(rhs, HermitianSolveMode.FullRank, budget).toOption.get
assert(solved.relativeResidual < 1e-12)
val inverse = factor.pseudoInverse(budget).toOption.get
assert(math.abs(inverse(0, 1).im + 1.0 / 3.0) < 1e-12)
```

The input is normalized by its largest complex entry magnitude. Cyclic unitary
Jacobi rotations act directly on the complex matrix; there is no real block
embedding or doubled dimension. The first contract admits dimensions 1..128.
Nonconvergence, cancellation, failed eigen residual/orthogonality certificates,
nonfinite arithmetic and indefinite input return typed refusals. Tolerances and
maximum sweeps are explicit. There is no hidden ridge, eigenvalue flooring or
replacement of an unconverged iterate with a result.

Input discrepancies within the declared Hermitian tolerance are projected to
their Hermitian mean, with the maximum relative correction reported. Strict
zero tolerance is available. Eigenvalues below the explicit negative-roundoff
allowance are refused; accepted negative roundoff remains in the rank record.
Retention uses strict `eigenvalue > threshold`, with a caller-selected absolute
or relative cutoff. Zero rank is representable. `psdFactor` returns the retained
factor F such that F F^H equals the retained PSD part. `pseudoInverse` does not
restore discarded directions. `FullRank` solves refuse deficient rank;
`RetainedLeastSquares` reports both actual residual and discarded RHS content,
and checks its residual against the retained projection.

Matrix constructors copy and public array accessors return copies. Factors own
their eigenvectors and retain an immutable input. Capacity preflight bounds
primitive payload before allocation: `Hermitian.numericalBytes` includes
normalization, iteration, sorting, certificates and output staging; individual
matrix operations separately report/admit output capacity. Immutable caller
inputs, object overhead, GC retention and RSS are excluded. A consumer must add
its own unit-conversion or read staging. The implementation is portable Scala
on JVM and Scala.js; accelerated routes and performance are unqualified.

The frozen independent NumPy fixtures include real and imaginary matrices,
full and deficient rank, repeated/zero spectra, solves and projectors. Broader
scientific interpretation, units, observation counts and source geometry belong
to the consumer.
