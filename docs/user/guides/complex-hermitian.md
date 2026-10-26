# Factor and solve a complex Hermitian matrix

Use `gale.complex` to store complex matrices and factor small Hermitian
positive-semidefinite (PSD) systems on JVM or Scala.js. A Hermitian matrix
equals its conjugate transpose. The factorization supports square matrices of
dimension 1–128 and returns eigenvalues, eigenvectors, and numerical-rank
information. General nonsymmetric complex decompositions are outside its scope.

## Construct and factor a matrix

`OwnedComplex.copying` takes separate real and imaginary arrays, both in
row-major order. The example represents `[[2, i], [-i, 2]]`, whose eigenvalues
are 1 and 3.

The `budget` argument is a byte limit for the operation's numerical payload.
Here 1 MiB is ample for the small example. It is not a total heap limit: object
overhead, caller-owned inputs, and garbage-collector retention are excluded.

```scala mdoc
import gale.complex.*

val budget = BigInt(1024 * 1024)
val matrix = OwnedComplex.copying(
  2, 2,
  Array(2.0, 0.0, 0.0, 2.0), // real entries
  Array(0.0, 1.0, -1.0, 0.0), // imaginary entries
  budget
).toOption.get
val policy = HermitianPolicy(HermitianCutoff.Relative(1e-10))
val factor = Hermitian.factor(matrix, policy, budget).toOption.get
assert(factor.record.rank == 2)
assert(math.abs(factor.record.eigenvalues.head - 1.0) < 1e-12)
assert(math.abs(factor.record.eigenvalues.last - 3.0) < 1e-12)
```

The cutoff determines which eigenvalues contribute to numerical rank and
rank-truncated operations. Both eigenvalues are retained here. The example
unwraps checked results because its input is fixed; application code should
handle typed failures, including insufficient capacity or indefinite input.

## Solve with the retained factor

The right-hand side below is the column vector `[1, i]`. `FullRank` requires
that every direction was retained by the factorization:

```scala mdoc
val rhs = OwnedComplex.copying(
  2, 1, Array(1.0, 0.0), Array(0.0, 1.0), budget
).toOption.get
val solved = factor.solve(rhs, HermitianSolveMode.FullRank, budget).toOption.get
assert(solved.relativeResidual < 1e-12)
```

Reuse `factor` for additional right-hand sides. If an explicit pseudoinverse is
needed, it uses the same retained eigenvalues:

```scala mdoc
val inverse = factor.pseudoInverse(budget).toOption.get
assert(math.abs(inverse(0, 1).im + 1.0 / 3.0) < 1e-12)
```

## Interpret rank and numerical checks

The input is normalized by its largest complex entry magnitude. Cyclic unitary
Jacobi rotations act directly on the complex matrix; there is no real block
embedding or doubled dimension.
Nonconvergence, cancellation, failed eigen residual/orthogonality certificates,
nonfinite arithmetic and indefinite input return typed failures. Tolerances and
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

## Budget and ownership

Matrix constructors copy and public array accessors return copies. Factors own
their eigenvectors and retain an immutable input. Capacity preflight bounds
primitive payload before allocation: `Hermitian.numericalBytes` includes
normalization, iteration, sorting, certificates and output staging; individual
matrix operations separately check their output capacity. Immutable caller
inputs, object overhead, GC retention and RSS are excluded. Account separately for any buffers your application uses to read or convert data. The implementation is portable Scala
on JVM and Scala.js; accelerated routes and performance are unqualified.

The frozen independent NumPy fixtures include real and imaginary matrices,
full and deficient rank, repeated/zero spectra, solves and projectors. Broader
scientific interpretation, units, observation counts and source geometry belong
to the consumer.
