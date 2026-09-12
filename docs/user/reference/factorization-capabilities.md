# Factorization capabilities

Gale exposes two small capabilities for code that consumes an already-computed
factor without depending on its concrete dense or provider-specific type:

- `ExactSolveFactor` solves square systems for vector or matrix right-hand
  sides. `LU`, `Cholesky` and `BandedCholesky` implement it.
- `LeastSquaresFactor` solves overdetermined systems for vector or matrix
  right-hand sides. `QR` implements it.

Both extend `FactorizationCapability`, which exposes the existing
`FactorizationDiagnostics` unchanged. Backend-produced `LU`, `Cholesky`, and
`QR` values implement the same capabilities because providers return Gale's
typed factor objects.

```scala
import gale.linalg.*

def solveSystem(factor: ExactSolveFactor, rhs: DMat) =
  factor.solve(rhs)

def fit(factor: LeastSquaresFactor, observations: DMat) =
  factor.solveLeastSquares(observations)
```

The distinction is mathematical, not an accuracy claim: "exact solve" means the
square problem `A X = B`, while least squares minimizes `||A X - B||`. Both use
floating-point arithmetic and return `Either[LinAlgError, ...]`.

## Deliberately separate operations

The common capabilities do not expose determinant, inverse, pivot/permutation,
orthogonal-factor, covariance, or residualization operations. Those stay on the
concrete factor that can state their correct contract:

- `LU.det` remains an LU operation.
- `QR.q`, `applyQ`, `applyQT`, `residualize`, and `normalizedCovariance` remain QR
  operations.
- Cholesky's lower factor and tolerance policy remain Cholesky-specific.

Concrete methods are unchanged, so `lu.solve(rhs)`, `cholesky.solve(rhs)`, and
`qr.solveLeastSquares(rhs)` remain the ergonomic surface. The capability types
are for generic algorithms and future sparse/provider factors that genuinely
satisfy the same semantics.

Workspace-aware solves intentionally remain on concrete `QR`. A
`DenseWorkspace` is an execution-resource and pure-kernel choice, not part of
the mathematical `LeastSquaresFactor` capability. Use
`qr.solveLeastSquaresWith(rhs, workspace)` when repeated transformed-right-hand-
side allocation is measured; generic code should keep using
`LeastSquaresFactor.solveLeastSquares`.

Likewise, `DMat.qr(options)` states a pivoting and rank policy and uses the
portable factorization, while the no-argument `DMat.qr` may use an eligible
backend provider. `DMat.qrWith(options, workspace)` is the explicit
allocation-controlled portable route. These names distinguish dispatch and
resource ownership; they are not interchangeable performance spellings.

For a one-off dense solve, `a.solve(rhs)` accepts either a `DVec` or a `DMat`.
The matrix overload factors `a` once and solves every column of the right-hand
side. Retain `a.lu` or `a.cholesky` only when separate calls reuse the same
coefficient matrix.

## Dimensions and failures

`ExactSolveFactor.size` is the square system dimension.
`LeastSquaresFactor.observationCount` and `coefficientCount` describe the design
shape. Mismatched right-hand sides, singular factors, non-positive pivots, rank
deficiency, and unsupported underdetermined solves remain typed `LinAlgError`
values; capability adaptation does not erase or reinterpret diagnostics.

## Packed banded SPD systems

`BandedCholesky.factorLower(bands)` factors a symmetric positive-definite
matrix using only its lower band. The input is `n` rows by `bandwidth + 1`
columns: `bands(i, d) = A(i, i - d)`, so column zero is the diagonal. Entries
with `d > i` are ignored padding and become zero in `factor.lowerBands`.
A nonempty matrix requires `0 <= bandwidth < n`; an empty system uses a
`0 × 1` band. A `gale.sparse.Banded` overload reads that format's lower band
using the same lower-triangle convention. It does not check upper symmetry.

```scala
import gale.linalg.*

val bands = Matrix.tabulate(100, 2)((_, d) => if d == 0 then 2.0 else -1.0)
val result = for
  factor <- BandedCholesky.factorLower(bands)
  solution <- factor.solve(Vec.fill(100)(1.0))
yield (solution, factor.logDet, factor.conditioning)
```

The pure constructor preserves logical input values, including strided views.
`DMatBuilder.consumeBandedCholesky()` transfers a builder's storage without
copying; the builder closes even when factorization fails. Neither route
expands to a dense square matrix. Factorization takes `O(n b²)` operations and
`O(n b)` storage, plus `O(n)` scratch for the matrix norm. Immutable factors
can be reused across independently owned destinations.

`solve` accepts vector and multiple-column matrix RHS values and preserves
them. `solveInPlace` accepts a `MutableDVec` or open `DMatBuilder`; each RHS
costs `O(n b)` operations, with no numerical scratch buffer. `solveTranspose`
and `solveTransposeInPlace` solve the same symmetric system. Use `solveLower`
and `solveLowerTranspose` (and their `InPlace` forms) for `L` and `Lᵀ`.
Dimension errors and nonfinite RHS entries are detected before writes.
Arithmetic overflow returns `LinAlgError.InvalidArgument` and can leave a
partially solved destination. A closed builder follows its existing throwing
ownership contract. Inputs and factors must not be mutated concurrently.

`CholeskyOptions.pivotTolerance` is an absolute threshold on squared diagonal
pivots, as for dense Cholesky; it is not a condition-number threshold. Active
nonfinite band entries and nonpositive pivots return typed failures. `logDet`
sums logarithms of pivots, avoiding determinant products.

The lazy `conditioning` diagnostic uses absolute-value triangular comparison
solves in `O(n b)` time and `O(n)` temporary storage. It reports the symmetric
matrix 1-norm, an inverse 1-norm upper bound, their condition-number upper
bound, and its reciprocal lower bound. These bounds hold in exact arithmetic;
the floating-point calculation is not a directed-rounding certificate.
Cancellation can make them conservative, and overflow gives infinity/zero.
The bound is exact for diagonal matrices. `pivotRatio` separately reports
minimum/maximum squared pivots and must not be read as reciprocal condition.
See [qualification evidence](https://github.com/canardlapin/gale/blob/main/docs/banded-cholesky-evidence.md) for analytic
checks, JVM/JS coverage, independent LAPACK parity and measured scaling.
