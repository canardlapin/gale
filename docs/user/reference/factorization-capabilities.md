# Factorization capabilities

Gale exposes two small capabilities for code that consumes an already-computed
factor without depending on its concrete dense or provider-specific type:

- `ExactSolveFactor` solves square systems for vector or matrix right-hand
  sides. `LU` and `Cholesky` implement it.
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
