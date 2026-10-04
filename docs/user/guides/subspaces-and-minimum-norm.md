# Choose a numerical rank and reuse its geometry

Use SVD for minimum-norm solves of rank-deficient or underdetermined systems.
Use one retained factor when solving many responses or computing both row and
column spaces under the same policy.

```scala mdoc:silent
import gale.linalg.*
import gale.spectral.*
```

## Select singular directions

`SvdCutoff.Relative(tau)` keeps exactly `sigma > tau * sigmaMax`.
`Absolute(value)` keeps `sigma > value`. `Default` uses
`max(rows, columns) * eps * sigmaMax`, preserving the existing `pinv` cutoff.
Explicit values must be finite and non-negative. Equality is discarded; zero
keeps all positive singular values. Zero matrices retain rank zero. Empty input
dimensions and non-finite input entries return typed errors.

```scala mdoc
val design = Matrix(2, 3)(1.0, 1.0, 0.0, 0.0, 0.0, 2.0)
val factor = design.truncatedSvd(SvdCutoff.Relative(1e-8)).orThrow
val fit = factor.solve(Vec(2.0, 4.0)).orThrow
(fit.solution.toSeq, fit.residual.toSeq, fit.rank, fit.cutoff)
```

`factor.solve` accepts vector and matrix right-hand sides; every matrix column
is a response. The residual is formed through the kept left singular vectors,
without multiplying the coefficients by the original design. `factor.pinv`
materializes the inverse. One-shot equivalents are `design.pinv(cutoff)` and
`design.minimumNormLeastSquares(rhs, cutoff)`.

These are the pseudo-inverse and minimum-norm least-squares solution of the
**truncated matrix** `A_k`. Discarding a nonzero singular value can change the
least-squares optimum for the original matrix. Moore-Penrose identities apply
to `A_k`. Retained rank is computed from this cutoff; it does not reuse
`design.svd`'s rank metadata or QR's rank estimate, which have different policies.

## Project without a square projector

```scala mdoc
val rows = factor.rowSpace
val columns = factor.columnSpace
val contrast = Vec(3.0, 1.0, 7.0)
(rows.project(contrast).orThrow.toSeq,
 rows.residual(contrast).orThrow.toSeq,
 rows.distance(contrast).orThrow)
```

`Subspace.rowSpace(matrix, cutoff)` and `columnSpace` are one-shot factories.
The basis has shape `dimension x rank` with orthonormal **columns**. Matrix
projection and residual operations act on each input column. `distance` is the
Euclidean distance of one vector; it makes no statistical estimability decision.
The source's `singularValues`, `sigmaMax`, and resolved `cutoff` are exposed so
callers can apply their own subspace uncertainty and covariance contracts.

Projection costs O(dimension * rank * responses), with no dimension-square
projector. Singular vectors can change sign or rotate within repeated singular
values across spectral backends; compare projections, not basis coordinates.

## Reduce rows and report null equations

```scala mdoc
val redundant = Matrix(1, 3)(1.0, 1.0, 1.0)
val echelon = RowReduction.reducedRowEchelon(redundant, SvdCutoff.Default).orThrow
val aliases = RowReduction.nullSpaceBasis(redundant).orThrow
(echelon.pivotColumns, aliases.basis.valuesRowMajor, aliases.nullity)
```

The cutoff overload first selects numerical row space by SVD, then reduces it
and pads zero rows to the input shape. Its rank matches the cutoff-selected
rank. The `Double` overload instead takes an **absolute elimination pivot
tolerance** in input units; it does not promise SVD rank agreement. Neither
overload prunes normalized coefficients using that input tolerance.

Null-space bases have shape `columns x (columns - rank)`, including the omitted
null directions of wide economy SVDs. `NullSpaceForm.Sparse` returns canonical
row-reduced null equations, transposed into columns: for `[1,1,1]` these are
`[1,0,-1]` and `[0,1,-1]`, the exact-arithmetic RREF convention of the null-space
projector. Sparse denotes interpretable elimination structure, not minimal
support or sparse storage. `Orthonormal` uses twice-reorthogonalized
Gram-Schmidt. No square projector is formed, although a full null basis can
itself require quadratic storage. Null vectors describe the truncated matrix;
their residual against the original matrix can be nonzero under truncation.

Elimination chooses columns in input order and the largest remaining absolute
row pivot, breaking ties by first row. Converting the retained orthonormal row
space uses a roundoff pivot floor `columns * eps`; a failure to retain its
selected rank is a typed error. Near a rank boundary, floating-point backend
differences can affect membership; the cutoff comparison itself is strict and
deterministic for the computed singular values.
