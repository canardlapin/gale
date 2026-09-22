# Understand Gale's core model

Five rules explain most of Gale's public API. Learn them once, then choose an
algorithm from the task guides.

## 1. Values are immutable-facing

`DVec` and `DMat` are dense `Double` values. Ordinary arithmetic, solves, and
factorizations return results that will not change behind the caller's back.

A transpose, row, column, or contiguous slice can be an `O(1)` view sharing
immutable storage. Gathers, exporters, and ordinary interop return copies.
Because the public core exposes no mutable alias to Gale-owned storage,
immutable sharing remains safe.

Builders and workspaces are a separate, single-owner tier. A builder transfers
its storage and closes. A workspace can be reused sequentially, while ordinary
results returned from workspace APIs still own their storage. An API that
deliberately exposes later mutation carries `unsafe` in its name.

## 2. Shape is checked where it matters

Every value records its dimensions. A matrix product, solve, factorization, or
operator application checks the shape required by that operation.

Primitive arithmetic methods assume their preconditions and throw a typed
`LinAlgError` when shapes do not match. Total numerical entry points return
`Either[LinAlgError, A]` when singularity, rank, definiteness, unsupported
selection, or another structural condition is part of the method's normal
failure model.

## 3. Success and convergence are different questions

A direct solve can return `Left(SingularMatrix(...))`. An iterative solve can
finish its legal iteration budget and still return the best available iterate
with `converged = false`. Partial eigensolvers may return only the pairs whose
residuals passed the requested tolerance.

Inspect the result's diagnostics when the application requires all requested
components, a residual threshold, orthogonality, or certification that a pair
belongs to the requested global spectral end. A successful return is not a
blanket accuracy claim.

## 4. Operators represent action, not storage

Use `DMat`, `CSR`, or `CSC` when entries and storage structure matter. Use
`DoubleLinearOperator` when code can compute `A * x` without materializing
`A`.

Algorithms state the extra evidence they need. LSQR and partial SVD need a
transpose action. A generalized symmetric eigensolver needs explicit symmetry
and positive-definite metric evidence. A preconditioner transforms a residual;
a `LinearSolveOperator` represents an actual system solve. Gale does not invent
one capability from another.

## 5. Execution policy is explicit

With only `gale-core`, Gale uses its portable implementation on the JVM and
Scala.js. An optional backend import can accelerate eligible JVM operations,
but it does not take over every shape or layout.

Explicit numerical policy and allocation control can pin a portable route. For
example, `a.qr(options)` preserves one pivot/rank policy across platforms, and
`a.qrWith(workspace)` promises caller-owned scratch rather than provider
routing. Result types, ownership, and documented numerical meaning remain the
same.

## What changes and what is preserved

| Operation | Values and shape | Storage or ownership | Failure or qualification |
| --- | --- | --- | --- |
| `a.t`, `a.row`, `a.col`, `a.slice` | logical view with derived shape | may share immutable storage | invalid index or slice throws `LinAlgError` |
| `a * b`, `a + b`, `a * scalar` | new dense result | owned result | shape precondition throws `LinAlgError` |
| `a.solve(rhs)` | solution with RHS column count | owned result | structural failure is `Left` |
| `a.qr(options)` | factor of `A(:, permutation)` when pivoted | factor owns storage | rank is diagnostic; least-squares can reject rank deficiency |
| iterative solve | last iterate | result owns iterate | inspect convergence and residual |
| workspace call | same mathematical result | scratch is borrowed; result is owned | workspace is sequential and grow-only |

## Next steps

Choose a task from [Guides](guides/index.md). Read the
[numerical contract](advanced/numerical-contract.md) before depending on exact
ordering, tolerance, sparse canonicalization, or backend behavior. Read
[Ownership](advanced/ownership.md) before retaining a borrowed view.
