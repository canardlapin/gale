# Core concepts

The same rules appear throughout Gale's API. Dense and sparse values keep their
shape. Ordinary results do not change after a method returns. Operations that
can fail report why. Optional backends may change how Gale computes an answer,
but they do not change what the answer means.

## Dense values and shapes

Use `DVec` for a dense vector of `Double` values and `DMat` for a dense matrix.
Use `CSR` or `CSC` when most entries are zero and you need compressed sparse
storage.

Each value records its dimensions. Gale checks those dimensions before it
multiplies matrices, solves a system, or applies an operator. A shape error is
reported at the operation that caused it; it does not emerge later as a bad
array index.

Matrix literals list their entries by row. That construction rule does not
expose Gale's storage. Code that needs to exchange data with another library
should use Gale's views and interop functions, not assume that a `DMat` contains
a particular JVM array or JavaScript typed array.

## Immutable results

Methods such as `a * b`, `a.solve(b)`, and `a.qr` return values that will not
change behind the caller's back. A transpose, row, column, or slice may share
storage with another immutable Gale value, but no public mutable handle can
alter that storage.

This is the right API for most programs. It is simple to retain an intermediate
result, pass it to another method, or use it from more than one thread.

Builders and workspaces provide a separate route for programs that spend too
much time allocating temporary storage. They are mutable and must have one
owner. A builder's `result()` closes the builder and returns an immutable
value. A workspace may be reused, but results returned from it remain stable.
Use these APIs after measurement shows that allocation matters.

## Failures and diagnostics

A solve can fail because the matrix is singular or because the right-hand side
has the wrong length. Cholesky can fail because the matrix is not positive
definite. An iterative method can stop before it reaches the requested
tolerance.

Gale reports such failures as `Either[LinAlgError, A]`. Library code can inspect
or transform the `Either` without throwing an exception. Tests and small
applications can call `.orThrow` when they want the program to stop
immediately.

Some successful calls also need qualification. Iterative solvers and partial
eigensolvers return diagnostics with residuals, iteration counts, and the
number of results that converged. Check those fields when your program requires
all requested results or a particular accuracy.

## Matrices and linear operators

Use a matrix when you have its entries. Use `DoubleLinearOperator` when you can
compute `A * x` without storing `A`.

The operator form is useful for a discretized differential equation, a large
graph, or a covariance calculation whose matrix would be expensive to build.
Gale's iterative solvers and partial eigensolvers accept operators when their
algorithms need only matrix-vector or matrix-block products.

Gale keeps related operations separate. A preconditioner transforms a residual;
it does not claim to solve a system. A `LinearSolveOperator` solves a system; it
is not merely another matrix-vector product. These types prevent an eigensolver
from inventing an inverse or factorization that the caller did not request.

## Optional backends

If a program imports only `gale-core`, Gale uses its portable implementation on
the JVM and Scala.js. JVM applications may import an optional backend for
Vector API or native BLAS/LAPACK support.

A backend does not take over every operation. Gale routes only the operations,
layouts, and problem sizes that the backend supports. Small or strided
operations may continue to use the portable code. With or without a backend,
Gale keeps the same result types, failures, ownership rules, and diagnostics.

## Next steps

Continue with [worked examples](guides/examples.md) for complete tasks. Read the
[numerical contract](advanced/numerical-contract.md) when you need the precise
rules for tolerances, ordering, sparse storage, or backend selection.
