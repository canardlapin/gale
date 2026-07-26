# Reusable primitive workspaces

Gale keeps ordinary one-shot APIs allocating and simple. Repeated numerical
pipelines can opt into `DenseWorkspace`, a single-owner, grow-only execution
resource that reuses primitive `Double` and index scratch on the JVM and
Scala.js. Workspaces never own result storage and never expose their arrays.

## Checked requirements

`ScratchRequirement` describes simultaneously live primitive cells:

```scala
import gale.linalg.*
import gale.spectral.*

val qr = DenseWorkspace.qrRequirement(rows = 200, cols = 40).toOption.get
val eigen = Eigen.symmetricScratchRequirement(
  order = 40,
  vectors = EigenVectors.ValuesOnly
).toOption.get
```

Counts are non-negative `Int`-addressable quantities. Construction and addition
use `Either[LinAlgError, ScratchRequirement]`, so negative counts and overflow
are reported before allocation. A zero-sized problem reports
`ScratchRequirement.empty`.

Composition names encode lifetimes:

- `a.simultaneous(b)` adds Double counts and index counts component-wise because
  both regions must coexist. Addition is checked for overflow.
- `a.alternative(b)` takes component-wise maxima because only one branch runs.

Both operations are commutative and associative for valid counts; the empty
requirement is their identity. These rules let a larger algorithm report one
requirement without exposing its internal array layout.

## Ownership and lifetime

```scala
val requirement = Eigen
  .symmetricScratchRequirement(128, EigenVectors.ValuesOnly)
  .toOption
  .get
val workspace = DenseWorkspace.forRequirement(requirement)

val first = Eigen.eigSymmetricWith(
  a,
  EigenSelection.All,
  EigenVectors.ValuesOnly,
  workspace
).toOption.get

val second = Eigen.eigSymmetricWith(
  b,
  EigenSelection.All,
  EigenVectors.ValuesOnly,
  workspace
).toOption.get
```

`reserve` grows the Double and index regions independently and never shrinks
them. A request already covered by the current capacities retains both backing
identities. Returned matrices, vectors, factors, and sparse structures always
own their storage: reusing a workspace cannot mutate `first` above.

A workspace is mutable, single-owner execution state. It is not thread-safe and
must not be used concurrently. Give each concurrent worker its own instance.
Keeping one workspace for a long-lived sequential pipeline is safe.

## Supported algorithms

The initial portable reuse tier covers three independently measured families:

- `DenseWorkspace.qrRequirement(rows, cols)` and `DMat.qrWith(workspace)` reuse
  Householder/panel scratch. The ordinary `qr` facade creates a suitable
  workspace internally.
- `Eigen.symmetricScratchRequirement(order, vectors)` and
  `Eigen.eigSymmetricWith(...)` reuse dense symmetric reduction scratch. This is
  deliberately the pure Gale route: an optional provider cannot promise to use
  the caller's workspace. Ordinary `Eigen.eigSymmetric(...)` retains provider
  routing and remains the default facade.
- `CSR.canonicalizeScratchRequirement` and `CSR.canonicalizeWith(workspace)`
  reuse per-row Double/index sorting buffers. The returned canonical CSR owns
  its exact-sized structure and values. Ordinary `canonicalize` creates its
  workspace internally.

The values-only symmetric eigensolver needs `n*n + n` Double cells: the
symmetrized reduction matrix and tridiagonal off-diagonal. When vectors are
requested, the `n*n` eigenvector matrix is result storage rather than scratch,
so the reported requirement is only `n` Double cells. CSR canonicalization
reports `nnz` cells of each primitive kind, the widest possible row.

Use a workspace only after profiling identifies a repeated allocation hot path.
It does not change numerical results, ordering, diagnostics, or error semantics.
