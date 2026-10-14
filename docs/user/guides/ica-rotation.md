# Bounded symmetric tanh rotation

`gale.optim.SymmetricTanhIca` supplies the fixed-point numerical iteration for
an already centered, population-whitened real matrix. Columns are observations;
rows are retained dimensions. Statistical preprocessing, rank selection, training
selection and component interpretation belong to the consuming library.

```scala mdoc
import gale.linalg.Matrix
import gale.optim.{SymmetricTanhIca, TanhIcaPolicy}

// A Cartesian product of two independent, centered, unit-variance sources.
val z = Matrix.tabulate(2, 192) { (row, sample) =>
  val a = (sample % 4 * 2 - 3) / math.sqrt(5.0)
  val b = ((sample / 4) % 3 - 1) * math.sqrt(1.5)
  if row == 0 then math.cos(0.37) * a + math.sin(0.37) * b
  else -math.sin(0.37) * a + math.cos(0.37) * b
}
val initial = Matrix(2, 2)(1.0, 0.1, 0.2, 1.3)
val policy = TanhIcaPolicy(tileSamples = 64)
val capacity = SymmetricTanhIca.inspect(z.rows, z.cols, policy).toOption.get
val rotation = SymmetricTanhIca.fitWhitened(z, initial, capacity, policy).toOption.get
assert(rotation.iterations == 3)
assert(rotation.orthogonalityError < 1e-12)
```

The update uses `tanh(W Z)` and its derivative, then symmetrically decorrelates
the new rotation with the portable symmetric eigensolver. A relative eigenvalue
cutoff of `1e-12` refuses a singular rotation. There is no regularization, rank
reduction, random initialization or hidden restart. Callers supply an immutable
initial matrix; this API performs no RNG work.

Stopping requires the maximum rowwise `abs(abs(dot(W_new, W_old)) - 1)` to meet
the requested tolerance after a complete iteration. Exhaustion returns
`TanhIcaFailure.NonConverged`, and cancellation returns `Cancelled`; neither
contains a reusable rotation. Arithmetic, whitening and orthogonality failures
are distinct from convergence. Convergence does not establish independence or
identify artifacts.

`inspect` is pure and can run before a consumer opens a source. The conservative
primitive-payload capacity includes owned result and portable eigensolver
staging, and excludes caller inputs, objects, allocator retention and RSS.
Projections/nonlinearities are tiled; they are never retained for the complete
observation matrix. Every scalar reduction visits observations in their original
order, independent of tile size. Cancellation is checked at tile boundaries and
before publication. The validation pass checks finite centered observations and
population identity covariance explicitly.

Only the portable numerical route is admitted by this capacity contract.
Provider acceleration and statistical model-family APIs are separate work.
