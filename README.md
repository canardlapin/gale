# gale

A cross-platform (JVM + Scala.js) linear algebra library for Scala 3, built
around a small, allocation-conscious kernel layer and an Either-first public API.

`gale` provides dense vectors and matrices, dense factorizations (LU, Cholesky,
QR), standalone triangular solves, iterative and least-squares solvers, and a
family of sparse matrix formats — all sharing one set of strided `Double` kernels
that run identically on the JVM (`Array[Double]`) and in the browser
(`Float64Array`).

## Modules

| Module | Path | What it is |
| --- | --- | --- |
| `gale-core` | `core/` | The library. `crossProject` (JVM + JS), `CrossType.Full`. Shared source in `core/shared`, platform storage/interop in `core/jvm` and `core/js`. |
| `gale-laws` | `laws/` | Reusable, munit/ScalaCheck-backed **law bundles** (`VecLaws`, `MatrixLaws`, `SparseLaws`, `SolverLaws`) expressed against the public API. munit and ScalaCheck are *main* dependencies here — the bundles are library code you can call from your own tests. Depends on `gale-core`. |
| `gale-backend-jvm-vector` | `backend-jvm-vector/` | Optional JDK Vector API GEMM backend with adaptive, measured dispatch. |
| `gale-backend-jvm-native` | `backend-jvm-native/` | Optional JDK 22+ `NativeDMat` storage over FFM `MemorySegment`. |
| `gale-backend-jvm-blas-ffm` | `backend-jvm-blas-ffm/` | Optional JDK 22+ runtime-discovered BLAS/LAPACK backend (Accelerate/OpenBLAS/reference/MKL). |
| benchmarks | `benchmarks/jvm`, `benchmarks/js` | JMH (JVM) and a Scala.js smoke runner. Compile-checked in CI; not published. |

See the [backend dashboard](benchmarks/dashboard.md) for the current conformance,
dispatch, and platform-specific performance evidence.

The [v1 compatibility and artifact policy](docs/release-policy.md) defines the
published module set, supported runtimes, and the boundary of the 1.x stability
promise.

The [numerical, sparse, and backend contract](docs/numerical-contract.md) states
accuracy/determinism guarantees, the sparse v1 boundary, and how to choose an
accelerator without widening the Breeze-equivalence claim.

The [immutable vector ownership contract](docs/immutable-vector-ownership.md)
defines the mutable-to-immutable boundary, the explicitly unsafe workspace view,
and the audited ownership status of every public `DVec` return path derived from
mutable storage.

The current [v1 acceptance audit](docs/v1-acceptance-audit.md) and
[release evidence](docs/release-evidence.md) distinguish locally verified code
readiness from the binary-publication and remote-CI steps still required for a
public release. Gale is licensed under [Apache-2.0](LICENSE).

### Package tour (`gale-core`)

- `gale.linalg` — `Vec`/`DVec`, `MutableVec`/`MutableDVec`, `Matrix`/`DMat`,
  `LinearOperator`, dense factorizations (`DenseDecompositions`, `LU`, `Cholesky`,
  `QR`), `TriangularSolve`, `LinAlgError`.
- `gale.kernel` — `DoubleKernels`, the strided BLAS-style inner loops
  (`private[gale]`).
- `gale.sparse` — `COO`, `CSR`, `CSC`, `Banded`, `Diagonal`, `Identity`, `Zero`,
  `Permutation`, and Matrix Market I/O.
- `gale.solvers` — `cg`, `bicgstab`, `gmres`, `cgnr`, `lsqr`, preconditioners.
- `gale.platform` — platform array abstractions (`DoubleArray`, `IndexArray`),
  distinct per platform.

## Build, test, benchmark

Core and the Vector backend require JDK 21; the finalized FFM modules require
JDK 22+. sbt is pinned by `project/build.properties`. The build targets Scala
3.3.8.

```sh
sbt compileAll        # compile core + laws, JVM and JS
sbt testAll           # test core + laws, JVM and JS
sbt testAllFull       # testAll, then a full-optimizing Scala.js link of the JS
                      # test bundles (a stricter check than fastLink)
sbt benchCompile      # compile the JMH and Scala.js benchmarks
sbt nativeBackendTest blasFfmBackendTest benchFfmCompile  # JDK 22+ native gates
sbt benchSmokeJS      # run the Scala.js benchmark smoke runner (fastOpt)
sbt benchSmokeJSFull  # the same under fullOpt
```

Per-platform, per-module tasks also work directly, e.g. `sbt coreJVM/test`,
`sbt lawsJS/test`, `sbt "coreJS/testOnly gale.sparse.*"`.

## Browser demo

`demo/` is a small Scala.js page that runs a full PCA — seeded synthetic 5-D
data, centering, scatter matrix, `Eigen.eigSymmetric`, top-2 projection — live
in the browser on gale's ordinary public API, with per-run timing and a
"Re-sample" button. Build and open it:

```sh
sbt demo/fastLinkJS   # or the alias: sbt demoBuild
open demo/index.html
```

The demo links with `ModuleKind.NoModule`, so the emitted script
(`demo/target/scala-3.3.8/gale-demo-fastopt/main.js`) loads from a plain
`<script>` tag directly off `file://` — no local server or bundler needed.

## Optional acceleration backends

No import uses Gale's pure JVM/JS kernels. On JDK 21 or 22, opt into the Vector
backend with:

```scala
import gale.backend.jvm.vector.given
val c = a * b
```

On JDK 22+, opt into runtime-discovered native BLAS/LAPACK with:

```scala
import gale.backend.jvm.blas.given
val c = a * b
```

The FFM loader checks `-Dgale.blas.library=/absolute/path` and
`GALE_BLAS_LIBRARY` first, then probes OpenBLAS, the platform BLAS (including
Accelerate), and MKL. Launch user applications with
`--enable-native-access=ALL-UNNAMED`. Only known optimized library families
dispatch automatically; generic/reference BLAS remains direct-callable but
default-disabled until measured. When the selected library exposes the required
Fortran LAPACK symbols, the same backend also advertises `NativeLapack` and
provides typed LU, Cholesky, QR, and symmetric-eigen operations. Defaults are
library-family-specific: the measured Accelerate route enables square GEMM at
`n >= 512` and LU at `n >= 128`, while GEMV, QR, and Cholesky remain explicit
opt-ins because their copy-inclusive behavior lost or was non-monotone.
OpenBLAS/MKL automatic routes remain disabled until equivalent platform sweeps
exist; those libraries are still loadable and direct-callable.

## Continuous integration

`.github/workflows/ci.yml` runs on pushes and PRs to `main`:

- **test** — JVM and JS test suites on Scala 3.3.8 (Temurin 21), as a matrix.
- **scala-next** — `testAll` on the latest Scala 3.8.x, as a non-failing-fast
  (`continue-on-error`) job so a compiler-next regression surfaces without gating
  merges.
- **bench** — `benchCompile`.
- **vector-backend** — Vector tests and JMH compilation on JDK 21 and 22.
- **ffm-blas-backend** — native storage and BLAS/LAPACK conformance plus JMH
  compilation on JDK 22 + OpenBLAS.
- **wasm** — an experimental WebAssembly link check, allow-failure (see below).

## Experimental WebAssembly (Scala.js)

The Scala.js output can target WebAssembly instead of JavaScript. It is **off by
default** — a plain `sbt testAll` produces exactly the JavaScript build. Set the
`GALE_WASM` environment variable to opt in; the build configures ES2022 modules
and Node's required exception-reference flag:

```sh
GALE_WASM=1 sbt coreJS/test              # link and execute core tests as Wasm
GALE_WASM=1 sbt benchSmokeJSFull         # execute the Wasm kernel profile
```

The current Node 24 profile is a correctness success but a performance failure:
Wasm is 23–43x slower than optimized JavaScript on the selected dense kernels.
It therefore remains experimental and default-off. The CI `wasm` job executes
both the core tests and the profile, but remains allow-failure while the Scala.js
backend and V8 support mature. See the
[Wasm profile receipt](benchmarks/results/2026-07-17-wasm-profile.md).

## Breeze replacement scope

Gale targets the dense real-`Double`, sparse-matvec, and selected spectral slice
that numerical Scala libraries commonly used from Breeze. It is not a Breeze
fork or a source-compatible replacement. Cross-library tests cover dense BLAS
operations, solves and factorizations, rank/condition overlap, symmetric eigen,
partial SVD, sparse/banded matvec, conversions, aliasing, and typed failure cases.

See [`docs/breeze-equivalence.md`](docs/breeze-equivalence.md) for the capability
matrix, measured performance statement, migration examples, and explicit
non-equivalence boundary.

## `vendor/breeze`

`vendor/breeze` is an **unreferenced reference checkout** of the Breeze library,
kept for consulting its numerics while developing `gale`. It is not part of any
sbt source root, is not compiled, and is excluded from packaging and publishing —
no `gale` artifact depends on it. It can be deleted without affecting the build.
