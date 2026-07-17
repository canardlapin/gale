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
`GALE_WASM` environment variable to opt in:

```sh
GALE_WASM=1 sbt coreJS/Test/fastLinkJS   # link core's JS tests to Wasm
```

Enabling the toggle switches the Scala.js linker to
`withExperimentalUseWebAssembly(true)` and ES-module output. **Executing** the
Wasm tests (not just linking them) additionally needs a recent Node.js
(20+) whose V8 supports the exception-handling proposal, launched with the
`--experimental-wasm-exnref` flag — configure it on the Scala.js `jsEnv`'s Node
arguments (e.g. via `jsEnvInput` / `NodeJSEnv` args) before `coreJS/test`.

The CI `wasm` job only performs the link check, since running requires that
runtime flag; it is marked allow-failure so the experimental backend never gates
the pipeline.

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
