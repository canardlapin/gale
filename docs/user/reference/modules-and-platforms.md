# Choose a module and platform

Gale is pre-release and source-only. This page describes the live build; it does
not claim that the listed artifact names are available from Maven Central.

## Public library modules

| Build project | Intended artifact | Platform | Choose it when... |
| --- | --- | --- | --- |
| `coreJVM`, `coreJS` | `gale-core` | JVM, Scala.js | you need Gale's dense, sparse, solver, optimization, spectral, sized, and public backend contracts |
| `lawsJVM`, `lawsJS` | `gale-laws` | JVM, Scala.js | a downstream implementation wants Gale's reusable MUnit/ScalaCheck laws |
| `interopBreeze` | `gale-interop-breeze` | JVM | a migration boundary still accepts or returns Breeze values |
| `vectorBackend` | `gale-backend-jvm-vector` | JVM | a JDK 21+ application opts into measured Vector API dispatch |
| `nativeBackend` | `gale-backend-jvm-native` | JVM | a JDK 22+ application explicitly owns FFM-native matrix storage |
| `blasFfmBackend` | `gale-backend-jvm-blas-ffm` | JVM | a JDK 22+ application opts into runtime-discovered BLAS/LAPACK |
| `interopRavelJVM`, `interopRavelJS` | `gale-interop-ravel` | JVM, Scala.js | an explicitly opt-in development boundary needs copy conversion between Ravel and Gale |

`gale-interop-ravel` is excluded from the 0.1 milestone artifact set because
its current `ravel-core` dependency is a development snapshot. Breeze interop
and the JVM acceleration backends are tested but provisional: they are not in
the first Central bundle. The dedicated `interopRavelTest` route remains
available when a matching Ravel checkout or published snapshot is present; do
not infer publication from the project definition alone. See the
[release manifest on GitHub](https://github.com/canardlapin/gale/blob/main/docs/release-manifest.md)
for the complete admitted, provisional, and excluded set.

`parity`, benchmark projects, the documentation site, browser demo, Scala-next
consumer probe, and the root aggregator are build or evidence projects and set
`publish / skip := true` where applicable.

## Runtime and compiler boundary

| Route | Current requirement or evidence boundary |
| --- | --- |
| core JVM and Vector backend | JDK 21 minimum; Vector backend also needs the incubator module |
| FFM storage and BLAS/LAPACK | JDK 22 minimum plus explicit native access |
| Scala.js JavaScript | Node 22 is the required CI runtime; evergreen-browser support is intended but not a browser CI release gate |
| Scala.js WebAssembly | experimental, default-off, and allow-failure |
| Scala compiler | published artifacts are intended to be produced with Scala 3.7.4; consumers need a compatible 3.7.4-or-newer compiler |

The advisory Scala-next lane currently tests 3.8.4 without making it a release
gate. Consult the [0.1 versioning plan](https://github.com/canardlapin/gale/blob/main/docs/release-policy.md) before
turning a current build observation into a future compatibility claim.

## Local API reference

Run:

```sh
sbt docsCheck
```

This builds JVM and Scala.js Scaladoc and executes/renders the guide. Gale does
not currently claim a stable hosted Scaladoc URL. Use exact source links or the
locally generated API until publication configures one.
