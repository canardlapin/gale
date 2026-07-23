# Gale v1 compatibility and artifact policy

This policy begins with the first `1.0.0` release. The current build version is
`1.0.0-SNAPSHOT`; a snapshot is not a published release or a compatibility
promise.

## Coordinates and published modules

The organization is `io.gale` and the Scala binary suffix is `_3`. The intended
v1 artifact set is:

| Artifact | Platforms | Purpose |
| --- | --- | --- |
| `gale-core` | JVM, Scala.js | dense, sparse, solver, spectral, kernel, and typed public APIs |
| `gale-laws` | JVM, Scala.js | reusable MUnit/ScalaCheck law and backend-conformance bundles |
| `gale-interop-breeze` | JVM | explicit Breeze 2.1 conversion and migration helpers |
| `gale-backend-jvm-vector` | JVM | opt-in JDK Vector API acceleration |
| `gale-backend-jvm-native` | JVM | explicit JDK FFM native matrix storage |
| `gale-backend-jvm-blas-ffm` | JVM | opt-in runtime-discovered BLAS/LAPACK provider |

`gale-parity`, all benchmark projects, and the root aggregator are
`publish / skip := true`. The packages `gale.kernel`, `gale.linalg`,
`gale.sparse`, `gale.solvers`, and `gale.spectral` are intentionally shipped in
`gale-core`; v1 does not create empty fine-grained artifacts for those packages.
`gale-laws` is a normal published module because downstream libraries are meant
to extend its suites, not a classifier-only test jar.

The JVM/Scala.js modules are released from the JDK 21 build. Native storage and
BLAS/LAPACK artifacts require a separate JDK 22 publication pass. A release is
incomplete unless both passes publish the same version.

## Compatibility promise

For the `1.x` line, Gale follows semantic versioning:

- patch releases preserve public source and binary compatibility and only make
  compatible additions or fixes;
- minor releases preserve public binary compatibility and may add APIs;
- removals, incompatible signature changes, or semantic contract breaks require
  a new major version;
- types or members documented `private[gale]`, `private`, experimental, or
  internal are outside the compatibility promise;
- experimental Wasm configuration and incubating JDK Vector implementation
  details are not stable APIs, but the public `Backend` contract and explicit
  backend import points are.

Scala 3 TASTy and compiler compatibility still constrain consumers. Gale
publishes for the Scala 3 binary line (`_3`) from Scala 3.3.8 and CI checks the
current Scala Next line as advisory evidence. Cross-building with Scala Next is
not a promise that every newer compiler can consume every older TASTy artifact.

There is no MiMa gate for `1.0.0`: it has no earlier stable baseline. Before the
first `1.1.0` release, the build must add an automated binary-compatibility check
against the latest `1.0.x` artifact. Until that gate exists, maintainers must
review exported API diffs manually and must not describe a candidate as
binary-verified.

## Supported runtimes

| Module/route | Minimum supported runtime | Notes |
| --- | --- | --- |
| `gale-core` JVM, laws, interop | JDK 21 | primary required CI line |
| Vector backend | JDK 21 | tested on JDK 21 and 22; requires `jdk.incubator.vector` at compile/run time |
| Native and BLAS/LAPACK FFM | JDK 22 | finalized FFM API; applications enable native access |
| Scala.js JavaScript | Node 22 and current evergreen browsers | optimized JS is the supported browser performance route |
| Scala.js Wasm | Node 22+ experimental profile | explicit, allow-failure, and not covered by v1 compatibility/performance promises |

“Current evergreen browsers” means the latest two stable major releases of
Chrome, Firefox, Safari, and Edge at the time a Gale release is cut. Browser CI
is not yet a release gate, so v1 release notes must state that Node is the tested
Scala.js runtime and must not overstate browser certification.

## Release blockers and provenance

A public release requires all of the following in addition to green tests:

- a repository license chosen by the owner and represented in both a root
  `LICENSE` file and published POM metadata;
- project homepage, SCM, and developer metadata in the POM;
- a configured publishing destination and credentials outside the repository;
- clean source, binary, and documentation artifacts for every intended module;
- the acceptance and release-evidence records for the exact release commit.

The owner selected Apache-2.0 and `https://github.com/canardlapin/gale` as the
canonical SCM repository on 2026-07-19. The root `LICENSE` and generated POMs
carry that provenance. Binary publication still requires a configured
destination, credentials, signing policy, and remote CI on the exact release
commit. Local `publishLocal` is useful packaging evidence but is not a public
release.
